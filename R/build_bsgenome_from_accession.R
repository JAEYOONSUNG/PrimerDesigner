#' Build BSgenome Package from GenBank Accession
#'
#' Builds a BSgenome package for a given GenBank accession number using a previously downloaded .fna file.
#'
#' @param genbank_accession A character string specifying the GenBank accession number (e.g., "GCF_030376745.1").
#' @param dir A character string specifying the working directory (default: current working directory).
#' @return A character string with the name of the generated BSgenome package.
#' @examples
#' \dontrun{
#' pkg_name <- build_bsgenome_from_accession(
#'   genbank_accession = "GCF_030376745.1",
#'   dir = getwd()
#' )
#' }
#' @export

build_bsgenome_from_accession <- function(genbank_accession, dir = base::getwd()) {
  # Normalize directory path
  dir <- base::normalizePath(dir, mustWork = TRUE)
  base::cat("Building BSgenome package for accession:", genbank_accession, "\n")

  # Define paths
  seqs_dir <- base::file.path(dir, "seqs_srcdir")
  fasta_file <- base::file.path(seqs_dir, base::paste0(genbank_accession, ".fna"))

  # Check if FASTA file exists
  if (!base::file.exists(fasta_file)) {
    base::stop("FASTA file not found: ", fasta_file, "\nPlease run download_genbank_fna first.")
  }
  base::cat("Using existing FASTA file:", fasta_file, "\n")

  # Check if BSgenome package already exists
  clean_accession <- base::gsub("[^a-zA-Z0-9]", "", genbank_accession)
  pkg_name <- base::paste0("BSgenome.", clean_accession)
  if (requireNamespace(pkg_name, quietly = TRUE)) {
    base::cat("BSgenome package already exists:", pkg_name, "\n")
    return(pkg_name)
  }

  # Get metadata (avoid redundant download)
  metadata <- base::list(
    fna_file = fasta_file,
    organism_name = "Unknown",
    release_date = base::format(base::Sys.Date(), "%Y/%m/%d"),
    seqs_dir = seqs_dir
  )
  base::tryCatch({
    search_result <- rentrez::entrez_search(
      db = "assembly",
      term = genbank_accession,
      retmode = "xml",
      retmax = 1
    )
    if (base::length(search_result$ids) > 0) {
      summary <- rentrez::entrez_summary(db = "assembly", id = search_result$ids[1])
      metadata$organism_name <- base::ifelse(base::is.null(summary$speciesname), "Unknown", summary$speciesname)
      metadata$release_date <- base::ifelse(base::is.null(summary$lastupdatedate), base::format(base::Sys.Date(), "%Y/%m/%d"), summary$lastupdatedate)
    }
  }, error = function(e) {
    base::cat("Warning: Failed to fetch metadata:", base::conditionMessage(e), "\nUsing default metadata.\n")
  })

  # ===== UCSC Suppression (MUST be before any Bioconductor operations) =====
  # Non-model organisms (NCBI accessions) trigger GenomeInfoDb UCSC validation
  # that fails. Suppress BEFORE forge_BSgenome and devtools::install.
  base::cat("Suppressing UCSC lookups for non-model organism...\n")
  utils::flush.console()
  gi_ns <- base::asNamespace("GenomeInfoDb")

  .safe_ucsc_fn <- function(genome,
                             goldenPath.url = base::getOption("UCSC.goldenPath.url"),
                             assembled.molecules.only = FALSE,
                             map.NCBI = FALSE,
                             add.ensembl.col = FALSE,
                             recache = FALSE,
                             as.Seqinfo = FALSE) {
    base::data.frame(
      chrom = base::character(0), size = base::integer(0),
      assembled = base::logical(0), circular = base::logical(0),
      stringsAsFactors = FALSE
    )
  }

  # Save originals for restoration
  .orig_getChromInfo <- NULL
  .orig_fetchExtended <- NULL

  # Override getChromInfoFromUCSC
  base::tryCatch({
    .orig_getChromInfo <- base::get("getChromInfoFromUCSC", envir = gi_ns)
    unlockBinding("getChromInfoFromUCSC", gi_ns)
    base::assign("getChromInfoFromUCSC", .safe_ucsc_fn, envir = gi_ns)
    lockBinding("getChromInfoFromUCSC", gi_ns)
    base::cat("  Overrode getChromInfoFromUCSC\n")
  }, error = function(e) {
    base::cat("  unlockBinding failed:", base::conditionMessage(e), "\n")
  })

  # Override fetchExtendedChromInfoFromUCSC (may not exist in all Bioc versions)
  base::tryCatch({
    .orig_fetchExtended <- base::get("fetchExtendedChromInfoFromUCSC", envir = gi_ns)
    unlockBinding("fetchExtendedChromInfoFromUCSC", gi_ns)
    base::assign("fetchExtendedChromInfoFromUCSC", function(genome, ...) {
      base::data.frame(
        UCSC_seqlevel = base::character(0), NCBI_seqlevel = base::character(0),
        stringsAsFactors = FALSE
      )
    }, envir = gi_ns)
    lockBinding("fetchExtendedChromInfoFromUCSC", gi_ns)
    base::cat("  Overrode fetchExtendedChromInfoFromUCSC\n")
  }, error = function(e) NULL)

  utils::flush.console()

  # Register cleanup to restore original functions on exit
  base::on.exit({
    base::tryCatch({
      if (!base::is.null(.orig_getChromInfo)) {
        unlockBinding("getChromInfoFromUCSC", gi_ns)
        base::assign("getChromInfoFromUCSC", .orig_getChromInfo, envir = gi_ns)
        lockBinding("getChromInfoFromUCSC", gi_ns)
      }
    }, error = function(e) NULL)
    base::tryCatch({
      if (!base::is.null(.orig_fetchExtended)) {
        unlockBinding("fetchExtendedChromInfoFromUCSC", gi_ns)
        base::assign("fetchExtendedChromInfoFromUCSC", .orig_fetchExtended, envir = gi_ns)
        lockBinding("fetchExtendedChromInfoFromUCSC", gi_ns)
      }
    }, error = function(e) NULL)
  }, add = TRUE)

  # Step 1: Create BSgenome package directory (UCSC already suppressed)
  base::cat("Creating BSgenome package for accession:", genbank_accession, "\n")
  utils::flush.console()
  base::tryCatch({
    pkg_name <- forge_BSgenome(dir = dir, accession = genbank_accession, metadata = metadata)
    base::cat("BSgenome package directory created:", pkg_name, "\n")
    utils::flush.console()
  }, error = function(e) {
    base::cat("Error creating BSgenome package:", base::conditionMessage(e), "\n")
    base::stop("Failed to create BSgenome package for accession: ", genbank_accession)
  })

  # ===== Step 2: Load BSgenome into current R session =====
  # Strategy: try multiple approaches in order of preference
  # 2a) devtools::load_all() -- runs in CURRENT session (UCSC suppression active)
  # 2b) devtools::install()  -- spawns NEW R process (UCSC suppression NOT active)
  # 2c) Direct BSgenome creation -- manual object construction

  pkg_dir <- base::file.path(dir, pkg_name)
  genome_obj_name <- clean_accession

  # --- Step 2a: Try devtools::load_all() (preferred for non-model organisms) ---
  base::cat("Step 2a: Attempting load_all() in current session...\n")
  utils::flush.console()
  load_all_ok <- FALSE
  base::tryCatch({
    devtools::load_all(pkg_dir, quiet = TRUE)
    load_all_ok <- TRUE
    base::cat("load_all() succeeded for:", pkg_name, "\n")
  }, error = function(e) {
    base::cat("load_all() failed:", base::conditionMessage(e), "\n")
  })
  utils::flush.console()

  if (load_all_ok) {
    # Try to get the BSgenome object from the loaded namespace
    bsg <- base::tryCatch({
      ns <- base::asNamespace(pkg_name)
      base::get(genome_obj_name, envir = ns)
    }, error = function(e) {
      base::cat("  Could not retrieve BSgenome from namespace:", base::conditionMessage(e), "\n")
      NULL
    })
    if (!base::is.null(bsg) && methods::is(bsg, "BSgenome")) {
      base::cat("BSgenome object loaded via load_all():", pkg_name, "\n")
      utils::flush.console()
      return(bsg)
    }
  }

  # --- Step 2b: Try devtools::install() (spawns new R process) ---
  base::cat("Step 2b: Attempting package installation...\n")
  utils::flush.console()
  install_ok <- FALSE
  base::tryCatch({
    devtools::install(pkg_dir, quiet = TRUE)
    install_ok <- TRUE
    base::cat("BSgenome package installed successfully:", pkg_name, "\n")
  }, error = function(e) {
    base::cat("Package install failed (expected for non-model organisms):",
              base::conditionMessage(e), "\n")
  })
  utils::flush.console()

  if (install_ok) {
    return(pkg_name)
  }

  # --- Step 2c: Direct BSgenome object creation (deepest fallback) ---
  # UCSC suppression is active in THIS process.
  base::cat("Step 2c: Creating BSgenome object directly (no package install)...\n")
  utils::flush.console()

  # Read sequences from FASTA
  seqs <- Biostrings::readDNAStringSet(fasta_file)
  seq_names_clean <- base::sub(" .*", "", base::names(seqs))
  base::names(seqs) <- seq_names_clean

  # Create Seqinfo with genome=NA (no UCSC validation)
  si <- methods::new("Seqinfo",
            seqnames = seq_names_clean,
            seqlengths = base::as.integer(Biostrings::width(seqs)),
            is_circular = base::rep(FALSE, base::length(seqs)),
            genome = base::rep(NA_character_, base::length(seqs)))

  bsg <- NULL
  extdata_dir <- base::file.path(dir, pkg_name, "inst", "extdata")

  # Ensure extdata has standard-named sequence files for OnDiskNamedSequences
  # BSgenome::getSeq() uses single_sequences directly, which requires
  # files named single_sequences.2bit or single_sequences.fa in extdata
  std_2bit <- base::file.path(extdata_dir, "single_sequences.2bit")
  std_fasta <- base::file.path(extdata_dir, "single_sequences.fa")
  std_fai <- base::file.path(extdata_dir, "single_sequences.fa.fai")

  if (!base::file.exists(std_2bit) && !base::file.exists(std_fasta)) {
    base::cat("  Ensuring standard-named sequence files in extdata...\n")
    # Check for any existing files with old naming and rename them
    old_2bit <- base::list.files(extdata_dir, pattern = "\\.2bit$", full.names = TRUE)
    old_fasta <- base::list.files(extdata_dir, pattern = "\\.fa$", full.names = TRUE)
    if (base::length(old_2bit) > 0) {
      base::file.copy(old_2bit[1], std_2bit, overwrite = TRUE)
      base::cat("  Renamed 2bit to single_sequences.2bit\n")
    } else if (base::length(old_fasta) > 0) {
      base::file.copy(old_fasta[1], std_fasta, overwrite = TRUE)
      old_fai <- base::paste0(old_fasta[1], ".fai")
      if (base::file.exists(old_fai)) {
        base::file.copy(old_fai, std_fai, overwrite = TRUE)
      } else {
        base::tryCatch({
          Rsamtools::indexFa(std_fasta)
        }, error = function(e) NULL)
      }
      base::cat("  Renamed FASTA to single_sequences.fa\n")
    } else {
      # Write a clean FASTA as single_sequences.fa
      Biostrings::writeXStringSet(seqs, std_fasta)
      base::tryCatch({
        Rsamtools::indexFa(std_fasta)
      }, error = function(e) NULL)
      base::cat("  Created single_sequences.fa in extdata\n")
    }
  }
  utils::flush.console()

  # Approach A: BSgenome() constructor (proper API, handles OnDiskNamedSequences)
  # With standard filenames, the constructor should find the sequence files
  base::tryCatch({
    bsg <- BSgenome::BSgenome(
      organism = metadata$organism_name,
      common_name = base::gsub(" .*", "", metadata$organism_name),
      genome = NA,
      provider = "NCBI",
      provider_version = genbank_accession,
      release_date = metadata$release_date,
      release_name = NA_character_,
      source_url = "",
      seqnames = seq_names_clean,
      circ_seqs = NA,
      mseqnames = base::character(0),
      seqs_pkgname = pkg_name,
      seqs_dirpath = extdata_dir
    )
    # CRITICAL: constructor may not set seqlengths from FASTA.
    # Without seqlengths, crisprBowtie drops ALL alignments (n0=0).
    cur_sl <- GenomeInfoDb::seqlengths(bsg)
    if (base::any(base::is.na(cur_sl))) {
      GenomeInfoDb::seqinfo(bsg) <- si
      base::cat("  seqlengths forced from Seqinfo\n")
    }
    base::cat("  BSgenome created via BSgenome() constructor\n")
  }, error = function(e) {
    base::cat("  BSgenome() constructor failed:", base::conditionMessage(e), "\n")
    bsg <<- NULL
  })
  utils::flush.console()

  # Approach B: manual BSgenome with proper OnDiskNamedSequences
  if (base::is.null(bsg)) {
    base::tryCatch({
      bsg <- methods::new("BSgenome")
      bsg@pkgname <- pkg_name
      bsg@seqinfo <- si
      bsg@user_seqnames <- base::setNames(seq_names_clean, seq_names_clean)
      # Create proper OnDiskNamedSequences (critical for getSeq to work)
      ondisk <- BSgenome:::OnDiskNamedSequences(extdata_dir, seqnames = seq_names_clean)
      methods::slot(bsg, "single_sequences", check = FALSE) <- ondisk
      # Also populate cache for [[ access
      cache_env <- bsg@.seqs_cache
      for (i in base::seq_along(seq_names_clean)) {
        base::assign(seq_names_clean[i], seqs[[i]], envir = cache_env)
      }
      base::cat("  BSgenome created with OnDiskNamedSequences + cache (",
                base::length(seqs), " sequences)\n")
    }, error = function(e) {
      base::cat("  OnDiskNamedSequences approach failed:", base::conditionMessage(e), "\n")
      bsg <<- NULL
    })
  }
  utils::flush.console()

  # Approach C: cache-only fallback (last resort -- getSeq may not work)
  if (base::is.null(bsg)) {
    base::tryCatch({
      bsg <- methods::new("BSgenome")
      bsg@pkgname <- pkg_name
      bsg@seqinfo <- si
      bsg@user_seqnames <- base::setNames(seq_names_clean, seq_names_clean)
      cache_env <- bsg@.seqs_cache
      for (i in base::seq_along(seq_names_clean)) {
        base::assign(seq_names_clean[i], seqs[[i]], envir = cache_env)
      }
      base::cat("  BSgenome created with cache only (", base::length(seqs),
                " sequences) -- WARNING: getSeq may not work\n")
    }, error = function(e) {
      base::cat("  Cache-only creation failed:", base::conditionMessage(e), "\n")
      bsg <<- NULL
    })
  }
  utils::flush.console()

  if (base::is.null(bsg)) {
    base::stop("All BSgenome creation methods failed for accession: ", genbank_accession,
               "\nPlease check your Bioconductor/BSgenome version.")
  }

  return(bsg)
}
