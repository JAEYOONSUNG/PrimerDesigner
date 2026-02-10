#' Forge BSgenome Package
#'
#' Creates a BSgenome package from a downloaded .fna file and metadata.
#' For non-model organisms, creates the BSgenome package manually to avoid
#' UCSC validation errors.
#'
#' @param dir A character string specifying the directory to save the BSgenome package.
#' @param accession A character string specifying the GenBank accession number (e.g., "GCF_030376745.1").
#' @param metadata A list containing metadata from the download_genbank_fna function (fna_file, organism_name, release_date, seqs_dir).
#' @return A character string specifying the name of the generated BSgenome package.
#' @examples
#' \dontrun{
#' metadata <- download_genbank_fna(dir = getwd(), accession = "GCF_030376745.1")
#' pkg_name <- forge_BSgenome(dir = getwd(), accession = "GCF_030376745.1", metadata = metadata)
#' print(pkg_name)
#' }
#' @export

forge_BSgenome <- function(dir, accession, metadata) {
  # Normalize paths to handle spaces
  dir <- base::normalizePath(dir, mustWork = TRUE)
  seqs_dir <- base::normalizePath(metadata$seqs_dir, mustWork = TRUE)
  clean_accession <- base::gsub("[^a-zA-Z0-9]", "", accession)
  pkg_name <- base::paste0("BSgenome.", clean_accession)
  pkg_dir <- base::file.path(dir, pkg_name)

  # Step 1: Read FASTA file
  base::cat("[forge] Reading FASTA file:", metadata$fna_file, "\n")
  utils::flush.console()
  seqs <- Biostrings::readDNAStringSet(metadata$fna_file)
  seq_names_clean <- base::sub(" .*", "", base::names(seqs))
  seq_lengths <- Biostrings::width(seqs)
  base::names(seq_lengths) <- seq_names_clean
  base::cat("[forge] Read", base::length(seqs), "sequences:",
            base::paste(seq_names_clean, collapse = ", "), "\n")
  utils::flush.console()

  # Step 2: Create 2bit file (reuse existing if available)
  seq_2bit_out <- base::file.path(seqs_dir, base::paste0(accession, "_out.2bit"))
  if (base::file.exists(seq_2bit_out)) {
    base::cat("[forge] 2bit file already exists, reusing:", seq_2bit_out, "\n")
  } else {
    base::cat("[forge] Creating 2bit file:", seq_2bit_out, "\n")
    utils::flush.console()
    # Strip all metadata from seqs to avoid any UCSC triggers
    plain_seqs <- seqs
    base::names(plain_seqs) <- seq_names_clean  # Clean names only, no descriptions
    base::tryCatch({
      rtracklayer::export.2bit(plain_seqs, seq_2bit_out)
      base::cat("[forge] 2bit file created successfully\n")
    }, error = function(e) {
      base::cat("[forge] rtracklayer::export.2bit failed:", base::conditionMessage(e), "\n")
      base::cat("[forge] Trying alternative 2bit creation...\n")
      utils::flush.console()
      # Alternative: write FASTA, then convert using rtracklayer::export
      fasta_tmp <- base::tempfile(fileext = ".fa")
      Biostrings::writeXStringSet(plain_seqs, fasta_tmp)
      base::tryCatch({
        tmp_seqs <- Biostrings::readDNAStringSet(fasta_tmp)
        rtracklayer::export(tmp_seqs, seq_2bit_out, format = "2bit")
        base::cat("[forge] Alternative 2bit creation succeeded\n")
      }, error = function(e2) {
        base::cat("[forge] Alternative also failed:", base::conditionMessage(e2), "\n")
        base::cat("[forge] Will use FASTA file directly instead of 2bit\n")
        # Copy FASTA as fallback -- zzz.R will handle this
        fasta_dest <- base::file.path(seqs_dir, base::paste0(accession, "_genome.fa"))
        Biostrings::writeXStringSet(plain_seqs, fasta_dest)
      })
      base::file.remove(fasta_tmp)
    })
  }
  utils::flush.console()

  # Step 3: Check and remove existing package directory
  if (base::dir.exists(pkg_dir)) {
    base::cat("[forge] Removing existing package directory:", pkg_dir, "\n")
    base::unlink(pkg_dir, recursive = TRUE, force = TRUE)
    max_attempts <- 5
    attempt <- 1
    while (base::dir.exists(pkg_dir) && attempt <= max_attempts) {
      Sys.sleep(1)
      base::unlink(pkg_dir, recursive = TRUE, force = TRUE)
      attempt <- attempt + 1
    }
    if (base::dir.exists(pkg_dir)) {
      base::stop("Failed to remove existing package directory: ", pkg_dir)
    }
  }

  # Step 4: Create BSgenome package manually
  base::cat("[forge] Creating BSgenome package structure...\n")
  utils::flush.console()

  # Create directory structure
  base::dir.create(base::file.path(pkg_dir, "R"), recursive = TRUE, showWarnings = FALSE)
  base::dir.create(base::file.path(pkg_dir, "inst", "extdata"), recursive = TRUE, showWarnings = FALSE)

  # Determine which sequence file to use (2bit or FASTA)
  # CRITICAL: BSgenome's OnDiskNamedSequences() constructor looks for EXACT filenames:
  #   single_sequences.2bit  -> TwobitNamedSequences
  #   single_sequences.fa    -> FastaNamedSequences
  # Using non-standard names causes the constructor to fail, breaking getSeq().
  extdata_dir <- base::file.path(pkg_dir, "inst", "extdata")
  use_twobit <- base::file.exists(seq_2bit_out)
  if (use_twobit) {
    # Copy as standard name: single_sequences.2bit
    seq_file_name <- "single_sequences.2bit"
    twobit_dest <- base::file.path(extdata_dir, seq_file_name)
    base::file.copy(seq_2bit_out, twobit_dest, overwrite = TRUE)
    base::cat("[forge] Copied 2bit as single_sequences.2bit\n")
  } else {
    # Fallback: use FASTA file + create index for random access
    fasta_fallback <- base::file.path(seqs_dir, base::paste0(accession, "_genome.fa"))
    if (!base::file.exists(fasta_fallback)) {
      plain_seqs <- seqs
      base::names(plain_seqs) <- seq_names_clean
      Biostrings::writeXStringSet(plain_seqs, fasta_fallback)
    }
    # Index the FASTA for efficient random access
    fai_file <- base::paste0(fasta_fallback, ".fai")
    if (!base::file.exists(fai_file)) {
      base::tryCatch({
        Rsamtools::indexFa(fasta_fallback)
        base::cat("[forge] Indexed FASTA:", fai_file, "\n")
      }, error = function(e) {
        base::cat("[forge] FASTA indexing failed:", base::conditionMessage(e), "\n")
      })
    }
    # Copy as standard name: single_sequences.fa (+ .fai index)
    seq_file_name <- "single_sequences.fa"
    fasta_dest <- base::file.path(extdata_dir, seq_file_name)
    base::file.copy(fasta_fallback, fasta_dest, overwrite = TRUE)
    # Also copy the .fai index with matching name
    fai_dest <- base::file.path(extdata_dir, "single_sequences.fa.fai")
    if (base::file.exists(fai_file)) {
      base::file.copy(fai_file, fai_dest, overwrite = TRUE)
      base::cat("[forge] Copied FASTA as single_sequences.fa + .fai index\n")
    } else {
      # Try to index the destination file directly
      base::tryCatch({
        Rsamtools::indexFa(fasta_dest)
        base::cat("[forge] Indexed single_sequences.fa in extdata\n")
      }, error = function(e) {
        base::cat("[forge] WARNING: FASTA index creation failed. getSeq may not work.\n")
      })
    }
  }
  utils::flush.console()

  # Write DESCRIPTION
  desc_lines <- base::c(
    base::paste0("Package: ", pkg_name),
    "Type: Package",
    base::paste0("Title: Full genome sequence for ", accession),
    "Version: 1.0.0",
    base::paste0("Description: Full genome sequence of ", accession, "."),
    "Author: PrimerDesigner",
    "Maintainer: PrimerDesigner <auto@auto.com>",
    "License: Artistic-2.0",
    "Depends: R (>= 4.0), BSgenome (>= 1.56.0)",
    "Imports: BSgenome, rtracklayer, GenomeInfoDb, Biostrings, Rsamtools, IRanges, methods",
    base::paste0("organism: ", metadata$organism_name),
    base::paste0("common_name: ", base::gsub(" .*", "", metadata$organism_name)),
    "genome: NA",
    "provider: NCBI",
    base::paste0("provider_version: ", accession),
    base::paste0("release_date: ", metadata$release_date),
    base::paste0("source_url: https://www.ncbi.nlm.nih.gov/assembly/", accession),
    base::paste0("BSgenomeObjname: ", clean_accession),
    base::paste0("seqfile_name: ", seq_file_name),
    base::paste0("seqnames: ", base::paste(seq_names_clean, collapse = ", ")),
    "circ_seqs: character(0)",
    base::paste0("organism_biocview: ", base::gsub(" ", "_", metadata$organism_name))
  )
  base::writeLines(desc_lines, base::file.path(pkg_dir, "DESCRIPTION"))
  base::cat("[forge] Wrote DESCRIPTION\n")

  # Write NAMESPACE
  ns_lines <- base::c(
    "import(BSgenome)",
    "import(methods)",
    base::paste0("export(", clean_accession, ")")
  )
  base::writeLines(ns_lines, base::file.path(pkg_dir, "NAMESPACE"))
  base::cat("[forge] Wrote NAMESPACE\n")

  # Write bsgenome_loader.R -- .onLoad() hook that creates the BSgenome object
  # Uses BSgenome() constructor (proper API) with OnDiskNamedSequences fallback
  # BSgenome() constructor handles OnDiskNamedSequences creation internally
  # Metadata (organism, provider, etc.) stored via Annotated, not as S4 slots
  loader_seqnames <- base::paste0('"', seq_names_clean, '"', collapse = ", ")
  loader_seqlengths <- base::paste0(seq_lengths, "L", collapse = ", ")
  loader_circulars <- base::paste0(base::rep("FALSE", base::length(seq_names_clean)), collapse = ", ")

  # Escape organism name for R string
  org_name_safe <- base::gsub('"', '\\\\"', metadata$organism_name)
  common_name_safe <- base::gsub('"', '\\\\"', base::gsub(" .*", "", metadata$organism_name))

  loader_lines <- base::c(
    base::paste0('.pkgname <- "', pkg_name, '"'),
    base::paste0('.objname <- "', clean_accession, '"'),
    base::paste0('.seqfile <- "', seq_file_name, '"'),
    "",
    ".onLoad <- function(libname, pkgname) {",
    "  extdata_dir <- system.file(\"extdata\", package = pkgname,",
    "                              lib.loc = libname, mustWork = TRUE)",
    "  seq_path <- file.path(extdata_dir, .seqfile)",
    "",
    "  # Embedded sequence metadata (avoids any UCSC lookups)",
    base::paste0("  sn <- c(", loader_seqnames, ")"),
    base::paste0("  sl <- c(", loader_seqlengths, ")"),
    "  names(sl) <- sn",
    base::paste0("  ic <- c(", loader_circulars, ")"),
    "  names(ic) <- sn",
    "",
    "  # Build Seqinfo with embedded metadata (used by both paths)",
    "  si <- new(\"Seqinfo\",",
    "            seqnames = sn,",
    "            seqlengths = sl,",
    "            is_circular = ic,",
    "            genome = rep(NA_character_, length(sn)))",
    "",
    "  # Try BSgenome() constructor first (proper API)",
    "  bsgenome <- tryCatch({",
    "    obj <- BSgenome::BSgenome(",
    base::paste0("      organism = \"", org_name_safe, "\","),
    base::paste0("      common_name = \"", common_name_safe, "\","),
    "      genome = NA,",
    "      provider = \"NCBI\",",
    base::paste0("      provider_version = \"", accession, "\","),
    base::paste0("      release_date = \"", metadata$release_date, "\","),
    "      release_name = NA_character_,",
    "      source_url = \"\",",
    base::paste0("      seqnames = c(", loader_seqnames, "),"),
    "      circ_seqs = NA,",
    "      mseqnames = character(0),",
    "      seqs_pkgname = pkgname,",
    "      seqs_dirpath = extdata_dir",
    "    )",
    "    # CRITICAL: constructor may not set seqlengths from FASTA.",
    "    # Without seqlengths, crisprBowtie drops ALL alignments (n0=0).",
    "    # Force-set Seqinfo with correct seqlengths.",
    "    cur_sl <- GenomeInfoDb::seqlengths(obj)",
    "    if (any(is.na(cur_sl))) {",
    "      GenomeInfoDb::seqinfo(obj) <- si",
    "    }",
    "    obj",
    "  }, error = function(e) NULL)",
    "",
    "  # Fallback: manual BSgenome with proper OnDiskNamedSequences",
    "  if (is.null(bsgenome)) {",
    "    cat(\"[BSgenome] Constructor failed, using fallback...\\n\")",
    "    bsgenome <- new(\"BSgenome\")",
    "    bsgenome@pkgname <- pkgname",
    "    bsgenome@seqinfo <- si",
    "    bsgenome@user_seqnames <- setNames(sn, sn)",
    "",
    "    # Try to create proper OnDiskNamedSequences from extdata_dir",
    "    # BSgenome::getSeq() bypasses cache and uses single_sequences directly,",
    "    # so this MUST be a valid OnDiskNamedSequences (not DNAStringSet!)",
    "    tryCatch({",
    "      ondisk <- BSgenome:::OnDiskNamedSequences(extdata_dir, seqnames = sn)",
    "      methods::slot(bsgenome, \"single_sequences\", check = FALSE) <- ondisk",
    "      cat(\"[BSgenome] OnDiskNamedSequences created from extdata\\n\")",
    "    }, error = function(e) {",
    "      cat(\"[BSgenome] OnDiskNamedSequences failed:\", conditionMessage(e), \"\\n\")",
    "      cat(\"[BSgenome] Falling back to cache-only (getSeq may not work)\\n\")",
    "    })",
    "",
    "    # Also populate cache for [[ access",
    "    tryCatch({",
    "      seqs_data <- Biostrings::readDNAStringSet(seq_path)",
    "      cache_env <- bsgenome@.seqs_cache",
    "      for (i in seq_along(sn)) {",
    "        assign(sn[i], seqs_data[[i]], envir = cache_env)",
    "      }",
    "    }, error = function(e) {",
    "      cat(\"[BSgenome] Cache population failed:\", conditionMessage(e), \"\\n\")",
    "    })",
    "  }",
    "",
    "  ns <- asNamespace(pkgname)",
    "  assign(.objname, bsgenome, envir = ns)",
    "  namespaceExport(ns, .objname)",
    "}"
  )
  base::writeLines(loader_lines, base::file.path(pkg_dir, "R", "bsgenome_loader.R"))
  base::cat("[forge] Wrote R/bsgenome_loader.R\n")

  # Verify package directory creation
  if (!base::dir.exists(pkg_dir)) {
    base::stop("BSgenome package directory not created: ", pkg_dir)
  }
  base::cat("[forge] BSgenome package ready:", pkg_dir, "\n")
  utils::flush.console()

  return(pkg_name)
}
