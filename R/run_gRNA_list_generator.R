
# ==============================================================================
# Internal: Extract FASTA from GenBank (.gbff) file
# Handles multi-contig genomes (multiple LOCUS/ORIGIN blocks)
# ==============================================================================
.extract_fasta_from_gbff <- function(gbff_file, fasta_output) {
  lines <- base::readLines(gbff_file, warn = FALSE)
  n <- base::length(lines)
  if (n == 0) base::stop("GBFF file is empty: ", gbff_file)

  # Parse all LOCUS..ORIGIN..// blocks
  fasta_entries <- base::list()
  current_locus <- NULL
  in_origin <- FALSE
  seq_parts <- base::character(0)

  for (i in base::seq_len(n)) {
    line <- lines[i]

    if (base::startsWith(line, "LOCUS")) {
      # Extract locus name (second field)
      parts <- base::strsplit(base::trimws(line), "\\s+")[[1]]
      if (base::length(parts) >= 2) {
        current_locus <- parts[2]
      } else {
        current_locus <- base::paste0("contig_", base::length(fasta_entries) + 1L)
      }
      in_origin <- FALSE
      seq_parts <- base::character(0)
      next
    }

    if (base::startsWith(line, "ORIGIN")) {
      in_origin <- TRUE
      next
    }

    if (base::startsWith(line, "//")) {
      if (in_origin && base::length(seq_parts) > 0) {
        seq_str <- base::toupper(base::gsub("[^a-zA-Z]", "",
                                             base::paste0(seq_parts, collapse = "")))
        if (base::nchar(seq_str) > 0) {
          fasta_entries[[base::length(fasta_entries) + 1L]] <- base::list(
            name = current_locus,
            seq  = seq_str
          )
        }
      }
      in_origin <- FALSE
      seq_parts <- base::character(0)
      current_locus <- NULL
      next
    }

    if (in_origin) {
      cleaned <- base::gsub("[^a-zA-Z]", "", line)
      if (base::nchar(cleaned) > 0) {
        seq_parts <- base::c(seq_parts, cleaned)
      }
    }
  }

  if (base::length(fasta_entries) == 0) {
    base::stop("No sequences found in GBFF file: ", gbff_file)
  }

  # Write FASTA
  out_dir <- base::dirname(fasta_output)
  if (!base::dir.exists(out_dir)) {
    base::dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  con <- base::file(fasta_output, "w")
  base::on.exit(base::close(con), add = TRUE)
  for (entry in fasta_entries) {
    base::writeLines(base::paste0(">", entry$name), con)
    # Write sequence in 70-char lines
    seq_str <- entry$seq
    seq_len <- base::nchar(seq_str)
    pos <- 1L
    while (pos <= seq_len) {
      end_pos <- base::min(pos + 69L, seq_len)
      base::writeLines(base::substring(seq_str, pos, end_pos), con)
      pos <- end_pos + 1L
    }
  }

  base::cat("  Extracted ", base::length(fasta_entries), " contig(s) from GBFF -> ",
            fasta_output, "\n", sep = "")
  base::invisible(fasta_output)
}


# Check if a sequence matches the locus_tag format
is_locus_tag <- function(seq) {
  # Pattern: underscore-separated, right part starts with RS or 5-7 digits
  pattern <- "^[A-Za-z0-9]+_[RS0-9]{5,7}$"
  return(base::grepl(pattern, seq))
}

# Validate if a sequence is a valid DNA sequence
is_valid_dna <- function(seq) {
  valid_pattern <- "^[ATGCNatgcnRYMKSWHBVDrymkswbdvh-]+$"
  return(base::grepl(valid_pattern, seq))
}

# Log invalid characters in a sequence
log_invalid_sequence <- function(seq, locus_tag) {
  invalid_chars <- base::unique(base::strsplit(base::gsub("[ATGCNatgcnRYMKSWHBVDrymkswbdvh-]", "", seq), "")[[1]])
  if (base::length(invalid_chars) > 0) {
    base::cat("Invalid characters found in sequence for locus_tag:", locus_tag, "- Invalid chars:", base::paste(invalid_chars, collapse = ", "), "\n")
  }
}

#' Run Complete gRNA Library Generation Workflow
#'
#' Executes the complete workflow for generating a genome-wide gRNA library.
#' Checks for and downloads the GenBank .fna and .gbff files (if not present), builds a BSgenome package,
#' creates a Bowtie index using Rbowtie, generates a gRNA list based on the target sequence,
#' retrieves a GenBank table (if not provided), builds a gRNA library, and saves the results to an Excel file.
#' The target sequence can be "all", a single locus tag (e.g., "XXXX_RSXXXXX"), a vector of locus tags, or a nucleotide sequence,
#' and the function automatically routes to the appropriate processing function.
#'
#' @param genbank_accession A character string specifying the GenBank accession number (e.g., "GCF_030376745.1").
#' @param target_sequence A character string or vector specifying the target(s): "all" to generate gRNA library for all genes, a single locus tag (e.g., "QT234_RS00005"), a vector of locus tags, or a nucleotide sequence (default: NULL).
#' @param genbank_table A data frame containing GenBank table data with columns including locus_tag and nt_seq (default: NULL, auto-retrieved if not provided).
#' @param nuclease A character string specifying the nuclease to use ("GeoCas9", "FnCas12a", or "FisCasI_B"; default: "GeoCas9").
#' @param bowtie_index A character string specifying the path to the Bowtie index (default: NULL, auto-generated).
#' @param genome A BSgenome object specifying the reference genome (default: NULL, auto-loaded).
#' @param output_file A character string specifying the output Excel file name (default: "gRNA_library.xlsx").
#' @param dir A character string specifying the working directory (default: current working directory).
#' @param position_range A numeric vector of length 2 specifying the position range (as percentage) for gRNA filtering (default: c(10, 90)).
#' @param tm_range A numeric vector of length 2 specifying the Tm range for gRNA filtering (default: c(50, 70)).
#' @param strand A character string specifying the strand to filter ("both", "5", or "3"; default: "both").
#' @param aligner A character string specifying the aligner to use ("bowtie", "bwa", or "biostrings"; default: "bowtie").
#' @param name_prefix A character string to prepend to gRNA names (default: NULL, auto-detected from nuclease). Examples: "Geo", "Sp", "ProjectX".
#' @param methylation_patterns Character vector of methylation recognition sequences
#'   using IUPAC ambiguity codes (default: NULL = no methylation filtering).
#'   Examples: \code{c("GCCAT", "CCANNNNNTTG", "GATC")}.
#'   When provided, gRNAs whose spacer + PAM region overlaps with any of these
#'   methylation motifs (on either strand) will be removed from the library.
#'   See \code{\link{filter_methylation_sites}} for details.
#' @return A data frame containing the final gRNA library with gRNA_name, gRNA_rank, composite_score, and merged GenBank table data.
#' @examples
#' \dontrun{
#' # Generate gRNA library for a single DNA sequence
#' accession <- "GCF_030376745.1"
#' target_seq <- "ttggcattgacaggtacagaccgcgtcaaacgcggcatggcggaaatgcaaaaaggcggcgtcattatggacgtcgtcaatgcagagcaagcgaagattgctgaggcggcaggggctgtcgcagtcatggcgctcgagcgtgtcccggcagacattcgcgccgctggcggtgtcgcgcgcatggctgatccgacgatcattgaagaagtgatgaacgccgtatcaatcccagtcatggcgaaagtgcgcatcgggcattatgtggaagcgcgtgttttagaggcgctcggcvtcgactatattgacgaaagtgaagtattgacgccggctgatgaagagttccatattgacaaacggcagtttacggtcccatttgtatgcggttgccgcgacttaggagaggccgcccgccgcattgctgaaggggcatcgatgttgcggacaaaaggggagccagggacaggaaacatcgtcgaggccgttcgccatatgcgcaaagtcaacgcgcaaatccgcaaagttgtcagcatgagcgaagacgaacttgtcgccgaggcgaaacagcttggggctccggttgaagtgctgcgtgaaatcaaacggcttggccgcctcccggtcgtcaacttcgccgccggcggtgtcgcgacaccagctgacgccgcgctcatgatgcacttaggcgccgacggtgtctfcgtcggttcgggcatttttaaatcggaaaatccggaaaaatacgctcgtgcgatcgttgaagcgacgactcattatgaagactatgagctgatcgcccatctatcgaaagggctgggcggcgcaatgcgcggcatcgatvtcgcgacactgctgccggagcatcggatgcaagaacgaggctggtaa"
#' gRNA_lib <- run_gRNA_list_generator(
#'   genbank_accession = accession,
#'   target_sequence = target_seq,
#'   genbank_table = NULL,
#'   bowtie_index = NULL,
#'   nuclease = "GeoCas9",
#'   position_range = c(20, 80),
#'   tm_range = c(55, 65),
#'   output_file = "gRNA_library.xlsx"
#' )
#'
#' # Generate gRNA library for all genes
#' gRNA_lib <- run_gRNA_list_generator(
#'   genbank_accession = accession,
#'   target_sequence = "all",
#'   genbank_table = NULL,
#'   bowtie_index = NULL,
#'   nuclease = "GeoCas9",
#'   output_file = "gRNA_library.xlsx"
#' )
#'
#' # Generate gRNA library for a specific locus tag
#' gRNA_lib <- run_gRNA_list_generator(
#'   genbank_accession = accession,
#'   target_sequence = "QT234_RS00005",
#'   genbank_table = NULL,
#'   bowtie_index = NULL,
#'   nuclease = "GeoCas9",
#'   output_file = "gRNA_library.xlsx"
#' )
#' print(gRNA_lib)
#' }
#' @export


run_gRNA_list_generator <- function(
    genbank_accession,
    target_sequence = NULL,
    genbank_table = NULL,
    nuclease = "GeoCas9",
    bowtie_index = NULL,
    genome = NULL,
    output_file = "gRNA_library.xlsx",
    dir = base::getwd(),
    position_range = c(10, 90),
    tm_range = c(50, 70),
    strand = "both",
    aligner = "bowtie",
    name_prefix = NULL,
    methylation_patterns = NULL
) {
  # Load plyr before dplyr to avoid conflicts
  if (!requireNamespace("plyr", quietly = TRUE)) {
    base::install.packages("plyr")
  }
  library(plyr)
  library(dplyr)

  # Step 0: Check for required R packages
  base::cat("Step 0: Checking for required R packages...\n")
  required_packages <- base::c("Rbowtie", "Biostrings", "rentrez", "dplyr", "xlsx", "BSgenome", "genbankr", "GenomicRanges", "devtools", "rtracklayer", "glue")
  missing_packages <- required_packages[!base::sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (base::length(missing_packages) > 0) {
    base::stop("The following R packages are missing: ", base::paste(missing_packages, collapse = ", "), ". Please install them using:\ninstall.packages(c('dplyr', 'xlsx', 'devtools', 'glue'))\nBiocManager::install(c('Rbowtie', 'Biostrings', 'rentrez', 'BSgenome', 'genbankr', 'GenomicRanges', 'rtracklayer'))")
  }
  base::cat("All required packages are available.\n")

  # Store original working directory
  original_dir <- base::getwd()

  # Step 1: Ensure .fna and .gbff files exist
  base::cat("Step 1: Checking genome files for accession:", genbank_accession, "\n")
  seqs_dir <- base::file.path(dir, "seqs_srcdir")
  base::dir.create(seqs_dir, recursive = TRUE, showWarnings = FALSE)
  fasta_file <- base::file.path(seqs_dir, base::paste0(genbank_accession, ".fna"))
  gbff_file  <- base::file.path(seqs_dir, base::paste0(genbank_accession, ".gbff"))

  # --- Ensure .gbff exists (check local first, then download) ---
  if (!base::file.exists(gbff_file)) {
    base::cat("  GBFF not found locally. Downloading from NCBI...\n")
    base::tryCatch({
      metadata <- download_genbank_gbff(dir = dir, accession = genbank_accession)
      gbff_file <- metadata$gbff_file
      if (!base::file.exists(gbff_file)) {
        base::stop("Failed to download GBFF file")
      }
      base::cat("  Downloaded GBFF:", gbff_file, "\n")
    }, error = function(e) {
      base::stop("Failed to obtain GenBank .gbff for accession: ", genbank_accession,
                 "\n  ", base::conditionMessage(e))
    })
  } else {
    base::cat("  Using existing GBFF:", gbff_file, "\n")
  }

  # --- Ensure .fna exists (extract from .gbff first, download as fallback) ---
  if (!base::file.exists(fasta_file)) {
    # Try extracting from .gbff (handles multi-contig genomes)
    if (base::file.exists(gbff_file)) {
      base::cat("  FASTA not found. Extracting from GBFF...\n")
      base::tryCatch({
        .extract_fasta_from_gbff(gbff_file, fasta_file)
      }, error = function(e) {
        base::cat("  GBFF extraction failed: ", base::conditionMessage(e), "\n")
      })
    }

    # Fallback: download from NCBI if extraction failed
    if (!base::file.exists(fasta_file)) {
      base::cat("  Downloading .fna from NCBI...\n")
      base::tryCatch({
        metadata <- download_genbank_fna(dir = dir, accession = genbank_accession)
        fasta_file <- metadata$fna_file
        if (!base::file.exists(fasta_file)) {
          base::stop("Failed to download FASTA file")
        }
        base::cat("  Downloaded FASTA:", fasta_file, "\n")
      }, error = function(e) {
        base::stop("Failed to obtain .fna for accession: ", genbank_accession,
                   "\n  ", base::conditionMessage(e))
      })
    }
  } else {
    base::cat("  Using existing FASTA:", fasta_file, "\n")
  }

  # Step 2: Clean accession for package name
  clean_accession <- base::gsub("[^a-zA-Z0-9]", "", genbank_accession)
  base::cat("Step 2: Cleaned accession for package name:", clean_accession, "\n")

  # Step 3: Set default bowtie_index
  if (base::is.null(bowtie_index)) {
    bowtie_index <- base::file.path(dir, base::paste0(genbank_accession, "_bowtie"), genbank_accession)
    base::cat("Using default bowtie_index:", bowtie_index, "\n")
  } else {
    base::cat("Using provided bowtie_index:", bowtie_index, "\n")
  }

  # Step 4: Bowtie index setup
  bowtie_dir <- base::file.path(dir, base::paste0(genbank_accession, "_bowtie"))
  index_prefix <- base::file.path(bowtie_dir, genbank_accession)
  ebwt_files <- base::list.files(
    bowtie_dir,
    pattern = base::paste0("^", genbank_accession, "\\..*\\.ebwt$"),
    full.names = TRUE
  )
  if (base::length(ebwt_files) == 0) {
    base::cat("Bowtie index files not found at:", bowtie_dir, "\nBuilding Bowtie index...\n")
    base::tryCatch({
      run_bowtie_build(directory = dir, accession = genbank_accession)
      ebwt_files <- base::list.files(
        bowtie_dir,
        pattern = base::paste0("^", genbank_accession, "\\..*\\.ebwt$"),
        full.names = TRUE
      )
      if (base::length(ebwt_files) == 0) {
        base::stop("Failed to build Bowtie index files at:", bowtie_dir)
      }
    }, error = function(e) {
      base::cat("Error building Bowtie index:", base::conditionMessage(e), "\n")
      base::stop("Failed to build Bowtie index at:", bowtie_dir)
    })
  }
  base::cat("Bowtie index files found:", ebwt_files, "\n")

  # Step 5: Build/load BSgenome
  base::cat("Step 5: Loading BSgenome for accession:", genbank_accession, "\n")
  pkg_name <- base::paste0("BSgenome.", clean_accession)

  .load_bsgenome_from_ns <- function(pkg_name) {
    obj_name <- base::gsub(".*\\.", "", pkg_name)
    base::get(obj_name, envir = base::asNamespace(pkg_name))
  }

  .validate_bsgenome_getSeq <- function(genome_obj) {
    # Validate that getSeq() works AND seqlengths are set
    # Missing seqlengths causes crisprBowtie to drop ALL alignments (n0=0)
    bsg_sn <- GenomeInfoDb::seqnames(genome_obj)
    if (base::length(bsg_sn) == 0) return(FALSE)
    # Check seqlengths are non-NA (critical for crisprBowtie)
    sl <- GenomeInfoDb::seqlengths(genome_obj)
    if (base::any(base::is.na(sl))) {
      base::cat("  seqlengths missing (NA) -- will rebuild\n")
      return(FALSE)
    }
    base::tryCatch({
      test_result <- BSgenome::getSeq(genome_obj,
        names = bsg_sn[1], start = 1L, end = 30L,
        strand = "+", as.character = TRUE)
      base::nchar(test_result) == 30L
    }, error = function(e) {
      base::cat("  getSeq validation failed:", base::conditionMessage(e), "\n")
      FALSE
    })
  }

  if (base::is.null(genome)) {
    # Try loading existing BSgenome first
    loaded_ok <- FALSE
    if (requireNamespace(pkg_name, quietly = TRUE)) {
      base::cat("BSgenome package found:", pkg_name, "\n")
      genome <- .load_bsgenome_from_ns(pkg_name)
      # Validate getSeq works (catches old BSgenome with broken single_sequences)
      if (.validate_bsgenome_getSeq(genome)) {
        base::cat("BSgenome getSeq validation: OK\n")
        loaded_ok <- TRUE
      } else {
        base::cat("BSgenome getSeq validation: FAILED -- will rebuild\n")
        genome <- NULL
      }
    }

    if (!loaded_ok) {
      # Build/create BSgenome (either first time or rebuild due to validation failure)
      base::cat("Building BSgenome for accession:", genbank_accession, "\n")
      # Force clean rebuild: unload old namespace if present
      if (requireNamespace(pkg_name, quietly = TRUE)) {
        base::tryCatch({
          base::unloadNamespace(pkg_name)
          base::cat("  Unloaded old BSgenome namespace\n")
        }, error = function(e) NULL)
      }
      # Delete old package directory to force rebuild with new file naming
      old_pkg_dir <- base::file.path(dir, pkg_name)
      if (base::dir.exists(old_pkg_dir)) {
        base::unlink(old_pkg_dir, recursive = TRUE, force = TRUE)
        base::cat("  Removed old package directory\n")
      }

      base::tryCatch({
        build_result <- build_bsgenome_from_accession(genbank_accession = genbank_accession, dir = dir)
        if (base::inherits(build_result, "BSgenome")) {
          genome <- build_result
          base::cat("Using in-memory BSgenome object\n")
        } else {
          genome <- .load_bsgenome_from_ns(pkg_name)
          base::cat("Using installed BSgenome package:", pkg_name, "\n")
        }
      }, error = function(e) {
        base::cat("Error building BSgenome:", base::conditionMessage(e), "\n")
        base::stop("Failed to build BSgenome for accession:", genbank_accession)
      })
    }
  } else {
    base::cat("Using provided genome:", base::deparse(base::substitute(genome)), "\n")
  }

  # Step 5.1: Suppress UCSC lookups for non-model organisms
  # Non-UCSC genome tags (e.g., NCBI accessions like GCF_030376765.1) trigger
  # GenomeInfoDb UCSC validation that fails for non-model organisms.
  # Use unlockBinding + assign to directly replace the UCSC function in the namespace.
  .ucsc_orig_fns <- base::list()
  .ucsc_gi_ns <- base::asNamespace("GenomeInfoDb")
  base::tryCatch({
    bsg_genome_tag <- genome@seqinfo@genome
    if (!base::all(base::is.na(bsg_genome_tag))) {
      base::cat("Step 5.1: Non-model organism detected (genome tag: ",
                base::paste(bsg_genome_tag, collapse = ", "), ")\n")
      base::cat("  Suppressing UCSC lookups...\n")

      .safe_getChromInfo <- function(genome,
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

      # Method: unlockBinding + assign (direct namespace modification)
      base::tryCatch({
        .ucsc_orig_fns[["getChromInfoFromUCSC"]] <- base::get("getChromInfoFromUCSC", envir = .ucsc_gi_ns)
        unlockBinding("getChromInfoFromUCSC", .ucsc_gi_ns)
        base::assign("getChromInfoFromUCSC", .safe_getChromInfo, envir = .ucsc_gi_ns)
        lockBinding("getChromInfoFromUCSC", .ucsc_gi_ns)
        base::cat("  Replaced getChromInfoFromUCSC (unlockBinding)\n")
      }, error = function(e) {
        base::cat("  unlockBinding failed:", base::conditionMessage(e), "\n")
      })

      # Also override fetchExtendedChromInfoFromUCSC (may not exist)
      base::tryCatch({
        .ucsc_orig_fns[["fetchExtendedChromInfoFromUCSC"]] <- base::get("fetchExtendedChromInfoFromUCSC", envir = .ucsc_gi_ns)
        unlockBinding("fetchExtendedChromInfoFromUCSC", .ucsc_gi_ns)
        base::assign("fetchExtendedChromInfoFromUCSC", function(genome, ...) {
          base::data.frame(UCSC_seqlevel = base::character(0), NCBI_seqlevel = base::character(0),
                           stringsAsFactors = FALSE)
        }, envir = .ucsc_gi_ns)
        lockBinding("fetchExtendedChromInfoFromUCSC", .ucsc_gi_ns)
        base::cat("  Replaced fetchExtendedChromInfoFromUCSC\n")
      }, error = function(e) NULL)

      base::cat("  UCSC lookups suppressed (",
                base::length(.ucsc_orig_fns), " functions replaced)\n")
    } else {
      base::cat("Step 5.1: BSgenome genome tag is NA -- no UCSC suppression needed\n")
    }
  }, error = function(e) {
    base::cat("  Note: Could not check/suppress UCSC lookups: ",
              base::conditionMessage(e), "\n")
  })

  # Cleanup: restore original UCSC functions when this function exits
  base::on.exit({
    for (.fn_name in base::names(.ucsc_orig_fns)) {
      base::tryCatch({
        unlockBinding(.fn_name, .ucsc_gi_ns)
        base::assign(.fn_name, .ucsc_orig_fns[[.fn_name]], envir = .ucsc_gi_ns)
        lockBinding(.fn_name, .ucsc_gi_ns)
      }, error = function(e) NULL)
    }
  }, add = TRUE)

  # Step 5.5: Validate .fna seqnames match BSgenome seqnames
  # Read ONLY the headers (not full sequences) for fast validation
  bsgenome_seqnames <- GenomeInfoDb::seqnames(genome)
  fna_names <- NULL
  base::tryCatch({
    # Fast header-only read: use fai index if available, otherwise scan headers
    fai_path <- base::paste0(fasta_file, ".fai")
    if (base::file.exists(fai_path)) {
      fai_data <- utils::read.delim(fai_path, header = FALSE, stringsAsFactors = FALSE)
      fna_names <- fai_data$V1
    } else {
      # Read headers without loading sequences (much faster)
      con <- base::file(fasta_file, "r")
      lines <- base::readLines(con, n = 500)  # enough to get all headers
      base::close(con)
      header_lines <- lines[base::grepl("^>", lines)]
      fna_names <- base::sub("^>([^ ]+).*", "\\1", header_lines)
    }
  }, error = function(e) {
    base::cat("Step 5.5: Could not read .fna headers:", base::conditionMessage(e), "\n")
  })

  if (!base::is.null(fna_names) && !base::all(fna_names %in% bsgenome_seqnames)) {
    base::cat("Step 5.5: Seqname mismatch detected!\n")
    base::cat("  .fna headers: ", base::paste(fna_names, collapse = ", "), "\n")
    base::cat("  BSgenome:     ", base::paste(bsgenome_seqnames, collapse = ", "), "\n")

    if (base::length(fna_names) == base::length(bsgenome_seqnames)) {
      base::cat("  Rewriting .fna headers to match BSgenome...\n")
      # Only read full sequences when we need to rewrite
      fna_seqs <- Biostrings::readDNAStringSet(fasta_file)
      base::names(fna_seqs) <- bsgenome_seqnames
      Biostrings::writeXStringSet(fna_seqs, fasta_file)
      base::cat("  .fna headers updated\n")

      # Rebuild Bowtie index with corrected .fna
      base::cat("  Rebuilding Bowtie index with corrected seqnames...\n")
      base::unlink(bowtie_dir, recursive = TRUE)
      run_bowtie_build(directory = dir, accession = genbank_accession)
      bowtie_index <- base::file.path(bowtie_dir, genbank_accession)
      base::cat("  Bowtie index rebuilt successfully\n")
    } else {
      base::warning("Cannot auto-fix seqname mismatch: .fna has ",
                     base::length(fna_names), " sequences but BSgenome has ",
                     base::length(bsgenome_seqnames),
                     ". Consider deleting the BSgenome package and re-running.")
    }
  } else {
    base::cat("Step 5.5: .fna seqnames match BSgenome -- OK\n")
  }

  # Step 6: Validate target_sequence
  if (base::is.null(target_sequence)) {
    base::stop("Error: 'target_sequence' must be provided.")
  }

  # Step 7: Get GenBank table if not provided
  if (base::is.null(genbank_table)) {
    genbank_table_file <- base::file.path(dir, base::paste0(genbank_accession, "_genbank_table.rds"))
    if (base::exists("genbank_table", envir = base::globalenv())) {
      base::cat("Step 7: Found genbank_table in global environment. Using existing object.\n")
      genbank_table <- base::get("genbank_table", envir = base::globalenv())
    } else if (base::file.exists(genbank_table_file)) {
      base::cat("Step 7: Loading cached GenBank table from:", genbank_table_file, "\n")
      base::tryCatch({
        genbank_table <- base::readRDS(genbank_table_file)
        if (!base::is.data.frame(genbank_table) || base::nrow(genbank_table) == 0) {
          base::cat("Warning: Loaded genbank_table is invalid or empty. Re-retrieving GenBank table.\n")
          # Change to seqs_srcdir for DNMB::run_DNMB()
          base::setwd(seqs_dir)
          warnings_log <- utils::capture.output({
            genbank_table <- DNMB::run_DNMB()
          }, type = "message")
          # Restore original directory
          base::setwd(original_dir)
          # Log warnings
          if (base::length(warnings_log) > 0) {
            base::cat("Warnings from DNMB::run_DNMB():\n", base::paste(warnings_log, collapse = "\n"), "\n")
          }
          # Check global environment if return value is NULL
          if (base::is.null(genbank_table) && base::exists("genbank_table", envir = base::globalenv())) {
            base::cat("DNMB::run_DNMB() returned NULL, retrieving genbank_table from global environment.\n")
            genbank_table <- base::get("genbank_table", envir = base::globalenv())
          }
        }
      }, error = function(e) {
        # Restore original directory
        base::setwd(original_dir)
        base::cat("Error loading cached GenBank table:", base::conditionMessage(e), "\n")
        base::cat("Re-retrieving GenBank table.\n")
        # Change to seqs_srcdir for DNMB::run_DNMB()
        base::setwd(seqs_dir)
        warnings_log <- utils::capture.output({
          genbank_table <- DNMB::run_DNMB()
        }, type = "message")
        # Restore original directory
        base::setwd(original_dir)
        # Log warnings
        if (base::length(warnings_log) > 0) {
          base::cat("Warnings from DNMB::run_DNMB():\n", base::paste(warnings_log, collapse = "\n"), "\n")
        }
        # Check global environment if return value is NULL
        if (base::is.null(genbank_table) && base::exists("genbank_table", envir = base::globalenv())) {
          base::cat("DNMB::run_DNMB() returned NULL, retrieving genbank_table from global environment.\n")
          genbank_table <- base::get("genbank_table", envir = base::globalenv())
        }
        if (base::is.null(genbank_table)) {
          base::stop("Failed to retrieve valid genbank_table from DNMB::run_DNMB(). Please check the .gbff file or DNMB::run_DNMB() implementation.")
        }
      })
    } else {
      base::cat("Step 7: Retrieving GenBank table from .gbff file in:", seqs_dir, "\n")
      base::tryCatch({
        # Change to seqs_srcdir for DNMB::run_DNMB()
        base::setwd(seqs_dir)
        warnings_log <- utils::capture.output({
          genbank_table <- DNMB::run_DNMB()
        }, type = "message")
        # Restore original directory
        base::setwd(original_dir)
        # Log warnings
        if (base::length(warnings_log) > 0) {
          base::cat("Warnings from DNMB::run_DNMB():\n", base::paste(warnings_log, collapse = "\n"), "\n")
        }
        # Check global environment if return value is NULL
        if (base::is.null(genbank_table) && base::exists("genbank_table", envir = base::globalenv())) {
          base::cat("DNMB::run_DNMB() returned NULL, retrieving genbank_table from global environment.\n")
          genbank_table <- base::get("genbank_table", envir = base::globalenv())
        }
        if (base::is.null(genbank_table)) {
          base::stop("Failed to retrieve valid genbank_table from DNMB::run_DNMB(). Please check the .gbff file or DNMB::run_DNMB() implementation.")
        }
      }, error = function(e) {
        # Restore original directory
        base::setwd(original_dir)
        base::cat("Error retrieving GenBank table:", base::conditionMessage(e), "\n")
        base::stop("Failed to retrieve GenBank table for accession: ", genbank_accession)
      })
    }

    # Debug: Inspect genbank_table
    base::cat("Inspecting genbank_table:\n")
    base::cat("Is data frame:", base::is.data.frame(genbank_table), "\n")
    base::cat("Dimensions:", base::paste(base::dim(genbank_table), collapse = " x "), "\n")
    base::cat("Column names:", base::paste(base::colnames(genbank_table), collapse = ", "), "\n")
    base::cat("Has locus_tag and nt_seq:", base::all(base::c("locus_tag", "nt_seq") %in% base::colnames(genbank_table)), "\n")

    if (!base::is.data.frame(genbank_table) || base::nrow(genbank_table) == 0) {
      base::stop("Error: genbank_table is invalid or empty. Please check DNMB::run_DNMB().")
    }
    if (!base::all(base::c("locus_tag", "nt_seq") %in% base::colnames(genbank_table))) {
      base::stop("genbank_table must contain 'locus_tag' and 'nt_seq' columns. Available columns: ", base::paste(base::colnames(genbank_table), collapse = ", "))
    }

    base::tryCatch({
      base::saveRDS(genbank_table, genbank_table_file)
      base::cat("GenBank table saved to:", genbank_table_file, "\n")
    }, error = function(e) {
      base::cat("Error saving GenBank table to RDS:", base::conditionMessage(e), "\n")
      base::stop("Failed to save genbank_table to RDS. Check file permissions or disk space.")
    })
  } else {
    base::cat("Step 7: Using provided genbank_table\n")
    if (!base::is.data.frame(genbank_table) || base::nrow(genbank_table) == 0) {
      base::stop("Error: Provided genbank_table is invalid or empty.")
    }
    if (!base::all(base::c("locus_tag", "nt_seq") %in% base::colnames(genbank_table))) {
      base::stop("genbank_table must contain 'locus_tag' and 'nt_seq' columns. Available columns: ", base::paste(base::colnames(genbank_table), collapse = ", "))
    }
  }

  base::cat("Columns in genbank_table:", base::paste(base::colnames(genbank_table), collapse = ", "), "\n")

  # Step 8: Clean nt_seq data (only if locus tags will be processed)
  if (base::length(target_sequence) > 1 || (base::length(target_sequence) == 1 && !is_valid_dna(target_sequence))) {
    base::cat("Step 8: Cleaning nt_seq data\n")
    genbank_table$nt_seq <- base::sapply(base::seq_along(genbank_table$nt_seq), function(i) {
      seq <- genbank_table$nt_seq[i]
      locus_tag <- genbank_table$locus_tag[i]
      if (!is_valid_dna(seq)) {
        log_invalid_sequence(seq, locus_tag)
      }
      seq <- base::gsub("[^ATGCNatgcnRYMKSWHBVDrymkswbdvh-]", "N", seq)
      return(seq)
    })
  } else {
    base::cat("Step 8: Skipping nt_seq cleaning as target_sequence appears to be a DNA sequence\n")
  }

  # Step 9: Process target_sequence
  gRNA_list <- NULL
  if (base::length(target_sequence) > 1) {
    base::cat("Step 9: Processing as a list of locus tags\n")
    if (!base::all(base::sapply(target_sequence, is_locus_tag))) {
      base::stop("Error: All elements in target_sequence must be valid locus tags when providing a vector.")
    }
    gRNA_list <- generate_gRNA_for_locus_tags(
      locus_tags = target_sequence,
      genbank_table = genbank_table,
      genbank_accession = genbank_accession,
      nuclease = nuclease,
      bowtie_index = bowtie_index,
      genome = genome,
      position_range = position_range,
      tm_range = tm_range,
      strand = strand,
      dir = dir,
      output_file = NULL,
      name_prefix = name_prefix
    )
  } else {
    target_seq_lower <- tolower(target_sequence)
    cat("Step 9: target_seq_lower (first 20 chars) =", substr(target_seq_lower, 1, 20), "...\n")
    cat("Step 9: is_locus_tag(target_sequence) =", is_locus_tag(target_sequence), "\n")
    is_all_target <- target_seq_lower == "all" || startsWith(target_seq_lower, "all") || endsWith(target_seq_lower, "all")
    is_single_locus_tag <- is_locus_tag(target_sequence)

    if (is_all_target) {
      cat("Step 9: Processing all locus tags (target_sequence = 'all')\n")
      genbank_filtered <- genbank_table %>%
        dplyr::filter(!is.na(protein_id) & protein_id != "")
      locus_tags <- genbank_filtered$locus_tag
      gRNA_list <- generate_gRNA_for_locus_tags(
        locus_tags = locus_tags,
        genbank_table = genbank_filtered,
        genbank_accession = genbank_accession,
        nuclease = nuclease,
        bowtie_index = bowtie_index,
        genome = genome,
        position_range = position_range,
        tm_range = tm_range,
        strand = strand,
        dir = dir,
        output_file = NULL,
        name_prefix = name_prefix
      )
    } else if (is_single_locus_tag) {
      cat("Step 9: Processing as a single locus tag\n")
      gRNA_list <- generate_gRNA_for_locus_tags(
        locus_tags = target_sequence,
        genbank_table = genbank_table,
        genbank_accession = genbank_accession,
        nuclease = nuclease,
        bowtie_index = bowtie_index,
        genome = genome,
        position_range = position_range,
        tm_range = tm_range,
        strand = strand,
        dir = dir,
        output_file = NULL,
        name_prefix = name_prefix
      )
    } else {
      cat("Step 9: Processing as a single DNA sequence\n")
      cat("Step 9: Passing aligner =", aligner, "to generate_gRNA_for_sequence\n")
      gRNA_list <- generate_gRNA_for_sequence(
        target_sequence = target_sequence,
        nuclease = nuclease,
        bowtie_index = bowtie_index,
        genome = genome,
        position_range = position_range,
        tm_range = tm_range,
        strand = strand,
        aligner = aligner  # Pass the aligner parameter
      )
    }
  }

  # Step 10: Build gRNA library per locus_tag
  is_locus_tag_processed <- (length(target_sequence) > 1 || is_all_target || is_single_locus_tag)
  if (is_locus_tag_processed) {
    cat("Step 10: Using gRNA_list from Step 9 as the final library\n")
    gRNA_library <- gRNA_list
  } else {
    cat("Step 10: Using gRNA_list from Step 9 for single DNA sequence\n")
    gRNA_list$locus_tag <- "custom_target"
    if ("locus_tag" %in% colnames(gRNA_list)) {
      gRNA_library <- gRNA_list %>% dplyr::select(locus_tag, dplyr::everything())
    } else {
      warning("locus_tag column not found in gRNA_list; proceeding without reordering columns.")
      gRNA_library <- gRNA_list
    }
  }

  # Step 10.5: Apply gRNA naming for custom DNA sequence targets
  if (!is_locus_tag_processed && nrow(gRNA_library) > 0) {
    cat("Step 10.5: Generating gRNA names for custom sequence target\n")
    gRNA_library <- generate_gRNA_names(
      gRNA_df = gRNA_library,
      prefix = name_prefix,
      nuclease = nuclease,
      rank_by = "composite_score"
    )
  }

  # Step 10.7: Methylation site filtering (optional)
  if (!base::is.null(methylation_patterns) && base::length(methylation_patterns) > 0 &&
      base::nrow(gRNA_library) > 0) {
    base::cat("Step 10.7: Filtering methylation site conflicts\n")
    nuclease_name <- base::ifelse(base::is.character(nuclease), nuclease, "GeoCas9")
    n_before <- base::nrow(gRNA_library)
    gRNA_library <- filter_methylation_sites(
      gRNA_df = gRNA_library,
      methylation_patterns = methylation_patterns,
      nuclease = nuclease_name,
      exclude = TRUE
    )
    n_after <- base::nrow(gRNA_library)
    base::cat("  Methylation filter: ", n_before, " → ", n_after,
              " gRNAs (removed ", n_before - n_after, ")\n", sep = "")
  }

  # Step 11: Merge with GenBank table and save to Excel
  cat("Step 11: Merging gRNA library with GenBank table and saving to Excel\n")
  if (!("custom_target" %in% gRNA_library$locus_tag)) {
    gRNA_library <- merge(
      x = gRNA_library,
      y = genbank_table %>% dplyr::select(locus_tag, gene, product, direction, protein_id, start, end, amino_acid),
      by = "locus_tag",
      suffix = c("_by_gene", "")
    )
    if ("locus_tag" %in% colnames(gRNA_library)) {
      gRNA_library <- gRNA_library %>% dplyr::select(locus_tag, dplyr::everything())
    } else {
      warning("locus_tag column not found in gRNA_library after merging; proceeding without reordering columns.")
    }
  }

  # Step 11.5: Reorder columns for clean output
  priority_cols <- c("gRNA_name", "gRNA_rank", "locus_tag", "gene", "product",
                     "protospacer", "pam", "strand", "composite_score",
                     "percentGC", "Tm", "position_percent", "n0", "n1")
  existing_priority <- intersect(priority_cols, colnames(gRNA_library))
  remaining_cols <- setdiff(colnames(gRNA_library), existing_priority)
  # Remove non-serializable and redundant columns
  remove_cols <- c("alignments", "score")
  remaining_cols <- setdiff(remaining_cols, remove_cols)
  gRNA_library <- gRNA_library[, c(existing_priority, remaining_cols), drop = FALSE]

  cat("Dimensions of gRNA_library:", dim(gRNA_library), "\n")
  cat("Column names of gRNA_library:", paste(colnames(gRNA_library), collapse = ", "), "\n")

  full_output_path <- file.path(dir, output_file)
  tryCatch({
    output_dir <- dir
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      cat("Created output directory:", output_dir, "\n")
    }
    xlsx::write.xlsx(gRNA_library, file = full_output_path, row.names = FALSE)
    cat("Successfully saved gRNA library to:", full_output_path, "\n")
  }, error = function(e) {
    cat("Error saving gRNA library to Excel:", conditionMessage(e), "\n")
    cat("gRNA_library structure:\n")
    str(gRNA_library)
    stop("Failed to save gRNA library to ", full_output_path, ". Please check file path, permissions, or xlsx package installation.")
  })

  cat("gRNA library generation completed. Output saved to:", full_output_path, "\n")
  return(gRNA_library)
}
