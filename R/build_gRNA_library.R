#' Build gRNA Library from GenBank Table
#'
#' Builds a gRNA library by generating gRNA lists for each locus tag in a GenBank table.
#' Uses nuclease_with_parameter to design gRNAs for each target sequence.
#'
#' @param genbank_table A data frame containing GenBank table data with columns including locus_tag and nt_seq.
#' @param accession A character string specifying the GenBank accession number (e.g., "GCF_030376745.1").
#' @param nuclease A character string specifying the nuclease to use ("GeoCas9", "FnCas12a", or "FisCasI_B"; default: "GeoCas9").
#' @param bowtie_index A character string specifying the path to the Bowtie index.
#' @param genome A BSgenome object specifying the reference genome (default: NULL, auto-loaded using accession).
#' @param position_range A numeric vector of length 2 specifying the position range (as percentage) for gRNA filtering (default: c(10, 90)).
#' @param tm_range A numeric vector of length 2 specifying the Tm range for gRNA filtering (default: c(50, 70)).
#' @param strand A character string specifying the strand to filter ("both", "5", or "3"; default: "both").
#' @param dir A character string specifying the working directory (default: current working directory).
#' @param name_prefix A character string to prepend to gRNA names (default: NULL, auto-detected from nuclease).
#' @return A data frame containing the gRNA library with gRNA_name, gRNA_rank, locus_tag, and gRNA details.
#' @examples
#' \dontrun{
#' genbank_table <- data.frame(locus_tag = "tag1", nt_seq = "atcg")
#' gRNA_lib <- build_gRNA_library(
#'   genbank_table = genbank_table,
#'   accession = "GCF_030376745.1",
#'   bowtie_index = "path/to/bowtie_index",
#'   nuclease = "GeoCas9",
#'   position_range = c(20, 80),
#'   tm_range = c(55, 65)
#' )
#' print(gRNA_lib)
#' }
#' @export

build_gRNA_library <- function(
    genbank_table,
    accession,
    nuclease = "GeoCas9",
    bowtie_index,
    genome = NULL,
    position_range = c(10, 90),
    tm_range = c(50, 70),
    strand = "both",
    dir = base::getwd(),
    name_prefix = NULL
) {
  # Validate inputs
  if (!base::is.data.frame(genbank_table)) {
    base::stop("genbank_table must be a data frame")
  }
  if (base::nrow(genbank_table) == 0) {
    base::stop("genbank_table is empty")
  }
  if (!base::all(base::c("locus_tag", "nt_seq") %in% base::colnames(genbank_table))) {
    base::stop("genbank_table must contain 'locus_tag' and 'nt_seq' columns")
  }
  if (base::missing(accession) || !base::is.character(accession) || base::nchar(accession) == 0) {
    base::stop("accession must be a non-empty character string")
  }
  if (base::missing(bowtie_index) || !base::is.character(bowtie_index) || base::nchar(bowtie_index) == 0) {
    base::stop("bowtie_index must be a non-empty character string")
  }
  if (base::length(position_range) != 2 || !base::is.numeric(position_range) || base::any(position_range < 0) || base::any(position_range > 100)) {
    base::stop("position_range must be a numeric vector of length 2 with values between 0 and 100")
  }
  if (base::length(tm_range) != 2 || !base::is.numeric(tm_range) || base::any(tm_range < 0)) {
    base::stop("tm_range must be a numeric vector of length 2 with non-negative values")
  }
  if (!strand %in% base::c("both", "5", "3")) {
    base::stop("strand must be one of 'both', '5', or '3'")
  }

  # Set default genome if NULL
  if (base::is.null(genome)) {
    clean_accession <- base::gsub("[^a-zA-Z0-9]", "", accession)
    pkg_name <- base::paste0("BSgenome.", clean_accession)
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      base::cat("Building BSgenome package for accession:", accession, "\n")
      build_bsgenome_from_accession(genbank_accession = accession, dir = dir)
    }
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      base::stop("BSgenome package not found:", pkg_name, "\nPlease ensure it is installed.")
    }
    genome <- base::get(base::gsub(".*\\.", "", pkg_name), envir = base::asNamespace(pkg_name))
    base::cat("Using default genome package:", pkg_name, "\n")
  }

  # Generate gRNA for each locus_tag
  failed_locus_tags <- base::c()
  gRNA_list <- base::lapply(1:base::nrow(genbank_table), function(i) {
    locus_tag <- genbank_table$locus_tag[i]
    target_seq <- genbank_table$nt_seq[i]

    # Validate sequence
    if (!is_valid_dna(target_seq)) {
      log_invalid_sequence(target_seq, locus_tag)
      failed_locus_tags <<- base::c(failed_locus_tags, locus_tag)
      return(NULL)  # Skip invalid sequences
    }

    base::cat("Generating gRNA for locus_tag:", locus_tag, "\n")
    base::tryCatch({
      gRNA_df <- nuclease_with_parameter(
        target_sequence = target_seq,
        strand = strand,
        bowtie_index = bowtie_index,
        genome = genome,
        nuclease = nuclease,
        position_range = position_range,
        tm_range = tm_range
      )
      gRNA_df$locus_tag <- locus_tag
      # Add gene name for downstream naming
      if ("gene" %in% base::colnames(genbank_table)) {
        gene_name <- genbank_table$gene[i]
        gRNA_df$gene <- base::ifelse(base::is.na(gene_name) || gene_name == "", NA_character_, gene_name)
      }
      return(gRNA_df)
    }, error = function(e) {
      base::cat("Error for locus_tag", locus_tag, ":", base::conditionMessage(e), "\n")
      failed_locus_tags <<- base::c(failed_locus_tags, locus_tag)
      return(NULL)
    })
  })

  # Remove NULL entries and combine
  gRNA_list <- base::Filter(base::Negate(base::is.null), gRNA_list)
  if (base::length(gRNA_list) == 0) {
    base::stop("No gRNAs generated for any locus_tag")
  }

  # Log summary of failed locus_tags
  if (base::length(failed_locus_tags) > 0) {
    base::cat("Summary: Failed to generate gRNA for", base::length(failed_locus_tags), "out of", base::nrow(genbank_table), "locus_tags\n")
    base::cat("Failed locus_tags:", base::paste(failed_locus_tags, collapse = ", "), "\n")
  }

  gRNA_library <- base::do.call(rbind, gRNA_list)

  # Generate systematic gRNA names
  base::cat("Generating systematic gRNA names for library\n")
  gRNA_library <- generate_gRNA_names(
    gRNA_df = gRNA_library,
    prefix = name_prefix,
    nuclease = nuclease,
    rank_by = "composite_score"
  )

  return(gRNA_library)
}
