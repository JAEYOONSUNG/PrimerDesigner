#' Build Bowtie Index for a Genome
#'
#' Builds a Bowtie index for a genome specified by a GenBank accession number.
#' The index is created from a previously downloaded .fna file.
#'
#' @param directory A character string specifying the directory containing the .fna file (default: NULL, uses current working directory).
#' @param accession A character string specifying the GenBank accession number (e.g., "GCF_030376745.1").
#' @return None. Creates Bowtie index files in the specified directory.
#' @examples
#' \dontrun{
#' run_bowtie_build(accession = "GCF_030376745.1")
#' }
#' @export

run_bowtie_build <- function(directory = NULL, accession) {
  # Step 1: Handle default directory
  if (base::is.null(directory)) {
    directory <- base::getwd()
    base::message("No directory specified. Using current working directory: ", directory)
  } else {
    base::message("Using specified directory: ", directory)
  }

  # Step 2: Define file paths
  sub_dir <- "seqs_srcdir"
  fa_file_name <- base::paste0(accession, ".fna")
  genome_fasta_file <- base::file.path(directory, sub_dir, fa_file_name)
  bowtie_out_dir <- base::file.path(directory, base::paste0(accession, "_bowtie"))
  base::cat("Bowtie output directory:", bowtie_out_dir, "\n")

  # Step 3: Check if FASTA file exists
  if (!base::file.exists(genome_fasta_file)) {
    base::stop("FASTA file not found: ", genome_fasta_file)
  }

  # Step 4: Execute Bowtie index build
  Rbowtie::bowtie_build(
    genome_fasta_file,
    outdir = bowtie_out_dir,
    force = TRUE,
    prefix = accession
  )
}
