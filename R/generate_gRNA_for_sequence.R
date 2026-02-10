#' Generate gRNA for a Single Sequence
#'
#' Generates a gRNA list for a single target DNA sequence using the specified nuclease.
#'
#' @param target_sequence A character string specifying the target DNA sequence.
#' @param nuclease A character string specifying the nuclease to use ("GeoCas9", "FnCas12a", or "FisCasI_B"; default: "GeoCas9").
#' @param bowtie_index A character string specifying the path to the Bowtie index.
#' @param genome A BSgenome object specifying the reference genome.
#' @param position_range A numeric vector of length 2 specifying the position range (as percentage) for gRNA filtering (default: c(10, 90)).
#' @param tm_range A numeric vector of length 2 specifying the Tm range for gRNA filtering (default: c(50, 70)).
#' @param strand A character string specifying the strand to filter ("both", "5", or "3"; default: "both").
#' @param aligner A character string specifying the aligner to use ("bowtie", "bwa", or "biostrings"; default: "bowtie").
#' @return A data frame containing the gRNA list for the target sequence.
#' @examples
#' \dontrun{
#' target_seq <- "ttggcattgacaggtacagaccgcgtcaaacgcggcatggcggaaatgcaaaaaggcggcgtcattatggacgtcgtcaatgcagagcaagcgaagattgctgaggcggcaggggctgtcgcagtcatggcgctcgagcgtgtcccggcagacattcgcgccgctggcggtgtcgcgcgcatggctgatccgacgatcattgaagaagtgatgaacgccgtatcaatcccagtcatggcgaaagtgcgcatcgggcattatgtggaagcgcgtgttttagaggcgctcggcvtcgactatattgacgaaagtgaagtattgacgccggctgatgaagagttccatattgacaaacggcagtttacggtcccatttgtatgcggttgccgcgacttaggagaggccgcccgccgcattgctgaaggggcatcgatgttgcggacaaaaggggagccagggacaggaaacatcgtcgaggccgttcgccatatgcgcaaagtcaacgcgcaaatccgcaaagttgtcagcatgagcgaagacgaacttgtcgccgaggcgaaacagcttggggctccggttgaagtgctgcgtgaaatcaaacggcttggccgcctcccggtcgtcaacttcgccgccggcggtgtcgcgacaccagctgacgccgcgctcatgatgcacttaggcgccgacggtgtctfcgtcggttcgggcatttttaaatcggaaaatccggaaaaatacgctcgtgcgatcgttgaagcgacgactcattatgaagactatgagctgatcgcccatctatcgaaagggctgggcggcgcaatgcgcggcatcgatvtcgcgacactgctgccggagcatcggatgcaagaacgaggctggtaa"
#' gRNA_list <- generate_gRNA_for_sequence(
#'   target_sequence = target_seq,
#'   nuclease = "GeoCas9",
#'   bowtie_index = "path/to/bowtie_index",
#'   genome = NULL,
#'   position_range = c(10, 90),
#'   tm_range = c(50, 70),
#'   aligner = "bowtie"
#' )
#' print(gRNA_list)
#' }
#' @export
# Modified generate_gRNA_for_sequence

generate_gRNA_for_sequence <- function(
    target_sequence,
    nuclease = "GeoCas9",
    bowtie_index,
    genome,
    position_range = c(10, 90),
    tm_range = c(50, 70),
    strand = "both",
    aligner = "bowtie"
) {
  base::cat("Generating gRNA list for target sequence\n")
  base::cat("Input target_sequence class:", base::class(target_sequence), "\n")

  # Clean target sequence (minimize warnings if already clean)
  target_sequence_clean <- base::toupper(target_sequence)
  if (!base::grepl("^[ATGCN]+$", target_sequence_clean)) {
    base::warning("Target sequence contains invalid characters (only ATGCN allowed). Replacing invalid characters with N.")
    log_invalid_sequence(target_sequence, "custom_target")
    target_sequence_clean <- base::gsub("[^ATGCN]", "N", target_sequence_clean)
  }

  # Validate cleaned sequence
  if (!is_valid_dna(target_sequence_clean)) {
    base::stop("Error: Invalid target sequence provided even after cleaning")
  }

  # BSgenome verification (skip here, delegated to nuclease_with_parameter)
  base::cat("Skipping BSgenome verification in generate_gRNA_for_sequence; handled in nuclease_with_parameter\n")

  # Pass cleaned sequence to nuclease_with_parameter
  base::cat("Passing target_sequence_clean (class:", base::class(target_sequence_clean), ") to nuclease_with_parameter\n")
  gRNA_list <- nuclease_with_parameter(
    target_sequence = target_sequence_clean,
    strand = strand,
    bowtie_index = bowtie_index,
    genome = genome,
    nuclease = nuclease,
    position_range = position_range,
    tm_range = tm_range,
    aligner = aligner
  )

  return(gRNA_list)
}
