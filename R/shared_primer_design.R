#' End-to-end shared-primer CRISPR deletion design
#'
#' Canonical top-level entry point for the shared-primer workflow. Wraps
#' \code{\link{design_shared_grna_and_deletion}} with the full set of knobs
#' — target resolution, shared UF/UR/DF/DR core selection, gRNA + cloning
#' primer design, Gibson overhang assembly, combined construct GenBank
#' export, colony-PCR (\emph{check}) primer pair design, and a 3-sheet
#' Excel workbook (primer_order / check_primers / final_construct_groups).
#'
#' This function is the recommended user-facing entry point. The underlying
#' \code{design_shared_grna_and_deletion()} is retained as the lower-level
#' implementation and remains exported for back-compat.
#'
#' @param target_table A data frame with at least \code{genome_id},
#'   \code{genbank_file}, \code{locus_tag}. Typically produced by
#'   \code{\link{find_target_across_genomes}}.
#' @param ... Any argument accepted by
#'   \code{\link{design_shared_grna_and_deletion}}. Common ones:
#'   \describe{
#'     \item{\code{overlap_policy}}{\code{"strict"} (default) preserves all
#'        CDS overlaps at the deletion boundary and enforces inner primers
#'        (UR, DF) exactly at the stop / start codon; \code{"flexible"}
#'        permits widening for shared primer coverage.}
#'     \item{\code{min_arm_bp}, \code{max_arm_bp}}{Homology arm bounds.}
#'     \item{\code{upstream_bp}, \code{downstream_bp}}{Preferred arm
#'        lengths used by the main UF/DR primer search.}
#'     \item{\code{check_outer_pad}, \code{check_search_window},
#'           \code{check_search_window_max}}{Colony-PCR primer placement:
#'        starting window just outside the effective homology arm (default
#'        50 / 800 bp) with progressive widening up to
#'        \code{check_search_window_max} (default 6000 bp) when a primer
#'        shared across every target genome cannot be found within the
#'        initial window.}
#'     \item{\code{check_tm_target}, \code{check_tm_tolerance},
#'           \code{check_pair_dtm_max}}{Check primer Tm constraints (55 ± 3
#'        °C, pair ΔTm ≤ 2 °C by default).}
#'     \item{\code{output_file}}{Destination .xlsx path. When set, the
#'        3-sheet workbook is written.}
#'     \item{\code{combined_construct_output_dir}}{Destination directory
#'        for the per-group combined construct GenBank files.}
#'   }
#' @return A list returned by
#'   \code{\link{design_shared_grna_and_deletion}}, additionally carrying
#'   \code{check_primer_cores} (per-cluster cF / cR) and
#'   \code{check_primer_pairs} (per-genome diagnostic PCR band sizes).
#' @seealso \code{\link{visualize_shared_design}},
#'   \code{\link{find_target_across_genomes}}.
#' @examples
#' \dontrun{
#' targets <- find_target_across_genomes(
#'   genbank_dirs = "path/to/gbk",
#'   query        = "TIGR02679",
#'   query_type   = "product"
#' )
#' res <- shared_primer_design(
#'   target_table   = targets[, c("genome_id", "genbank_file", "locus_tag")],
#'   nuclease       = "GeoCas9",
#'   overlap_policy = "strict",
#'   upstream_bp    = 500,
#'   downstream_bp  = 500,
#'   min_arm_bp     = 150L,
#'   max_arm_bp     = 3000L,
#'   grna_vector_file  = "pJET.gb",
#'   grna_start        = 7823,
#'   grna_end          = 7850,
#'   combined_vector_file       = "pJET.gb",
#'   combined_grna_start        = 7823,
#'   combined_grna_end          = 7850,
#'   combined_deletion_start    = 3977,
#'   combined_deletion_end      = 4986,
#'   combined_construct_output_dir = "out/constructs",
#'   output_file                   = "out/design.xlsx"
#' )
#' visualize_shared_design(
#'   result           = res,
#'   genbank_dir      = "path/to/gbk",
#'   construct_gbk_dir = "out/constructs",
#'   output_file      = "out/summary.pdf"
#' )
#' }
#' @export
shared_primer_design <- function(target_table, ...) {
  design_shared_grna_and_deletion(target_table = target_table, ...)
}
