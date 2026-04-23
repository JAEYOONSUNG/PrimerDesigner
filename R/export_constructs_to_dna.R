#' Convert every `.gbk` in a directory to SnapGene `.dna` format
#'
#' Runs the SnapGene CLI once per file, killing the SnapGene GUI ONCE up
#' front so the CLI can acquire the files (SnapGene holds an exclusive
#' lock while open). Use this as a post-processing step when you want
#' native SnapGene `.dna` copies of the combined construct files produced
#' by \code{\link{shared_primer_design}} — the Primers panel, custom
#' colors, and other SnapGene-specific metadata round-trip cleanly, which
#' they do not when SnapGene imports a GenBank file.
#'
#' @param dir_path Directory containing one or more `.gbk` / `.gb` files.
#' @param kill_snapgene Logical. If \code{TRUE} (default) and SnapGene is
#'   running, it is closed first so the CLI can convert.
#' @return Character vector of `.dna` paths that were produced.
#' @examples
#' \dontrun{
#' export_constructs_to_dna("out/constructs")
#' }
#' @export
export_constructs_to_dna <- function(dir_path, kill_snapgene = TRUE) {
  if (!base::dir.exists(dir_path)) {
    base::stop("directory not found: ", dir_path)
  }
  .snapgene_export_dna_dir(dir_path, kill_snapgene = kill_snapgene)
}
