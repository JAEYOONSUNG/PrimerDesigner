# End-to-end smoke test for the shared-design workflow.
#
# What it does:
#   1. Source every R/ file so syntax errors surface immediately.
#   2. Run find_target_across_genomes() in non-interactive mode against a
#      small GenBank directory.
#   3. Run visualize_shared_design() against a previously saved pipeline
#      result (so we don't need to run the full crisprVerse pipeline in the
#      test, which requires BSgenome index setup).
#
# Usage:
#   Rscript inst/tests/test_shared_design_e2e.R \
#       /path/to/genome_gbk /path/to/design_result.rds /path/to/combined_gbk
#
# The test prints PASS/FAIL lines and exits with non-zero status on any
# failure, so it can be wired into CI later.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript test_shared_design_e2e.R <gbk_dir> <result.rds> <construct_gbk_dir>")
}
gbk_dir          <- args[1]
result_rds       <- args[2]
construct_gbk_dir <- args[3]

# Locate package root: prefer script path, fall back to current directory.
self_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
pkg_root <- if (!is.null(self_path)) {
  normalizePath(file.path(dirname(self_path), "..", ".."), mustWork = FALSE)
} else {
  getwd()
}
if (!dir.exists(file.path(pkg_root, "R"))) pkg_root <- getwd()

report <- function(ok, msg) {
  cat(if (ok) "PASS" else "FAIL", msg, "\n")
  if (!ok) quit(status = 1)
}

# --- 1) Source every R/ file --------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(ggtext); library(ggdendro)
})
r_files <- list.files(file.path(pkg_root, "R"),
                      pattern = "\\.R$", full.names = TRUE)
for (rf in r_files) {
  ok <- tryCatch({ source(rf, local = FALSE); TRUE },
                 error = function(e) { cat("    ", conditionMessage(e), "\n"); FALSE })
  report(ok, paste0("source ", basename(rf)))
}

# --- 2) find_target_across_genomes -------------------------------------------
tt <- tryCatch(
  find_target_across_genomes(
    genbank_dirs = gbk_dir,
    query = "TIGR02680",
    query_type = "product",
    interactive = FALSE,
    one_per_genome = FALSE),
  error = function(e) { cat("    ", conditionMessage(e), "\n"); NULL })
report(!is.null(tt) && nrow(tt) > 0,
       sprintf("find_target_across_genomes returned %d rows",
               if (is.null(tt)) 0 else nrow(tt)))

# --- 3) visualize_shared_design -----------------------------------------------
if (file.exists(result_rds) && dir.exists(construct_gbk_dir)) {
  res <- readRDS(result_rds)
  out_pdf <- tempfile(fileext = ".pdf")
  ok <- tryCatch({
    visualize_shared_design(
      result = res,
      genbank_dir = gbk_dir,
      construct_gbk_dir = construct_gbk_dir,
      output_file = out_pdf)
    file.exists(out_pdf) && file.info(out_pdf)$size > 1000
  }, error = function(e) { cat("    ", conditionMessage(e), "\n"); FALSE })
  report(ok, sprintf("visualize_shared_design wrote %s (%d bytes)",
                     out_pdf,
                     if (file.exists(out_pdf)) as.integer(file.info(out_pdf)$size) else 0L))
  if (ok) file.remove(out_pdf)
} else {
  cat("SKIP visualize_shared_design (no saved result fixture provided)\n")
}

cat("\nAll checks completed.\n")
