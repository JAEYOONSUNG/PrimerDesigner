#' Filter gRNAs by Methylation Site Conflicts
#'
#' Checks whether methylation recognition motifs overlap with the gRNA
#' target region (spacer + PAM context). gRNAs with methylation conflicts
#' are flagged and optionally removed.
#'
#' Methyltransferases recognize specific DNA motifs and methylate bases
#' within them. If a methylation motif overlaps the PAM or spacer region,
#' the methylated base(s) may interfere with Cas nuclease binding or
#' cutting. This function identifies such conflicts.
#'
#' Both strands are checked: the motif is searched on the sense strand
#' and its reverse complement is searched on the antisense strand.
#'
#' @param gRNA_df A data frame containing gRNA data. Must have
#'   \code{protospacer} and \code{pam} columns.
#' @param methylation_patterns Character vector of methylation recognition
#'   sequences using IUPAC ambiguity codes. Examples:
#'   \code{c("GCCAT", "CCANNNNNTTG", "GATC")}.
#'   IUPAC codes: N = any, R = A/G, Y = C/T, S = G/C, W = A/T,
#'   K = G/T, M = A/C, B = C/G/T, D = A/G/T, H = A/C/T, V = A/C/G.
#' @param nuclease Character. Nuclease name to determine PAM side.
#'   Supported: "GeoCas9", "SpCas9" (PAM 3' of spacer),
#'   "FnCas12a", "FisCasI_B" (PAM 5' of spacer).
#'   Default: "GeoCas9".
#' @param exclude Logical. If TRUE (default), remove gRNAs with
#'   methylation conflicts. If FALSE, only add flag columns without
#'   removing any rows.
#' @return The input data frame with added columns:
#'   \itemize{
#'     \item \code{methylation_conflict}: Logical. TRUE if any motif was found.
#'     \item \code{methylation_motifs}: Character. Comma-separated list of
#'       matched motif(s) and their strand (e.g., "GCCAT(+), CCANNNNNTTG(-)").
#'   }
#'   If \code{exclude = TRUE}, rows with conflicts are removed.
#' @examples
#' \dontrun{
#' # Geobacillus stearothermophilus methylation patterns
#' gRNA_filtered <- filter_methylation_sites(
#'   gRNA_df = gRNA_lib,
#'   methylation_patterns = c("GCCAT", "CCANNNNNTTG"),
#'   nuclease = "GeoCas9",
#'   exclude = TRUE
#' )
#'
#' # Just flag without removing (for inspection)
#' gRNA_flagged <- filter_methylation_sites(
#'   gRNA_df = gRNA_lib,
#'   methylation_patterns = c("GCCAT", "CCANNNNNTTG"),
#'   nuclease = "GeoCas9",
#'   exclude = FALSE
#' )
#' table(gRNA_flagged$methylation_conflict)
#' }
#' @export
filter_methylation_sites <- function(
    gRNA_df,
    methylation_patterns,
    nuclease = "GeoCas9",
    exclude = TRUE
) {
  # --- Input validation ---
  if (!base::is.data.frame(gRNA_df) || base::nrow(gRNA_df) == 0) {
    base::warning("Empty gRNA data frame. Returning as-is.")
    return(gRNA_df)
  }

  required_cols <- base::c("protospacer", "pam")
  missing_cols <- base::setdiff(required_cols, base::colnames(gRNA_df))
  if (base::length(missing_cols) > 0) {
    base::stop("gRNA_df missing required columns: ",
               base::paste(missing_cols, collapse = ", "),
               ". Both 'protospacer' and 'pam' are needed for methylation checking.")
  }

  if (base::is.null(methylation_patterns) ||
      base::length(methylation_patterns) == 0) {
    gRNA_df$methylation_conflict <- FALSE
    gRNA_df$methylation_motifs <- ""
    return(gRNA_df)
  }

  # --- Determine PAM side ---
  pam_side <- .get_pam_side(nuclease)

  base::cat("=== Methylation Site Filter ===\n")
  base::cat("  Nuclease:", nuclease, "(PAM", pam_side, "of spacer)\n")
  base::cat("  Patterns:", base::paste(methylation_patterns, collapse = ", "), "\n")
  base::cat("  Checking", base::nrow(gRNA_df), "gRNAs...\n")

  # --- Build regex patterns (forward + reverse complement) ---
  pattern_info <- base::list()
  for (pat in methylation_patterns) {
    pat_upper <- base::toupper(pat)

    # Forward motif regex
    fwd_regex <- .iupac_to_regex(pat_upper)

    # Reverse complement motif regex
    rc_seq <- base::as.character(
      Biostrings::reverseComplement(Biostrings::DNAString(pat_upper))
    )
    rc_regex <- .iupac_to_regex(rc_seq)

    pattern_info[[base::length(pattern_info) + 1]] <- base::list(
      original = pat_upper,
      fwd_regex = fwd_regex,
      rc_regex = rc_regex,
      rc_seq = rc_seq
    )
  }

  # --- Check each gRNA ---
  n <- base::nrow(gRNA_df)
  conflict <- base::logical(n)
  conflict_detail <- base::character(n)

  for (i in base::seq_len(n)) {
    pam_seq <- base::toupper(gRNA_df$pam[i])
    spacer_seq <- base::toupper(gRNA_df$protospacer[i])

    # Skip if pam or spacer is NA/empty
    if (base::is.na(pam_seq) || base::is.na(spacer_seq) ||
        base::nchar(pam_seq) == 0 || base::nchar(spacer_seq) == 0) {
      conflict[i] <- FALSE
      conflict_detail[i] <- ""
      next
    }

    # Build genomic context based on PAM side
    # 3prime: 5'-spacer-PAM-3'  (GeoCas9, SpCas9)
    # 5prime: 5'-PAM-spacer-3'  (FnCas12a, FisCasI_B)
    if (pam_side == "3prime") {
      context <- base::paste0(spacer_seq, pam_seq)
    } else {
      context <- base::paste0(pam_seq, spacer_seq)
    }

    # Check each methylation motif
    matched <- base::character(0)

    for (pinfo in pattern_info) {
      # Check forward strand
      if (base::grepl(pinfo$fwd_regex, context, ignore.case = FALSE, perl = TRUE)) {
        matched <- base::c(matched, base::paste0(pinfo$original, "(+)"))
      }
      # Check reverse complement (antisense strand)
      if (pinfo$fwd_regex != pinfo$rc_regex) {
        if (base::grepl(pinfo$rc_regex, context, ignore.case = FALSE, perl = TRUE)) {
          matched <- base::c(matched, base::paste0(pinfo$original, "(-)"))
        }
      }
    }

    conflict[i] <- base::length(matched) > 0
    conflict_detail[i] <- base::paste(matched, collapse = ", ")
  }

  # --- Add columns ---
  gRNA_df$methylation_conflict <- conflict
  gRNA_df$methylation_motifs <- conflict_detail

  n_conflict <- base::sum(conflict)
  base::cat("  Result:", n_conflict, "of", n, "gRNAs have methylation conflicts\n")

  if (n_conflict > 0) {
    # Show a few examples
    conflict_idx <- base::which(conflict)
    n_show <- base::min(5, n_conflict)
    for (j in base::seq_len(n_show)) {
      idx <- conflict_idx[j]
      tag_str <- ""
      if ("locus_tag" %in% base::colnames(gRNA_df)) {
        tag_str <- base::paste0(gRNA_df$locus_tag[idx], " ")
      }
      base::cat("    ", tag_str, gRNA_df$protospacer[idx],
                " | PAM=", gRNA_df$pam[idx],
                " | ", conflict_detail[idx], "\n", sep = "")
    }
    if (n_conflict > n_show) {
      base::cat("    ... and", n_conflict - n_show, "more\n")
    }
  }

  # --- Exclude if requested ---
  if (exclude && n_conflict > 0) {
    gRNA_df <- gRNA_df[!conflict, , drop = FALSE]
    base::rownames(gRNA_df) <- NULL
    base::cat("  Removed", n_conflict, "gRNAs →", base::nrow(gRNA_df), "remaining\n")
  }

  return(gRNA_df)
}


# ============================================================
# Internal helpers
# ============================================================

#' Convert IUPAC ambiguity sequence to regex pattern
#'
#' @param iupac_seq Character. DNA sequence with IUPAC codes.
#' @return Character. Regex pattern string.
#' @keywords internal
.iupac_to_regex <- function(iupac_seq) {
  iupac_map <- base::c(
    "A" = "A", "C" = "C", "G" = "G", "T" = "T",
    "R" = "[AG]",  "Y" = "[CT]",  "S" = "[GC]",  "W" = "[AT]",
    "K" = "[GT]",  "M" = "[AC]",
    "B" = "[CGT]", "D" = "[AGT]", "H" = "[ACT]", "V" = "[ACG]",
    "N" = "[ACGT]"
  )

  chars <- base::strsplit(base::toupper(iupac_seq), "")[[1]]
  regex_parts <- base::character(base::length(chars))

  for (i in base::seq_along(chars)) {
    ch <- chars[i]
    if (ch %in% base::names(iupac_map)) {
      regex_parts[i] <- iupac_map[ch]
    } else {
      # Unknown character — escape and keep literal
      regex_parts[i] <- ch
    }
  }

  return(base::paste0(regex_parts, collapse = ""))
}


#' Get PAM side from nuclease name
#'
#' @param nuclease Character. Nuclease name.
#' @return Character. "3prime" or "5prime".
#' @keywords internal
.get_pam_side <- function(nuclease) {
  pam_side_map <- base::c(
    "GeoCas9"   = "3prime",
    "SpCas9"    = "3prime",
    "GtCas9"    = "3prime",
    "FnCas12a"  = "5prime",
    "FisCasI_B" = "5prime"
  )

  if (nuclease %in% base::names(pam_side_map)) {
    return(pam_side_map[nuclease])
  }

  base::warning("Unknown nuclease '", nuclease, "'. Defaulting to pam_side = '3prime'.")
  return("3prime")
}
