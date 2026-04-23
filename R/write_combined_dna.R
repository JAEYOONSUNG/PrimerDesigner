#' Write a combined-construct `.dna` by directly editing the source plasmid
#'
#' Preserves every SnapGene-native annotation (Primers panel, Features
#' colors, StrandColors, history, etc.) that is dropped when round-tripping
#' through GenBank. Uses the bundled Python editor at
#' `inst/python/dna_editor.py` via reticulate.
#'
#' Parameters:
#'   plasmid_path      source `.dna` file
#'   output_path       destination `.dna`
#'   edits             list of lists: each with start, end, replacement
#'   new_features      list of lists: name, type, start, end, color, directionality
#'   new_primers       list of lists: name, sequence, start, end, strand,
#'                     annealed_bases, description
#'
#' `edits` are applied in order. For combined constructs with both a gRNA
#' stuffer swap and a donor stuffer swap, pass the higher-coordinate edit
#' first so the lower-coordinate edit's coordinates don't need re-anchoring.
#' @keywords internal
.shared_write_combined_dna <- function(
    plasmid_path, output_path,
    edits = base::list(),
    new_features = base::list(),
    new_primers = base::list()
) {
  if (!base::requireNamespace("reticulate", quietly = TRUE)) {
    base::stop("reticulate is required for .dna writing; install it first.")
  }
  # Tell reticulate the Python deps we need. py_require is a no-op if the
  # package is already present in the active env.
  base::tryCatch(reticulate::py_require("xmltodict"),
                 error = function(e) NULL)

  py_src <- base::system.file("python", "dna_editor.py",
                                package = "PrimerDesigner")
  if (!base::nzchar(py_src) || !base::file.exists(py_src)) {
    # devtools::load_all fallback — resolve against the source tree.
    root <- base::tryCatch(base::find.package("PrimerDesigner"),
                            error = function(e) NULL)
    for (r in base::c(root, base::dirname(root))) {
      cand <- base::file.path(r, "inst", "python", "dna_editor.py")
      if (!base::is.null(r) && base::file.exists(cand)) { py_src <- cand; break }
      cand <- base::file.path(r, "python", "dna_editor.py")
      if (!base::is.null(r) && base::file.exists(cand)) { py_src <- cand; break }
    }
  }
  if (!base::nzchar(py_src) || !base::file.exists(py_src)) {
    base::stop("bundled dna_editor.py not found in PrimerDesigner package")
  }
  ed <- reticulate::py_run_file(py_src, convert = TRUE)
  ed$build_combined_construct(plasmid_path, output_path,
                                edits, new_features, new_primers)
  base::invisible(output_path)
}

# --------------------------------------------------------------------------
# Helpers to build the edits / features / primers specs from a design result
# --------------------------------------------------------------------------

#' Reverse-complement a DNA string.
#' @keywords internal
.shared_rc <- function(s) {
  base::chartr("ACGTacgtNnKkMmRrYySsWwBbDdHhVv",
                "TGCAtgcaNnMmKkYyRrSsWwVvHhDdBb",
                base::paste(base::rev(base::strsplit(s, "")[[1]]),
                             collapse = ""))
}

#' Find first occurrence of `oligo` on forward or reverse strand, falling
#' back to progressively shorter 3' anchors. Returns list(start, end,
#' strand, anchor) with 1-based inclusive coords; NULL if no match.
#' @keywords internal
.shared_find_primer_binding <- function(seq, oligo, min_anchor = 15L) {
  seq_u <- base::toupper(base::as.character(seq))
  oligo_u <- base::toupper(base::as.character(oligo))
  n <- base::nchar(oligo_u)
  for (k in base::seq.int(n, min_anchor, by = -1L)) {
    anchor <- base::substr(oligo_u, n - k + 1L, n)
    pos <- base::regexpr(anchor, seq_u, fixed = TRUE)[1]
    if (pos > 0) return(base::list(start = pos, end = pos + k - 1L,
                                     strand = "+", anchor = anchor))
    anchor_rc <- .shared_rc(anchor)
    pos <- base::regexpr(anchor_rc, seq_u, fixed = TRUE)[1]
    if (pos > 0) return(base::list(start = pos, end = pos + k - 1L,
                                     strand = "-", anchor = anchor))
  }
  NULL
}

#' Build spec dicts for one combined construct and return the list of
#' (edits, new_features, new_primers) ready to pass to
#' `.shared_write_combined_dna()`.
#'
#' Inputs:
#'   plasmid_seq  character(1) — original plasmid sequence (for Tm not needed here)
#'   grna_start, grna_end          — 1-based stuffer region for the gRNA spacer
#'   donor_start, donor_end        — 1-based stuffer region for the donor insert
#'   protospacer                    — string to insert at the gRNA stuffer
#'   donor_insert                   — string to insert at the donor stuffer
#'   arm_upstream_bp, arm_downstream_bp — lengths within donor_insert
#'   locus_tag                       — used in primer names
#'   method_abbrev                   — "inv" (Gibson) or "OA" (Golden Gate)
#'   primer_cluster_nums             — named list with UF, UR, DF, DR, gRNA integer cluster numbers
#'   arm_primer_cores                — named list UF/UR/DF/DR -> full primer sequences
#'   grna_primer_pair                — list(primer_F, primer_R)
#' @keywords internal
.shared_build_combined_dna_spec <- function(
    plasmid_seq, grna_start, grna_end, donor_start, donor_end,
    protospacer, donor_insert,
    arm_upstream_bp, arm_downstream_bp,
    locus_tag, method_abbrev, primer_cluster_nums,
    arm_primer_cores, grna_primer_pair
) {
  # Edits: apply HIGHER-coordinate edit first so lower-coord edit keeps
  # its original positions.
  if (grna_start > donor_end) {
    edits <- base::list(
      base::list(start = base::as.integer(grna_start),
                  end   = base::as.integer(grna_end),
                  replacement = protospacer),
      base::list(start = base::as.integer(donor_start),
                  end   = base::as.integer(donor_end),
                  replacement = donor_insert)
    )
  } else {
    edits <- base::list(
      base::list(start = base::as.integer(donor_start),
                  end   = base::as.integer(donor_end),
                  replacement = donor_insert),
      base::list(start = base::as.integer(grna_start),
                  end   = base::as.integer(grna_end),
                  replacement = protospacer)
    )
  }

  # Final-construct coordinates for annotation (after all edits).
  len_diff_donor <- base::nchar(donor_insert) - (donor_end - donor_start + 1L)
  len_diff_grna  <- base::nchar(protospacer)  - (grna_end  - grna_start  + 1L)
  # Donor coords in final construct (unchanged start, new end).
  donor_final_start <- donor_start
  donor_final_end   <- donor_start + base::nchar(donor_insert) - 1L
  arm_up_s <- donor_final_start
  arm_up_e <- donor_final_start + arm_upstream_bp - 1L
  arm_dn_s <- arm_up_e + 1L
  arm_dn_e <- donor_final_end
  # gRNA coords in final construct. If gRNA edit came after donor edit in
  # coordinate space, its start shifts by len_diff_donor.
  grna_final_start <- grna_start + len_diff_donor
  grna_final_end   <- grna_final_start + base::nchar(protospacer) - 1L

  new_features <- base::list(
    base::list(name = "upstream_arm",  type = "misc_feature",
                start = arm_up_s, end = arm_up_e,
                color = "#aecd87", directionality = 0L),
    base::list(name = "downstream_arm", type = "misc_feature",
                start = arm_dn_s, end = arm_dn_e,
                color = "#ffad73", directionality = 0L),
    base::list(name = "protospacer", type = "misc_feature",
                start = grna_final_start, end = grna_final_end,
                color = "#3cb371", directionality = 0L)
  )

  # Build the final-construct sequence in R so we can locate each primer.
  fs <- plasmid_seq
  if (grna_start > donor_end) {
    fs <- base::paste0(base::substr(fs, 1L, grna_start - 1L),
                        protospacer,
                        base::substr(fs, grna_end + 1L, base::nchar(fs)))
    fs <- base::paste0(base::substr(fs, 1L, donor_start - 1L),
                        donor_insert,
                        base::substr(fs, donor_end + 1L, base::nchar(fs)))
  } else {
    fs <- base::paste0(base::substr(fs, 1L, donor_start - 1L),
                        donor_insert,
                        base::substr(fs, donor_end + 1L, base::nchar(fs)))
    new_grna_start <- grna_start + len_diff_donor
    new_grna_end   <- grna_end   + len_diff_donor
    fs <- base::paste0(base::substr(fs, 1L, new_grna_start - 1L),
                        protospacer,
                        base::substr(fs, new_grna_end + 1L, base::nchar(fs)))
  }

  arm_name <- function(role, n) {
    fr <- if (role %in% c("UF", "DF")) "F" else "R"
    base::sprintf("%s_%s_clstr%d_%s", locus_tag, role, n, fr)
  }
  grna_name <- function(fr) base::sprintf("sgRNA_%s_clstr%d_%s_%s",
    locus_tag, primer_cluster_nums$gRNA, method_abbrev, fr)

  primer_specs <- base::list(
    UF = base::list(name = arm_name("UF", primer_cluster_nums$UF),
                     seq = arm_primer_cores$UF),
    UR = base::list(name = arm_name("UR", primer_cluster_nums$UR),
                     seq = arm_primer_cores$UR),
    DF = base::list(name = arm_name("DF", primer_cluster_nums$DF),
                     seq = arm_primer_cores$DF),
    DR = base::list(name = arm_name("DR", primer_cluster_nums$DR),
                     seq = arm_primer_cores$DR),
    gRNA_F = base::list(name = grna_name("F"),
                         seq = grna_primer_pair$primer_F),
    gRNA_R = base::list(name = grna_name("R"),
                         seq = grna_primer_pair$primer_R)
  )

  new_primers <- base::list()
  for (role in base::names(primer_specs)) {
    sp <- primer_specs[[role]]
    hit <- .shared_find_primer_binding(fs, sp$seq, min_anchor = 15L)
    if (base::is.null(hit)) next
    new_primers[[base::length(new_primers) + 1L]] <- base::list(
      name            = sp$name,
      sequence        = sp$seq,
      start           = base::as.integer(hit$start),
      end             = base::as.integer(hit$end),
      strand          = hit$strand,
      annealed_bases  = hit$anchor,
      description     = base::sprintf("designed %s", role)
    )
  }

  base::list(edits = edits,
              new_features = new_features,
              new_primers = new_primers)
}
