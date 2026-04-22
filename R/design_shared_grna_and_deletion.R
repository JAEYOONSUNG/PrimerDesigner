#' Design Shared gRNA and Deletion Primer Strategy Across Multiple Genomes
#'
#' High-level multi-genome workflow that jointly designs:
#' \itemize{
#'   \item shared or subgroup-specific gRNA candidates,
#'   \item shared or subgroup-specific deletion primer cores (\code{UF/UR/DF/DR}),
#'   \item grouped donor constructs when the same left+right arm insert is identical,
#'   \item grouped final combined constructs when both donor insert and gRNA are shared.
#' }
#'
#' This function is intended for the experimental case where the same functional
#' target is edited across several related genomes and the goal is to minimize
#' synthesis burden while keeping the design exact and auditable.
#'
#' The grouping logic is deliberately strict:
#' \enumerate{
#'   \item prefer one design shared across all genomes,
#'   \item if impossible, fall back to subgroup-shared designs,
#'   \item finally fall back to strain-specific designs,
#'   \item enforce exact uniqueness in each intended target genome when requested.
#' }
#'
#' gRNA selection uses exact shared spacer+PAM candidates extracted directly from
#' the annotated target CDS sequences. This is intentionally separate from the
#' existing single-genome CrisprVerse workflow. The output here is suitable for
#' multi-genome planning and can then be converted into full cloning primers via
#' \code{design_cloning_primers()}.
#'
#' @param target_table Data frame with columns \code{genome_id},
#'   \code{genbank_file}, and \code{locus_tag}.
#' @param nuclease Character. Nuclease name. Currently supports
#'   \code{GeoCas9}, \code{SpCas9}, \code{GtCas9}, \code{FnCas12a},
#'   \code{FisCasI_B}.
#' @param top_n Integer. Number of gRNA groups to keep per shared/subgroup
#'   search pass. Default 1.
#' @param position_range Numeric length-2 vector. Desired position window on the
#'   target CDS as percentage.
#' @param tm_range Numeric length-2 vector. Desired spacer Tm range.
#' @param strand Character. \code{"both"}, \code{"5"}, or \code{"3"}.
#' @param require_unique Logical. If \code{TRUE}, both gRNA context candidates
#'   and deletion primer cores must be exact 1-hit matches in each intended
#'   target genome.
#' @param n_mismatches Integer. Off-target mismatch depth for shared gRNA
#'   evaluation. Default 0 for fast exact-match-only mode.
#' @param max_n1 Optional integer. If provided, reject shared gRNA candidates
#'   whose worst intended-target 1-mismatch off-target count exceeds this value.
#' @param max_n1_per_primer Optional integer. If provided, any deletion primer
#'   whose worst-case 1-mismatch genome hit count exceeds this value triggers
#'   a warning (per-genome values are still surfaced in \code{primer_cores}).
#' @param overlap_policy Either \code{"strict"} (default) or \code{"flexible"}.
#'   Controls how overlaps between the target CDS and its neighbors are
#'   handled at the deletion boundary:
#'   \itemize{
#'     \item \code{"strict"}: the deletion is trimmed so no bp of a
#'       neighboring CDS is removed. If an upstream gene extends INTO the
#'       target, deletion starts after that gene ends; same for downstream
#'       overlap. This is the safe default when the user says "only delete
#'       the target, keep everything else intact".
#'     \item \code{"flexible"}: the deletion spans the full target
#'       coordinates, which may clip up to a few bp of an overlapping
#'       neighbor. A warning lists how many bp are clipped per genome.
#'   }
#' @param preserve_upstream_gene Logical or \code{"auto"}. Legacy knob used
#'   to push the deletion 5' end even further (stop at the full upstream CDS
#'   end, not just the overlap point). Usually left at the default now that
#'   \code{overlap_policy} handles the common case.
#' @param preserve_upstream_buffer_bp Integer. Safety buffer retained from the
#'   upstream CDS end when \code{preserve_upstream_gene} is active. Default 0.
#' @param auto_preserve_gap_bp Integer. Intergenic gap threshold (bp) used by
#'   the \code{"auto"} preservation heuristic. Default 100.
#' @param upstream_bp Integer. Preferred left-arm length in bp (soft target
#'   used for scoring). Default 500.
#' @param downstream_bp Integer. Preferred right-arm length in bp. Default 500.
#' @param min_arm_bp Integer. Minimum effective homology arm length. Default 200.
#' @param max_arm_bp Integer. Maximum effective homology arm length; also the
#'   length of the extended arm extracted for flexible shared-primer discovery.
#'   Default 1000.
#' @param primer_search_window Integer or NULL. Legacy boundary search window.
#'   NULL (default) scans the entire extended arm.
#' @param primer_min_length Integer. Minimum core length for deletion primers.
#' @param primer_max_length Integer. Maximum core length for deletion primers.
#' @param grna_vector_file Optional vector for gRNA cloning primer generation.
#' @param grna_start,grna_end Optional gRNA stuffer coordinates for
#'   \code{design_cloning_primers()}.
#' @param grna_cloning_method Character. Passed to \code{design_cloning_primers()}.
#' @param grna_enzyme Character. Golden Gate enzyme when applicable.
#' @param grna_custom_overhangs Optional named list of custom overhangs.
#' @param grna_name_prefix Optional name prefix for gRNA cluster labels.
#' @param donor_vector_file Optional donor vector path for donor-only grouped
#'   construct maps.
#' @param donor_start,donor_end Optional donor insertion coordinates.
#' @param donor_construct_output_dir Optional directory for donor construct .gbk
#'   files grouped by identical donor inserts.
#' @param combined_vector_file Optional all-in-one vector path for final grouped
#'   combined construct maps.
#' @param combined_grna_start,combined_grna_end Optional gRNA insertion site on
#'   the combined vector.
#' @param combined_deletion_start,combined_deletion_end Optional deletion insert
#'   site on the combined vector.
#' @param combined_construct_output_dir Optional directory for grouped final
#'   combined construct .gbk files.
#' @param output_file Optional Excel output path.
#' @return A list with components:
#'   \describe{
#'     \item{grna_groups}{Selected shared / subgroup / strain-specific gRNA groups.}
#'     \item{primer_cores}{Selected shared / subgroup / strain-specific deletion cores.}
#'     \item{target_context}{Per-genome extracted target and arm context.}
#'     \item{construct_groups}{Donor insert grouping table.}
#'     \item{final_construct_groups}{Grouping of final combined constructs by shared donor + shared gRNA.}
#'     \item{construct_deletion_primers}{Construct-group-specific full Gibson deletion primer table.}
#'     \item{grna_cloning_primers}{Optional cloning primer table for gRNA groups.}
#'   }
#' @export
design_shared_grna_and_deletion <- function(
    target_table,
    reference_genome_id = NULL,
    reference_locus_tag = NULL,
    auto_resolve_targets = FALSE,
    nuclease = "GeoCas9",
    top_n = 1,
    position_range = c(5, 95),
    tm_range = c(50, 70),
    strand = "both",
    require_unique = TRUE,
    n_mismatches = 0,
    max_n1 = NULL,
    max_n1_per_primer = NULL,
    audit_n_mismatches = 0L,
    verbose = TRUE,
    offtarget_candidate_cap = 100,
    preserve_upstream_gene = FALSE,
    preserve_upstream_buffer_bp = 0L,
    auto_preserve_gap_bp = 100L,
    overlap_policy = c("strict", "flexible"),
    upstream_bp = 500,
    downstream_bp = 500,
    min_arm_bp = 200L,
    max_arm_bp = 1000L,
    primer_search_window = NULL,
    primer_min_length = 20,
    primer_max_length = 35,
    grna_vector_file = NULL,
    grna_start = NULL,
    grna_end = NULL,
    grna_cloning_method = "gibson",
    grna_enzyme = "BbsI",
    grna_custom_overhangs = NULL,
    grna_name_prefix = NULL,
    donor_vector_file = NULL,
    donor_start = NULL,
    donor_end = NULL,
    donor_construct_output_dir = NULL,
    combined_vector_file = NULL,
    combined_grna_start = NULL,
    combined_grna_end = NULL,
    combined_deletion_start = NULL,
    combined_deletion_end = NULL,
    combined_construct_output_dir = NULL,
    output_file = NULL
) {
  overlap_policy <- base::match.arg(overlap_policy)
  if (base::isTRUE(verbose)) {
    base::message(sprintf(
      "[overlap] policy = '%s' (strict = preserve all CDS overlaps, flexible = delete the full target range).",
      overlap_policy))
  }
  target_table <- .shared_prepare_target_table(
    target_table = target_table,
    reference_genome_id = reference_genome_id,
    reference_locus_tag = reference_locus_tag,
    auto_resolve_targets = auto_resolve_targets
  )

  deletion_res <- design_shared_deletion_primer_cores(
    target_table = target_table,
    upstream_bp = upstream_bp,
    downstream_bp = downstream_bp,
    min_arm_bp = min_arm_bp,
    max_arm_bp = max_arm_bp,
    preserve_upstream_gene = preserve_upstream_gene,
    preserve_upstream_buffer_bp = preserve_upstream_buffer_bp,
    auto_preserve_gap_bp = auto_preserve_gap_bp,
    overlap_policy = overlap_policy,
    search_window = primer_search_window,
    min_primer_length = primer_min_length,
    max_primer_length = primer_max_length,
    tm_target = 60,
    require_unique = require_unique,
    audit_n_mismatches = audit_n_mismatches,
    max_n1_per_primer = max_n1_per_primer,
    verbose = verbose,
    vector_file = donor_vector_file,
    start = donor_start,
    end = donor_end,
    construct_output_dir = donor_construct_output_dir,
    output_file = NULL
  )

  target_context <- deletion_res$target_context
  grna_groups <- .shared_design_grna_groups(
    target_context = target_context,
    whole_genome_map = deletion_res$genome_sequences,
    nuclease = nuclease,
    top_n = top_n,
    position_range = position_range,
    tm_range = tm_range,
    strand = strand,
    require_unique = require_unique,
    n_mismatches = n_mismatches,
    max_n1 = max_n1,
    offtarget_candidate_cap = offtarget_candidate_cap
  )

  grna_cloning_primers <- NULL
  if (!base::is.null(grna_vector_file) &&
      !base::is.null(grna_start) &&
      !base::is.null(grna_end) &&
      base::nrow(grna_groups) > 0) {
    grna_df <- dplyr::transmute(
      grna_groups,
      locus_tag = cluster_id,
      gene = cluster_id,
      protospacer = protospacer,
      pam = pam,
      spacer_context = spacer_context,
      spacer_length = base::nchar(protospacer),
      strand = strand,
      composite_score = composite_score,
      percentGC = percentGC,
      Tm = Tm,
      position_percent = position_percent,
      n0 = worst_n0,
      n1 = worst_n1
    )
    grna_df <- generate_gRNA_names(
      gRNA_df = grna_df,
      prefix = grna_name_prefix,
      nuclease = nuclease,
      rank_by = "composite_score"
    )
    grna_cloning_primers <- design_cloning_primers(
      gRNA_df = grna_df,
      vector_file = grna_vector_file,
      start = grna_start,
      end = grna_end,
      cloning_method = grna_cloning_method,
      enzyme = grna_enzyme,
      custom_overhangs = grna_custom_overhangs,
      output_file = NULL
    )
  }

  final_construct_groups <- .shared_build_final_construct_groups(
    construct_groups = deletion_res$construct_groups,
    grna_groups = grna_groups
  )
  construct_deletion_primers <- NULL
  if (!base::is.null(combined_vector_file) &&
      !base::is.null(combined_deletion_start) &&
      !base::is.null(combined_deletion_end) &&
      base::nrow(deletion_res$construct_groups) > 0) {
    construct_deletion_primers <- .shared_design_construct_deletion_primers(
      construct_groups = deletion_res$construct_groups,
      target_context = target_context,
      primer_cores = deletion_res$primer_cores,
      vector_file = combined_vector_file,
      start = combined_deletion_start,
      end = combined_deletion_end
    )
  }

  if (!base::is.null(combined_construct_output_dir)) {
    if (base::is.null(combined_vector_file) ||
        base::is.null(combined_grna_start) ||
        base::is.null(combined_grna_end) ||
        base::is.null(combined_deletion_start) ||
        base::is.null(combined_deletion_end)) {
      base::stop("Combined construct writing requires combined_vector_file, combined_grna_start/end, and combined_deletion_start/end.")
    }
    if (!base::dir.exists(combined_construct_output_dir)) {
      base::dir.create(combined_construct_output_dir, recursive = TRUE)
    }
    .shared_write_final_construct_groups(
      final_construct_groups = final_construct_groups,
      target_context = target_context,
      grna_groups = grna_groups,
      construct_deletion_primers = construct_deletion_primers,
      vector_file = combined_vector_file,
      grna_start = combined_grna_start,
      grna_end = combined_grna_end,
      deletion_start = combined_deletion_start,
      deletion_end = combined_deletion_end,
      upstream_bp = upstream_bp,
      downstream_bp = downstream_bp,
      output_dir = combined_construct_output_dir
    )
  }

  if (!base::is.null(output_file)) {
    sheets <- base::list(
      grna_groups = grna_groups,
      primer_cores = deletion_res$primer_cores,
      target_context = target_context,
      construct_groups = deletion_res$construct_groups,
      final_construct_groups = final_construct_groups
    )
    if (!base::is.null(grna_cloning_primers)) {
      sheets$grna_cloning_primers <- grna_cloning_primers
    }
    if (base::requireNamespace("openxlsx", quietly = TRUE)) {
      .shared_write_design_workbook(
        output_file = output_file,
        grna_groups = grna_groups,
        primer_cores = deletion_res$primer_cores,
        target_context = target_context,
        construct_groups = deletion_res$construct_groups,
        final_construct_groups = final_construct_groups,
        grna_cloning_primers = grna_cloning_primers,
        construct_deletion_primers = construct_deletion_primers
      )
    } else {
      writexl::write_xlsx(sheets, path = output_file)
    }
  }

  base::list(
    grna_groups = grna_groups,
    primer_cores = deletion_res$primer_cores,
    target_context = target_context,
    construct_groups = deletion_res$construct_groups,
    final_construct_groups = final_construct_groups,
    grna_cloning_primers = grna_cloning_primers,
    construct_deletion_primers = construct_deletion_primers
  )
}

.shared_write_design_workbook <- function(
    output_file,
    grna_groups,
    primer_cores,
    target_context,
    construct_groups,
    final_construct_groups,
    grna_cloning_primers = NULL,
    construct_deletion_primers = NULL
) {
  wb <- openxlsx::createWorkbook()

  primer_order <- .shared_build_primer_order_sheet(primer_cores)
  actual_primer_order <- .shared_build_actual_primer_order_sheet(
    construct_deletion_primers, primer_cores = primer_cores
  )
  target_context_light <- .shared_target_context_light(target_context)
  construct_requirements <- .shared_build_construct_requirements(
    primer_cores = primer_cores,
    construct_groups = construct_groups,
    final_construct_groups = final_construct_groups,
    grna_groups = grna_groups
  )
  summary_df <- .shared_build_shared_design_summary(
    grna_groups = grna_groups,
    primer_cores = primer_cores,
    primer_order = primer_order,
    actual_primer_order = actual_primer_order,
    construct_groups = construct_groups,
    final_construct_groups = final_construct_groups
  )

  sheets <- base::list(
    summary = summary_df,
    primer_order = primer_order,
    full_primer_order = actual_primer_order,
    grna_groups = grna_groups,
    primer_cores = primer_cores,
    target_context = target_context_light,
    construct_groups = construct_groups,
    final_construct_groups = final_construct_groups,
    construct_requirements = construct_requirements
  )
  if (!base::is.null(grna_cloning_primers)) {
    sheets$grna_cloning_primers <- grna_cloning_primers
  }
  if (!base::is.null(construct_deletion_primers) && base::nrow(construct_deletion_primers) > 0) {
    sheets$construct_deletion_primers <- construct_deletion_primers
  }

  for (sn in base::names(sheets)) {
    openxlsx::addWorksheet(wb, sn)
    openxlsx::writeData(wb, sn, sheets[[sn]])
    .shared_style_basic_sheet(wb, sn, sheets[[sn]])
  }

  .shared_style_cluster_sheet(wb, "primer_order", primer_order, "color_key")
  .shared_style_cluster_sheet(wb, "full_primer_order", actual_primer_order, "genomes_shared")
  .shared_style_cluster_sheet(wb, "grna_groups", grna_groups, "used_by")
  .shared_style_cluster_sheet(wb, "construct_groups", construct_groups, "used_by")
  .shared_style_cluster_sheet(wb, "final_construct_groups", final_construct_groups, "members")
  .shared_style_cluster_sheet(wb, "construct_requirements", construct_requirements, "construct_label")
  if (!base::is.null(grna_cloning_primers)) {
    .shared_style_cluster_sheet(wb, "grna_cloning_primers", grna_cloning_primers, "locus_tag")
  }

  openxlsx::saveWorkbook(wb, file = output_file, overwrite = TRUE)
}

.shared_target_context_light <- function(target_context) {
  drop_cols <- c(
    "left_arm_seq", "right_arm_seq", "insert_seq", "target_nt_seq", "genome_seq",
    "effective_left_arm_seq", "effective_right_arm_seq", "effective_insert_seq"
  )
  keep <- base::setdiff(base::colnames(target_context), drop_cols)
  target_context[, keep, drop = FALSE]
}

.shared_build_shared_design_summary <- function(
    grna_groups,
    primer_cores,
    primer_order,
    actual_primer_order,
    construct_groups,
    final_construct_groups
) {
  n_final_unique_oligos <- if (!base::is.null(actual_primer_order) &&
                               base::nrow(actual_primer_order) > 0) {
    base::nrow(actual_primer_order)
  } else {
    NA_integer_
  }
  n_constructs <- base::nrow(construct_groups)
  # Naive "one set of 4 Gibson primers per construct" baseline for cost
  # comparison — the sharing is worthwhile only if final_unique_oligos is
  # substantially smaller than this.
  baseline_oligos <- if (n_constructs > 0) 4L * n_constructs else NA_integer_
  oligos_saved <- if (!base::is.na(n_final_unique_oligos) && !base::is.na(baseline_oligos)) {
    base::max(0L, baseline_oligos - n_final_unique_oligos)
  } else {
    NA_integer_
  }
  worst_n1 <- if (!base::is.null(primer_cores) && "worst_n1" %in% base::colnames(primer_cores) && base::nrow(primer_cores) > 0) {
    base::suppressWarnings(base::max(primer_cores$worst_n1, na.rm = TRUE))
  } else {
    NA_integer_
  }
  base::data.frame(
    metric = c(
      "gRNA_groups",
      "primer_core_rows",
      "distinct_core_sequences",
      "shared_core_sets",
      "construct_groups",
      "final_construct_groups",
      "baseline_oligos_4x_constructs",
      "final_unique_oligos_to_order",
      "oligos_saved_vs_baseline",
      "worst_per_primer_1mismatch_hits"
    ),
    value = c(
      base::nrow(grna_groups),
      base::nrow(primer_cores),
      base::length(base::unique(primer_order$primer_core_5to3)),
      base::length(base::unique(primer_order$used_by)),
      n_constructs,
      base::nrow(final_construct_groups),
      baseline_oligos,
      n_final_unique_oligos,
      oligos_saved,
      worst_n1
    ),
    stringsAsFactors = FALSE
  )
}

.shared_build_primer_order_sheet <- function(primer_cores) {
  if (base::nrow(primer_cores) == 0) return(primer_cores)

  df <- primer_cores
  df$used_by <- vapply(strsplit(df$used_by, ", ", fixed = TRUE), function(x) base::paste(sort(x), collapse = ", "), character(1))
  df$n_targets <- vapply(strsplit(df$used_by, ", ", fixed = TRUE), base::length, integer(1))
  df$order_name <- base::paste0(df$role, "_", df$cluster_id)
  df$color_key <- df$used_by
  keep_cols <- c(
    "order_name", "role", "scope", "used_by", "n_targets",
    "primer_core_5to3", "primer_length", "tm_core", "gc_percent",
    "boundary_distance_bp", "exact_unique_in_targets",
    "worst_n0", "worst_n1", "per_genome_n0", "per_genome_n1",
    "color_key"
  )
  keep_cols <- base::intersect(keep_cols, base::colnames(df))
  df <- df[, keep_cols, drop = FALSE]
  df <- df[base::order(-df$n_targets, df$role, df$order_name), , drop = FALSE]
  rownames(df) <- NULL
  df
}

.shared_build_actual_primer_order_sheet <- function(construct_deletion_primers,
                                                    primer_cores = NULL) {
  if (base::is.null(construct_deletion_primers) || base::nrow(construct_deletion_primers) == 0) {
    return(base::data.frame())
  }
  # Role-keyed column maps pointing at per-construct primer records.
  role_map <- base::list(
    UF = c("upstream_forward_name", "upstream_forward_primer", "upstream_forward_tm_full"),
    UR = c("upstream_reverse_name", "upstream_reverse_primer", "upstream_reverse_tm_full"),
    DF = c("downstream_forward_name", "downstream_forward_primer", "downstream_forward_tm_full"),
    DR = c("downstream_reverse_name", "downstream_reverse_primer", "downstream_reverse_tm_full")
  )

  # Per-cluster audit columns from primer_cores (worst-case n0 / n1) that we
  # want to carry onto the final dedup sheet so the user can see, for every
  # ordered oligo, whether it stays unique in every genome that uses it.
  audit_lookup <- base::list()
  if (!base::is.null(primer_cores) && base::nrow(primer_cores) > 0 &&
      base::all(c("cluster_id", "role") %in% base::colnames(primer_cores))) {
    for (ri in base::seq_len(base::nrow(primer_cores))) {
      key <- base::paste(primer_cores$role[ri], primer_cores$cluster_id[ri], sep = "::")
      audit_lookup[[key]] <- base::list(
        worst_n0 = if ("worst_n0" %in% base::colnames(primer_cores)) primer_cores$worst_n0[ri] else NA_integer_,
        worst_n1 = if ("worst_n1" %in% base::colnames(primer_cores)) primer_cores$worst_n1[ri] else NA_integer_,
        per_genome_n0 = if ("per_genome_n0" %in% base::colnames(primer_cores)) primer_cores$per_genome_n0[ri] else NA_character_,
        per_genome_n1 = if ("per_genome_n1" %in% base::colnames(primer_cores)) primer_cores$per_genome_n1[ri] else NA_character_
      )
    }
  }

  rows <- base::list()
  idx <- 1L
  for (i in base::seq_len(base::nrow(construct_deletion_primers))) {
    row <- construct_deletion_primers[i, , drop = FALSE]
    for (role in base::names(role_map)) {
      nm_col <- role_map[[role]][1]
      seq_col <- role_map[[role]][2]
      tm_col <- role_map[[role]][3]
      cluster_name <- base::as.character(row[[nm_col]][1])
      audit_key <- base::paste(role, cluster_name, sep = "::")
      audit <- if (!base::is.null(audit_lookup[[audit_key]])) audit_lookup[[audit_key]] else base::list(
        worst_n0 = NA_integer_, worst_n1 = NA_integer_,
        per_genome_n0 = NA_character_, per_genome_n1 = NA_character_
      )
      rows[[idx]] <- base::data.frame(
        cluster_id = cluster_name,
        construct_label = row$construct_label[1],
        genomes_in_construct = row$used_by[1],
        role = role,
        sequence_to_order_5to3 = row[[seq_col]][1],
        full_length = base::nchar(row[[seq_col]][1]),
        tm_full = row[[tm_col]][1],
        worst_n0 = audit$worst_n0,
        worst_n1 = audit$worst_n1,
        per_genome_n0 = audit$per_genome_n0,
        per_genome_n1 = audit$per_genome_n1,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  df <- dplyr::bind_rows(rows)
  if (base::nrow(df) == 0) return(df)

  # Round-3 sharing recomputation: the FINAL full primer sequence (core +
  # vector overhang) is what actually gets ordered. Any two constructs that
  # yield identical sequences after overhang assembly are now collapsed into a
  # single row so the user only pays for one oligo.
  grouped <- base::split(df, df$sequence_to_order_5to3)
  out <- lapply(grouped, function(g) {
    genomes_all <- base::sort(base::unique(base::unlist(
      base::strsplit(g$genomes_in_construct, ", ", fixed = TRUE)
    )))
    cluster_ids <- base::sort(base::unique(g$cluster_id))
    construct_labels <- base::unique(g$construct_label)
    # Canonical order name: if exactly one cluster yields this full sequence
    # use its cluster_id; otherwise compose a dedup tag that keeps the mapping
    # traceable back to the source clusters.
    order_name <- if (base::length(cluster_ids) == 1L) {
      cluster_ids[1]
    } else {
      base::paste(cluster_ids, collapse = "+")
    }
    base::data.frame(
      order_name = order_name,
      role = g$role[1],
      n_genomes_shared = base::length(genomes_all),
      n_constructs_shared = base::length(construct_labels),
      genomes_shared = base::paste(genomes_all, collapse = ", "),
      construct_labels = base::paste(construct_labels, collapse = " | "),
      cluster_ids = base::paste(cluster_ids, collapse = ", "),
      sequence_to_order_5to3 = g$sequence_to_order_5to3[1],
      full_length = g$full_length[1],
      tm_full = g$tm_full[1],
      worst_n0 = base::suppressWarnings(base::max(g$worst_n0, na.rm = TRUE)),
      worst_n1 = base::suppressWarnings(base::max(g$worst_n1, na.rm = TRUE)),
      per_genome_n0 = g$per_genome_n0[1],
      per_genome_n1 = g$per_genome_n1[1],
      stringsAsFactors = FALSE
    )
  })
  out <- dplyr::bind_rows(out)
  out[base::order(out$role, -out$n_genomes_shared, out$order_name), , drop = FALSE]
}

.shared_build_construct_requirements <- function(primer_cores, construct_groups, final_construct_groups, grna_groups) {
  if (base::nrow(final_construct_groups) == 0) return(base::data.frame())
  rows <- base::list()
  idx <- 1L
  for (i in base::seq_len(base::nrow(final_construct_groups))) {
    fg <- final_construct_groups[i, , drop = FALSE]
    members <- sort(base::strsplit(fg$used_by[1], ", ", fixed = TRUE)[[1]])
    usable_primers <- primer_cores[vapply(strsplit(primer_cores$used_by, ", ", fixed = TRUE), function(x) all(members %in% x), logical(1)), , drop = FALSE]
    rows[[idx]] <- base::data.frame(
      construct_label = fg$construct_label[1],
      members = base::paste(members, collapse = ", "),
      donor_group_id = fg$donor_group_id[1],
      grna_cluster_id = fg$grna_group_id[1],
      grna_protospacer = grna_groups$protospacer[grna_groups$cluster_id == fg$grna_group_id[1]][1],
      primer_sets = base::paste(base::paste0(usable_primers$role, ":", usable_primers$cluster_id), collapse = " | "),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
  dplyr::bind_rows(rows)
}

.shared_design_construct_deletion_primers <- function(
    construct_groups,
    target_context,
    primer_cores,
    vector_file,
    start,
    end
) {
  if (base::nrow(construct_groups) == 0) return(base::data.frame())
  core_map <- stats::setNames(base::split(primer_cores, primer_cores$role), unique(primer_cores$role))
  rows <- base::list()
  for (i in base::seq_len(base::nrow(construct_groups))) {
    rep_id <- construct_groups$representative_genome[i]
    label <- construct_groups$construct_label[i]
    ctx <- target_context[target_context$genome_id == rep_id, , drop = FALSE]
    if (base::nrow(ctx) != 1L) next
    ur_row <- core_map$UR[0, , drop = FALSE]
    df_row <- core_map$DF[0, , drop = FALSE]
    dr_row <- core_map$DR[0, , drop = FALSE]
    uf_row <- core_map$UF[0, , drop = FALSE]
    for (role in c("UF", "UR", "DF", "DR")) {
      cand <- primer_cores[primer_cores$role == role, , drop = FALSE]
      keep <- vapply(strsplit(cand$used_by, ", ", fixed = TRUE), function(x) rep_id %in% x, logical(1))
      picked <- cand[keep, , drop = FALSE]
      if (base::nrow(picked) != 1L) {
        base::stop("Expected exactly one ", role, " core for construct representative ", rep_id)
      }
      if (role == "UF") uf_row <- picked
      if (role == "UR") ur_row <- picked
      if (role == "DF") df_row <- picked
      if (role == "DR") dr_row <- picked
    }
    primer_res <- design_deletion_primers(
      vector_file = vector_file,
      insert_seq = ctx$effective_insert_seq[1],
      start = start,
      end = end,
      upstream_bp = ctx$effective_upstream_bp[1],
      downstream_bp = ctx$effective_downstream_bp[1],
      locus_tag = label,
      upstream_forward_target = base::as.character(uf_row$shared_core[1]),
      upstream_reverse_target = base::as.character(ur_row$shared_core[1]),
      downstream_forward_target = base::as.character(df_row$shared_core[1]),
      downstream_reverse_target = base::as.character(dr_row$shared_core[1]),
      kill_snapgene = FALSE
    )
    # Canonical primer names are keyed on cluster_id so identical sequences
    # from different construct groups collapse into a single oligo order.
    # Construct-specific annotation names are kept separately for .gbk display.
    uf_order_name <- base::as.character(uf_row$cluster_id[1])
    ur_order_name <- base::as.character(ur_row$cluster_id[1])
    df_order_name <- base::as.character(df_row$cluster_id[1])
    dr_order_name <- base::as.character(dr_row$cluster_id[1])
    uf_anno_name <- base::paste0(label, "::", uf_order_name)
    ur_anno_name <- base::paste0(label, "::", ur_order_name)
    df_anno_name <- base::paste0(label, "::", df_order_name)
    dr_anno_name <- base::paste0(label, "::", dr_order_name)

    rows[[i]] <- base::data.frame(
      construct_label = label,
      representative_genome = rep_id,
      used_by = construct_groups$used_by[i],
      effective_upstream_bp = ctx$effective_upstream_bp[1],
      effective_downstream_bp = ctx$effective_downstream_bp[1],
      effective_deletion_start = ctx$effective_deletion_start[1],
      effective_deletion_end = ctx$effective_deletion_end[1],
      upstream_forward_name = uf_order_name,
      upstream_forward_annotation_name = uf_anno_name,
      upstream_forward_primer = primer_res$upstream_forward_primer,
      upstream_forward_tm_target = primer_res$upstream_forward_tm_target,
      upstream_forward_tm_full = primer_res$upstream_forward_tm_full,
      upstream_forward_start = primer_res$upstream_forward_start,
      upstream_forward_end = primer_res$upstream_forward_end,
      upstream_reverse_name = ur_order_name,
      upstream_reverse_annotation_name = ur_anno_name,
      upstream_reverse_primer = primer_res$upstream_reverse_primer,
      upstream_reverse_tm_target = primer_res$upstream_reverse_tm_target,
      upstream_reverse_tm_full = primer_res$upstream_reverse_tm_full,
      upstream_reverse_start = primer_res$upstream_reverse_start,
      upstream_reverse_end = primer_res$upstream_reverse_end,
      downstream_forward_name = df_order_name,
      downstream_forward_annotation_name = df_anno_name,
      downstream_forward_primer = primer_res$downstream_forward_primer,
      downstream_forward_tm_target = primer_res$downstream_forward_tm_target,
      downstream_forward_tm_full = primer_res$downstream_forward_tm_full,
      downstream_forward_start = primer_res$downstream_forward_start,
      downstream_forward_end = primer_res$downstream_forward_end,
      downstream_reverse_name = dr_order_name,
      downstream_reverse_annotation_name = dr_anno_name,
      downstream_reverse_primer = primer_res$downstream_reverse_primer,
      downstream_reverse_tm_target = primer_res$downstream_reverse_tm_target,
      downstream_reverse_tm_full = primer_res$downstream_reverse_tm_full,
      downstream_reverse_start = primer_res$downstream_reverse_start,
      downstream_reverse_end = primer_res$downstream_reverse_end,
      overlap_tm = primer_res$overlap_tm,
      overlap_length = primer_res$overlap_length,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

.shared_style_basic_sheet <- function(wb, sheet_name, df) {
  if (base::nrow(df) == 0) return(invisible(NULL))
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#D9E2F3",
    border = "Bottom"
  )
  openxlsx::addStyle(
    wb, sheet = sheet_name, style = header_style,
    rows = 1, cols = base::seq_len(base::ncol(df)), gridExpand = TRUE, stack = TRUE
  )
  openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  base::tryCatch(
    openxlsx::setColWidths(wb, sheet = sheet_name, cols = base::seq_len(base::ncol(df)), widths = "auto"),
    error = function(e) NULL
  )
}

.shared_style_cluster_sheet <- function(wb, sheet_name, df, cluster_col) {
  if (base::nrow(df) == 0) return(invisible(NULL))
  if (!cluster_col %in% base::colnames(df)) return(invisible(NULL))

  palette <- c(
    "#DCEBFA", "#DFF3E4", "#FFF0D9", "#FBE0E0", "#E9E1FB",
    "#FCECCB", "#E4F7F5", "#F6E6F5", "#E8EEF9", "#FBEFD8"
  )
  groups <- base::unique(df[[cluster_col]])
  idx <- ((base::seq_along(groups) - 1L) %% base::length(palette)) + 1L
  group_map <- stats::setNames(palette[idx], groups)
  for (g in groups) {
    rows <- base::which(df[[cluster_col]] == g) + 1L
    style <- openxlsx::createStyle(fgFill = group_map[[g]])
    openxlsx::addStyle(
      wb, sheet = sheet_name, style = style,
      rows = rows, cols = base::seq_len(base::ncol(df)), gridExpand = TRUE, stack = TRUE
    )
  }
}

.shared_design_grna_groups <- function(
    target_context,
    whole_genome_map,
    nuclease,
    top_n,
    position_range,
    tm_range,
    strand,
    require_unique,
    n_mismatches,
    max_n1,
    offtarget_candidate_cap
) {
  genome_ids <- target_context$genome_id
  hit_count_cache <- new.env(parent = emptyenv())

  candidate_map <- base::lapply(base::seq_len(base::nrow(target_context)), function(i) {
    row <- target_context[i, , drop = FALSE]
    .shared_find_guides_in_target(
      target_seq = row$target_nt_seq[1],
      genome_id = row$genome_id[1],
      locus_tag = row$locus_tag[1],
      nuclease = nuclease,
      position_range = position_range,
      tm_range = tm_range,
      strand = strand
    )
  })
  names(candidate_map) <- genome_ids
  candidate_pool <- .shared_build_grna_candidate_pool(
    candidate_map = candidate_map,
    whole_genome_map = whole_genome_map,
    tm_target = base::mean(tm_range),
    require_unique = require_unique,
    n_mismatches = n_mismatches,
    max_n1 = max_n1,
    offtarget_candidate_cap = offtarget_candidate_cap,
    hit_count_cache = hit_count_cache
  )
  .shared_select_grna_groups_from_pool(candidate_pool, genome_ids, top_n = top_n)
}

.shared_build_grna_candidate_pool <- function(
    candidate_map,
    whole_genome_map,
    tm_target,
    require_unique,
    n_mismatches,
    max_n1,
    offtarget_candidate_cap,
    hit_count_cache
) {
  all_candidates <- dplyr::bind_rows(candidate_map)
  if (base::nrow(all_candidates) == 0) {
    return(base::list())
  }

  split_keys <- base::split(all_candidates, all_candidates$candidate_key)
  pre_meta <- lapply(split_keys, function(df) {
    base::list(
      key = df$candidate_key[1],
      genomes = unique(df$genome_id),
      protospacer = df$protospacer[1],
      pam = df$pam[1],
      spacer_context = df$spacer_context[1],
      strand = df$strand[1],
      local_min_score = base::min(df$composite_score),
      local_mean_score = base::mean(df$composite_score),
      local_max_pos_dev = base::max(base::abs(df$position_percent - 50)),
      tm = df$Tm[1],
      percentGC = df$percentGC[1],
      per_genome_score = tapply(df$composite_score, df$genome_id, max),
      per_genome_pos = tapply(df$position_percent, df$genome_id, function(x) x[which.min(abs(x - 50))[1]]),
      duplicate_hits_in_target = any(duplicated(df$genome_id))
    )
  })

  pre_order <- base::order(
    vapply(pre_meta, function(x) -base::length(x$genomes), numeric(1)),
    vapply(pre_meta, function(x) -x$local_min_score, numeric(1)),
    vapply(pre_meta, function(x) -x$local_mean_score, numeric(1)),
    vapply(pre_meta, function(x) x$local_max_pos_dev, numeric(1)),
    vapply(pre_meta, function(x) base::abs(x$tm - tm_target), numeric(1))
  )
  if (!base::is.null(offtarget_candidate_cap) && base::length(pre_order) > offtarget_candidate_cap) {
    pre_meta <- pre_meta[pre_order[base::seq_len(offtarget_candidate_cap)]]
  } else {
    pre_meta <- pre_meta[pre_order]
  }

  out <- base::list()
  idx <- 1L
  for (cand in pre_meta) {
    offtarget_list <- lapply(cand$genomes, function(id) {
      .shared_spacer_offtarget_counts(
        genome_id = id,
        protospacer = cand$protospacer,
        whole_genome_map = whole_genome_map,
        hit_count_cache = hit_count_cache,
        n_mismatches = n_mismatches
      )
    })
    names(offtarget_list) <- cand$genomes
    n0_vals <- vapply(offtarget_list, function(x) x$n0, numeric(1))
    n1_vals <- vapply(offtarget_list, function(x) x$n1, numeric(1))
    if (require_unique && !base::all(n0_vals == 1L)) next
    if (!base::is.null(max_n1) && base::any(n1_vals > max_n1)) next

    out[[idx]] <- base::list(
      candidate_key = cand$key,
      genomes = cand$genomes,
      protospacer = cand$protospacer,
      pam = cand$pam,
      spacer_context = cand$spacer_context,
      strand = cand$strand,
      percentGC = cand$percentGC,
      Tm = cand$tm,
      per_genome_score = cand$per_genome_score,
      per_genome_pos = cand$per_genome_pos,
      per_genome_n0 = n0_vals,
      per_genome_n1 = n1_vals
    )
    idx <- idx + 1L
  }
  out
}

.shared_select_grna_groups_from_pool <- function(candidate_pool, genome_ids, top_n = 1L) {
  if (base::length(candidate_pool) == 0) {
    base::stop("No shared gRNA candidates survived uniqueness/off-target filtering.")
  }

  remaining <- genome_ids
  out <- base::list()
  cluster_index <- 1L
  while (base::length(remaining) > 0) {
    best <- NULL
    best_score <- NULL
    for (cand in candidate_pool) {
      current_members <- base::intersect(cand$genomes, remaining)
      if (base::length(current_members) == 0) next
      score_vec <- c(
        base::length(current_members),
        -base::max(cand$per_genome_n0[current_members]),
        -base::max(cand$per_genome_n1[current_members]),
        -base::sum(cand$per_genome_n1[current_members]),
        base::min(cand$per_genome_score[current_members]),
        base::mean(cand$per_genome_score[current_members]),
        -base::max(base::abs(cand$per_genome_pos[current_members] - 50)),
        -base::abs(cand$Tm - 60)
      )
      if (base::is.null(best_score) || .shared_numeric_score_better(score_vec, best_score)) {
        best <- base::list(
          candidate = cand,
          genomes = current_members,
          worst_n0 = base::max(cand$per_genome_n0[current_members]),
          worst_n1 = base::max(cand$per_genome_n1[current_members]),
          sum_n1 = base::sum(cand$per_genome_n1[current_members]),
          composite_score = base::mean(cand$per_genome_score[current_members]),
          min_score = base::min(cand$per_genome_score[current_members]),
          position_percent = base::mean(cand$per_genome_pos[current_members])
        )
        best_score <- score_vec
      }
    }
    if (base::is.null(best)) {
      base::stop("No shared gRNA candidate found for remaining genomes: ", base::paste(remaining, collapse = ", "))
    }

    cluster_label <- if (base::length(best$genomes) == base::length(genome_ids)) {
      "all"
    } else if (base::length(best$genomes) == 1L) {
      base::paste0(best$genomes, " only")
    } else {
      base::paste(best$genomes, collapse = ", ")
    }
    cand <- best$candidate
    out[[base::length(out) + 1L]] <- base::data.frame(
      cluster_id = base::paste0("gRNA_cluster_", cluster_index),
      scope = cluster_label,
      used_by = base::paste(best$genomes, collapse = ", "),
      protospacer = cand$protospacer,
      pam = cand$pam,
      spacer_context = cand$spacer_context,
      strand = cand$strand,
      percentGC = base::round(cand$percentGC, 2),
      Tm = base::round(cand$Tm, 2),
      position_percent = base::round(best$position_percent, 2),
      exact_unique_in_targets = best$worst_n0 == 1L,
      worst_n0 = best$worst_n0,
      worst_n1 = best$worst_n1,
      sum_n1 = best$sum_n1,
      composite_score = base::round(best$composite_score, 2),
      n_genomes = base::length(best$genomes),
      stringsAsFactors = FALSE
    )
    remaining <- base::setdiff(remaining, best$genomes)
    cluster_index <- cluster_index + 1L
    if (!base::is.null(top_n) && base::length(out) >= top_n) break
  }

  dplyr::bind_rows(out)
}

.shared_numeric_score_better <- function(a, b) {
  for (i in base::seq_along(a)) {
    if (a[i] > b[i]) return(TRUE)
    if (a[i] < b[i]) return(FALSE)
  }
  FALSE
}

.shared_best_grna_for_subset <- function(
    genome_ids,
    candidate_map,
    whole_genome_map,
    tm_target,
    require_unique,
    n_mismatches,
    max_n1,
    offtarget_candidate_cap,
    hit_count_cache
) {
  shared_keys <- NULL
  for (id in genome_ids) {
    keys <- candidate_map[[id]]$candidate_key
    shared_keys <- if (base::is.null(shared_keys)) keys else base::intersect(shared_keys, keys)
  }
  if (base::length(shared_keys) == 0) return(NULL)

  # Fast pre-ranking using only target-local scores before genome-wide
  # off-target checks. This reduces the number of expensive whole-genome
  # searches while preserving the highest-quality candidates.
  pre_rank <- lapply(shared_keys, function(key) {
    rows <- lapply(genome_ids, function(id) {
      candidate_map[[id]][candidate_map[[id]]$candidate_key == key, , drop = FALSE]
    })
    meta <- rows[[1]][1, , drop = FALSE]
    base::list(
      key = key,
      min_score = base::min(vapply(rows, function(x) x$composite_score[1], numeric(1))),
      mean_score = base::mean(vapply(rows, function(x) x$composite_score[1], numeric(1))),
      max_pos_dev = base::max(vapply(rows, function(x) base::abs(x$position_percent[1] - 50), numeric(1))),
      tm = meta$Tm[1]
    )
  })
  pre_order <- base::order(
    vapply(pre_rank, function(x) -x$min_score, numeric(1)),
    vapply(pre_rank, function(x) -x$mean_score, numeric(1)),
    vapply(pre_rank, function(x) x$max_pos_dev, numeric(1)),
    vapply(pre_rank, function(x) base::abs(x$tm - tm_target), numeric(1))
  )
  if (!base::is.null(offtarget_candidate_cap) && base::length(pre_order) > offtarget_candidate_cap) {
    shared_keys <- vapply(pre_rank[pre_order[base::seq_len(offtarget_candidate_cap)]], function(x) x$key, character(1))
  }

  candidates <- base::list()
  for (key in shared_keys) {
    rows <- lapply(genome_ids, function(id) {
      candidate_map[[id]][candidate_map[[id]]$candidate_key == key, , drop = FALSE]
    })
    meta <- rows[[1]][1, , drop = FALSE]
    offtarget_list <- lapply(genome_ids, function(id) {
      .shared_spacer_offtarget_counts(
        genome_id = id,
        protospacer = meta$protospacer[1],
        whole_genome_map = whole_genome_map,
        hit_count_cache = hit_count_cache,
        n_mismatches = n_mismatches
      )
    })
    n0_vals <- vapply(offtarget_list, function(x) x$n0, numeric(1))
    n1_vals <- vapply(offtarget_list, function(x) x$n1, numeric(1))
    unique_flags <- if (require_unique) n0_vals == 1L else rep(TRUE, length(genome_ids))
    if (require_unique && !base::all(unique_flags)) next
    if (!base::is.null(max_n1) && base::any(n1_vals > max_n1)) next

    min_score <- base::min(vapply(rows, function(x) x$composite_score[1], numeric(1)))
    mean_score <- base::mean(vapply(rows, function(x) x$composite_score[1], numeric(1)))
    max_pos_dev <- base::max(vapply(rows, function(x) base::abs(x$position_percent[1] - 50), numeric(1)))
    candidates[[base::length(candidates) + 1L]] <- base::list(
      genome_ids = genome_ids,
      protospacer = meta$protospacer[1],
      pam = meta$pam[1],
      spacer_context = meta$spacer_context[1],
      strand = meta$strand[1],
      percentGC = meta$percentGC[1],
      Tm = meta$Tm[1],
      position_percent = base::mean(vapply(rows, function(x) x$position_percent[1], numeric(1))),
      all_unique = base::all(unique_flags),
      worst_n0 = base::max(n0_vals),
      worst_n1 = base::max(n1_vals),
      sum_n1 = base::sum(n1_vals),
      composite_score = mean_score,
      min_score = min_score,
      max_pos_dev = max_pos_dev
    )
  }
  if (base::length(candidates) == 0) return(NULL)
  ranked <- base::order(
    vapply(candidates, function(x) x$worst_n0, numeric(1)),
    vapply(candidates, function(x) x$worst_n1, numeric(1)),
    vapply(candidates, function(x) x$sum_n1, numeric(1)),
    vapply(candidates, function(x) -x$min_score, numeric(1)),
    vapply(candidates, function(x) -x$composite_score, numeric(1)),
    vapply(candidates, function(x) x$max_pos_dev, numeric(1)),
    vapply(candidates, function(x) base::abs(x$Tm - tm_target), numeric(1))
  )
  candidates[[ranked[1]]]
}

.shared_prepare_target_table <- function(
    target_table,
    reference_genome_id = NULL,
    reference_locus_tag = NULL,
    auto_resolve_targets = FALSE
) {
  if (!base::is.data.frame(target_table) || base::nrow(target_table) == 0) {
    base::stop("target_table must be a non-empty data frame.")
  }
  required_cols <- c("genome_id", "genbank_file")
  missing_cols <- base::setdiff(required_cols, base::colnames(target_table))
  if (base::length(missing_cols) > 0) {
    base::stop("target_table must contain columns: ", base::paste(required_cols, collapse = ", "))
  }

  if (!"locus_tag" %in% base::colnames(target_table)) {
    target_table$locus_tag <- NA_character_
  }

  if (!auto_resolve_targets && base::any(base::is.na(target_table$locus_tag) | target_table$locus_tag == "")) {
    base::stop("target_table contains missing locus_tag values. Set auto_resolve_targets = TRUE and provide reference_genome_id + reference_locus_tag.")
  }

  if (!auto_resolve_targets) return(target_table)
  if (base::is.null(reference_genome_id) || base::is.null(reference_locus_tag)) {
    base::stop("auto_resolve_targets = TRUE requires reference_genome_id and reference_locus_tag.")
  }

  .shared_resolve_target_locus_tags(
    target_table = target_table,
    reference_genome_id = reference_genome_id,
    reference_locus_tag = reference_locus_tag
  )
}

.shared_resolve_target_locus_tags <- function(target_table, reference_genome_id, reference_locus_tag) {
  ref_row <- target_table[target_table$genome_id == reference_genome_id, , drop = FALSE]
  if (base::nrow(ref_row) != 1) {
    base::stop("reference_genome_id must match exactly one row in target_table.")
  }

    parsed_ref <- .shared_read_genome_cached(base::as.character(ref_row$genbank_file[1]))
  ref_tbl <- parsed_ref$genbank_table
  ref_hit <- ref_tbl[ref_tbl$locus_tag == reference_locus_tag, , drop = FALSE]
  if (base::nrow(ref_hit) != 1) {
    base::stop("reference_locus_tag not found in reference genome: ", reference_locus_tag)
  }
  ref_gene <- if ("gene" %in% base::colnames(ref_hit)) base::as.character(ref_hit$gene[1]) else NA_character_
  ref_product <- if ("product" %in% base::colnames(ref_hit)) base::as.character(ref_hit$product[1]) else NA_character_
  ref_aa <- .shared_translate_nt(ref_hit$nt_seq[1])
  ref_len <- base::nchar(ref_aa)

  resolved <- target_table
  resolved$resolved_by <- if ("resolved_by" %in% base::colnames(resolved)) resolved$resolved_by else NA_character_
  resolved$mapping_score <- if ("mapping_score" %in% base::colnames(resolved)) resolved$mapping_score else NA_real_
  resolved$mapping_note <- if ("mapping_note" %in% base::colnames(resolved)) resolved$mapping_note else NA_character_

  for (i in base::seq_len(base::nrow(resolved))) {
    if (!base::is.na(resolved$locus_tag[i]) && resolved$locus_tag[i] != "") next

    parsed <- .shared_read_genome_cached(base::as.character(resolved$genbank_file[i]))
    tbl <- parsed$genbank_table

    # direct gene-name shortcut first
    if (!base::is.na(ref_gene) && ref_gene != "" && "gene" %in% base::colnames(tbl)) {
      gene_hits <- tbl[!base::is.na(tbl$gene) & tbl$gene == ref_gene, , drop = FALSE]
      if (base::nrow(gene_hits) == 1) {
        resolved$locus_tag[i] <- gene_hits$locus_tag[1]
        resolved$resolved_by[i] <- "gene_exact"
        resolved$mapping_score[i] <- 1
        resolved$mapping_note[i] <- ref_gene
        next
      }
    }

    tbl$aa_seq <- vapply(tbl$nt_seq, .shared_translate_nt, character(1))
    tbl$aa_len <- base::nchar(tbl$aa_seq)
    tbl$len_ratio <- 1 - base::abs(tbl$aa_len - ref_len) / base::pmax(ref_len, 1)
    tbl$product_sim <- vapply(tbl$product, function(x) .shared_token_jaccard(ref_product, x), numeric(1))
    tbl$gene_match <- if ("gene" %in% base::colnames(tbl)) {
      base::ifelse(!base::is.na(ref_gene) & ref_gene != "" & !base::is.na(tbl$gene) & tbl$gene == ref_gene, 1, 0)
    } else 0
    tbl$kmer_sim <- vapply(tbl$aa_seq, function(x) .shared_kmer_jaccard(ref_aa, x, k = 3), numeric(1))
    tbl$prefilter_score <- 4 * tbl$gene_match + 2 * tbl$product_sim + 3 * tbl$kmer_sim + 1 * tbl$len_ratio

    top_n <- base::min(50, base::nrow(tbl))
    top_tbl <- tbl[base::order(tbl$prefilter_score, decreasing = TRUE)[base::seq_len(top_n)], , drop = FALSE]
    # Avoid expensive full dynamic-programming alignment across thousands of CDSs.
    # For shared-target resolution, token/product agreement + protein k-mer overlap
    # + length consistency is usually enough to find the ortholog in closely related genomes.
    top_tbl$align_score <- top_tbl$kmer_sim
    top_tbl$final_score <- 4 * top_tbl$gene_match + 2 * top_tbl$product_sim + 4 * top_tbl$kmer_sim + 1 * top_tbl$len_ratio

    best <- top_tbl[base::order(top_tbl$final_score, decreasing = TRUE), , drop = FALSE][1, , drop = FALSE]
    resolved$locus_tag[i] <- best$locus_tag[1]
    resolved$resolved_by[i] <- "protein_similarity"
    resolved$mapping_score[i] <- best$final_score[1]
    resolved$mapping_note[i] <- base::sprintf("gene=%s; product_sim=%.3f; kmer_sim=%.3f; align=%.3f",
                                              base::ifelse("gene" %in% base::colnames(best), best$gene[1], ""),
                                              best$product_sim[1], best$kmer_sim[1], best$align_score[1])
  }

  resolved
}

.shared_translate_nt <- function(nt_seq) {
  nt_seq <- base::toupper(base::gsub("[^ACGTN]", "N", base::as.character(nt_seq)))
  aa <- base::as.character(Biostrings::translate(Biostrings::DNAString(nt_seq), if.fuzzy.codon = "X"))
  aa <- base::sub("\\*$", "", aa)
  aa
}

.shared_token_jaccard <- function(a, b) {
  if (base::is.na(a) || base::is.na(b) || a == "" || b == "") return(0)
  tok <- function(x) unique(stringr::str_extract_all(base::tolower(x), "[a-z0-9]{3,}")[[1]])
  ta <- tok(a); tb <- tok(b)
  if (base::length(ta) == 0 || base::length(tb) == 0) return(0)
  base::length(base::intersect(ta, tb)) / base::length(base::union(ta, tb))
}

.shared_kmer_jaccard <- function(a, b, k = 3) {
  if (base::nchar(a) < k || base::nchar(b) < k) return(0)
  kmers <- function(x) unique(vapply(base::seq_len(base::nchar(x) - k + 1), function(i) base::substring(x, i, i + k - 1), character(1)))
  ka <- kmers(a); kb <- kmers(b)
  base::length(base::intersect(ka, kb)) / base::length(base::union(ka, kb))
}


.shared_is_better_grna_candidate <- function(a, b) {
  key_a <- c(base::length(a$genome_ids), -a$worst_n0, -a$worst_n1, -a$sum_n1, a$min_score, a$composite_score, -a$max_pos_dev, -base::abs(a$Tm - 60))
  key_b <- c(base::length(b$genome_ids), -b$worst_n0, -b$worst_n1, -b$sum_n1, b$min_score, b$composite_score, -b$max_pos_dev, -base::abs(b$Tm - 60))
  for (i in base::seq_along(key_a)) {
    if (key_a[i] > key_b[i]) return(TRUE)
    if (key_a[i] < key_b[i]) return(FALSE)
  }
  FALSE
}

.shared_find_guides_in_target <- function(
    target_seq,
    genome_id,
    locus_tag,
    nuclease,
    position_range,
    tm_range,
    strand
) {
  nuc <- .shared_get_nuclease_definition(nuclease)
  target_seq <- base::toupper(target_seq)
  seq_len <- base::nchar(target_seq)
  scans <- base::list(
    `+` = target_seq,
    `-` = .rc(target_seq)
  )
  if (strand == "3") scans <- scans["+"] else if (strand == "5") scans <- scans["-"]

  out <- base::list()
  idx <- 1L
  for (strand_label in base::names(scans)) {
    seq_i <- scans[[strand_label]]
    for (pam in nuc$pams) {
      pam_len <- base::nchar(pam)
      context_len <- nuc$spacer_length + pam_len
      if (base::nchar(seq_i) < context_len) next
      for (i in base::seq_len(base::nchar(seq_i) - context_len + 1L)) {
        ctx <- base::substring(seq_i, i, i + context_len - 1L)
        if (nuc$pam_side == "3prime") {
          spacer <- base::substring(ctx, 1L, nuc$spacer_length)
          pam_seq <- base::substring(ctx, nuc$spacer_length + 1L, context_len)
        } else {
          pam_seq <- base::substring(ctx, 1L, pam_len)
          spacer <- base::substring(ctx, pam_len + 1L, context_len)
        }
        if (!.shared_match_iupac(pam_seq, pam)) next
        gc <- 100 * (stringr::str_count(spacer, "[GC]") / base::nchar(spacer))
        tm <- .safe_tm(spacer)
        if (base::is.na(tm) || tm < tm_range[1] || tm > tm_range[2]) next

        if (strand_label == "+") {
          spacer_start <- if (nuc$pam_side == "3prime") i else i + pam_len
        } else {
          if (nuc$pam_side == "3prime") {
            spacer_start <- seq_len - (i + context_len - 1L) + 1L
          } else {
            spacer_start <- seq_len - (i + context_len - 1L) + pam_len + 1L
          }
        }
        pos_pct <- 100 * (spacer_start / seq_len)
        if (pos_pct < position_range[1] || pos_pct > position_range[2]) next

        tmp <- base::data.frame(
          genome_id = genome_id,
          locus_tag = locus_tag,
          protospacer = spacer,
          pam = pam_seq,
          spacer_context = if (nuc$pam_side == "3prime") base::paste0(spacer, pam_seq) else base::paste0(pam_seq, spacer),
          strand = strand_label,
          percentGC = gc,
          Tm = tm,
          position_percent = pos_pct,
          n0 = 1,
          n1 = 0,
          stringsAsFactors = FALSE
        )
        tmp <- calculate_composite_score(tmp)
        tmp$candidate_key <- base::paste(tmp$protospacer[1], tmp$pam[1], tmp$strand[1], sep = "::")
        out[[idx]] <- tmp
        idx <- idx + 1L
      }
    }
  }
  if (base::length(out) == 0L) {
    return(base::data.frame(
      genome_id = character(0), locus_tag = character(0), protospacer = character(0),
      pam = character(0), spacer_context = character(0), strand = character(0),
      percentGC = numeric(0), Tm = numeric(0), position_percent = numeric(0),
      n0 = numeric(0), n1 = numeric(0), composite_score = numeric(0), candidate_key = character(0),
      stringsAsFactors = FALSE
    ))
  }
  df <- dplyr::bind_rows(out)
  df <- df[base::order(df$composite_score, decreasing = TRUE), , drop = FALSE]
  df
}

.shared_get_nuclease_definition <- function(nuclease) {
  defs <- base::list(
    GeoCas9 = base::list(pams = c("NNNNCAAA"), pam_side = "3prime", spacer_length = 21L),
    GtCas9 = base::list(pams = c("NNNNCAAA"), pam_side = "3prime", spacer_length = 21L),
    SpCas9 = base::list(pams = c("NGG"), pam_side = "3prime", spacer_length = 20L),
    FnCas12a = base::list(pams = c("TTV"), pam_side = "5prime", spacer_length = 18L),
    FisCasI_B = base::list(pams = c("TTCA"), pam_side = "5prime", spacer_length = 28L)
  )
  if (!nuclease %in% base::names(defs)) {
    base::stop("Unsupported nuclease for shared design: ", nuclease)
  }
  defs[[nuclease]]
}

.shared_match_iupac <- function(seq, motif) {
  seq <- base::toupper(seq)
  motif <- base::toupper(motif)
  map <- base::list(
    A = c("A"), C = c("C"), G = c("G"), T = c("T"),
    R = c("A","G"), Y = c("C","T"), S = c("G","C"), W = c("A","T"),
    K = c("G","T"), M = c("A","C"), B = c("C","G","T"), D = c("A","G","T"),
    H = c("A","C","T"), V = c("A","C","G"), N = c("A","C","G","T")
  )
  if (base::nchar(seq) != base::nchar(motif)) return(FALSE)
  seq_chars <- base::strsplit(seq, "", fixed = TRUE)[[1]]
  motif_chars <- base::strsplit(motif, "", fixed = TRUE)[[1]]
  for (i in base::seq_along(seq_chars)) {
    allowed <- map[[motif_chars[i]]]
    if (base::is.null(allowed) || !seq_chars[i] %in% allowed) return(FALSE)
  }
  TRUE
}

.shared_exact_hit_count_cached_seq <- function(genome_id, primer_seq, whole_genome_map, hit_count_cache) {
  key <- base::paste(genome_id, primer_seq, sep = "::")
  if (base::exists(key, envir = hit_count_cache, inherits = FALSE)) {
    return(base::get(key, envir = hit_count_cache, inherits = FALSE))
  }
  genome_seq <- base::toupper(whole_genome_map[[genome_id]])
  primer_seq <- base::toupper(primer_seq)
  rc_seq <- .rc(primer_seq)
  hits_fwd <- gregexpr(primer_seq, genome_seq, fixed = TRUE)[[1]]
  hits_rev <- gregexpr(rc_seq, genome_seq, fixed = TRUE)[[1]]
  count_fwd <- if (hits_fwd[1] < 0) 0L else base::length(hits_fwd)
  count_rev <- if (hits_rev[1] < 0) 0L else base::length(hits_rev)
  total <- count_fwd + count_rev
  base::assign(key, total, envir = hit_count_cache)
  total
}

.shared_spacer_offtarget_counts <- function(genome_id, protospacer, whole_genome_map, hit_count_cache, n_mismatches = 1) {
  key <- base::paste(genome_id, protospacer, n_mismatches, sep = "::")
  if (base::exists(key, envir = hit_count_cache, inherits = FALSE)) {
    return(base::get(key, envir = hit_count_cache, inherits = FALSE))
  }
  genome_seq <- Biostrings::DNAString(base::toupper(whole_genome_map[[genome_id]]))
  spacer <- Biostrings::DNAString(base::toupper(protospacer))
  spacer_rc <- Biostrings::reverseComplement(spacer)

  fwd_0 <- Biostrings::matchPattern(spacer, genome_seq, max.mismatch = 0)
  rev_0 <- Biostrings::matchPattern(spacer_rc, genome_seq, max.mismatch = 0)
  n0 <- base::length(fwd_0) + base::length(rev_0)

  n1 <- 0L
  if (n_mismatches >= 1) {
    fwd_1 <- Biostrings::matchPattern(spacer, genome_seq, max.mismatch = 1)
    rev_1 <- Biostrings::matchPattern(spacer_rc, genome_seq, max.mismatch = 1)
    n1 <- base::length(fwd_1) + base::length(rev_1) - n0
  }

  val <- base::list(n0 = as.integer(n0), n1 = as.integer(n1))
  base::assign(key, val, envir = hit_count_cache)
  val
}

.shared_build_final_construct_groups <- function(construct_groups, grna_groups) {
  if (base::nrow(construct_groups) == 0 || base::nrow(grna_groups) == 0) {
    return(base::data.frame())
  }
  out <- base::list()
  idx <- 1L
  for (i in base::seq_len(base::nrow(construct_groups))) {
    donor_members <- base::strsplit(construct_groups$used_by[i], ", ", fixed = TRUE)[[1]]
    for (j in base::seq_len(base::nrow(grna_groups))) {
      grna_members <- base::strsplit(grna_groups$used_by[j], ", ", fixed = TRUE)[[1]]
      members <- base::intersect(donor_members, grna_members)
      if (base::length(members) == 0) next
      label_core <- if (base::length(members) == 1L) {
        members[1]
      } else {
        base::paste0(base::paste(members, collapse = "_"), "_common")
      }
      label <- base::paste0(label_core, "__", grna_groups$cluster_id[j])
      out[[idx]] <- base::data.frame(
        final_group_id = base::paste0("final_group_", idx),
        construct_label = label,
        used_by = base::paste(members, collapse = ", "),
        donor_group_id = construct_groups$construct_group_id[i],
        grna_group_id = grna_groups$cluster_id[j],
        representative_genome = members[1],
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  dplyr::bind_rows(out)
}

.shared_write_final_construct_groups <- function(
    final_construct_groups,
    target_context,
    grna_groups,
    construct_deletion_primers,
    vector_file,
    grna_start,
    grna_end,
    deletion_start,
    deletion_end,
    upstream_bp,
    downstream_bp,
    output_dir
) {
  if (base::nrow(final_construct_groups) == 0) return(invisible(NULL))
  for (i in base::seq_len(base::nrow(final_construct_groups))) {
    rep_id <- final_construct_groups$representative_genome[i]
    label <- final_construct_groups$construct_label[i]
    # construct_deletion_primers$construct_label is the donor label (e.g.
    # "..._common"), while `label` above includes the "__<grna_cluster>"
    # suffix. Match on the donor portion so the primer annotations get
    # attached to the combined .gbk output.
    donor_label <- base::sub("__.*$", "", label)
    ctx <- target_context[target_context$genome_id == rep_id, , drop = FALSE]
    grna_row <- grna_groups[grna_groups$cluster_id == final_construct_groups$grna_group_id[i], , drop = FALSE]
    dp_row <- if (!base::is.null(construct_deletion_primers) && base::nrow(construct_deletion_primers) > 0) {
      construct_deletion_primers[construct_deletion_primers$construct_label == donor_label, , drop = FALSE]
    } else {
      base::data.frame()
    }
    out_path <- base::file.path(output_dir, base::paste0(label, "_combined_construct.gbk"))
    .write_combined_construct_genbank(
      vector_file = vector_file,
      spacer_seq = grna_row$protospacer[1],
      grna_start = grna_start,
      grna_end = grna_end,
      insert_seq = if ("effective_insert_seq" %in% base::colnames(ctx)) ctx$effective_insert_seq[1] else ctx$insert_seq[1],
      deletion_start = deletion_start,
      deletion_end = deletion_end,
      locus_tag = label,
      gene = ctx$target_gene[1],
      gRNA_name = grna_row$cluster_id[1],
      nuclease = "shared_design",
      grna_primer_info = NULL,
      deletion_primers = if (base::nrow(dp_row) == 1L) {
        pick_name <- function(row, anno_col, order_col) {
          if (anno_col %in% base::colnames(row) && !base::is.na(row[[anno_col]][1])) {
            row[[anno_col]][1]
          } else {
            row[[order_col]][1]
          }
        }
        get_col <- function(row, col) {
          if (col %in% base::colnames(row)) row[[col]][1] else NA_real_
        }
        base::list(
          upstream_forward_name = pick_name(dp_row, "upstream_forward_annotation_name", "upstream_forward_name"),
          upstream_forward_primer = dp_row$upstream_forward_primer[1],
          upstream_forward_tm_target = get_col(dp_row, "upstream_forward_tm_target"),
          upstream_forward_tm_full = get_col(dp_row, "upstream_forward_tm_full"),
          upstream_forward_start = dp_row$upstream_forward_start[1],
          upstream_forward_end = dp_row$upstream_forward_end[1],
          upstream_reverse_name = pick_name(dp_row, "upstream_reverse_annotation_name", "upstream_reverse_name"),
          upstream_reverse_primer = dp_row$upstream_reverse_primer[1],
          upstream_reverse_tm_target = get_col(dp_row, "upstream_reverse_tm_target"),
          upstream_reverse_tm_full = get_col(dp_row, "upstream_reverse_tm_full"),
          upstream_reverse_start = dp_row$upstream_reverse_start[1],
          upstream_reverse_end = dp_row$upstream_reverse_end[1],
          downstream_forward_name = pick_name(dp_row, "downstream_forward_annotation_name", "downstream_forward_name"),
          downstream_forward_primer = dp_row$downstream_forward_primer[1],
          downstream_forward_tm_target = get_col(dp_row, "downstream_forward_tm_target"),
          downstream_forward_tm_full = get_col(dp_row, "downstream_forward_tm_full"),
          downstream_forward_start = dp_row$downstream_forward_start[1],
          downstream_forward_end = dp_row$downstream_forward_end[1],
          downstream_reverse_name = pick_name(dp_row, "downstream_reverse_annotation_name", "downstream_reverse_name"),
          downstream_reverse_primer = dp_row$downstream_reverse_primer[1],
          downstream_reverse_tm_target = get_col(dp_row, "downstream_reverse_tm_target"),
          downstream_reverse_tm_full = get_col(dp_row, "downstream_reverse_tm_full"),
          downstream_reverse_start = dp_row$downstream_reverse_start[1],
          downstream_reverse_end = dp_row$downstream_reverse_end[1],
          overlap_tm = dp_row$overlap_tm[1],
          overlap_length = dp_row$overlap_length[1]
        )
      } else NULL,
      upstream_bp = if ("effective_upstream_bp" %in% base::colnames(ctx)) ctx$effective_upstream_bp[1] else upstream_bp,
      downstream_bp = if ("effective_downstream_bp" %in% base::colnames(ctx)) ctx$effective_downstream_bp[1] else downstream_bp,
      output_path = out_path,
      kill_snapgene = FALSE
    )
  }
  invisible(NULL)
}
