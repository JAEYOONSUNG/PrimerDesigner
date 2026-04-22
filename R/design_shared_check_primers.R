# Shared colony-PCR ("check") primer design
#
# Places a forward primer (cF) outside the UF homology arm and a reverse
# primer (cR) outside the DR homology arm, shared across genomes where
# possible. Product length on the wild-type genome and after deletion is
# computed per (genome, construct_group) pair so the user can see both
# expected bands on a colony PCR gel.
#
# Reuses `.shared_assign_role_groups()` to get the same subgroup-splitting
# behaviour as the main primer design, running against a synthetic
# "upstream_check_seq" / "downstream_check_seq" arm field layout.

.shared_check_primer_heterodimer <- function(a, b, threshold = 6L) {
  if (base::is.na(a) || base::is.na(b) ||
      base::nchar(a) < threshold || base::nchar(b) < threshold) return(FALSE)
  a <- base::toupper(a); b <- base::toupper(b)
  rc_b <- .rc(b)
  n <- base::nchar(a) - threshold + 1L
  for (i in base::seq_len(n)) {
    w <- base::substring(a, i, i + threshold - 1L)
    if (base::grepl(w, rc_b, fixed = TRUE)) return(TRUE)
  }
  FALSE
}

.shared_find_oligo_genomic_pos <- function(oligo, genome_seq) {
  if (base::is.na(oligo) || base::nchar(oligo) == 0) return(NULL)
  oligo_u <- base::toupper(oligo)
  gs_u <- base::toupper(genome_seq)
  rc <- .rc(oligo_u)
  fwd <- base::unlist(base::gregexpr(oligo_u, gs_u, fixed = TRUE))
  rev <- base::unlist(base::gregexpr(rc,      gs_u, fixed = TRUE))
  fwd <- fwd[fwd > 0]
  rev <- rev[rev > 0]
  if (base::length(fwd) == 1L && base::length(rev) == 0L) {
    return(base::list(pos = fwd, strand = "+"))
  }
  if (base::length(rev) == 1L && base::length(fwd) == 0L) {
    return(base::list(pos = rev, strand = "-"))
  }
  # ambiguous / multi-hit -> let caller report worst_n0 > 1.
  NULL
}

#' Internal: design shared colony-PCR ("check") primers outside the arms.
#'
#' @param target_context data frame with per-genome target / arm coordinates.
#' @param genome_sequences named character vector of whole-genome sequences.
#' @param dnastr_cache named list of Biostrings::DNAString objects (used by
#'   the pool builder for fast uniqueness checks).
#' @param forbidden_oligos character vector of main-design primer oligos
#'   that check primers must not revcomp-collide with.
#' @return list with $check_primer_cores (primer_cores-shaped data frame for
#'   roles cF and cR), $check_pairs (per-genome pairing + PCR product sizes),
#'   $warnings character vector.
#' @keywords internal
.shared_design_check_primers <- function(
    target_context,
    genome_sequences,
    dnastr_cache = NULL,
    check_outer_pad = 50L,
    check_search_window = 800L,
    check_search_window_max = 6000L,
    tm_target = 55,
    tm_tolerance = 3,
    pair_dtm_max = 2,
    primer_min_length = 18L,
    primer_max_length = 25L,
    require_unique = TRUE,
    forbidden_oligos = character(0),
    dimer_threshold = 6L,
    verbose = TRUE
) {
  check_outer_pad <- base::as.integer(check_outer_pad)
  check_search_window <- base::as.integer(check_search_window)
  check_search_window_max <- base::as.integer(check_search_window_max)
  primer_min_length <- base::as.integer(primer_min_length)
  primer_max_length <- base::as.integer(primer_max_length)

  genome_ids <- base::as.character(target_context$genome_id)
  warnings_acc <- base::character(0)

  # Builder for the per-genome check context, parameterised on window size so
  # progressive widening can call it with larger values until all genomes
  # collapse into a single cluster.
  build_context <- function(window_bp) {
    window_bp <- base::as.integer(window_bp)
    check_ctx_list <- base::vector("list", base::length(genome_ids))
    for (i in base::seq_along(genome_ids)) {
      gid <- genome_ids[i]
      gs <- base::as.character(genome_sequences[[gid]])
      glen <- base::nchar(gs)
      del_start <- base::as.integer(target_context$deletion_start[i])
      del_end   <- base::as.integer(target_context$deletion_end[i])
      pref_up   <- base::as.integer(target_context$preferred_upstream_bp[i])
      pref_dn   <- base::as.integer(target_context$preferred_downstream_bp[i])

      la_outer <- del_start - pref_up
      ra_outer <- del_end   + pref_dn

      up_end   <- la_outer - check_outer_pad - 1L
      up_start <- up_end - window_bp + 1L
      dn_start <- ra_outer + check_outer_pad + 1L
      dn_end   <- dn_start + window_bp - 1L

      up_start <- base::max(1L, up_start)
      dn_end   <- base::min(glen, dn_end)

      up_seq <- if (up_end   >= up_start) base::substring(gs, up_start, up_end) else ""
      dn_seq <- if (dn_end   >= dn_start) base::substring(gs, dn_start, dn_end) else ""

      check_ctx_list[[i]] <- base::data.frame(
        genome_id = gid,
        upstream_check_seq   = up_seq,
        downstream_check_seq = dn_seq,
        upstream_check_start = up_start,
        upstream_check_end   = up_end,
        downstream_check_start = dn_start,
        downstream_check_end   = dn_end,
        preferred_cf_bp = base::nchar(up_seq),
        preferred_cr_bp = base::nchar(dn_seq),
        stringsAsFactors = FALSE
      )
    }
    check_ctx <- base::do.call(base::rbind, check_ctx_list)
    stats::setNames(
      base::lapply(base::seq_along(genome_ids), function(i) {
        base::list(
          genome_id = genome_ids[i],
          genome_seq = base::as.character(genome_sequences[[genome_ids[i]]]),
          context = check_ctx[i, , drop = FALSE]
        )
      }), genome_ids)
  }

  hit_count_cache <- base::new.env(parent = base::emptyenv())

  # 2) Call the shared pool builder + selector for cF and cR.
  # max_arm_bp is set well beyond check_search_window so no candidate is
  # dropped by the implied_arm filter (scoring still uses preferred_bp).
  tm_window <- base::c(tm_target - tm_tolerance, tm_target + tm_tolerance)

  run_role <- function(context_map, role, arm_field, region, reverse,
                        preferred_bp_field, extra_forbidden, window_bp) {
    full_forbidden <- base::union(base::as.character(forbidden_oligos),
                                   base::as.character(extra_forbidden))
    pool <- .shared_build_primer_candidate_pool(
      genome_ids = genome_ids,
      context_map = context_map,
      dnastr_cache = dnastr_cache,
      role = role,
      arm_field = arm_field,
      region = region,
      reverse = reverse,
      boundary = "outer",
      preferred_bp_field = preferred_bp_field,
      min_arm_bp = 1L,
      max_arm_bp = window_bp * 4L,
      search_window = window_bp,
      min_primer_length = primer_min_length,
      max_primer_length = primer_max_length,
      tm_target = tm_target,
      require_unique = require_unique,
      hit_count_cache = hit_count_cache,
      candidate_cap = NA_integer_,
      forbidden_oligos = full_forbidden,
      inner_boundary_tolerance_bp = NULL,
      tm_window = tm_window
    )
    .shared_select_primer_groups_from_pool(pool, genome_ids, role, reverse)
  }

  # Progressive widening loop: the user expects check primers shared across
  # every genome whenever biologically possible. Start with the requested
  # window, and if the role came back split into multiple clusters, double
  # the window until all genomes fall into one cluster OR the window hits
  # check_search_window_max. Any remaining split at that point is a genuine
  # biological divergence (e.g. mobile elements with different integration
  # sites per strain) and we fall back to subgroup-split.
  run_role_widen <- function(role, arm_field, region, reverse,
                              preferred_bp_field, extra_forbidden) {
    window_bp <- check_search_window
    best <- NULL
    while (TRUE) {
      attempt <- base::tryCatch(
        run_role(build_context(window_bp), role, arm_field, region, reverse,
                  preferred_bp_field, extra_forbidden, window_bp),
        error = function(e) e
      )
      if (!inherits(attempt, "error")) {
        best <- attempt
        if (base::nrow(attempt) == 1L) {
          if (verbose && window_bp != check_search_window) {
            base::message(sprintf(
              "[check] role %s unified across all genomes at window = %d bp",
              role, window_bp))
          }
          return(attempt)
        }
      }
      if (window_bp >= check_search_window_max) break
      next_bp <- base::min(check_search_window_max, window_bp * 2L)
      if (next_bp == window_bp) break
      if (verbose) {
        base::message(sprintf(
          "[check] role %s: %d cluster(s) at window = %d bp; widening to %d",
          role,
          if (base::is.null(best)) NA_integer_ else base::nrow(best),
          window_bp, next_bp))
      }
      window_bp <- next_bp
    }
    if (base::is.null(best)) base::stop("role ", role, ": no candidate after widening")
    best
  }

  cf_rows <- base::tryCatch(
    run_role_widen("cF", "upstream_check_seq",   "end",   FALSE,
                    "preferred_cf_bp", base::character(0)),
    error = function(e) {
      warnings_acc <<- base::c(warnings_acc,
        base::paste0("check primer cF design failed: ",
                     base::conditionMessage(e)))
      NULL
    }
  )
  if (base::is.null(cf_rows)) return(base::list(
    check_primer_cores = .shared_empty_check_primer_cores(),
    check_pairs = .shared_empty_check_pairs(),
    warnings = warnings_acc
  ))

  cr_rows <- base::tryCatch(
    run_role_widen("cR", "downstream_check_seq", "start", TRUE,
                    "preferred_cr_bp",
                    base::as.character(cf_rows$primer_core_5to3)),
    error = function(e) {
      warnings_acc <<- base::c(warnings_acc,
        base::paste0("check primer cR design failed: ",
                     base::conditionMessage(e)))
      NULL
    }
  )
  if (base::is.null(cr_rows)) return(base::list(
    check_primer_cores = cf_rows,
    check_pairs = .shared_empty_check_pairs(),
    warnings = warnings_acc
  ))

  check_primer_cores <- dplyr::bind_rows(cf_rows, cr_rows)

  # 3) Tm-window and self-dimer sanity pass (the pool builder already filters
  #    these but we annotate explicit warnings for transparency).
  for (i in base::seq_len(base::nrow(check_primer_cores))) {
    oli <- check_primer_cores$primer_core_5to3[i]
    tm  <- check_primer_cores$tm_core[i]
    if (!base::is.na(tm) && base::abs(tm - tm_target) > tm_tolerance) {
      warnings_acc <- base::c(warnings_acc, base::sprintf(
        "%s Tm %.1f falls outside target %g ± %g",
        check_primer_cores$cluster_id[i], tm, tm_target, tm_tolerance))
    }
  }

  # 4) Pair cF and cR per genome, reject dimer pairs, compute product sizes.
  cf_rows <- check_primer_cores[check_primer_cores$role == "cF", , drop = FALSE]
  cr_rows <- check_primer_cores[check_primer_cores$role == "cR", , drop = FALSE]

  lookup_cluster_for_genome <- function(rows, gid) {
    m <- base::vapply(base::strsplit(rows$used_by, ",\\s*"),
                       function(x) gid %in% base::trimws(x), logical(1))
    if (!base::any(m)) return(NULL)
    rows[base::which(m)[1], , drop = FALSE]
  }

  pair_rows <- base::list()
  for (gi in base::seq_along(genome_ids)) {
    gid <- genome_ids[gi]
    cf <- lookup_cluster_for_genome(cf_rows, gid)
    cr <- lookup_cluster_for_genome(cr_rows, gid)
    if (base::is.null(cf) || base::is.null(cr)) next

    # Heterodimer check.
    has_dimer <- .shared_check_primer_heterodimer(
      cf$primer_core_5to3, cr$primer_core_5to3, threshold = dimer_threshold)
    dtm <- if (base::is.na(cf$tm_core) || base::is.na(cr$tm_core)) NA_real_ else
      base::abs(cf$tm_core - cr$tm_core)
    pair_dtm_ok <- !base::is.na(dtm) && dtm <= pair_dtm_max

    # Genomic positions of each primer in the WT genome.
    gs <- base::as.character(genome_sequences[[gid]])
    cf_hit <- .shared_find_oligo_genomic_pos(cf$primer_core_5to3, gs)
    cr_hit <- .shared_find_oligo_genomic_pos(cr$primer_core_5to3, gs)

    del_len <- base::as.integer(target_context$deletion_end[gi] -
                                 target_context$deletion_start[gi] + 1L)

    if (base::is.null(cf_hit) || base::is.null(cr_hit)) {
      warnings_acc <- base::c(warnings_acc, base::sprintf(
        "%s: cF/cR did not resolve to a unique genomic position (skipping product-size calc).",
        gid))
      product_wt <- NA_integer_
      product_del <- NA_integer_
    } else {
      # cF is on forward strand (5'->3'), starts at cf_hit$pos.
      # cR is reverse primer, so its RC matches forward strand at cr_hit$pos;
      # the WT amplicon extends from cF start to cR RC-match + k - 1.
      k_cf <- base::as.integer(cf$primer_length)
      k_cr <- base::as.integer(cr$primer_length)
      cF_fwd_start <- base::as.integer(cf_hit$pos)
      cR_fwd_start <- base::as.integer(cr_hit$pos)
      product_wt <- cR_fwd_start + k_cr - 1L - cF_fwd_start + 1L
      product_del <- product_wt - del_len
    }

    if (has_dimer) {
      warnings_acc <- base::c(warnings_acc, base::sprintf(
        "%s: cF %s and cR %s may form a primer dimer (>=%d bp 3' complementarity).",
        gid, cf$cluster_id, cr$cluster_id, dimer_threshold))
    }
    if (!pair_dtm_ok && !base::is.na(dtm)) {
      warnings_acc <- base::c(warnings_acc, base::sprintf(
        "%s: cF/cR pair ΔTm %.1f > %g",
        gid, dtm, pair_dtm_max))
    }

    pair_rows[[base::length(pair_rows) + 1L]] <- base::data.frame(
      genome_id = gid,
      cF_cluster = cf$cluster_id,
      cR_cluster = cr$cluster_id,
      cF_sequence = cf$primer_core_5to3,
      cR_sequence = cr$primer_core_5to3,
      cF_length = base::as.integer(cf$primer_length),
      cR_length = base::as.integer(cr$primer_length),
      cF_tm = cf$tm_core,
      cR_tm = cr$tm_core,
      pair_dtm = dtm,
      pair_dimer = has_dimer,
      product_wt_bp = product_wt,
      product_deletion_bp = product_del,
      delta_bp = if (base::is.na(product_wt)) NA_integer_ else product_wt - product_del,
      stringsAsFactors = FALSE
    )
  }
  check_pairs <- if (base::length(pair_rows) == 0L)
    .shared_empty_check_pairs()
  else
    dplyr::bind_rows(pair_rows)

  base::list(
    check_primer_cores = check_primer_cores,
    check_pairs = check_pairs,
    warnings = base::unique(warnings_acc)
  )
}

.shared_empty_check_primer_cores <- function() {
  base::data.frame(
    role = character(0), cluster_id = character(0), scope = character(0),
    used_by = character(0), primer_core_5to3 = character(0),
    shared_core = character(0), reverse_primer = logical(0),
    primer_length = integer(0), gc_percent = numeric(0), tm_core = numeric(0),
    worst_arm_dev_bp = integer(0), mean_arm_dev_bp = numeric(0),
    boundary_distance_bp = integer(0), exact_unique_in_targets = logical(0),
    stringsAsFactors = FALSE
  )
}

.shared_empty_check_pairs <- function() {
  base::data.frame(
    genome_id = character(0), cF_cluster = character(0), cR_cluster = character(0),
    cF_sequence = character(0), cR_sequence = character(0),
    cF_length = integer(0), cR_length = integer(0),
    cF_tm = numeric(0), cR_tm = numeric(0),
    pair_dtm = numeric(0), pair_dimer = logical(0),
    product_wt_bp = integer(0), product_deletion_bp = integer(0),
    delta_bp = integer(0),
    stringsAsFactors = FALSE
  )
}
