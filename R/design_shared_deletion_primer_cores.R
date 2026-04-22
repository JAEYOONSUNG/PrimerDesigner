#' Design Shared Deletion Primer Cores Across Multiple Genomes
#'
#' Selects genome-annealing core primers for deletion-style Gibson assembly
#' across multiple related genomes. This function is designed for the case where
#' the same functional target is edited in several strains and the user wants to
#' minimize the number of ordered primers by preferring shared cores, while still
#' enforcing exact uniqueness inside each intended target genome.
#'
#' This function does **not** generate full Gibson oligos with vector tails.
#' Instead, it chooses the core annealing segments (\code{UF}, \code{UR},
#' \code{DF}, \code{DR}) that can later be passed into a vector-tail workflow.
#' The existing \code{design_deletion_primers()} remains the right tool for
#' per-genome full oligo generation once a final vector insertion site has been
#' chosen.
#'
#' Supported workflow:
#' \enumerate{
#'   \item Parse multiple genome GenBank files.
#'   \item Extract deletion arms for each target locus.
#'   \item Optionally preserve the full upstream CDS (useful when the upstream
#'         gene overlaps the deletion target).
#'   \item Search boundary windows for exact shared primer cores.
#'   \item Enforce exact uniqueness in every intended target genome.
#'   \item Greedily fall back from all-genome sharing to subgroup sharing to
#'         strain-specific primers when required.
#' }
#'
#' Expected input table columns:
#' \describe{
#'   \item{genome_id}{Unique strain identifier.}
#'   \item{genbank_file}{Genome GenBank path.}
#'   \item{locus_tag}{Target locus tag to delete.}
#' }
#'
#' @param target_table Data frame with at least \code{genome_id},
#'   \code{genbank_file}, and \code{locus_tag}.
#' @param upstream_bp Integer. Preferred left arm size in bp (soft target used
#'   for scoring). Actual arm length may flex within
#'   \code{[min_arm_bp, max_arm_bp]} so that a better shared primer can be found.
#' @param downstream_bp Integer. Preferred right arm size in bp (soft target).
#' @param min_arm_bp Integer. Hard lower bound on the effective homology arm
#'   length. Default 200.
#' @param max_arm_bp Integer. Hard upper bound on the effective homology arm
#'   length; also used as the length of the extended arm extracted for primer
#'   search. Default 1000.
#' @param preserve_upstream_gene Either a logical or \code{"auto"}. If
#'   \code{TRUE}, the left arm ends at the end of the immediately upstream CDS
#'   instead of \code{target_start - 1}. If \code{"auto"}, preservation is
#'   enabled when the upstream CDS overlaps the target or sits within
#'   \code{auto_preserve_gap_bp} of it.
#' @param preserve_upstream_buffer_bp Integer. Safety buffer in bp kept between
#'   the upstream CDS end and the deletion start when \code{preserve_upstream_gene}
#'   is active. Default 0.
#' @param auto_preserve_gap_bp Integer. Maximum intergenic gap (bp) that still
#'   triggers preservation in \code{"auto"} mode. Default 100.
#' @param search_window Integer or NULL. Legacy boundary window (bp) used if
#'   provided. When NULL (default) the entire extracted arm is scanned so that
#'   shared primers can be discovered anywhere in the flexible arm range.
#' @param min_primer_length Integer. Minimum primer core length.
#' @param max_primer_length Integer. Maximum primer core length.
#' @param tm_target Numeric. Preferred Tm target for ranking.
#' @param require_unique Logical. Require exact 1-hit uniqueness in each target
#'   genome for a candidate to be accepted.
#' @param max_n1_per_primer Optional integer. After core selection each primer
#'   is re-audited with \code{Biostrings::matchPattern} on every assigned
#'   genome. If any primer's worst 1-mismatch hit count exceeds this value a
#'   warning is emitted. Per-primer, per-genome \code{n0}/\code{n1} counts are
#'   always added to \code{primer_cores}.
#' @param vector_file Optional reference vector path. When supplied together with
#'   \code{start}, \code{end}, and \code{construct_output_dir}, representative
#'   donor construct GenBank files are generated for each unique insert
#'   sequence group.
#' @param start Integer. Vector insertion site start position (1-based).
#' @param end Integer. Vector insertion site end position (1-based).
#' @param construct_output_dir Optional directory for representative donor
#'   construct GenBank files. Identical donor inserts across genomes are written
#'   once as a shared construct, e.g. \code{strainA_strainB_common.gbk}.
#' @param output_file Optional path to write an Excel file with sheets
#'   \code{primer_cores}, \code{target_context}, and \code{construct_groups}.
#' @return A list with:
#'   \describe{
#'     \item{primer_cores}{Data frame of chosen shared / subgroup / strain-specific cores.}
#'     \item{target_context}{Data frame summarizing extracted arm and deletion coordinates.}
#'     \item{construct_groups}{Data frame grouping genomes by identical donor insert sequence.}
#'     \item{genome_sequences}{Named list of full genome sequences used for uniqueness checks.}
#'   }
#' @examples
#' \dontrun{
#' targets <- data.frame(
#'   genome_id = c("strainA", "strainB"),
#'   genbank_file = c("strainA.gbk", "strainB.gbk"),
#'   locus_tag = c("GENE_001", "GENE_104"),
#'   stringsAsFactors = FALSE
#' )
#'
#' res <- design_shared_deletion_primer_cores(
#'   target_table = targets,
#'   upstream_bp = 500,
#'   downstream_bp = 500,
#'   preserve_upstream_gene = TRUE
#' )
#' }
#' @export
design_shared_deletion_primer_cores <- function(
    target_table,
    upstream_bp = 500,
    downstream_bp = 500,
    min_arm_bp = 200L,
    max_arm_bp = 1000L,
    preserve_upstream_gene = FALSE,
    preserve_upstream_buffer_bp = 0L,
    auto_preserve_gap_bp = 100L,
    overlap_policy = c("strict", "flexible"),
    search_window = 500L,
    min_primer_length = 20,
    max_primer_length = 35,
    tm_target = 60,
    require_unique = TRUE,
    audit_n_mismatches = 0L,
    max_n1_per_primer = NULL,
    verbose = TRUE,
    vector_file = NULL,
    start = NULL,
    end = NULL,
    construct_output_dir = NULL,
    output_file = NULL,
    inner_boundary_tolerance_bp = 50L,
    design_check_primers = TRUE,
    check_outer_pad = 50L,
    check_search_window = 800L,
    check_search_window_max = 6000L,
    check_tm_target = 55,
    check_tm_tolerance = 3,
    check_pair_dtm_max = 2,
    check_primer_min_length = 18L,
    check_primer_max_length = 25L
) {
  min_arm_bp <- base::as.integer(min_arm_bp)
  max_arm_bp <- base::as.integer(max_arm_bp)
  overlap_policy <- base::match.arg(overlap_policy)
  if (max_arm_bp < min_arm_bp) {
    base::stop("max_arm_bp must be >= min_arm_bp.")
  }
  if (upstream_bp < min_arm_bp || upstream_bp > max_arm_bp) {
    base::stop("upstream_bp (preferred) must lie within [min_arm_bp, max_arm_bp].")
  }
  if (downstream_bp < min_arm_bp || downstream_bp > max_arm_bp) {
    base::stop("downstream_bp (preferred) must lie within [min_arm_bp, max_arm_bp].")
  }
  required_cols <- c("genome_id", "genbank_file", "locus_tag")
  missing_cols <- base::setdiff(required_cols, base::colnames(target_table))
  if (!base::is.data.frame(target_table) || base::nrow(target_table) == 0) {
    base::stop("target_table must be a non-empty data frame.")
  }
  if (base::length(missing_cols) > 0) {
    base::stop(
      "target_table is missing required columns: ",
      base::paste(missing_cols, collapse = ", ")
    )
  }

  parsed_targets <- lapply(base::seq_len(base::nrow(target_table)), function(i) {
    genome_id <- base::as.character(target_table$genome_id[i])
    locus_tag <- base::as.character(target_table$locus_tag[i])
    source_file <- base::as.character(target_table$genbank_file[i])
    parsed <- .shared_read_genome_cached(source_file)
    .shared_extract_target_context(
      genome_id = genome_id,
      locus_tag = locus_tag,
      genome_seq = parsed$genome_seq,
      genbank_table = parsed$genbank_table,
      preferred_upstream_bp = upstream_bp,
      preferred_downstream_bp = downstream_bp,
      max_arm_bp = max_arm_bp,
      min_arm_bp = min_arm_bp,
      preserve_upstream_gene = preserve_upstream_gene,
      preserve_upstream_buffer_bp = base::as.integer(preserve_upstream_buffer_bp),
      auto_preserve_gap_bp = base::as.integer(auto_preserve_gap_bp),
      overlap_policy = overlap_policy
    )
  })

  target_context <- dplyr::bind_rows(lapply(parsed_targets, `[[`, "context"))
  context_map <- stats::setNames(parsed_targets, vapply(parsed_targets, function(x) x$genome_id, character(1)))
  genome_ids <- target_context$genome_id
  hit_count_cache <- new.env(parent = emptyenv())
  genome_sequences <- stats::setNames(lapply(parsed_targets, `[[`, "genome_seq"), genome_ids)
  # Pre-build Biostrings DNAString objects once per genome; used for fast
  # C-level matchPattern uniqueness checks instead of R's gregexpr.
  dnastr_cache <- if (base::requireNamespace("Biostrings", quietly = TRUE)) {
    stats::setNames(
      base::lapply(genome_sequences, function(seq) {
        Biostrings::DNAString(base::toupper(seq))
      }),
      genome_ids
    )
  } else {
    NULL
  }

  role_defs <- base::list(
    UF = base::list(arm = "left_arm_seq", region = "start", reverse = FALSE,
                    boundary = "outer", preferred_bp_field = "preferred_upstream_bp"),
    UR = base::list(arm = "left_arm_seq", region = "end", reverse = TRUE,
                    boundary = "inner", preferred_bp_field = "preferred_upstream_bp"),
    DF = base::list(arm = "right_arm_seq", region = "start", reverse = FALSE,
                    boundary = "inner", preferred_bp_field = "preferred_downstream_bp"),
    DR = base::list(arm = "right_arm_seq", region = "end", reverse = TRUE,
                    boundary = "outer", preferred_bp_field = "preferred_downstream_bp")
  )

  # Cross-role exclusion: when selecting the reverse-strand role on an arm
  # (UR after UF, DR after DF), forbid any candidate whose reverse complement
  # matches a primer_core already chosen for the forward-strand role on the
  # same arm. Otherwise the pair can end up as RC-of-each-other, which binds
  # the same strand and reduces the effective arm to the primer length.
  primer_rows <- base::list()
  pair_of <- base::list(UF = NULL, UR = "UF", DF = NULL, DR = "DF")
  for (role in base::names(role_defs)) {
    forbidden <- base::character(0)
    partner <- pair_of[[role]]
    if (!base::is.null(partner) && !base::is.null(primer_rows[[partner]])) {
      forbidden <- base::as.character(primer_rows[[partner]]$primer_core_5to3)
      forbidden <- forbidden[!base::is.na(forbidden) & base::nchar(forbidden) > 0]
    }
    if (verbose) base::message(sprintf("[shared] selecting role %s ...", role))
    t_role <- base::Sys.time()
    primer_rows[[role]] <- .shared_assign_role_groups(
      genome_ids = genome_ids,
      context_map = context_map,
      dnastr_cache = dnastr_cache,
      role = role,
      arm_field = role_defs[[role]]$arm,
      region = role_defs[[role]]$region,
      reverse = role_defs[[role]]$reverse,
      boundary = role_defs[[role]]$boundary,
      preferred_bp_field = role_defs[[role]]$preferred_bp_field,
      min_arm_bp = min_arm_bp,
      max_arm_bp = max_arm_bp,
      search_window = search_window,
      min_primer_length = min_primer_length,
      max_primer_length = max_primer_length,
      tm_target = tm_target,
      require_unique = require_unique,
      hit_count_cache = hit_count_cache,
      forbidden_oligos = forbidden,
      inner_boundary_tolerance_bp = inner_boundary_tolerance_bp,
      overlap_policy = overlap_policy
    )
    if (verbose) {
      base::message(sprintf("[shared]   %s done in %.1fs (%d cluster(s))", role,
        base::as.numeric(base::difftime(base::Sys.time(), t_role, units = "secs")),
        base::nrow(primer_rows[[role]])))
    }
  }

  primer_cores <- dplyr::bind_rows(primer_rows)
  primer_cores <- dplyr::mutate(
    primer_cores,
    preserve_upstream_gene = base::as.character(preserve_upstream_gene),
    preferred_upstream_bp = upstream_bp,
    preferred_downstream_bp = downstream_bp,
    min_arm_bp = min_arm_bp,
    max_arm_bp = max_arm_bp
  )
  # Secondary verification: the primer must only bind at its intended arm
  # position and nowhere else in every genome it is assigned to. We annotate
  # exact (n0) hit counts always; 1-mismatch (n1) counts only when
  # `audit_n_mismatches >= 1` OR `max_n1_per_primer` is supplied (otherwise
  # the matchPattern(mismatch=1) scan is skipped — it is the slowest step
  # on multi-megabase bacterial genomes).
  if (verbose) base::message("[shared] running primer uniqueness audit ...")
  t_audit <- base::Sys.time()
  effective_mm <- base::as.integer(audit_n_mismatches)
  if (!base::is.null(max_n1_per_primer)) effective_mm <- base::max(effective_mm, 1L)
  primer_cores <- .shared_annotate_primer_uniqueness(
    primer_cores = primer_cores,
    genome_sequences = genome_sequences,
    max_n1_per_primer = max_n1_per_primer,
    audit_n_mismatches = effective_mm
  )
  if (verbose) {
    base::message(sprintf("[shared] audit done in %.1fs",
      base::as.numeric(base::difftime(base::Sys.time(), t_audit, units = "secs"))))
  }
  target_context <- .shared_apply_primer_geometry(
    target_context, primer_cores,
    min_arm_bp = min_arm_bp,
    max_arm_bp = max_arm_bp
  )
  construct_groups <- .shared_build_construct_groups(target_context)

  check_primer_cores <- .shared_empty_check_primer_cores()
  check_pairs <- .shared_empty_check_pairs()
  if (base::isTRUE(design_check_primers)) {
    if (verbose) base::message("[shared] designing colony-PCR check primers ...")
    main_oligos <- base::as.character(primer_cores$primer_core_5to3)
    main_oligos <- main_oligos[!base::is.na(main_oligos) & base::nchar(main_oligos) > 0]
    chk <- .shared_design_check_primers(
      target_context = target_context,
      genome_sequences = genome_sequences,
      dnastr_cache = dnastr_cache,
      check_outer_pad = check_outer_pad,
      check_search_window = check_search_window,
      check_search_window_max = check_search_window_max,
      tm_target = check_tm_target,
      tm_tolerance = check_tm_tolerance,
      pair_dtm_max = check_pair_dtm_max,
      primer_min_length = check_primer_min_length,
      primer_max_length = check_primer_max_length,
      require_unique = require_unique,
      forbidden_oligos = main_oligos,
      verbose = verbose
    )
    check_primer_cores <- chk$check_primer_cores
    check_pairs <- chk$check_pairs
    for (w in chk$warnings) base::warning(w, call. = FALSE)
  }

  if (!base::is.null(construct_output_dir)) {
    if (base::is.null(vector_file) || base::is.null(start) || base::is.null(end)) {
      base::stop("vector_file, start, and end are required when construct_output_dir is set.")
    }
    if (!base::dir.exists(construct_output_dir)) {
      base::dir.create(construct_output_dir, recursive = TRUE)
    }
    .shared_write_construct_groups(
      construct_groups = construct_groups,
      target_context = target_context,
      vector_file = vector_file,
      start = start,
      end = end,
      fallback_upstream_bp = upstream_bp,
      fallback_downstream_bp = downstream_bp,
      construct_output_dir = construct_output_dir
    )
  }

  if (!base::is.null(output_file)) {
    writexl::write_xlsx(
      x = base::list(
        primer_cores = primer_cores,
        target_context = target_context,
        construct_groups = construct_groups
      ),
      path = output_file
    )
  }

  return(base::list(
    primer_cores = primer_cores,
    target_context = target_context,
    construct_groups = construct_groups,
    check_primer_cores = check_primer_cores,
    check_primer_pairs = check_pairs,
    genome_sequences = genome_sequences
  ))
}

.shared_read_genome_cached <- function(genbank_file) {
  cache_root <- tools::R_user_dir("PrimerDesigner", which = "cache")
  cache_dir <- base::file.path(cache_root, "parsed_genomes")
  if (!base::dir.exists(cache_dir)) {
    base::dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  fi <- base::file.info(genbank_file)
  ext <- base::tolower(tools::file_ext(genbank_file))
  key_raw <- base::paste(
    normalizePath(genbank_file, winslash = "/", mustWork = FALSE),
    ext,
    fi$size,
    as.numeric(fi$mtime),
    sep = "__"
  )
  key <- gsub("[^A-Za-z0-9._-]", "_", key_raw)
  if (base::nchar(key) > 180) key <- base::substring(key, 1, 180)
  cache_file <- base::file.path(cache_dir, base::paste0(key, ".rds"))

  if (base::file.exists(cache_file)) {
    return(base::readRDS(cache_file))
  }

  if (ext == "dna") {
    gb_lines <- .get_genbank_from_dna(genbank_file, kill_snapgene = FALSE)
    tmp_gbk <- tempfile(fileext = ".gbk")
    base::on.exit(if (base::file.exists(tmp_gbk)) base::unlink(tmp_gbk), add = TRUE)
    base::writeLines(gb_lines, tmp_gbk)
    parsed <- read_genome_genbank(tmp_gbk)
  } else {
    parsed <- read_genome_genbank(genbank_file)
  }
  base::saveRDS(parsed, cache_file)
  parsed
}

.shared_extract_target_context <- function(
    genome_id,
    locus_tag,
    genome_seq,
    genbank_table,
    preferred_upstream_bp,
    preferred_downstream_bp,
    max_arm_bp,
    min_arm_bp,
    preserve_upstream_gene,
    preserve_upstream_buffer_bp = 0L,
    auto_preserve_gap_bp = 100L,
    overlap_policy = c("strict", "flexible")
) {
  overlap_policy <- base::match.arg(overlap_policy)
  cds_tbl <- genbank_table[!base::is.na(genbank_table$locus_tag), , drop = FALSE]
  target_row <- cds_tbl[cds_tbl$locus_tag == locus_tag, , drop = FALSE]
  if (base::nrow(target_row) != 1) {
    base::stop("Expected exactly one target row for ", genome_id, " / ", locus_tag)
  }

  target_start <- base::as.integer(target_row$start[1])
  target_end <- base::as.integer(target_row$end[1])
  target_gene <- if ("gene" %in% base::colnames(target_row)) base::as.character(target_row$gene[1]) else NA_character_
  target_product <- if ("product" %in% base::colnames(target_row)) base::as.character(target_row$product[1]) else NA_character_
  target_nt_seq <- if ("nt_seq" %in% base::colnames(target_row)) {
    base::toupper(base::as.character(target_row$nt_seq[1]))
  } else {
    .shared_subseq_circular(genome_seq, target_start, target_end)
  }
  target_strand <- if ("strand" %in% base::colnames(target_row)) {
    base::as.character(target_row$strand[1])
  } else {
    "+"
  }

  # Find the CDS immediately upstream of the target (regardless of overlap)
  # so that downstream reports (upstream_locus_tag / upstream_gene) stay
  # populated for annotation purposes.
  upstream_candidates <- cds_tbl[
    cds_tbl$end <= target_end & cds_tbl$locus_tag != locus_tag,
    , drop = FALSE
  ]
  upstream_candidates <- upstream_candidates[order(upstream_candidates$end, decreasing = TRUE), , drop = FALSE]
  if (base::nrow(upstream_candidates) == 0) {
    base::stop("No upstream CDS found for ", genome_id, " / ", locus_tag)
  }
  upstream_row <- upstream_candidates[1, , drop = FALSE]
  upstream_end_val <- base::as.integer(upstream_row$end[1])

  # ----- Deletion window driven by overlap policy -----
  # "strict"   : deletion must never clip another CDS. If an upstream gene
  #              extends INTO the target, we shrink the deletion to start
  #              AFTER that gene ends; similarly for downstream genes.
  # "flexible" : deletion spans the full target coordinates regardless of
  #              overlapping CDS (user accepts possible clipping of neighbor
  #              genes to preserve the natural target boundary).
  overlapping <- cds_tbl[
    cds_tbl$locus_tag != locus_tag &
      cds_tbl$end >= target_start &
      cds_tbl$start <= target_end,
    , drop = FALSE
  ]
  upstream_overlap_bp <- 0L
  downstream_overlap_bp <- 0L
  if (base::nrow(overlapping) > 0) {
    up_ov <- overlapping[overlapping$start < target_start, , drop = FALSE]
    if (base::nrow(up_ov) > 0) {
      max_up_end <- base::min(base::max(base::as.integer(up_ov$end)), target_end)
      upstream_overlap_bp <- max_up_end - target_start + 1L
    }
    dn_ov <- overlapping[overlapping$end > target_end, , drop = FALSE]
    if (base::nrow(dn_ov) > 0) {
      min_dn_start <- base::max(base::min(base::as.integer(dn_ov$start)), target_start)
      downstream_overlap_bp <- target_end - min_dn_start + 1L
    }
    nested <- overlapping[overlapping$start >= target_start &
                            overlapping$end <= target_end, , drop = FALSE]
    if (base::nrow(nested) > 0) {
      base::warning(sprintf(
        "[overlap] %s / %s: %d CDS nested inside target; their sequence cannot be preserved by a 4-primer Gibson deletion.",
        genome_id, locus_tag, base::nrow(nested)), call. = FALSE)
    }
  }

  if (overlap_policy == "strict") {
    deletion_start <- target_start + base::as.integer(upstream_overlap_bp)
    deletion_end <- target_end - base::as.integer(downstream_overlap_bp)
    if (upstream_overlap_bp > 0L || downstream_overlap_bp > 0L) {
      base::message(sprintf(
        "[overlap/strict] %s / %s: preserving %d bp upstream and %d bp downstream neighbor-gene overlap (deletion trimmed to %d-%d).",
        genome_id, locus_tag, upstream_overlap_bp, downstream_overlap_bp,
        deletion_start, deletion_end))
    }
  } else {
    # flexible: keep the natural target boundary
    deletion_start <- target_start
    deletion_end <- target_end
    if (upstream_overlap_bp > 0L || downstream_overlap_bp > 0L) {
      base::warning(sprintf(
        "[overlap/flexible] %s / %s: deletion will clip %d bp of an upstream neighbor and %d bp of a downstream neighbor.",
        genome_id, locus_tag, upstream_overlap_bp, downstream_overlap_bp), call. = FALSE)
    }
  }

  if (deletion_end < deletion_start) {
    base::stop(sprintf(
      "%s / %s: overlapping CDS consume the entire target (up=%d, dn=%d, target=%d-%d). Nothing left to delete.",
      genome_id, locus_tag, upstream_overlap_bp, downstream_overlap_bp,
      target_start, target_end))
  }

  # Back-compat knob: `preserve_upstream_gene = TRUE` still lets the user
  # push the deletion even further in the 5' direction so it starts right at
  # the end of the upstream CDS (instead of the target CDS start). Rarely
  # needed once overlap_policy handles the common case.
  force_preserve_upstream <- if (base::is.character(preserve_upstream_gene) &&
                                 base::identical(base::tolower(preserve_upstream_gene), "auto")) {
    FALSE
  } else {
    base::isTRUE(preserve_upstream_gene)
  }
  if (force_preserve_upstream) {
    deletion_start <- base::max(deletion_start,
                                upstream_end_val + base::as.integer(preserve_upstream_buffer_bp) + 1L)
  }

  left_arm_end <- deletion_start - 1L
  right_arm_start <- deletion_end + 1L

  # Extract extended arms (max_arm_bp) to give the shared-primer search room
  # to flex; effective arm bounds are enforced downstream in candidate
  # filtering.
  left_arm_start <- left_arm_end - base::as.integer(max_arm_bp) + 1L
  right_arm_end <- right_arm_start + base::as.integer(max_arm_bp) - 1L

  genome_len <- base::nchar(genome_seq)
  left_arm_seq <- .shared_subseq_circular(genome_seq, left_arm_start, left_arm_end)
  right_arm_seq <- .shared_subseq_circular(genome_seq, right_arm_start, right_arm_end)

  context <- base::data.frame(
    genome_id = genome_id,
    locus_tag = locus_tag,
    target_start = target_start,
    target_end = target_end,
    target_strand = target_strand,
    target_gene = target_gene,
    target_product = target_product,
    target_nt_seq = target_nt_seq,
    upstream_locus_tag = base::as.character(upstream_row$locus_tag[1]),
    upstream_gene = if ("gene" %in% base::colnames(upstream_row)) base::as.character(upstream_row$gene[1]) else NA_character_,
    upstream_product = if ("product" %in% base::colnames(upstream_row)) base::as.character(upstream_row$product[1]) else NA_character_,
    upstream_start = base::as.integer(upstream_row$start[1]),
    upstream_end = upstream_end_val,
    upstream_overlap_bp = base::as.integer(upstream_overlap_bp),
    downstream_overlap_bp = base::as.integer(downstream_overlap_bp),
    overlap_policy = overlap_policy,
    preferred_upstream_bp = base::as.integer(preferred_upstream_bp),
    preferred_downstream_bp = base::as.integer(preferred_downstream_bp),
    min_arm_bp = base::as.integer(min_arm_bp),
    max_arm_bp = base::as.integer(max_arm_bp),
    extracted_arm_bp = base::as.integer(max_arm_bp),
    left_arm_start = left_arm_start,
    left_arm_end = left_arm_end,
    right_arm_start = right_arm_start,
    right_arm_end = right_arm_end,
    deletion_start = deletion_start,
    deletion_end = deletion_end,
    genome_length = genome_len,
    left_arm_seq = left_arm_seq,
    right_arm_seq = right_arm_seq,
    insert_seq = base::paste0(left_arm_seq, right_arm_seq),
    stringsAsFactors = FALSE
  )

  return(base::list(
    genome_id = genome_id,
    genome_seq = genome_seq,
    context = context
  ))
}

.shared_build_construct_groups <- function(target_context) {
  seq_col <- if ("effective_insert_seq" %in% base::colnames(target_context)) {
    "effective_insert_seq"
  } else {
    "insert_seq"
  }
  seq_keys <- target_context[[seq_col]]
  split_ids <- base::split(target_context$genome_id, seq_keys)
  split_indices <- base::split(base::seq_len(base::nrow(target_context)), seq_keys)

  groups <- base::lapply(base::seq_along(split_ids), function(i) {
    idx <- split_indices[[i]]
    members <- target_context$genome_id[idx]
    label <- if (base::length(members) == 1L) {
      members[1]
    } else {
      base::paste0(base::paste(members, collapse = "_"), "_common")
    }
    base::data.frame(
      construct_group_id = base::paste0("construct_group_", i),
      construct_label = label,
      used_by = base::paste(members, collapse = ", "),
      n_genomes = base::length(members),
      representative_genome = members[1],
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(groups)
}

.shared_write_construct_groups <- function(
    construct_groups,
    target_context,
    vector_file,
    start,
    end,
    fallback_upstream_bp,
    fallback_downstream_bp,
    construct_output_dir
) {
  for (i in base::seq_len(base::nrow(construct_groups))) {
    rep_id <- construct_groups$representative_genome[i]
    label <- construct_groups$construct_label[i]
    ctx <- target_context[target_context$genome_id == rep_id, , drop = FALSE]
    out_path <- base::file.path(
      construct_output_dir,
      base::paste0(label, "_deletion_construct.gbk")
    )
    write_deletion_genbank(
      vector_file = vector_file,
      insert_seq = if ("effective_insert_seq" %in% base::colnames(ctx)) ctx$effective_insert_seq[1] else ctx$insert_seq[1],
      start = start,
      end = end,
      locus_tag = label,
      primers = NULL,
      upstream_bp = if ("effective_upstream_bp" %in% base::colnames(ctx)) ctx$effective_upstream_bp[1] else fallback_upstream_bp,
      downstream_bp = if ("effective_downstream_bp" %in% base::colnames(ctx)) ctx$effective_downstream_bp[1] else fallback_downstream_bp,
      output_path = out_path,
      kill_snapgene = FALSE
    )
  }
}

.shared_apply_primer_geometry <- function(target_context, primer_cores, min_arm_bp = NULL, max_arm_bp = NULL) {
  if (base::nrow(target_context) == 0 || base::nrow(primer_cores) == 0) return(target_context)

  # preferred_anchor_fn maps (arm_len, preferred_bp, k) to the ideal k-mer
  # position (start-coordinate) inside the arm so that the chosen occurrence
  # matches the position this candidate was scored against during pool
  # construction. Without this, .shared_locate_in_arm would blindly pick the
  # earliest/latest arm match and could place UF beyond UR (invalid geometry).
  role_cfg <- base::list(
    UF = base::list(arm = "left_arm_seq", region = "start", reverse = FALSE,
                    preferred_anchor = function(arm_len, pref, k) arm_len - pref + 1L),
    UR = base::list(arm = "left_arm_seq", region = "end", reverse = TRUE,
                    preferred_anchor = function(arm_len, pref, k) arm_len - k + 1L),
    DF = base::list(arm = "right_arm_seq", region = "start", reverse = FALSE,
                    preferred_anchor = function(arm_len, pref, k) 1L),
    DR = base::list(arm = "right_arm_seq", region = "end", reverse = TRUE,
                    preferred_anchor = function(arm_len, pref, k) pref - k + 1L)
  )
  role_preferred_bp_field <- base::list(
    UF = "preferred_upstream_bp", UR = "preferred_upstream_bp",
    DF = "preferred_downstream_bp", DR = "preferred_downstream_bp"
  )

  out_rows <- base::vector("list", base::nrow(target_context))
  for (i in base::seq_len(base::nrow(target_context))) {
    row <- target_context[i, , drop = FALSE]
    genome_id <- base::as.character(row$genome_id[1])
    role_hits <- lapply(base::names(role_cfg), function(role) {
      cand <- primer_cores[primer_cores$role == role, , drop = FALSE]
      keep <- vapply(strsplit(cand$used_by, ", ", fixed = TRUE), function(x) genome_id %in% x, logical(1))
      cand <- cand[keep, , drop = FALSE]
      if (base::nrow(cand) != 1L) {
        base::stop("Expected exactly one ", role, " primer assignment for genome ", genome_id)
      }
      cand
    })
    base::names(role_hits) <- base::names(role_cfg)

    locs <- lapply(base::names(role_cfg), function(role) {
      cfg <- role_cfg[[role]]
      arm_seq <- base::as.character(row[[cfg$arm]][1])
      cand <- role_hits[[role]]
      target_seq <- if ("shared_core" %in% base::colnames(cand) && !base::is.na(cand$shared_core[1]) && base::nchar(cand$shared_core[1]) > 0) {
        base::as.character(cand$shared_core[1])
      } else if (cfg$reverse) {
        .rc(base::as.character(cand$primer_core_5to3[1]))
      } else {
        base::as.character(cand$primer_core_5to3[1])
      }
      arm_len_i <- base::nchar(arm_seq)
      k_i <- base::nchar(target_seq)
      pref_bp <- if (role_preferred_bp_field[[role]] %in% base::colnames(row)) {
        base::as.integer(row[[role_preferred_bp_field[[role]]]][1])
      } else {
        arm_len_i
      }
      anchor_start <- base::as.integer(cfg$preferred_anchor(arm_len_i, pref_bp, k_i))
      .shared_locate_in_arm(
        arm_seq = arm_seq,
        target_seq = target_seq,
        region = cfg$region,
        preferred_start = anchor_start
      )
    })
    base::names(locs) <- base::names(role_cfg)

    left_start_off <- locs$UF$start
    left_end_off <- locs$UR$end
    right_start_off <- locs$DF$start
    right_end_off <- locs$DR$end

    # Geometry repair: when UF/UR (or DF/DR) are out of order for this genome,
    # re-search the entire arm for both shared cores and try to pick an
    # occurrence pair that yields a valid (UF.start < UR.end) frame with
    # effective arm in [min_arm_bp, max_arm_bp].
    fix_pair <- function(arm_seq, outer_seq, inner_seq, arm_len_i,
                          strict_min = TRUE) {
      outer_hits <- base::unlist(gregexpr(outer_seq, arm_seq, fixed = TRUE), use.names = FALSE)
      inner_hits <- base::unlist(gregexpr(inner_seq, arm_seq, fixed = TRUE), use.names = FALSE)
      outer_hits <- outer_hits[outer_hits > 0]
      inner_hits <- inner_hits[inner_hits > 0]
      if (base::length(outer_hits) == 0 || base::length(inner_hits) == 0) {
        return(NULL)
      }
      best <- NULL
      best_dev <- Inf
      pref_arm <- if (!base::is.null(min_arm_bp) && !base::is.null(max_arm_bp)) {
        base::as.integer((min_arm_bp + max_arm_bp) / 2L)
      } else {
        base::as.integer(arm_len_i / 2L)
      }
      for (o in outer_hits) {
        for (ih in inner_hits) {
          inner_end <- ih + base::nchar(inner_seq) - 1L
          if (inner_end <= o) next
          eff_arm <- inner_end - o + 1L
          if (strict_min && !base::is.null(min_arm_bp) && eff_arm < min_arm_bp) next
          if (!base::is.null(max_arm_bp) && eff_arm > max_arm_bp) next
          dev <- base::abs(eff_arm - pref_arm)
          if (dev < best_dev) {
            best <- base::list(outer_start = o, inner_start = ih, inner_end = inner_end,
                                eff_arm = eff_arm)
            best_dev <- dev
          }
        }
      }
      best
    }
    # Two-stage: first try strict min_arm_bp; if that yields nothing, relax the
    # minimum to pick the largest achievable arm (so the design still produces
    # a result, with a loud warning).
    fix_pair_or_best <- function(arm_seq, outer_seq, inner_seq, arm_len_i) {
      r <- fix_pair(arm_seq, outer_seq, inner_seq, arm_len_i, strict_min = TRUE)
      if (!base::is.null(r)) return(r)
      fix_pair(arm_seq, outer_seq, inner_seq, arm_len_i, strict_min = FALSE)
    }

    # Trigger repair when either anchors are out of order OR the resulting
    # effective arm is below min_arm_bp. When primers have multiple binding
    # sites in the arm (e.g. transposon repeats), fix_pair picks the
    # occurrence pair whose effective arm stays within [min_arm_bp, max_arm_bp]
    # and is closest to the preferred arm length.
    left_eff_bp <- if (!base::is.null(min_arm_bp) &&
                      left_start_off <= left_end_off)
      left_end_off - left_start_off + 1L else NA_integer_
    trigger_left <- left_start_off > left_end_off ||
      (!base::is.null(min_arm_bp) && !base::is.na(left_eff_bp) &&
       left_eff_bp < base::as.integer(min_arm_bp))
    if (trigger_left) {
      base::warning(sprintf(
        "[geom] %s left-arm anchors %s (UF.start=%d, UR.end=%d%s); searching full arm for a valid pair.",
        genome_id,
        if (left_start_off > left_end_off) "out of order" else "give arm below min_arm_bp",
        left_start_off, left_end_off,
        if (!base::is.na(left_eff_bp)) base::sprintf(", eff=%dbp", left_eff_bp) else ""),
        call. = FALSE)
      uf_seq <- base::as.character(role_hits$UF$shared_core[1])
      if (base::is.na(uf_seq) || base::nchar(uf_seq) == 0) uf_seq <- base::as.character(role_hits$UF$primer_core_5to3[1])
      ur_seq <- base::as.character(role_hits$UR$shared_core[1])
      if (base::is.na(ur_seq) || base::nchar(ur_seq) == 0) ur_seq <- .rc(base::as.character(role_hits$UR$primer_core_5to3[1]))
      arm_seq_l <- base::as.character(row$left_arm_seq[1])
      fix <- fix_pair_or_best(arm_seq_l, uf_seq, ur_seq, base::nchar(arm_seq_l))
      if (base::is.null(fix)) {
        base::warning(sprintf(
          "[geom] %s: no valid UF/UR placement in left arm; skipping this genome (strain-specific design recommended).",
          genome_id), call. = FALSE)
        out_rows[[i]] <- NULL
        next
      }
      if (!base::is.null(min_arm_bp) && fix$eff_arm < base::as.integer(min_arm_bp)) {
        base::warning(sprintf(
          "[geom] %s left arm: best achievable effective arm is %d bp (< min_arm_bp=%s). Kept for design but verify manually.",
          genome_id, fix$eff_arm, base::format(min_arm_bp)), call. = FALSE)
      }
      left_start_off <- fix$outer_start
      left_end_off <- fix$inner_end
      locs$UF$start <- fix$outer_start
      locs$UF$end <- fix$outer_start + base::nchar(uf_seq) - 1L
      locs$UR$start <- fix$inner_start
      locs$UR$end <- fix$inner_end
    }

    right_eff_bp <- if (!base::is.null(min_arm_bp) &&
                       right_start_off <= right_end_off)
      right_end_off - right_start_off + 1L else NA_integer_
    trigger_right <- right_start_off > right_end_off ||
      (!base::is.null(min_arm_bp) && !base::is.na(right_eff_bp) &&
       right_eff_bp < base::as.integer(min_arm_bp))
    if (trigger_right) {
      base::warning(sprintf(
        "[geom] %s right-arm anchors %s (DF.start=%d, DR.end=%d%s); searching full arm for a valid pair.",
        genome_id,
        if (right_start_off > right_end_off) "out of order" else "give arm below min_arm_bp",
        right_start_off, right_end_off,
        if (!base::is.na(right_eff_bp)) base::sprintf(", eff=%dbp", right_eff_bp) else ""),
        call. = FALSE)
      df_seq <- base::as.character(role_hits$DF$shared_core[1])
      if (base::is.na(df_seq) || base::nchar(df_seq) == 0) df_seq <- base::as.character(role_hits$DF$primer_core_5to3[1])
      dr_seq <- base::as.character(role_hits$DR$shared_core[1])
      if (base::is.na(dr_seq) || base::nchar(dr_seq) == 0) dr_seq <- .rc(base::as.character(role_hits$DR$primer_core_5to3[1]))
      arm_seq_r <- base::as.character(row$right_arm_seq[1])
      fix <- fix_pair_or_best(arm_seq_r, df_seq, dr_seq, base::nchar(arm_seq_r))
      if (base::is.null(fix)) {
        base::warning(sprintf(
          "[geom] %s: no valid DF/DR placement in right arm; skipping this genome.",
          genome_id), call. = FALSE)
        out_rows[[i]] <- NULL
        next
      }
      if (!base::is.null(min_arm_bp) && fix$eff_arm < base::as.integer(min_arm_bp)) {
        base::warning(sprintf(
          "[geom] %s right arm: best achievable effective arm is %d bp (< min_arm_bp=%s). Kept for design but verify manually.",
          genome_id, fix$eff_arm, base::format(min_arm_bp)), call. = FALSE)
      }
      right_start_off <- fix$outer_start
      right_end_off <- fix$inner_end
      locs$DF$start <- fix$outer_start
      locs$DF$end <- fix$outer_start + base::nchar(df_seq) - 1L
      locs$DR$start <- fix$inner_start
      locs$DR$end <- fix$inner_end
    }

    left_eff <- base::substring(base::as.character(row$left_arm_seq[1]), left_start_off, left_end_off)
    right_eff <- base::substring(base::as.character(row$right_arm_seq[1]), right_start_off, right_end_off)
    genome_len <- base::as.integer(row$genome_length[1])

    row$UF_cluster_id <- role_hits$UF$cluster_id[1]
    row$UR_cluster_id <- role_hits$UR$cluster_id[1]
    row$DF_cluster_id <- role_hits$DF$cluster_id[1]
    row$DR_cluster_id <- role_hits$DR$cluster_id[1]

    row$effective_left_arm_start <- .shared_wrap_coord(base::as.integer(row$left_arm_start[1]) + left_start_off - 1L, genome_len)
    row$effective_left_arm_end <- .shared_wrap_coord(base::as.integer(row$left_arm_start[1]) + left_end_off - 1L, genome_len)
    row$effective_right_arm_start <- .shared_wrap_coord(base::as.integer(row$right_arm_start[1]) + right_start_off - 1L, genome_len)
    row$effective_right_arm_end <- .shared_wrap_coord(base::as.integer(row$right_arm_start[1]) + right_end_off - 1L, genome_len)
    row$effective_deletion_start <- .shared_wrap_coord(base::as.integer(row$effective_left_arm_end[1]) + 1L, genome_len)
    row$effective_deletion_end <- .shared_wrap_coord(base::as.integer(row$effective_right_arm_start[1]) - 1L, genome_len)
    row$effective_upstream_bp <- base::nchar(left_eff)
    row$effective_downstream_bp <- base::nchar(right_eff)
    row$effective_left_arm_seq <- left_eff
    row$effective_right_arm_seq <- right_eff
    row$effective_insert_seq <- base::paste0(left_eff, right_eff)
    row$uf_offset_bp <- left_start_off - 1L
    row$ur_offset_bp <- base::as.integer(row$left_arm_end[1]) - base::as.integer(row$effective_left_arm_end[1])
    row$df_offset_bp <- right_start_off - 1L
    row$dr_offset_bp <- base::as.integer(row$right_arm_end[1]) - base::as.integer(row$effective_right_arm_end[1])

    pref_up <- if ("preferred_upstream_bp" %in% base::colnames(row)) base::as.integer(row$preferred_upstream_bp[1]) else NA_integer_
    pref_dn <- if ("preferred_downstream_bp" %in% base::colnames(row)) base::as.integer(row$preferred_downstream_bp[1]) else NA_integer_
    row$upstream_arm_deviation_bp <- if (!base::is.na(pref_up)) base::abs(base::nchar(left_eff) - pref_up) else NA_integer_
    row$downstream_arm_deviation_bp <- if (!base::is.na(pref_dn)) base::abs(base::nchar(right_eff) - pref_dn) else NA_integer_

    if (!base::is.null(min_arm_bp) || !base::is.null(max_arm_bp)) {
      lo <- if (!base::is.null(min_arm_bp)) base::as.integer(min_arm_bp) else -Inf
      hi <- if (!base::is.null(max_arm_bp)) base::as.integer(max_arm_bp) else Inf
      if (base::nchar(left_eff) < lo || base::nchar(left_eff) > hi) {
        base::warning(sprintf(
          "Effective left arm (%d bp) for %s is outside [%s, %s] bp.",
          base::nchar(left_eff), genome_id,
          base::format(lo), base::format(hi)
        ), call. = FALSE)
      }
      if (base::nchar(right_eff) < lo || base::nchar(right_eff) > hi) {
        base::warning(sprintf(
          "Effective right arm (%d bp) for %s is outside [%s, %s] bp.",
          base::nchar(right_eff), genome_id,
          base::format(lo), base::format(hi)
        ), call. = FALSE)
      }
    }
    out_rows[[i]] <- row
  }
  out_rows <- Filter(base::Negate(base::is.null), out_rows)
  dplyr::bind_rows(out_rows)
}

.shared_locate_in_arm <- function(arm_seq, target_seq,
                                  region = c("start", "end"),
                                  preferred_start = NULL) {
  region <- base::match.arg(region)
  starts <- unlist(gregexpr(target_seq, arm_seq, fixed = TRUE), use.names = FALSE)
  starts <- starts[starts > 0]
  if (base::length(starts) == 0) {
    base::stop("Target sequence not found in arm: ", target_seq)
  }
  # When a preferred start position is supplied (computed from the role's
  # scored position), pick the occurrence closest to it. This prevents the
  # earliest/latest heuristic from placing a shared oligo at an arm location
  # different from the one the scorer picked, which can make UF.start > UR.end
  # and break Gibson geometry.
  start_pos <- if (!base::is.null(preferred_start)) {
    starts[base::which.min(base::abs(starts - base::as.integer(preferred_start)))]
  } else if (region == "start") {
    base::min(starts)
  } else {
    starts[base::which.max(starts + base::nchar(target_seq) - 1L)]
  }
  base::list(start = start_pos, end = start_pos + base::nchar(target_seq) - 1L)
}

.shared_wrap_coord <- function(coord, genome_len) {
  if (coord < 1L) return(genome_len + coord)
  if (coord > genome_len) return(coord - genome_len)
  coord
}

.shared_subseq_circular <- function(seq, start, end) {
  seq <- base::toupper(seq)
  n <- base::nchar(seq)
  if (start < 1L) {
    start_wrap <- n + start + 1L
    return(base::paste0(
      base::substring(seq, start_wrap, n),
      base::substring(seq, 1L, end)
    ))
  }
  if (end > n) {
    end_wrap <- end - n
    return(base::paste0(
      base::substring(seq, start, n),
      base::substring(seq, 1L, end_wrap)
    ))
  }
  return(base::substring(seq, start, end))
}

.shared_assign_role_groups <- function(
    genome_ids,
    context_map,
    dnastr_cache,
    role,
    arm_field,
    region,
    reverse,
    boundary,
    preferred_bp_field,
    min_arm_bp,
    max_arm_bp,
    search_window,
    min_primer_length,
    max_primer_length,
    tm_target,
    require_unique,
    hit_count_cache,
    forbidden_oligos = character(0),
    inner_boundary_tolerance_bp = 50L,
    overlap_policy = c("strict", "flexible"),
    tm_window = NULL
) {
  overlap_policy <- base::match.arg(overlap_policy)
  # Under strict, inner primers (UR, DF) must abut the deletion boundary
  # exactly (0 bp tolerance). Any non-zero value the caller supplied is
  # honoured only in flexible mode.
  if (overlap_policy == "strict" && boundary == "inner") {
    inner_boundary_tolerance_bp <- 0L
  }
  # Single strict-uniqueness pool with a generous cap (covers most roles in
  # one pass). Building the pool + Biostrings uniqueness dominates runtime,
  # so re-building several times for escalating caps was wasteful.
  candidate_pool <- .shared_build_primer_candidate_pool(
    genome_ids = genome_ids,
    context_map = context_map,
    dnastr_cache = dnastr_cache,
    role = role,
    arm_field = arm_field,
    region = region,
    reverse = reverse,
    boundary = boundary,
    preferred_bp_field = preferred_bp_field,
    min_arm_bp = min_arm_bp,
    max_arm_bp = max_arm_bp,
    search_window = search_window,
    min_primer_length = min_primer_length,
    max_primer_length = max_primer_length,
    tm_target = tm_target,
    require_unique = require_unique,
    hit_count_cache = hit_count_cache,
    candidate_cap = 150L,
    forbidden_oligos = forbidden_oligos,
    inner_boundary_tolerance_bp = inner_boundary_tolerance_bp,
    tm_window = tm_window
  )
  attempt <- base::tryCatch(
    .shared_select_primer_groups_from_pool(candidate_pool, genome_ids, role, reverse),
    error = function(e) e
  )
  if (!inherits(attempt, "error")) return(attempt)

  if (require_unique) {
    base::warning(sprintf(
      "Role %s: no unique primer covers every target genome in the search window; retrying with uniqueness relaxed for the residual genomes.",
      role), call. = FALSE)
    candidate_pool <- .shared_build_primer_candidate_pool(
      genome_ids = genome_ids,
      context_map = context_map,
      dnastr_cache = dnastr_cache,
      role = role,
      arm_field = arm_field,
      region = region,
      reverse = reverse,
      boundary = boundary,
      preferred_bp_field = preferred_bp_field,
      min_arm_bp = min_arm_bp,
      max_arm_bp = max_arm_bp,
      search_window = search_window,
      min_primer_length = min_primer_length,
      max_primer_length = max_primer_length,
      tm_target = tm_target,
      require_unique = FALSE,
      hit_count_cache = hit_count_cache,
      candidate_cap = 300L,
      forbidden_oligos = forbidden_oligos,
      inner_boundary_tolerance_bp = inner_boundary_tolerance_bp,
      tm_window = tm_window
    )
    attempt <- base::tryCatch(
      .shared_select_primer_groups_from_pool(candidate_pool, genome_ids, role, reverse),
      error = function(e) e
    )
    if (!inherits(attempt, "error")) return(attempt)
  }
  # Progressive fallback for inner roles (UR/DF): widen the inner-boundary
  # pin step by step. Many genomes at once may have no conserved sequence
  # close enough to the deletion edge to share a primer; widening lets a
  # candidate survive, which the subsequent selector may still subgroup-split.
  # Under overlap_policy == "strict", widening is FORBIDDEN: inner primers
  # must sit at the deletion boundary. If no shared candidate exists at the
  # boundary, the selector must subgroup-split (per-strain inner primers) or
  # the run must fail — never relax the pin.
  if (overlap_policy != "strict" &&
      boundary == "inner" &&
      !base::is.null(inner_boundary_tolerance_bp) &&
      base::is.finite(inner_boundary_tolerance_bp)) {
    current <- base::as.integer(inner_boundary_tolerance_bp)
    # Try 4x, then 16x, then unlimited (NULL).
    for (widened in base::list(current * 4L, current * 16L, NULL)) {
      base::warning(sprintf(
        "Role %s: widening inner_boundary_tolerance_bp from %d to %s as a last-resort to find any valid candidate.",
        role, current,
        if (base::is.null(widened)) "unlimited" else base::as.character(widened)),
        call. = FALSE)
      candidate_pool <- .shared_build_primer_candidate_pool(
        genome_ids = genome_ids,
        context_map = context_map,
        dnastr_cache = dnastr_cache,
        role = role,
        arm_field = arm_field,
        region = region,
        reverse = reverse,
        boundary = boundary,
        preferred_bp_field = preferred_bp_field,
        min_arm_bp = min_arm_bp,
        max_arm_bp = max_arm_bp,
        search_window = search_window,
        min_primer_length = min_primer_length,
        max_primer_length = max_primer_length,
        tm_target = tm_target,
        require_unique = FALSE,
        hit_count_cache = hit_count_cache,
        candidate_cap = 300L,
        forbidden_oligos = forbidden_oligos,
        inner_boundary_tolerance_bp = widened,
        tm_window = tm_window
      )
      attempt <- base::tryCatch(
        .shared_select_primer_groups_from_pool(candidate_pool, genome_ids, role, reverse),
        error = function(e) e
      )
      if (!inherits(attempt, "error")) return(attempt)
      if (!base::is.null(widened)) current <- base::as.integer(widened)
    }
  }
  if (overlap_policy == "strict" && boundary == "inner") {
    base::stop(sprintf(
      "Role %s: no primer candidate touches the deletion boundary (0 bp tolerance) for every genome, even with subgroup-splitting. Under overlap_policy='strict' the inner primer MUST abut the deletion boundary exactly. Options: (a) remove the offending genome(s) from the target_table, (b) switch to overlap_policy='flexible', or (c) widen primer_min_length/primer_max_length so some primer fits at the boundary. Upstream error: %s",
      role, base::conditionMessage(attempt)
    ), call. = FALSE)
  }
  base::stop(base::conditionMessage(attempt))
}

# Implied effective arm length if the OTHER role in the arm sits at its
# boundary of maximum flexibility — this is the permissive upper bound used
# for filtering with `[min_arm_bp, max_arm_bp]`. Separate from the scoring
# term (`.shared_role_arm_dev`) which biases toward the preferred position.
.shared_role_arm_length <- function(role, pos_start, k, arm_len, preferred_bp) {
  if (role == "UF") {
    arm_len - pos_start + 1L            # UR at arm_len
  } else if (role == "UR") {
    pos_start + k - 1L                  # UF at 1
  } else if (role == "DF") {
    arm_len - pos_start + 1L            # DR at arm_len
  } else if (role == "DR") {
    pos_start + k - 1L                  # DF at 1
  } else if (role == "cF") {
    # cF: colony-PCR forward primer, sits on the upstream-CHECK window
    # (further upstream of the UF homology arm). We want it close to the
    # right edge of the check window (i.e. near the upstream arm).
    pos_start + k - 1L
  } else if (role == "cR") {
    # cR: colony-PCR reverse primer, sits on the downstream-CHECK window
    # close to the left edge (i.e. near the downstream arm).
    arm_len - pos_start + 1L
  } else {
    NA_integer_
  }
}

# Scoring deviation from the IDEAL boundary for this role:
#   - Outer roles (UF, DR): we want effective arm ≈ preferred_bp, so the
#     deviation is |implied arm (paired at boundary) - preferred_bp|.
#   - Inner roles (UR, DF): we want the k-mer right up against the target,
#     so the deviation is simply distance from the inner boundary
#     (i.e. bp lost relative to `preferred_bp` when paired role is at
#     its preferred point).
.shared_role_arm_dev <- function(role, pos_start, k, arm_len, preferred_bp) {
  preferred_bp <- base::as.integer(preferred_bp)
  if (role == "UF") {
    base::abs((arm_len - pos_start + 1L) - preferred_bp)
  } else if (role == "UR") {
    arm_len - (pos_start + k - 1L)
  } else if (role == "DF") {
    pos_start - 1L
  } else if (role == "DR") {
    base::abs((pos_start + k - 1L) - preferred_bp)
  } else if (role == "cF") {
    # Deviation = distance from the right edge of the check window (how far
    # the primer's 3' end lies from being adjacent to the homology arm).
    arm_len - (pos_start + k - 1L)
  } else if (role == "cR") {
    # Deviation = distance from the left edge of the check window.
    pos_start - 1L
  } else {
    NA_integer_
  }
}

.shared_build_primer_candidate_pool <- function(
    genome_ids,
    context_map,
    dnastr_cache = NULL,
    role,
    arm_field,
    region,
    reverse,
    boundary,
    preferred_bp_field,
    min_arm_bp,
    max_arm_bp,
    search_window,
    min_primer_length,
    max_primer_length,
    tm_target,
    require_unique,
    hit_count_cache,
    candidate_cap = 20L,
    forbidden_oligos = character(0),
    inner_boundary_tolerance_bp = 50L,
    seed_prefilter_enabled = FALSE,
    tm_window = NULL
) {
  # --- Vectorized k-mer scan -------------------------------------------------
  # Build a long data.frame of (genome, oligo, position, k, implied_arm, ...)
  # in one shot per genome, using R's vectorized substring(). This avoids the
  # previous per-position env lookups (slow for 10^5 positions).
  per_genome_frames <- base::vector("list", base::length(genome_ids))
  for (gi in base::seq_along(genome_ids)) {
    id <- genome_ids[gi]
    ctx <- context_map[[id]]$context
    arm_seq <- base::as.character(ctx[[arm_field]][1])
    arm_len <- base::nchar(arm_seq)
    preferred_bp <- base::as.integer(ctx[[preferred_bp_field]][1])

    if (base::is.null(search_window)) {
      scan_start <- 1L
      scan_end <- arm_len
    } else {
      sw <- base::as.integer(search_window)
      if (boundary == "outer" && region == "start") {
        pref_pos <- arm_len - preferred_bp + 1L
        scan_start <- base::max(1L, pref_pos - sw %/% 2L)
        scan_end <- base::min(arm_len, pref_pos + sw %/% 2L)
      } else if (boundary == "outer" && region == "end") {
        pref_end <- preferred_bp
        scan_start <- base::max(1L, pref_end - sw %/% 2L)
        scan_end <- base::min(arm_len, pref_end + sw %/% 2L)
      } else if (region == "start") {
        scan_start <- 1L
        scan_end <- base::min(arm_len, sw)
      } else {
        scan_start <- base::max(1L, arm_len - sw + 1L)
        scan_end <- arm_len
      }
    }

    parts <- base::list()
    for (k in base::seq.int(min_primer_length, max_primer_length)) {
      if (arm_len < k) next
      eff_end <- base::min(scan_end, arm_len - k + 1L)
      if (eff_end < scan_start) next
      pos <- base::seq.int(scan_start, eff_end)
      core_vec <- base::substring(arm_seq, pos, pos + k - 1L)
      oligo_vec <- if (reverse) .rc_vec(core_vec) else core_vec
      implied_vec <- .shared_role_arm_length(role, pos, k, arm_len, preferred_bp)
      keep <- implied_vec >= min_arm_bp & implied_vec <= max_arm_bp
      arm_dev_vec <- .shared_role_arm_dev(role, pos, k, arm_len, preferred_bp)
      if (boundary == "inner") {
        boundary_dist_vec <- if (region == "end") arm_len - (pos + k - 1L) else pos - 1L
      } else {
        boundary_dist_vec <- if (region == "start") pos - 1L else arm_len - (pos + k - 1L)
      }
      # Pin inner-boundary roles (UR, DF) to the deletion boundary; only the
      # outer roles (UF, DR) should have positional flexibility. Without this
      # the outer-role search can drag the inner primer far from the start /
      # stop codon, giving an effectively-shifted homology arm.
      if (boundary == "inner" &&
          !base::is.null(inner_boundary_tolerance_bp) &&
          base::is.finite(inner_boundary_tolerance_bp)) {
        keep <- keep & (boundary_dist_vec <= base::as.integer(inner_boundary_tolerance_bp))
      }
      if (!base::any(keep)) next
      parts[[base::length(parts) + 1L]] <- base::data.frame(
        core = core_vec[keep],
        oligo = oligo_vec[keep],
        position = pos[keep],
        k = k,
        implied_arm = implied_vec[keep],
        arm_dev = arm_dev_vec[keep],
        boundary_dist = boundary_dist_vec[keep],
        stringsAsFactors = FALSE
      )
    }
    if (base::length(parts) == 0) next
    gframe <- base::do.call(base::rbind, parts)
    gframe$genome_id <- id

    # Outer-role seed pre-filter (gated by `seed_prefilter_enabled`): most of
    # the outer-role candidates sit in repetitive regions (transposons, rRNA
    # repeats) whose 3' 15-mer also appears elsewhere in the genome. Building
    # a PDict over the unique 15-mer seeds and running countPDict once against
    # the whole genome lets us drop non-unique-seed positions in a single
    # pass. Disabled by default because, for targets adjacent to mobile
    # elements, hard-filtering by seed uniqueness can wipe out the entire
    # outer-role candidate pool. Inner-role candidates are already pinned to
    # a small window so no prefilter needed.
    if (base::isTRUE(seed_prefilter_enabled) &&
        boundary == "outer" &&
        !base::is.null(dnastr_cache) &&
        !base::is.null(dnastr_cache[[id]]) &&
        base::requireNamespace("Biostrings", quietly = TRUE) &&
        base::nrow(gframe) > 0) {
      seed_len <- 15L
      gdna <- dnastr_cache[[id]]
      nc <- base::nchar(gframe$oligo)
      seed_ok <- nc >= seed_len
      if (base::any(seed_ok)) {
        seeds <- base::rep(NA_character_, base::nrow(gframe))
        seeds[seed_ok] <- base::substr(gframe$oligo[seed_ok],
                                         nc[seed_ok] - seed_len + 1L,
                                         nc[seed_ok])
        seeds_valid <- seed_ok &
          base::grepl("^[ACGTacgt]+$", base::ifelse(base::is.na(seeds), "", seeds))
        if (base::any(seeds_valid)) {
          useeds <- base::unique(base::toupper(seeds[seeds_valid]))
          rc_useeds <- .rc_vec(useeds)
          counts <- base::tryCatch({
            pd_fwd <- Biostrings::PDict(useeds)
            pd_rev <- Biostrings::PDict(rc_useeds)
            f <- Biostrings::countPDict(pd_fwd, gdna)
            r <- Biostrings::countPDict(pd_rev, gdna)
            stats::setNames(base::as.integer(f + r), useeds)
          }, error = function(e) NULL)
          if (!base::is.null(counts)) {
            cand_counts <- base::rep(NA_integer_, base::nrow(gframe))
            cand_counts[seeds_valid] <- counts[base::toupper(seeds[seeds_valid])]
            keep <- base::is.na(cand_counts) | cand_counts == 1L
            # Guard against over-filtering: if fewer than ~5 % of candidates
            # (or <50 absolute) survive, the arm is likely in a deeply
            # repetitive region (e.g. transposon cluster) and hard-filtering
            # by seed uniqueness would eliminate every viable primer. Keep
            # the original pool in that case; the later uniqueness audit will
            # still report n0 > 1 for the non-unique survivors.
            n_before <- base::nrow(gframe)
            n_after  <- base::sum(keep)
            if (n_after >= 50L && n_after >= 0.05 * n_before) {
              gframe <- gframe[keep, , drop = FALSE]
            }
          }
        }
      }
    }

    if (base::nrow(gframe) == 0) next
    per_genome_frames[[gi]] <- gframe
  }
  long <- base::do.call(base::rbind, per_genome_frames)
  if (base::is.null(long) || base::nrow(long) == 0) return(base::list())

  # --- Aggregate per oligo, pick the best position per (oligo, genome) -----
  # Best = smallest arm_dev (closest to preferred arm length).
  long <- long[base::order(long$oligo, long$genome_id, long$arm_dev), , drop = FALSE]
  first_per_genome <- !base::duplicated(base::paste(long$oligo, long$genome_id))
  long <- long[first_per_genome, , drop = FALSE]

  # Build the raw pool WITHOUT Tm / self-dimer yet -- these per-oligo
  # evaluations dominate wall time when arm length is extended (e.g. 30 s on
  # ~30 k candidates for Tm_NN alone). Both are needed only for the top-cap
  # candidates that survive the fast GC/arm-based pre-rank, so we defer them.
  split_oligos <- base::split(long, long$oligo)
  raw_pool <- base::lapply(split_oligos, function(df) {
    core <- df$core[1]
    oligo <- df$oligo[1]
    gc <- 100 * .count_gc(core) / base::nchar(core)
    base::list(
      core = core,
      oligo = oligo,
      genomes = df$genome_id,
      gc = gc,
      tm_core = NA_real_,  # filled after pre-rank cap
      positions = stats::setNames(df$position, df$genome_id),
      implied_arms = stats::setNames(df$implied_arm, df$genome_id),
      arm_devs = stats::setNames(df$arm_dev, df$genome_id),
      boundary_dists = stats::setNames(df$boundary_dist, df$genome_id)
    )
  })
  # Cheap filter first -- GC only. self-dimer check runs after cap.
  # Also drop any candidate whose reverse complement matches a partner-role
  # oligo already selected for the SAME arm (UF<->UR, DF<->DR); otherwise the
  # two primers bind the same arm segment and the effective arm collapses.
  if (base::length(forbidden_oligos) > 0) {
    forbidden_up <- base::toupper(base::as.character(forbidden_oligos))
    raw_pool <- Filter(function(x) {
      !(base::toupper(.rc(x$oligo)) %in% forbidden_up)
    }, raw_pool)
    if (base::length(raw_pool) == 0) return(base::list())
  }
  raw_pool <- Filter(function(x) {
    !base::is.null(x) && x$gc >= 40 && x$gc <= 70
  }, raw_pool)
  if (base::length(raw_pool) == 0) return(base::list())

  # Lexicographic pre-ranking (Tm intentionally omitted here; it's added as a
  # tiebreaker only for the top-cap survivors below).
  pre_order <- base::order(
    vapply(raw_pool, function(x) -base::length(x$genomes), numeric(1)),
    vapply(raw_pool, function(x) base::max(x$arm_devs), numeric(1)),
    vapply(raw_pool, function(x) base::mean(x$arm_devs), numeric(1)),
    vapply(raw_pool, function(x) base::max(x$boundary_dists), numeric(1)),
    vapply(raw_pool, function(x) base::abs(x$gc - 55), numeric(1)),
    vapply(raw_pool, function(x) base::nchar(x$oligo), numeric(1))
  )
  # Apply a generous pre-cap multiplier so we still have room for Tm-based
  # tiebreaking among the survivors without evaluating Tm on the full pool.
  pre_cap <- if (base::is.na(candidate_cap)) base::length(pre_order)
             else base::min(base::length(pre_order),
                              base::as.integer(candidate_cap * 4L))
  raw_pool <- raw_pool[pre_order[base::seq_len(pre_cap)]]

  # Post-cap: compute Tm once per survivor and drop self-dimers.
  raw_pool <- base::lapply(raw_pool, function(x) {
    x$tm_core <- .safe_tm(x$oligo); x
  })
  raw_pool <- Filter(function(x) !.check_self_dimerization(x$oligo), raw_pool)
  if (base::length(raw_pool) == 0) return(base::list())

  # Hard Tm-window filter (used for check primers where Tm must land within a
  # narrow band). For the main UF/UR/DF/DR design this is left NULL and Tm is
  # used only as a scoring tiebreaker below.
  if (!base::is.null(tm_window)) {
    tm_window <- base::as.numeric(tm_window)
    raw_pool <- Filter(function(x) {
      !base::is.na(x$tm_core) &&
        x$tm_core >= tm_window[1] &&
        x$tm_core <= tm_window[2]
    }, raw_pool)
    if (base::length(raw_pool) == 0) return(base::list())
  }

  # Final rank (Tm tiebreaker included now) and hard cap to candidate_cap.
  fine_order <- base::order(
    vapply(raw_pool, function(x) -base::length(x$genomes), numeric(1)),
    vapply(raw_pool, function(x) base::max(x$arm_devs), numeric(1)),
    vapply(raw_pool, function(x) base::mean(x$arm_devs), numeric(1)),
    vapply(raw_pool, function(x) base::max(x$boundary_dists), numeric(1)),
    vapply(raw_pool, function(x) base::abs(.na0(x$tm_core) - tm_target), numeric(1)),
    vapply(raw_pool, function(x) base::abs(x$gc - 55), numeric(1)),
    vapply(raw_pool, function(x) base::nchar(x$oligo), numeric(1))
  )
  if (!base::is.na(candidate_cap)) {
    candidate_cap <- base::min(base::length(fine_order), candidate_cap)
    raw_pool <- raw_pool[fine_order[base::seq_len(candidate_cap)]]
  } else {
    raw_pool <- raw_pool[fine_order]
  }

  # Batch-prefill the exact-hit-count cache across all (candidate, genome)
  # pairs using Biostrings::matchPDict. Populating this once (Aho-Corasick
  # over all candidates of the same length in a single genome pass) replaces
  # hundreds of per-candidate matchPattern calls and dominates the wall-clock
  # win of the role loop.
  if (require_unique) {
    .shared_prefill_exact_hits(
      candidates = raw_pool,
      context_map = context_map,
      hit_count_cache = hit_count_cache,
      dnastr_cache = dnastr_cache
    )
  }

  pool <- base::lapply(raw_pool, function(obj) {
    unique_flags <- vapply(obj$genomes, function(id) {
      if (!require_unique) return(TRUE)
      .shared_exact_hit_count_cached(
        genome_id = id,
        primer_seq = obj$oligo,
        context_map = context_map,
        hit_count_cache = hit_count_cache,
        dnastr_cache = dnastr_cache
      ) == 1L
    }, logical(1))
    # Partial-sharing: drop genomes where uniqueness fails, but keep the
    # candidate alive for the genomes where it does pass. Previously any
    # single uniqueness failure rejected the whole candidate, which caused
    # the multi-pass fallback (cap 20/100/500/NA) to blow up when a few
    # genomes had promiscuous matches.
    if (require_unique) {
      keep <- unique_flags
      if (!base::any(keep)) return(NULL)
      kept_ids <- obj$genomes[keep]
      obj$genomes <- kept_ids
      if (!base::is.null(obj$positions))       obj$positions <- obj$positions[kept_ids]
      if (!base::is.null(obj$implied_arms))    obj$implied_arms <- obj$implied_arms[kept_ids]
      if (!base::is.null(obj$arm_devs))        obj$arm_devs <- obj$arm_devs[kept_ids]
      if (!base::is.null(obj$boundary_dists))  obj$boundary_dists <- obj$boundary_dists[kept_ids]
      obj$all_unique <- base::all(unique_flags)
    } else {
      obj$all_unique <- TRUE
    }
    obj
  })
  Filter(base::Negate(base::is.null), pool)
}

.shared_select_primer_groups_from_pool <- function(candidate_pool, genome_ids, role, reverse) {
  if (base::length(candidate_pool) == 0) {
    base::stop("No valid candidate found for role ", role)
  }
  remaining <- genome_ids
  out <- base::list()
  group_index <- 1L
  while (base::length(remaining) > 0) {
    best <- NULL
    best_score <- NULL
    for (cand in candidate_pool) {
      current_members <- base::intersect(cand$genomes, remaining)
      if (base::length(current_members) == 0) next
      arm_devs_m <- cand$arm_devs[current_members]
      boundary_m <- cand$boundary_dists[current_members]
      tm_dev <- base::abs(.na0(cand$tm_core) - 60)
      score_vec <- c(
        base::length(current_members),               # 1. share as many genomes as possible
        -base::max(arm_devs_m),                      # 2. worst-case arm deviation from preferred
        -base::mean(arm_devs_m),                     # 3. mean arm deviation
        -base::max(boundary_m),                      # 4. keep primer near its ideal boundary
        -(2 * tm_dev),                               # 5. weighted Tm polishing (max term)
        -base::abs(cand$gc - 55),                    # 6. GC polishing
        -base::nchar(cand$oligo)                     # 7. shorter primer preferred
      )
      if (base::is.null(best_score) || .shared_numeric_score_better(score_vec, best_score)) {
        best <- list(
          candidate = cand,
          genomes = current_members,
          worst_arm_dev = base::max(arm_devs_m),
          mean_arm_dev = base::mean(arm_devs_m),
          worst_boundary = base::max(boundary_m)
        )
        best_score <- score_vec
      }
    }
    if (base::is.null(best)) {
      base::stop("No valid candidate found for role ", role, " among remaining genomes.")
    }
    scope <- if (base::length(best$genomes) == 1L) {
      base::paste0(best$genomes, " only")
    } else if (base::length(best$genomes) == base::length(genome_ids)) {
      "all"
    } else {
      base::paste(best$genomes, collapse = ", ")
    }
    cand <- best$candidate
    members <- best$genomes
    out[[base::length(out) + 1L]] <- base::data.frame(
      role = role,
      cluster_id = base::paste0(role, "_cluster_", group_index),
      scope = scope,
      used_by = base::paste(members, collapse = ", "),
      primer_core_5to3 = cand$oligo,
      shared_core = cand$core,
      reverse_primer = reverse,
      primer_length = base::nchar(cand$oligo),
      gc_percent = base::round(cand$gc, 1),
      tm_core = .safe_round(cand$tm_core),
      worst_arm_dev_bp = base::as.integer(best$worst_arm_dev),
      mean_arm_dev_bp = base::round(best$mean_arm_dev, 1),
      boundary_distance_bp = base::as.integer(best$worst_boundary),
      exact_unique_in_targets = cand$all_unique,
      stringsAsFactors = FALSE
    )
    remaining <- base::setdiff(remaining, members)
    group_index <- group_index + 1L
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

# Batch exact-hit-count prefill via Biostrings::countPDict.
#
# Populates `hit_count_cache` for every (genome_id, oligo) pair that the
# per-candidate loop would otherwise look up one call at a time. countPDict
# (vs matchPDict) avoids allocating the position-list MIndex; combined with
# a once-per-length PDict build it roughly doubles throughput on short primer
# patterns against ~3-4 Mbp genomes.
.shared_prefill_exact_hits <- function(candidates, context_map,
                                       hit_count_cache, dnastr_cache) {
  if (base::length(candidates) == 0) return(invisible(NULL))
  if (!base::requireNamespace("Biostrings", quietly = TRUE)) return(invisible(NULL))
  if (base::is.null(dnastr_cache)) return(invisible(NULL))

  # Collect (oligo, genomes-needed) pairs; drop ambiguity-code sequences that
  # PDict does not accept.
  oligo_vec <- base::vapply(candidates, `[[`, character(1), "oligo")
  genomes_vec <- base::lapply(candidates, `[[`, "genomes")
  valid <- base::grepl("^[ACGTacgt]+$", oligo_vec)
  if (!base::any(valid)) return(invisible(NULL))

  pair_oligo  <- base::character(0)
  pair_genome <- base::character(0)
  for (i in base::which(valid)) {
    oligo <- base::toupper(oligo_vec[i])
    gids <- base::unique(base::as.character(genomes_vec[[i]]))
    for (gid in gids) {
      if (base::is.null(dnastr_cache[[gid]])) next
      key <- base::paste(gid, oligo, sep = "::")
      if (!base::exists(key, envir = hit_count_cache, inherits = FALSE)) {
        pair_oligo  <- base::c(pair_oligo,  oligo)
        pair_genome <- base::c(pair_genome, gid)
      }
    }
  }
  if (base::length(pair_oligo) == 0) return(invisible(NULL))

  # Build one PDict (fwd + rc) per oligo length, then scan each genome.
  klen <- base::nchar(pair_oligo)
  for (k in base::unique(klen)) {
    at_k <- klen == k
    oligos_k <- base::unique(pair_oligo[at_k])
    rc_k <- base::vapply(oligos_k, .rc, character(1), USE.NAMES = FALSE)
    pd_fwd <- base::tryCatch(Biostrings::PDict(oligos_k), error = function(e) NULL)
    pd_rev <- base::tryCatch(Biostrings::PDict(rc_k),     error = function(e) NULL)
    gids_k <- base::unique(pair_genome[at_k])
    for (gid in gids_k) {
      gdna <- dnastr_cache[[gid]]
      if (base::is.null(gdna)) next
      if (!base::is.null(pd_fwd) && !base::is.null(pd_rev)) {
        fwd <- Biostrings::countPDict(pd_fwd, gdna)
        rev <- Biostrings::countPDict(pd_rev, gdna)
      } else {
        # Fallback: per-pattern countPattern when PDict is unavailable.
        fwd <- base::vapply(oligos_k, function(o)
          Biostrings::countPattern(Biostrings::DNAString(o), gdna, fixed = TRUE),
          integer(1))
        rev <- base::vapply(rc_k, function(o)
          Biostrings::countPattern(Biostrings::DNAString(o), gdna, fixed = TRUE),
          integer(1))
      }
      total <- base::as.integer(fwd + rev)
      for (j in base::seq_along(oligos_k)) {
        key <- base::paste(gid, oligos_k[j], sep = "::")
        base::assign(key, total[j], envir = hit_count_cache)
      }
    }
  }
  base::invisible(NULL)
}

.shared_exact_hit_count_cached <- function(genome_id, primer_seq, context_map,
                                           hit_count_cache, dnastr_cache = NULL) {
  key <- base::paste(genome_id, primer_seq, sep = "::")
  if (base::exists(key, envir = hit_count_cache, inherits = FALSE)) {
    return(base::get(key, envir = hit_count_cache, inherits = FALSE))
  }
  primer_seq <- base::toupper(primer_seq)
  total <- if (!base::is.null(dnastr_cache) && !base::is.null(dnastr_cache[[genome_id]])) {
    # Fast path: Biostrings countPattern on pre-built DNAString, O(n+m) in C.
    # countPattern avoids allocating the full match-position XStringViews
    # object that matchPattern would return, a ~1.7x saving for counts-only
    # uniqueness queries.
    pat <- Biostrings::DNAString(primer_seq)
    pat_rc <- Biostrings::reverseComplement(pat)
    fwd <- Biostrings::countPattern(pat, dnastr_cache[[genome_id]], fixed = TRUE)
    rev <- Biostrings::countPattern(pat_rc, dnastr_cache[[genome_id]], fixed = TRUE)
    base::as.integer(fwd + rev)
  } else {
    genome_seq <- context_map[[genome_id]]$genome_seq
    rc_seq <- .rc(primer_seq)
    hits_fwd <- gregexpr(primer_seq, genome_seq, fixed = TRUE)[[1]]
    hits_rev <- gregexpr(rc_seq, genome_seq, fixed = TRUE)[[1]]
    count_fwd <- if (hits_fwd[1] < 0) 0L else base::length(hits_fwd)
    count_rev <- if (hits_rev[1] < 0) 0L else base::length(hits_rev)
    base::as.integer(count_fwd + count_rev)
  }
  base::assign(key, total, envir = hit_count_cache)
  total
}

# Vectorized RC via Biostrings::DNAStringSet. Reverse-role candidate pool
# building calls this on tens of thousands of oligos per genome; per-element
# Biostrings::DNAString allocation (what .rc does) was ~370us/call, a batched
# DNAStringSet is under 1us/call (~265x faster on a 10k batch). Falls back
# to the per-element helper when Biostrings is unavailable or the batch
# contains ambiguity codes DNAStringSet would reject.
.rc_vec <- function(seqs) {
  if (base::length(seqs) == 0) return(character(0))
  if (base::requireNamespace("Biostrings", quietly = TRUE)) {
    out <- base::tryCatch(
      base::as.character(
        Biostrings::reverseComplement(
          Biostrings::DNAStringSet(base::toupper(seqs)))),
      error = function(e) NULL)
    if (!base::is.null(out)) return(out)
  }
  base::vapply(seqs, .rc, character(1), USE.NAMES = FALSE)
}
.count_gc <- function(seq_chr) {
  base::nchar(base::gsub("[^GCgc]", "", seq_chr))
}

# Whole-genome primer off-target profile using Biostrings. Returns exact (n0)
# and cumulative mismatch counts up to `max_mismatch`; both strands counted.
.shared_primer_offtarget_hits <- function(primer_seq, genome_seq, max_mismatch = 1L) {
  primer <- Biostrings::DNAString(base::toupper(primer_seq))
  primer_rc <- Biostrings::reverseComplement(primer)
  if (base::inherits(genome_seq, "DNAString")) {
    genome <- genome_seq
  } else {
    genome <- Biostrings::DNAString(base::toupper(genome_seq))
  }
  fwd0 <- Biostrings::matchPattern(primer, genome, max.mismatch = 0L)
  rev0 <- Biostrings::matchPattern(primer_rc, genome, max.mismatch = 0L)
  n0 <- base::as.integer(base::length(fwd0) + base::length(rev0))
  n1 <- 0L
  if (base::as.integer(max_mismatch) >= 1L) {
    fwd1 <- Biostrings::matchPattern(primer, genome, max.mismatch = 1L)
    rev1 <- Biostrings::matchPattern(primer_rc, genome, max.mismatch = 1L)
    n1 <- base::as.integer(base::length(fwd1) + base::length(rev1) - n0)
  }
  base::list(n0 = n0, n1 = n1)
}

# Per-primer per-genome audit of off-target binding. Adds the following
# columns to primer_cores: worst_n0, worst_n1, per_genome_n0, per_genome_n1.
# Emits a warning if any primer still has duplicate exact hits (worst_n0 > 1)
# or exceeds `max_n1_per_primer` one-mismatch hits in any assigned genome.
# `audit_n_mismatches` controls how deeply we scan: 0 = exact only (fast),
# 1 = also count 1-mismatch hits via `Biostrings::matchPattern` (slow on
# multi-megabase genomes).
.shared_annotate_primer_uniqueness <- function(primer_cores, genome_sequences,
                                               max_n1_per_primer = NULL,
                                               audit_n_mismatches = 0L) {
  if (base::nrow(primer_cores) == 0) return(primer_cores)
  if (!base::requireNamespace("Biostrings", quietly = TRUE)) {
    base::warning("Biostrings is not available; skipping primer uniqueness audit.",
                  call. = FALSE)
    return(primer_cores)
  }
  audit_n_mismatches <- base::as.integer(audit_n_mismatches)

  # Cache genome DNAString objects once per id; reused across primers.
  dnastr_cache <- base::lapply(genome_sequences, function(seq) {
    Biostrings::DNAString(base::toupper(seq))
  })

  worst_n0 <- base::integer(base::nrow(primer_cores))
  worst_n1 <- base::integer(base::nrow(primer_cores))
  per_n0_str <- base::character(base::nrow(primer_cores))
  per_n1_str <- base::character(base::nrow(primer_cores))

  for (i in base::seq_len(base::nrow(primer_cores))) {
    oligo <- base::as.character(primer_cores$primer_core_5to3[i])
    members <- base::strsplit(primer_cores$used_by[i], ", ", fixed = TRUE)[[1]]
    n0_vec <- base::integer(base::length(members))
    n1_vec <- base::integer(base::length(members))
    for (j in base::seq_along(members)) {
      id <- members[j]
      dna <- dnastr_cache[[id]]
      if (base::is.null(dna)) {
        n0_vec[j] <- NA_integer_
        n1_vec[j] <- NA_integer_
        next
      }
      hits <- .shared_primer_offtarget_hits(oligo, dna, max_mismatch = audit_n_mismatches)
      n0_vec[j] <- hits$n0
      n1_vec[j] <- hits$n1
    }
    worst_n0[i] <- base::suppressWarnings(base::max(n0_vec, na.rm = TRUE))
    worst_n1[i] <- base::suppressWarnings(base::max(n1_vec, na.rm = TRUE))
    per_n0_str[i] <- base::paste(members, n0_vec, sep = "=", collapse = ", ")
    per_n1_str[i] <- if (audit_n_mismatches >= 1L) {
      base::paste(members, n1_vec, sep = "=", collapse = ", ")
    } else {
      NA_character_
    }
  }

  primer_cores$worst_n0 <- worst_n0
  primer_cores$worst_n1 <- if (audit_n_mismatches >= 1L) worst_n1 else NA_integer_
  primer_cores$per_genome_n0 <- per_n0_str
  primer_cores$per_genome_n1 <- per_n1_str

  non_unique <- primer_cores$worst_n0 > 1L
  if (base::any(non_unique, na.rm = TRUE)) {
    base::warning(sprintf(
      "Primer(s) with >1 exact genome hit(s) survived selection (worst_n0 > 1): %s",
      base::paste(primer_cores$cluster_id[non_unique], collapse = ", ")
    ), call. = FALSE)
  }

  if (!base::is.null(max_n1_per_primer) && audit_n_mismatches >= 1L) {
    bad <- primer_cores$worst_n1 > base::as.integer(max_n1_per_primer)
    bad[base::is.na(bad)] <- FALSE
    if (base::any(bad)) {
      base::warning(sprintf(
        "Primer(s) exceed max_n1_per_primer=%d: %s",
        base::as.integer(max_n1_per_primer),
        base::paste(primer_cores$cluster_id[bad], collapse = ", ")
      ), call. = FALSE)
    }
  }

  primer_cores
}
