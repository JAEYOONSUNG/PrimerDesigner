#' Visualize a Shared gRNA / Deletion Primer Design
#'
#' Generates a single-page PDF that summarizes the output of
#' \code{design_shared_grna_and_deletion()}. Seven panels are produced:
#' \enumerate{
#'   \item A - Combined constructs (full plasmid, ~13 kb).
#'   \item B - Insert-site zoom with per-role primer lanes.
#'   \item C - Primer sharing matrix (role x genome).
#'   \item D - Unique oligos to order (one bar per distinct full primer).
#'   \item E - Upstream arm alignment heatmap + primer bands.
#'   \item F - Downstream arm alignment heatmap + primer bands.
#'   \item G - Flanking gene context (which CDSes occupy the arms).
#' }
#' Genomes are ordered by hierarchical clustering on the combined arm sequences
#' so similar arms sit next to each other; the dendrogram is drawn alongside
#' the heatmaps in E/F.
#'
#' @param result The list returned by \code{design_shared_grna_and_deletion()}.
#'   Must contain at least \code{target_context}, \code{primer_cores},
#'   \code{construct_deletion_primers}, and \code{construct_groups}.
#' @param genbank_dir Character. Directory containing the per-genome genome
#'   files used as inputs to the design (for Panel G flanking annotations).
#'   Each genome_id in \code{result$target_context} must have a matching
#'   \code{.gbk/.gb/.gbff/.genbank} or SnapGene \code{.dna} file whose
#'   basename (without extension) equals the genome_id.
#' @param construct_gbk_dir Character. Directory that holds the combined
#'   construct \code{.gbk} files produced by the pipeline (typically the same
#'   \code{combined_construct_output_dir} passed to
#'   \code{design_shared_grna_and_deletion()}). One file per construct group.
#' @param output_file Output PDF path. A sibling \code{.png} will also be
#'   written if \code{also_png = TRUE}.
#' @param width,height Numeric. PDF canvas dimensions (inches). Defaults tuned
#'   for an A3-like single-page summary.
#' @param also_png Logical. Write a matching PNG at 180 dpi. Default \code{FALSE}.
#' @return Invisibly returns the path of the written PDF.
#' @examples
#' \dontrun{
#' res <- design_shared_grna_and_deletion(
#'   target_table = tt, nuclease = "GeoCas9",
#'   combined_construct_output_dir = "out/constructs/combined",
#'   combined_vector_file = "pJET.dna",
#'   combined_grna_start = 100, combined_grna_end = 121,
#'   combined_deletion_start = 3900, combined_deletion_end = 5050)
#'
#' visualize_shared_design(
#'   result = res,
#'   genbank_dir = "data/gbk",
#'   construct_gbk_dir = "out/constructs/combined",
#'   output_file = "out/shared_design_summary.pdf")
#' }
#' @export
visualize_shared_design <- function(
    result,
    genbank_dir,
    construct_gbk_dir,
    output_file,
    width = 14,
    height = 17,
    also_png = FALSE
) {
  .need_pkg <- function(pkg) if (!base::requireNamespace(pkg, quietly = TRUE))
    base::stop("Package '", pkg, "' is required for visualize_shared_design().")
  base::lapply(c("ggplot2", "dplyr", "tidyr", "patchwork", "ggtext",
                 "ggdendro"), .need_pkg)

  if (!base::is.list(result) ||
      !base::all(c("target_context", "primer_cores",
                   "construct_deletion_primers") %in% base::names(result))) {
    base::stop("result does not look like a design_shared_grna_and_deletion() output.")
  }
  if (!base::dir.exists(genbank_dir))
    base::stop("genbank_dir not found: ", genbank_dir)
  if (!base::dir.exists(construct_gbk_dir))
    base::stop("construct_gbk_dir not found: ", construct_gbk_dir)

  tc <- result$target_context
  pc <- result$primer_cores
  cdp <- result$construct_deletion_primers

  # Resolve per-genome GenBank paths for Panel G
  gbk_paths <- .shared_viz_resolve_gbk(tc$genome_id, genbank_dir)

  # 1) Parse construct .gbk (A, B panels)
  construct_files <- base::list.files(construct_gbk_dir, pattern = "\\.gbk$",
                                      full.names = TRUE)
  if (base::length(construct_files) == 0) {
    base::stop("No .gbk files in construct_gbk_dir: ", construct_gbk_dir)
  }

  combo <- .shared_viz_build(
    tc = tc, pc = pc, cdp = cdp,
    construct_files = construct_files,
    gbk_paths = gbk_paths
  )

  # Save PDF (cairo_pdf preserves text as vector)
  ggplot2::ggsave(output_file, combo, width = width, height = height,
                  device = grDevices::cairo_pdf)
  if (also_png) {
    png_path <- base::sub("\\.pdf$", ".png", output_file, ignore.case = TRUE)
    ggplot2::ggsave(png_path, combo, width = width, height = height, dpi = 180)
  }
  base::cat("Wrote:", output_file, "\n")
  base::invisible(output_file)
}

# =============================================================================
# Internal: GenBank path resolution
# =============================================================================
.shared_viz_resolve_gbk <- function(genome_ids, dir) {
  exts <- c("gbk", "gb", "gbff", "genbank", "dna")
  out <- base::character(base::length(genome_ids))
  base::names(out) <- genome_ids
  for (g in genome_ids) {
    candidates <- base::file.path(dir, base::paste0(g, ".", exts))
    hit <- candidates[base::file.exists(candidates)]
    if (base::length(hit) == 0) {
      base::stop("Genome file for genome_id '", g, "' not found in ", dir,
                 ". Expected one of: ", base::paste(base::basename(candidates),
                                                    collapse = ", "))
    }
    picked <- hit[1]
    # .dna (SnapGene) is converted on the fly to a temp GenBank file.
    if (base::tolower(tools::file_ext(picked)) == "dna") {
      gb_lines <- .get_genbank_from_dna(picked, kill_snapgene = FALSE)
      tmp_gbk <- base::file.path(base::tempdir(), base::paste0(g, ".gbk"))
      base::writeLines(gb_lines, tmp_gbk)
      picked <- tmp_gbk
    }
    out[[g]] <- picked
  }
  out
}

# =============================================================================
# Internal: main plot assembly
# =============================================================================
.shared_viz_build <- function(tc, pc, cdp, construct_files, gbk_paths) {
  ggplot2 <- NULL; dplyr <- NULL  # avoid R CMD check NOTE
  `%>%` <- magrittr::`%>%`

  th <- ggplot2::theme_minimal(base_size = 9.5, base_family = "") +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 10.5,
                                                colour = "grey10",
                                                margin = ggplot2::margin(b = 3)),
      plot.subtitle    = ggplot2::element_text(colour = "grey45", size = 8.5,
                                                margin = ggplot2::margin(b = 4)),
      panel.grid       = ggplot2::element_blank(),
      axis.ticks.x     = ggplot2::element_line(colour = "grey75", linewidth = 0.25),
      axis.ticks.y     = ggplot2::element_blank(),
      axis.text        = ggplot2::element_text(colour = "grey30"),
      axis.title.x     = ggplot2::element_text(colour = "grey30", size = 8.5,
                                                margin = ggplot2::margin(t = 3)),
      legend.position  = "bottom",
      legend.box.spacing = grid::unit(1, "pt"),
      legend.margin    = ggplot2::margin(t = 2, b = 0),
      legend.key.size  = grid::unit(0.32, "cm"),
      legend.key.height = grid::unit(0.32, "cm"),
      legend.title     = ggplot2::element_text(size = 8.5, colour = "grey25"),
      legend.text      = ggplot2::element_text(size = 8),
      plot.margin      = ggplot2::margin(6, 10, 6, 10)
    )

  # Wong 2011 color-blind safe role palette (F=blue family, R=orange family)
  pal_role <- c(UF = "#0072B2", UR = "#D55E00", DF = "#56B4E9", DR = "#E69F00")
  pal_track <- c(
    "Gibson primer"      = "#D55E00",
    "gRNA spacer"        = "#CC79A7",
    "Upstream arm"       = "#009E73",
    "Downstream arm"     = "#B3E0C6",
    "CDS"                = "#B8B8B8",
    "Promoter"           = "#F0E442",
    "Terminator"         = "#7D7D7D",
    "Origin"             = "#8C6BB1",
    "Vector primer_bind" = "#DDDDDD",
    "Other"              = "#EFEFEF"
  )

  # ---- Panels A & B: construct overview + insert-site zoom ------------------
  ab <- .shared_viz_panels_AB(construct_files, th, pal_role, pal_track)
  p_overview <- ab$p_overview
  p_zoom     <- ab$p_zoom

  # ---- Panel C: sharing matrix ---------------------------------------------
  p_share <- .shared_viz_panel_C(pc, tc, th, pal_role)

  # ---- Panel D: unique oligos to order -------------------------------------
  p_order <- .shared_viz_panel_D(cdp, pc, pal_role)

  # ---- Panels E/F + dendrogram: arm alignments ------------------------------
  ef <- .shared_viz_panels_EF(tc, pc, th, pal_track)

  # ---- Panel G: flanking gene context --------------------------------------
  p_G <- .shared_viz_panel_G(tc, gbk_paths, th, pal_track,
                              ef$clust_order_wrapped)

  # ---- Assemble -------------------------------------------------------------
  mid    <- (p_share | p_order) + patchwork::plot_layout(widths = c(1.2, 1.0))
  dendro_w <- 28
  bottom <- (ef$p_dendro | ef$p_left | ef$p_right) +
    patchwork::plot_layout(
      widths = c(dendro_w, ef$min_l, ef$min_r) /
                base::sum(c(dendro_w, ef$min_l, ef$min_r)))

  combo <- p_overview / p_zoom / mid / bottom / p_G +
    patchwork::plot_layout(heights = c(0.80, 0.70, 0.70, 0.85, 0.70)) +
    patchwork::plot_annotation(
      title = "Shared-design summary",
      subtitle = base::sprintf(
        "%d genomes . %d construct groups . %d unique oligo records",
        base::nrow(tc),
        base::length(base::unique(
          tools::file_path_sans_ext(base::basename(construct_files)))),
        base::nrow(pc)),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13),
        plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9))
    )
  combo
}

# =============================================================================
# Internal: panels A + B (combined construct overview + zoom)
# =============================================================================
.shared_viz_panels_AB <- function(construct_files, th, pal_role, pal_track) {
  `%>%` <- magrittr::`%>%`
  tracks <- base::lapply(construct_files, function(f) {
    short <- base::sub("__gRNA_cluster_1_combined_construct$", "",
                       tools::file_path_sans_ext(base::basename(f)))
    p <- .shared_viz_parse_gbk(f)
    p$features$construct <- short
    p$features$total_len <- p$length
    p
  })
  feat_df <- dplyr::bind_rows(base::lapply(tracks, function(t) t$features))

  classify <- function(type, label) {
    lab <- base::ifelse(base::is.na(label), "", label)
    dplyr::case_when(
      base::grepl("::UF_cluster|::UR_cluster|::DF_cluster|::DR_cluster", lab) ~ "Gibson primer",
      base::grepl("^gRNA", lab) | base::grepl("spacer", lab, ignore.case = TRUE) ~ "gRNA spacer",
      base::grepl("^UP ", lab) ~ "Upstream arm",
      base::grepl("^DN ", lab) ~ "Downstream arm",
      type == "primer_bind" ~ "Vector primer_bind",
      type == "CDS" ~ "CDS",
      type == "promoter" ~ "Promoter",
      type == "terminator" ~ "Terminator",
      type == "rep_origin" ~ "Origin",
      type == "source" ~ "source",
      TRUE ~ "Other")
  }
  feat_df$category <- classify(feat_df$type, feat_df$label)
  feat_df <- feat_df %>% dplyr::filter(category != "source")
  construct_order <- base::unique(feat_df$construct)

  display_label <- function(x) {
    x <- base::gsub("_CP\\d+", "", x)
    x <- base::gsub("_NZ", "", x)
    out <- x
    common_idx <- base::grepl("_common$", x)
    if (base::any(common_idx)) {
      stripped <- base::sub("_common$", "", x[common_idx])
      out[common_idx] <- base::vapply(
        base::strsplit(stripped, "_", fixed = TRUE),
        function(t) base::paste0(base::paste(t, collapse = "\n"), "\n(common)"),
        character(1))
    }
    out
  }
  feat_df$construct_label <- display_label(base::as.character(feat_df$construct))
  lab_order <- display_label(construct_order)
  feat_df$construct_label <- base::factor(feat_df$construct_label,
                                          levels = base::rev(lab_order))

  design_cats <- c("Gibson primer", "gRNA spacer", "Upstream arm", "Downstream arm")
  design_df <- feat_df %>%
    dplyr::filter(category %in% design_cats) %>%
    dplyr::mutate(role = base::sub(".*::", "", label) %>% base::sub("_.*", "", .),
                  fill_cat = base::ifelse(category == "Gibson primer" &
                                            role %in% base::names(pal_role),
                                          role, category))
  vector_df <- feat_df %>%
    dplyr::filter(!category %in% design_cats) %>%
    dplyr::mutate(fill_cat = category)

  p_overview <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = dplyr::distinct(feat_df, construct_label, total_len),
      ggplot2::aes(x = 1, xend = total_len,
                   y = base::as.numeric(construct_label),
                   yend = base::as.numeric(construct_label)),
      colour = "grey75", linewidth = 0.4) +
    ggplot2::geom_rect(
      data = vector_df,
      ggplot2::aes(xmin = start, xmax = end,
                   ymin = base::as.numeric(construct_label) - 0.18,
                   ymax = base::as.numeric(construct_label) + 0.18,
                   fill = fill_cat), color = NA) +
    ggplot2::geom_rect(
      data = design_df,
      ggplot2::aes(xmin = start, xmax = end,
                   ymin = base::as.numeric(construct_label) - 0.18,
                   ymax = base::as.numeric(construct_label) + 0.18,
                   fill = fill_cat), color = NA) +
    ggplot2::scale_y_continuous(
      breaks = base::seq_along(base::levels(feat_df$construct_label)),
      labels = base::levels(feat_df$construct_label),
      expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_fill_manual(
      values = c(pal_track, pal_role),
      breaks = c("gRNA spacer", "Upstream arm", "Downstream arm",
                 "UF", "UR", "DF", "DR",
                 "CDS", "Promoter", "Terminator", "Origin",
                 "Vector primer_bind", "Other"),
      name = NULL,
      guide = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
    ggplot2::labs(title = "A . Combined constructs (full)",
                  subtitle = "vector backbone + arms + gRNA spacer + Gibson primers",
                  x = "Position (bp)", y = NULL) + th

  # Zoom around insert site — auto-detect from design features
  insert_x <- base::range(design_df$start, design_df$end, na.rm = TRUE)
  insert_pad <- 150
  xlim_zoom <- c(base::max(1, insert_x[1] - insert_pad),
                 insert_x[2] + insert_pad)
  zoom_df <- feat_df %>% dplyr::filter(end >= xlim_zoom[1], start <= xlim_zoom[2])
  zoom_design <- zoom_df %>%
    dplyr::filter(category %in% design_cats) %>%
    dplyr::mutate(role = base::sub(".*::", "", label) %>% base::sub("_.*", "", .),
                  fill_cat = base::ifelse(category == "Gibson primer" &
                                            role %in% base::names(pal_role),
                                          role, category))
  zoom_vector <- zoom_df %>%
    dplyr::filter(!category %in% design_cats) %>%
    dplyr::mutate(fill_cat = category)
  gibson_df <- zoom_design %>% dplyr::filter(category == "Gibson primer") %>%
    dplyr::mutate(short = base::sub(".*::", "", label),
                  role = base::sub("_.*", "", short))

  p_zoom <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = dplyr::distinct(zoom_df, construct_label),
      ggplot2::aes(x = xlim_zoom[1], xend = xlim_zoom[2],
                   y = base::as.numeric(construct_label),
                   yend = base::as.numeric(construct_label)),
      colour = "grey80", linewidth = 0.3) +
    ggplot2::geom_rect(data = zoom_vector,
      ggplot2::aes(xmin = start, xmax = end,
                   ymin = base::as.numeric(construct_label) - 0.08,
                   ymax = base::as.numeric(construct_label) + 0.08,
                   fill = fill_cat), color = NA) +
    ggplot2::geom_rect(
      data = zoom_design %>% dplyr::filter(category != "Gibson primer"),
      ggplot2::aes(xmin = start, xmax = end,
                   ymin = base::as.numeric(construct_label) - 0.08,
                   ymax = base::as.numeric(construct_label) + 0.08,
                   fill = fill_cat), color = NA) +
    ggplot2::geom_rect(
      data = gibson_df %>% dplyr::filter(role %in% c("UF", "UR")),
      ggplot2::aes(xmin = start, xmax = end,
                   ymin = base::as.numeric(construct_label) + 0.12,
                   ymax = base::as.numeric(construct_label) + 0.26,
                   fill = role), color = NA) +
    ggplot2::geom_rect(
      data = gibson_df %>% dplyr::filter(role %in% c("DF", "DR")),
      ggplot2::aes(xmin = start, xmax = end,
                   ymin = base::as.numeric(construct_label) - 0.26,
                   ymax = base::as.numeric(construct_label) - 0.12,
                   fill = role), color = NA) +
    ggplot2::geom_text(
      data = gibson_df,
      ggplot2::aes(x = (start + end) / 2,
                   y = base::as.numeric(construct_label) +
                       base::ifelse(role %in% c("UF","UR"), 0.33, -0.33),
                   label = role,
                   vjust = base::ifelse(role %in% c("UF","UR"), 0, 1)),
      size = 2.3, colour = "grey10", fontface = "bold") +
    ggplot2::scale_y_continuous(
      breaks = base::seq_along(base::levels(feat_df$construct_label)),
      labels = base::levels(feat_df$construct_label),
      expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_fill_manual(values = c(pal_track, pal_role),
                               breaks = c("UF", "UR", "DF", "DR"),
                               name = NULL,
                               guide = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::coord_cartesian(xlim = xlim_zoom) +
    ggplot2::labs(title = "B . Insert-site zoom",
                  subtitle = "per-role primer lanes",
                  x = "Position (bp)", y = NULL) + th

  base::list(p_overview = p_overview, p_zoom = p_zoom)
}

# Minimal GenBank feature parser (type, location, /label) for panel A/B.
.shared_viz_parse_gbk <- function(path) {
  lines <- base::readLines(path, warn = FALSE)
  seq_len <- base::as.integer(base::strsplit(base::trimws(
    lines[base::grep("^LOCUS", lines)[1]]), "\\s+")[[1]][3])
  feat_start <- base::grep("^FEATURES", lines)[1]
  origin <- base::grep("^ORIGIN", lines)[1]
  feats <- lines[(feat_start + 1):(origin - 1)]
  rows <- base::list(); cur <- NULL; in_qualifiers <- FALSE
  for (ln in feats) {
    is_new <- base::grepl("^     \\S", ln) && !base::grepl("^     /", ln)
    if (is_new) {
      if (!base::is.null(cur)) rows[[base::length(rows) + 1]] <- cur
      parts <- base::regmatches(ln, base::regexec("^     (\\S+)\\s+(.*)$", ln))[[1]]
      cur <- base::list(type = parts[2], loc_raw = parts[3], label = NA_character_)
      in_qualifiers <- FALSE
      next
    }
    if (!base::is.null(cur) && base::grepl("^                     /", ln)) {
      in_qualifiers <- TRUE
      if (base::is.na(cur$label) && base::grepl("^                     /label=", ln)) {
        cur$label <- base::sub('^.*?/label=["]?(.*?)["]?$', "\\1", ln)
      }
      next
    }
    if (!base::is.null(cur) && !in_qualifiers &&
        base::grepl("^                     ", ln) &&
        !base::grepl("^                     /", ln)) {
      cur$loc_raw <- base::paste0(cur$loc_raw, base::trimws(ln))
    }
  }
  if (!base::is.null(cur)) rows[[base::length(rows) + 1]] <- cur
  df <- base::do.call(base::rbind, base::lapply(rows, function(r) {
    strand <- if (base::grepl("^complement", r$loc_raw)) -1 else 1
    nums <- base::as.integer(base::unlist(base::regmatches(
      r$loc_raw, base::gregexpr("\\d+", r$loc_raw))))
    if (base::length(nums) < 2) return(NULL)
    base::data.frame(type = r$type, start = nums[1], end = utils::tail(nums, 1),
                     strand = strand, label = r$label, stringsAsFactors = FALSE)
  }))
  base::list(length = seq_len, features = df)
}

# =============================================================================
# Internal: panel C (sharing matrix)
# =============================================================================
.shared_viz_panel_C <- function(pc, tc, th, pal_role) {
  `%>%` <- magrittr::`%>%`
  wrap_genome <- function(x) base::sub("^([^_]+)_(.*)$", "\\1\n\\2", x)

  share_df <- pc %>% dplyr::select(role, cluster_id, used_by) %>%
    dplyr::mutate(genome = base::strsplit(used_by, ", ")) %>%
    tidyr::unnest(genome)
  share_grid <- base::expand.grid(
    role = c("UF", "UR", "DF", "DR"),
    genome = tc$genome_id, stringsAsFactors = FALSE) %>%
    dplyr::left_join(share_df, by = c("role", "genome")) %>%
    dplyr::mutate(cluster_id = base::ifelse(base::is.na(cluster_id), "-",
                                             cluster_id))
  share_grid$role <- base::factor(share_grid$role,
                                   levels = c("DR", "DF", "UR", "UF"))
  share_grid$role_factor <- base::as.character(share_grid$role)

  ggplot2::ggplot(share_grid,
    ggplot2::aes(x = genome, y = role,
                  fill = base::ifelse(base::is.na(cluster_id) | cluster_id == "-",
                                       NA_character_, role_factor))) +
    ggplot2::geom_tile(color = "white", linewidth = 0.7, height = 0.62) +
    ggplot2::geom_text(ggplot2::aes(label = cluster_id), size = 2.1,
                        colour = "white", fontface = "bold") +
    ggplot2::scale_fill_manual(values = pal_role, na.value = "grey90",
                                guide = "none") +
    ggplot2::scale_x_discrete(labels = function(x) wrap_genome(x)) +
    ggplot2::labs(title = "C . Primer cluster sharing",
                   subtitle = "role x genome . cell color = role . text = cluster ID",
                   x = NULL, y = NULL) +
    th + ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 0, hjust = 0.5, size = 7.2, lineheight = 0.85))
}

# =============================================================================
# Internal: panel D (unique oligos to order)
# =============================================================================
.shared_viz_panel_D <- function(cdp, pc, pal_role) {
  `%>%` <- magrittr::`%>%`
  if (base::is.null(cdp) || base::nrow(cdp) == 0)
    base::stop("construct_deletion_primers missing — rerun the design with combined construct output enabled.")

  role_long <- dplyr::bind_rows(
    cdp %>% dplyr::transmute(role = "UF", cluster_id = upstream_forward_name,
                              construct_label,
                              full_primer = upstream_forward_primer,
                              full_tm = upstream_forward_tm_full),
    cdp %>% dplyr::transmute(role = "UR", cluster_id = upstream_reverse_name,
                              construct_label,
                              full_primer = upstream_reverse_primer,
                              full_tm = upstream_reverse_tm_full),
    cdp %>% dplyr::transmute(role = "DF", cluster_id = downstream_forward_name,
                              construct_label,
                              full_primer = downstream_forward_primer,
                              full_tm = downstream_forward_tm_full),
    cdp %>% dplyr::transmute(role = "DR", cluster_id = downstream_reverse_name,
                              construct_label,
                              full_primer = downstream_reverse_primer,
                              full_tm = downstream_reverse_tm_full))
  n_per_genome <- pc %>% dplyr::select(role, cluster_id, used_by) %>%
    dplyr::mutate(n_genomes = base::vapply(
      base::strsplit(used_by, ", ", fixed = TRUE), base::length, integer(1)))
  full_dedup <- role_long %>%
    dplyr::left_join(n_per_genome, by = c("role", "cluster_id")) %>%
    dplyr::group_by(role, full_primer) %>%
    dplyr::summarise(
      cluster_id = base::paste(base::sort(base::unique(cluster_id)),
                                collapse = "+"),
      n_constructs = dplyr::n_distinct(construct_label),
      n_genomes = base::max(n_genomes, na.rm = TRUE),
      full_tm = dplyr::first(full_tm),
      .groups = "drop") %>%
    dplyr::mutate(role = base::factor(role, levels = c("UF","UR","DF","DR"))) %>%
    dplyr::arrange(role)

  full_dedup <- full_dedup %>%
    dplyr::left_join(pc %>% dplyr::select(role, cluster_id,
                                           core = primer_core_5to3,
                                           tm_core_val = tm_core),
                      by = c("role", "cluster_id")) %>%
    dplyr::mutate(
      core = dplyr::coalesce(core,
        base::vapply(cluster_id, function(k) {
          ks <- base::strsplit(k, "+", fixed = TRUE)[[1]]
          m <- pc$primer_core_5to3[pc$cluster_id %in% ks]
          if (base::length(m) == 0) NA_character_ else m[1]
        }, character(1))),
      tm_core_val = dplyr::coalesce(tm_core_val,
        base::vapply(cluster_id, function(k) {
          ks <- base::strsplit(k, "+", fixed = TRUE)[[1]]
          m <- pc$tm_core[pc$cluster_id %in% ks]
          if (base::length(m) == 0) NA_real_ else m[1]
        }, numeric(1))),
      full_len     = base::nchar(full_primer),
      core_len     = base::nchar(core),
      overhang_len = full_len - core_len,
      overhang_seq = base::substr(full_primer, 1, overhang_len),
      core_seq     = base::substr(full_primer, overhang_len + 1, full_len))

  full_dedup <- full_dedup %>% dplyr::arrange(role) %>%
    dplyr::mutate(row_y = base::rev(base::seq(from = 1, by = 2.2,
                                               length.out = dplyr::n())))
  full_dedup <- full_dedup %>%
    dplyr::mutate(y_label = base::sprintf(
      "%s<br><span style='color:#555555;font-size:6.5pt;'>n=%d genomes . %d constructs</span><br><span style='color:#111111;font-size:7pt;'>**core: %d bp . %.1f C**</span><br><span style='color:#555555;font-size:6.5pt;'>full: %d bp . %.1f C</span>",
      cluster_id, n_genomes, n_constructs,
      core_len, tm_core_val, full_len, full_tm))

  overhang_text_pal <- c(UF = "#99BFDC", UR = "#EDB699",
                         DF = "#B9D9F1", DR = "#F2CC80")
  full_dedup <- full_dedup %>%
    dplyr::mutate(seq_display = base::sprintf(
      "5'- <span style='color:%s;'>%s</span><span style='color:white;'>**%s**</span> -3'",
      overhang_text_pal[base::as.character(role)], overhang_seq, core_seq))

  bar_xmax <- 1
  ggplot2::ggplot(full_dedup) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = 0, xmax = bar_xmax,
                   ymin = row_y - 0.30, ymax = row_y + 0.30,
                   fill = base::as.character(role)),
      colour = "grey20", linewidth = 0.18) +
    ggtext::geom_richtext(
      ggplot2::aes(x = bar_xmax / 2, y = row_y, label = seq_display),
      family = "mono", size = 2.4, colour = "white",
      fill = NA, label.color = NA,
      label.padding = grid::unit(0, "pt")) +
    ggplot2::scale_fill_manual(values = pal_role, guide = "none") +
    ggplot2::scale_y_continuous(breaks = full_dedup$row_y,
                                 labels = full_dedup$y_label,
                                 expand = ggplot2::expansion(add = 0.45)) +
    ggplot2::scale_x_continuous(limits = c(0, bar_xmax),
                                 expand = c(0, 0)) +
    ggplot2::labs(
      title = base::sprintf("D . Unique oligos to order (%d total)",
                             base::nrow(full_dedup)),
      subtitle = "uniform bars . 5' to 3' . blended text = 5' overhang . white bold = genome-binding core",
      x = NULL, y = NULL) +
    # theme_grey base: element_markdown dispatches correctly on ggplot2 >= 4.0
    ggplot2::theme_grey(base_size = 9.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10.5,
                                          colour = "grey10",
                                          margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(colour = "grey45", size = 8.5,
                                             margin = ggplot2::margin(b = 4)),
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(6, 10, 6, 6),
      axis.text.y = ggtext::element_markdown(size = 7.5, lineheight = 1.25,
                                              colour = "grey25", hjust = 1))
}

# =============================================================================
# Internal: panels E/F (arm alignments + dendrogram)
# =============================================================================
.shared_viz_panels_EF <- function(tc, pc, th, pal_track) {
  `%>%` <- magrittr::`%>%`
  wrap_genome <- function(x) base::sub("^([^_]+)_(.*)$", "\\1\n\\2", x)

  genome_ids <- tc$genome_id
  left_arms  <- stats::setNames(tc$effective_left_arm_seq,  genome_ids)
  right_arms <- stats::setNames(tc$effective_right_arm_seq, genome_ids)
  min_l <- base::min(base::nchar(left_arms))
  min_r <- base::min(base::nchar(right_arms))
  left_aligned  <- base::sapply(left_arms,
    function(s) base::substr(s, base::nchar(s) - min_l + 1, base::nchar(s)))
  right_aligned <- base::sapply(right_arms,
    function(s) base::substr(s, 1, min_r))

  to_long <- function(aligned_vec, arm_tag) {
    m <- base::do.call(base::rbind, base::strsplit(aligned_vec, "", fixed = TRUE))
    base::rownames(m) <- base::names(aligned_vec)
    df <- base::as.data.frame(m, stringsAsFactors = FALSE)
    df$genome <- base::rownames(df)
    long <- tidyr::pivot_longer(df, cols = -genome, names_to = "pos",
                                 values_to = "base")
    long$pos <- base::as.integer(base::sub("^V", "", long$pos))
    long$arm <- arm_tag; long
  }
  base_long <- dplyr::bind_rows(to_long(left_aligned,  "Upstream arm"),
                                 to_long(right_aligned, "Downstream arm"))
  consensus <- base_long %>% dplyr::group_by(arm, pos) %>%
    dplyr::summarise(
      cons = base::names(base::sort(base::table(base), decreasing = TRUE))[1],
      .groups = "drop")
  base_long <- base_long %>% dplyr::left_join(consensus, by = c("arm", "pos")) %>%
    dplyr::mutate(match = base == cons)
  base_long$cell_fill <- base::ifelse(base_long$match, "match", base_long$base)

  # Hamming clustering
  hamming_mat <- function(mat) {
    n <- base::nrow(mat)
    d <- base::matrix(0, n, n,
                      dimnames = base::list(base::rownames(mat),
                                             base::rownames(mat)))
    for (i in base::seq_len(n - 1)) for (j in (i + 1):n) {
      diff <- base::sum(mat[i, ] != mat[j, ])
      d[i, j] <- diff; d[j, i] <- diff
    }
    stats::as.dist(d)
  }
  left_mat  <- base::do.call(base::rbind, base::strsplit(left_aligned, "", fixed = TRUE))
  right_mat <- base::do.call(base::rbind, base::strsplit(right_aligned, "", fixed = TRUE))
  base::rownames(left_mat)  <- base::names(left_aligned)
  base::rownames(right_mat) <- base::names(right_aligned)
  combined_mat <- base::cbind(left_mat, right_mat)
  hc <- stats::hclust(hamming_mat(combined_mat), method = "average")
  clust_order <- hc$labels[hc$order]
  clust_order_wrapped <- wrap_genome(clust_order)
  base_long$genome <- base::factor(wrap_genome(base_long$genome),
                                    levels = base::rev(clust_order_wrapped))

  cell_pal <- c(match = "#ECECEC",
                A = "#00A087", C = "#3C5488",
                G = "#F39B7F", `T` = "#E64B35",
                N = "#9AA0A6")

  pal_role <- c(UF = "#0072B2", UR = "#D55E00", DF = "#56B4E9", DR = "#E69F00")
  get_primer_positions <- function(arm_tag, aligned_len, arms_original, role_df) {
    result <- base::list()
    for (rid in base::seq_len(base::nrow(role_df))) {
      core <- role_df$shared_core[rid]
      search_seq <- if (!base::is.na(core) && base::nchar(core) > 0) core
                    else role_df$primer_core_5to3[rid]
      members <- base::strsplit(role_df$used_by[rid], ", ", fixed = TRUE)[[1]]
      for (gid in members) {
        arm <- arms_original[[gid]]; if (base::is.null(arm)) next
        hits <- base::unlist(base::gregexpr(search_seq, arm, fixed = TRUE))
        hits <- hits[hits > 0]
        for (h in hits) {
          ap <- if (arm_tag == "Upstream arm") h - (base::nchar(arm) - aligned_len) else h
          if (ap < 1 || ap + base::nchar(search_seq) - 1 > aligned_len) next
          result[[base::length(result) + 1]] <- base::data.frame(
            arm = arm_tag, genome = gid, role = role_df$role[rid],
            cluster_id = role_df$cluster_id[rid],
            start = ap, end = ap + base::nchar(search_seq) - 1,
            stringsAsFactors = FALSE)
        }
      }
    }
    base::do.call(base::rbind, result)
  }
  spans_left  <- get_primer_positions("Upstream arm",  min_l, left_arms,
                                       pc[pc$role %in% c("UF","UR"), ])
  spans_right <- get_primer_positions("Downstream arm", min_r, right_arms,
                                       pc[pc$role %in% c("DF","DR"), ])
  primer_spans <- dplyr::bind_rows(spans_left, spans_right) %>%
    dplyr::group_by(arm, role, cluster_id) %>%
    dplyr::summarise(start = base::min(start), end = base::max(end),
                      .groups = "drop")

  make_arm <- function(arm_tag, panel_tag, show_y = TRUE) {
    df <- base_long %>% dplyr::filter(arm == arm_tag)
    spans <- primer_spans %>% dplyr::filter(arm == arm_tag)
    nlev <- base::length(base::levels(df$genome))
    arm_col <- if (arm_tag == "Upstream arm") pal_track[["Upstream arm"]]
               else pal_track[["Downstream arm"]]
    p <- ggplot2::ggplot(df) +
      ggplot2::geom_tile(
        ggplot2::aes(x = pos, y = genome, fill = cell_fill), height = 0.62) +
      ggplot2::scale_fill_manual(values = cell_pal,
        breaks = c("match", "A", "C", "G", "T"),
        labels = c("match", "A", "C", "G", "T"),
        name = NULL)
    for (i in base::seq_len(base::nrow(spans))) {
      sp <- spans[i, ]
      role_col <- pal_role[[sp$role]]
      p <- p +
        ggplot2::annotate("rect",
          xmin = sp$start - 0.5, xmax = sp$end + 0.5,
          ymin = 0.5, ymax = nlev + 1.05,
          fill = role_col, alpha = 0.08) +
        ggplot2::annotate("rect",
          xmin = sp$start, xmax = sp$end,
          ymin = nlev + 0.55, ymax = nlev + 1.05,
          fill = role_col, colour = "grey25", linewidth = 0.15) +
        ggplot2::annotate("text",
          x = (sp$start + sp$end) / 2, y = nlev + 0.80,
          label = sp$role,
          size = 3.2, colour = "white", fontface = "bold")
    }
    p +
      ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.2, 1.15))) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.005)) +
      ggplot2::labs(
        title = base::sprintf("%s . %s alignment (%d bp)", panel_tag, arm_tag,
          if (arm_tag == "Upstream arm") min_l else min_r),
        subtitle = base::sprintf("grey = match . coloured = mismatch . tinted column = primer binding (%s)",
          if (arm_tag == "Upstream arm") "UF / UR" else "DF / DR"),
        x = base::sprintf("aligned bp - %s",
          if (arm_tag == "Upstream arm") "inner boundary at right edge"
          else "inner boundary at left edge"),
        y = NULL) +
      th +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(colour = arm_col, fill = NA,
                                              linewidth = 1.4),
        plot.title = ggplot2::element_text(face = "bold", size = 10.5,
                                            colour = arm_col,
                                            margin = ggplot2::margin(b = 3)),
        axis.text.y = if (show_y) ggplot2::element_text(
                        colour = "grey30", lineheight = 0.85, size = 7.5)
                       else ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        plot.margin = if (show_y) ggplot2::margin(6, 4, 6, 6)
                      else ggplot2::margin(6, 6, 6, 2))
  }
  p_left  <- make_arm("Upstream arm",   "E", show_y = TRUE)
  p_right <- make_arm("Downstream arm", "F", show_y = FALSE)

  # Dendrogram
  dd <- ggdendro::dendro_data(stats::as.dendrogram(hc), type = "rectangle")
  segs <- dd$segments
  n_leaves <- base::length(clust_order)
  dx_to_row <- function(x) n_leaves - x + 1
  segs$y1    <- dx_to_row(segs$x)
  segs$yend1 <- dx_to_row(segs$xend)
  bmax <- base::max(segs$y, segs$yend, 1)
  segs$x1    <- segs$y / bmax
  segs$xend1 <- segs$yend / bmax

  p_dendro <- ggplot2::ggplot(segs) +
    ggplot2::geom_segment(
      ggplot2::aes(x = x1, xend = xend1, y = y1, yend = yend1),
      colour = "grey25", linewidth = 0.55, lineend = "round") +
    ggplot2::scale_x_reverse(
      expand = ggplot2::expansion(mult = c(0.05, 0.02))) +
    ggplot2::scale_y_continuous(limits = c(0.5, n_leaves + 1.3),
                                 expand = c(0, 0)) +
    ggplot2::labs(title = NULL, x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 7.5) +
    ggplot2::theme(plot.margin = ggplot2::margin(34, 1, 48, 1))

  base::list(p_dendro = p_dendro, p_left = p_left, p_right = p_right,
             min_l = min_l, min_r = min_r,
             clust_order_wrapped = clust_order_wrapped)
}

# =============================================================================
# Internal: panel G (flanking gene context)
# =============================================================================
.shared_viz_panel_G <- function(tc, gbk_paths, th, pal_track,
                                 clust_order_wrapped) {
  `%>%` <- magrittr::`%>%`
  wrap_genome <- function(x) base::sub("^([^_]+)_(.*)$", "\\1\n\\2", x)

  flank_df <- .shared_viz_flanking(tc, gbk_paths)

  flank_df <- flank_df %>% dplyr::mutate(
    us_rel_start  = us_start  - tgt_start,
    us_rel_end    = us_end    - tgt_start,
    tgt_rel_end   = tgt_end   - tgt_start,
    ds_rel_start  = ds_start  - tgt_start,
    ds_rel_end    = ds_end    - tgt_start,
    la_rel_start  = eff_left_start  - tgt_start,
    la_rel_end    = eff_left_end    - tgt_start,
    ra_rel_start  = eff_right_start - tgt_start,
    ra_rel_end    = eff_right_end   - tgt_start)

  flank_df$genome <- base::factor(wrap_genome(flank_df$genome_id),
                                   levels = base::rev(clust_order_wrapped))

  short_prod <- function(x, n = 28) {
    x <- base::ifelse(base::is.na(x), "-", x)
    base::ifelse(base::nchar(x) > n,
                  base::paste0(base::substr(x, 1, n - 1), "..."), x)
  }
  xlim_g <- c(-2500, base::max(flank_df$ra_rel_end, flank_df$ds_rel_end,
                                 na.rm = TRUE) + 300)

  arm_rects <- dplyr::bind_rows(
    base::data.frame(genome = flank_df$genome,
      xmin = flank_df$la_rel_start, xmax = flank_df$la_rel_end,
      arm = "Upstream arm", stringsAsFactors = FALSE),
    base::data.frame(genome = flank_df$genome,
      xmin = flank_df$ra_rel_start, xmax = flank_df$ra_rel_end,
      arm = "Downstream arm", stringsAsFactors = FALSE))

  OV_WIDTH <- 80
  gene_segs <- .shared_viz_build_gene_segs(flank_df, OV_WIDTH, short_prod)
  gene_labels <- gene_segs %>% dplyr::filter(!is_ov) %>%
    dplyr::mutate(label_x_clip = base::pmax(xmin, xlim_g[1] + 20) +
      (base::pmin(xmax, xlim_g[2] - 20) - base::pmax(xmin, xlim_g[1] + 20)) / 2)

  gene_pal <- c(upstream = "#4A6FA5", target = "#C9364D", downstream = "#C47B3F")

  backbone_df <- base::data.frame(genome = base::levels(flank_df$genome),
                                   stringsAsFactors = FALSE)

  ggplot2::ggplot() +
    ggplot2::geom_rect(data = backbone_df,
      ggplot2::aes(xmin = xlim_g[1], xmax = xlim_g[2],
                   ymin = base::match(genome, base::levels(flank_df$genome)) - 0.18,
                   ymax = base::match(genome, base::levels(flank_df$genome)) + 0.18),
      fill = "#E8E8E8", colour = NA) +
    ggplot2::geom_rect(data = gene_segs,
      ggplot2::aes(xmin = xmin, xmax = xmax,
                   ymin = base::as.numeric(genome) + y_lo,
                   ymax = base::as.numeric(genome) + y_hi,
                   fill = role),
      colour = "grey25", linewidth = 0.12) +
    ggplot2::geom_text(data = gene_labels,
      ggplot2::aes(x = label_x_clip, y = base::as.numeric(genome), label = label),
      size = 1.95, colour = "white", fontface = "bold") +
    ggplot2::geom_rect(data = arm_rects,
      ggplot2::aes(xmin = xmin, xmax = xmax,
                   ymin = base::as.numeric(genome) - 0.36,
                   ymax = base::as.numeric(genome) - 0.22,
                   fill = arm), colour = NA) +
    ggplot2::geom_segment(
      data = dplyr::distinct(flank_df, genome, tgt_rel_end),
      ggplot2::aes(x = 0, xend = 0,
                   y = base::as.numeric(genome) - 0.36,
                   yend = base::as.numeric(genome) + 0.22),
      colour = "grey15", linewidth = 0.3, linetype = "dashed") +
    ggplot2::geom_segment(
      data = dplyr::distinct(flank_df, genome, tgt_rel_end),
      ggplot2::aes(x = tgt_rel_end, xend = tgt_rel_end,
                   y = base::as.numeric(genome) - 0.36,
                   yend = base::as.numeric(genome) + 0.22),
      colour = "grey15", linewidth = 0.3, linetype = "dashed") +
    ggplot2::scale_y_continuous(
      breaks = base::seq_along(base::levels(flank_df$genome)),
      labels = base::levels(flank_df$genome),
      expand = ggplot2::expansion(add = 0.4)) +
    ggplot2::scale_fill_manual(
      values = c(gene_pal,
                 "Upstream arm"   = pal_track[["Upstream arm"]],
                 "Downstream arm" = pal_track[["Downstream arm"]]),
      breaks = c("upstream", "target", "downstream",
                 "Upstream arm", "Downstream arm"),
      labels = c("upstream CDS", "target", "downstream CDS",
                 "upstream arm", "downstream arm"),
      name = NULL, guide = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::coord_cartesian(xlim = xlim_g) +
    ggplot2::labs(title = "G . Flanking gene context",
      subtitle = base::sprintf(
        "grey bar = continuous DNA . colored box = CDS . half-height step at junction = CDS overlap (drawn as %d bp for legibility) . thin strip = effective arm . dashed = deletion boundary",
        OV_WIDTH),
      x = "bp from target_start", y = NULL) +
    th + ggplot2::theme(axis.text.y = ggplot2::element_text(
      lineheight = 0.85, size = 7.5))
}

# ---- Flanking-gene extraction (pure R, per-genome CDS-only) -----------------
.shared_viz_flanking <- function(tc, gbk_paths) {
  parse_cds <- function(path) {
    lines <- base::readLines(path, warn = FALSE)
    feat_line <- base::grep("^FEATURES", lines)[1]
    origin_line <- base::grep("^ORIGIN|^CONTIG", lines)
    origin_line <- origin_line[origin_line > feat_line][1]
    if (base::is.na(origin_line)) origin_line <- base::length(lines) + 1L
    blk <- lines[(feat_line + 1):(origin_line - 1)]
    is_feat_header <- base::grepl("^     \\S", blk) & !base::grepl("^     /", blk)
    idxs <- base::which(is_feat_header)
    idxs2 <- base::c(idxs, base::length(blk) + 1L)
    rows <- base::list()
    for (k in base::seq_along(idxs)) {
      i_start <- idxs2[k]; i_end <- idxs2[k + 1] - 1L
      header <- blk[i_start]
      m <- base::regmatches(header,
        base::regexec("^     (\\S+)\\s+(.*)$", header))[[1]]
      type <- m[2]; loc <- m[3]
      if (base::is.na(type) || type != "CDS") next
      sub_lines <- blk[i_start:i_end]
      q_idx <- base::grep("^                     /", sub_lines)
      loc_block <- sub_lines[1]
      if (base::length(q_idx) > 0 && q_idx[1] > 2) {
        cont <- sub_lines[2:(q_idx[1] - 1)]
        loc_block <- base::paste0(loc_block,
                                   base::paste0(base::trimws(cont), collapse = ""))
      } else if (base::length(q_idx) == 0 && base::length(sub_lines) > 1) {
        loc_block <- base::paste0(loc_block,
          base::paste0(base::trimws(sub_lines[-1]), collapse = ""))
      }
      pick <- function(key) {
        q_idx2 <- base::c(q_idx, base::length(sub_lines) + 1L)
        for (qi in base::seq_along(q_idx)) {
          qs <- q_idx[qi]; qe <- q_idx2[qi + 1] - 1L
          qtext <- base::paste(base::trimws(sub_lines[qs:qe]), collapse = " ")
          m2 <- base::regmatches(qtext,
            base::regexec(base::sprintf('^/%s=(?:"([^"]*)"|(\\S+))', key),
                           qtext))[[1]]
          if (base::length(m2) > 0) {
            v <- if (base::nzchar(m2[2])) m2[2] else m2[3]
            if (!base::is.na(v) && base::nzchar(v)) return(v)
          }
        }
        NA_character_
      }
      locus_tag <- if (base::length(q_idx) > 0) pick("locus_tag") else NA_character_
      gene      <- if (base::length(q_idx) > 0) pick("gene")      else NA_character_
      product   <- if (base::length(q_idx) > 0) pick("product")   else NA_character_
      strand <- if (base::grepl("^complement", loc_block)) -1L else 1L
      nums <- base::as.integer(base::unlist(base::regmatches(
        loc_block, base::gregexpr("\\d+", loc_block))))
      if (base::length(nums) < 2) next
      rows[[base::length(rows) + 1]] <- base::data.frame(
        start = nums[1], end = utils::tail(nums, 1), strand = strand,
        locus_tag = locus_tag, gene = gene, product = product,
        stringsAsFactors = FALSE)
    }
    base::do.call(base::rbind, rows)
  }

  flank_rows <- base::list()
  for (gid in tc$genome_id) {
    cds <- parse_cds(gbk_paths[[gid]])
    tcr <- tc[tc$genome_id == gid, ]
    ds <- cds[cds$start > tcr$target_end, , drop = FALSE]
    ds <- ds[base::order(ds$start), , drop = FALSE]
    ds_row <- if (base::nrow(ds) > 0) ds[1, ] else
      base::data.frame(start = NA_integer_, end = NA_integer_, strand = 0L,
                        locus_tag = NA_character_, gene = NA_character_,
                        product = NA_character_)
    us <- cds[cds$locus_tag == tcr$upstream_locus_tag &
              !base::is.na(cds$locus_tag), , drop = FALSE]
    us_strand <- if (base::nrow(us) > 0) us$strand[1] else NA_integer_
    flank_rows[[gid]] <- base::data.frame(
      genome_id = gid,
      us_start = tcr$upstream_start, us_end = tcr$upstream_end,
      us_strand = us_strand,
      us_locus = tcr$upstream_locus_tag,
      us_gene = if (base::nrow(us) > 0) us$gene[1] else NA_character_,
      us_product = tcr$upstream_product,
      us_overlap = tcr$upstream_overlap_bp,
      tgt_start = tcr$target_start, tgt_end = tcr$target_end,
      tgt_strand = base::ifelse(tcr$target_strand == "+", 1L, -1L),
      tgt_locus = tcr$locus_tag, tgt_product = tcr$target_product,
      ds_start = ds_row$start, ds_end = ds_row$end, ds_strand = ds_row$strand,
      ds_locus = ds_row$locus_tag, ds_gene = ds_row$gene,
      ds_product = ds_row$product,
      eff_left_start = tcr$effective_left_arm_start,
      eff_left_end = tcr$effective_left_arm_end,
      eff_right_start = tcr$effective_right_arm_start,
      eff_right_end = tcr$effective_right_arm_end,
      deletion_start = tcr$deletion_start, deletion_end = tcr$deletion_end,
      stringsAsFactors = FALSE)
  }
  base::do.call(base::rbind, flank_rows)
}

# ---- half-height step builder (for CDS overlap visualisation) ---------------
.shared_viz_build_gene_segs <- function(fd, OV_WIDTH, short_prod) {
  out <- base::list()
  for (i in base::seq_len(base::nrow(fd))) {
    r <- fd[i, ]
    us_xmax  <- r$us_rel_end
    tgt_xmin <- 0
    tgt_xmax <- r$tgt_rel_end
    ds_xmin  <- r$ds_rel_start
    ds_xmax  <- r$ds_rel_end
    has_us_tgt <- !base::is.na(us_xmax) && us_xmax > tgt_xmin
    has_tgt_ds <- !base::is.na(ds_xmin) && !base::is.na(tgt_xmax) &&
                  ds_xmin < tgt_xmax
    if (has_us_tgt) {
      ov_mid <- (tgt_xmin + us_xmax) / 2
      ov_lo  <- ov_mid - OV_WIDTH / 2
      ov_hi  <- ov_mid + OV_WIDTH / 2
      out[[base::length(out)+1]] <- base::data.frame(
        genome = r$genome, role = "upstream",
        xmin = r$us_rel_start, xmax = ov_lo,
        y_lo = -0.18, y_hi = 0.18, is_ov = FALSE,
        label = base::sprintf("%s . %s", r$us_locus, short_prod(r$us_product)),
        label_x = (r$us_rel_start + ov_lo) / 2,
        stringsAsFactors = FALSE)
      out[[base::length(out)+1]] <- base::data.frame(
        genome = r$genome, role = "upstream",
        xmin = ov_lo, xmax = ov_hi, y_lo = -0.18, y_hi = 0, is_ov = TRUE,
        label = NA_character_, label_x = NA_real_, stringsAsFactors = FALSE)
      tgt_start_eff <- ov_hi
      tgt_ov_step <- base::data.frame(
        genome = r$genome, role = "target",
        xmin = ov_lo, xmax = ov_hi, y_lo = 0, y_hi = 0.18, is_ov = TRUE,
        label = NA_character_, label_x = NA_real_, stringsAsFactors = FALSE)
    } else {
      tgt_start_eff <- tgt_xmin; tgt_ov_step <- NULL
    }
    if (has_tgt_ds) {
      ov_mid <- (ds_xmin + tgt_xmax) / 2
      ov_lo  <- ov_mid - OV_WIDTH / 2
      ov_hi  <- ov_mid + OV_WIDTH / 2
      tgt_end_eff <- ov_lo
      tgt_ov_step2 <- base::data.frame(
        genome = r$genome, role = "target",
        xmin = ov_lo, xmax = ov_hi, y_lo = 0, y_hi = 0.18, is_ov = TRUE,
        label = NA_character_, label_x = NA_real_, stringsAsFactors = FALSE)
      ds_ov_step <- base::data.frame(
        genome = r$genome, role = "downstream",
        xmin = ov_lo, xmax = ov_hi, y_lo = -0.18, y_hi = 0, is_ov = TRUE,
        label = NA_character_, label_x = NA_real_, stringsAsFactors = FALSE)
      ds_start_eff <- ov_hi
    } else {
      tgt_end_eff <- tgt_xmax; tgt_ov_step2 <- NULL
      ds_ov_step <- NULL; ds_start_eff <- ds_xmin
    }
    out[[base::length(out)+1]] <- base::data.frame(
      genome = r$genome, role = "target",
      xmin = tgt_start_eff, xmax = tgt_end_eff,
      y_lo = -0.18, y_hi = 0.18, is_ov = FALSE,
      label = base::sprintf("%s . %s", r$tgt_locus, short_prod(r$tgt_product)),
      label_x = (tgt_start_eff + tgt_end_eff) / 2,
      stringsAsFactors = FALSE)
    if (!base::is.null(tgt_ov_step))  out[[base::length(out)+1]] <- tgt_ov_step
    if (!base::is.null(tgt_ov_step2)) out[[base::length(out)+1]] <- tgt_ov_step2
    if (!base::is.na(ds_start_eff) && !base::is.na(ds_xmax)) {
      out[[base::length(out)+1]] <- base::data.frame(
        genome = r$genome, role = "downstream",
        xmin = ds_start_eff, xmax = ds_xmax,
        y_lo = -0.18, y_hi = 0.18, is_ov = FALSE,
        label = base::sprintf("%s . %s", r$ds_locus, short_prod(r$ds_product)),
        label_x = (ds_start_eff + ds_xmax) / 2,
        stringsAsFactors = FALSE)
    }
    if (!base::is.null(ds_ov_step)) out[[base::length(out)+1]] <- ds_ov_step
  }
  base::do.call(base::rbind, out)
}

# Avoid R CMD check NOTE for unbound globals.
utils::globalVariables(c(
  "construct_label", "type", "label", "category", "fill_cat", "role",
  "start", "end", "short", "cluster_id", "n_genomes", "n_constructs",
  "role_factor", "used_by", "genome", "cons", "base", "match",
  "cell_fill", "pos", "arm", "x1", "xend1", "y1", "yend1",
  "upstream_forward_name", "upstream_forward_primer", "upstream_forward_tm_full",
  "upstream_reverse_name", "upstream_reverse_primer", "upstream_reverse_tm_full",
  "downstream_forward_name", "downstream_forward_primer", "downstream_forward_tm_full",
  "downstream_reverse_name", "downstream_reverse_primer", "downstream_reverse_tm_full",
  "full_primer", "full_tm", "primer_core_5to3", "tm_core", "tm_core_val",
  "full_len", "core_len", "overhang_len", "overhang_seq", "core_seq",
  "row_y", "y_label", "seq_display", "total_len", "us_rel_start", "us_rel_end",
  "tgt_rel_end", "ds_rel_start", "ds_rel_end", "la_rel_start", "la_rel_end",
  "ra_rel_start", "ra_rel_end", "xmin", "xmax", "y_lo", "y_hi", "is_ov",
  "label_x_clip", "us_start", "us_end", "ds_start", "ds_end",
  "eff_left_start", "eff_left_end", "eff_right_start", "eff_right_end",
  "tgt_start", "tgt_end", "target_end", "construct_files"
))
