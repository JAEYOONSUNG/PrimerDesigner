#' Calculate Composite Score for gRNAs
#'
#' Calculates a composite quality score for each gRNA based on multiple criteria
#' including GC content, off-target specificity, melting temperature, position,
#' and sequence structure quality (palindrome/hairpin avoidance).
#'
#' @param gRNA_df A data frame containing gRNA data with columns: percentGC, Tm, n0, n1, position_percent.
#'   Optionally includes palindrome_score and hairpin_score (from pre-filtering step).
#' @param weights A named numeric vector specifying the weight for each scoring component.
#'   Default: c(gc = 0.25, offtarget = 0.35, tm = 0.10, position = 0.10, structure = 0.20).
#' @param optimal_gc Numeric. The optimal GC content percentage (default: 50).
#' @param optimal_tm Numeric. The optimal melting temperature in Celsius (default: 60).
#' @param optimal_position Numeric. The optimal position percentage in the target (default: 50).
#' @return The input data frame with an added \code{composite_score} column (0-100 scale).
#' @examples
#' \dontrun{
#' gRNA_df <- data.frame(
#'   protospacer = c("ATCG...", "GCTA..."),
#'   percentGC = c(48, 62),
#'   Tm = c(59, 65),
#'   n0 = c(1, 2),
#'   n1 = c(0, 1),
#'   position_percent = c(45, 75)
#' )
#' scored_df <- calculate_composite_score(gRNA_df)
#' }
#' @export
calculate_composite_score <- function(
    gRNA_df,
    weights = base::c(gc = 0.25, offtarget = 0.35, tm = 0.10, position = 0.10,
                      structure = 0.20),
    optimal_gc = 50,
    optimal_tm = 60,
    optimal_position = 50
) {
  # Validate inputs
  if (!base::is.data.frame(gRNA_df) || base::nrow(gRNA_df) == 0) {
    base::warning("Empty or invalid gRNA data frame. Returning as-is.")
    gRNA_df$composite_score <- base::numeric(0)
    return(gRNA_df)
  }

  required_cols <- base::c("percentGC", "Tm", "n0", "n1", "position_percent")
  missing_cols <- base::setdiff(required_cols, base::colnames(gRNA_df))
  if (base::length(missing_cols) > 0) {
    base::stop("Missing required columns for composite scoring: ", base::paste(missing_cols, collapse = ", "))
  }

  # Normalize weights
  weights <- weights / base::sum(weights)

  # --- GC Score (0-1): closeness to optimal GC ---
  gc_score <- base::pmax(0, 1 - base::abs(gRNA_df$percentGC - optimal_gc) / 50)

  # --- Off-target Score (0-1) ---
  n0_score <- base::ifelse(gRNA_df$n0 == 1, 1.0,
              base::ifelse(gRNA_df$n0 == 0, 0.3,
              base::pmax(0, 1 - (gRNA_df$n0 - 1) * 0.25)))
  n1_score <- base::pmax(0, 1 - gRNA_df$n1 * 0.2)
  offtarget_score <- 0.6 * n0_score + 0.4 * n1_score

  # --- Tm Score (0-1): closeness to optimal Tm ---
  tm_score <- base::pmax(0, 1 - base::abs(gRNA_df$Tm - optimal_tm) / 20)

  # --- Position Score (0-1): closeness to center of target ---
  position_score <- base::pmax(0, 1 - base::abs(gRNA_df$position_percent - optimal_position) / 50)

  # --- Structure Score (0-1): penalty for palindrome & hairpin ---
  # Lower palindrome_score and hairpin_score = better (less self-complementary)
  if ("palindrome_score" %in% base::colnames(gRNA_df) &&
      "hairpin_score" %in% base::colnames(gRNA_df)) {
    # Invert: 0 palindrome/hairpin → structure_score = 1.0 (best)
    palin_penalty <- base::pmax(0, 1 - gRNA_df$palindrome_score * 2.5)
    hairpin_penalty <- base::pmax(0, 1 - gRNA_df$hairpin_score * 2.0)
    structure_score <- 0.5 * palin_penalty + 0.5 * hairpin_penalty
  } else {
    # If columns absent, assume no penalty (backward compatibility)
    structure_score <- base::rep(1, base::nrow(gRNA_df))
  }

  # --- Composite Score (0-100) ---
  # Ensure all weight names exist
  w_gc <- base::ifelse("gc" %in% base::names(weights), weights["gc"], 0.25)
  w_ot <- base::ifelse("offtarget" %in% base::names(weights), weights["offtarget"], 0.35)
  w_tm <- base::ifelse("tm" %in% base::names(weights), weights["tm"], 0.10)
  w_pos <- base::ifelse("position" %in% base::names(weights), weights["position"], 0.10)
  w_str <- base::ifelse("structure" %in% base::names(weights), weights["structure"], 0.20)

  composite <- (w_gc * gc_score +
                w_ot * offtarget_score +
                w_tm * tm_score +
                w_pos * position_score +
                w_str * structure_score) * 100

  gRNA_df$composite_score <- base::round(base::as.numeric(composite), 2)

  return(gRNA_df)
}


#' Generate Systematic Names for gRNAs
#'
#' Assigns human-readable, systematic names to each gRNA in a data frame.
#' Names follow the format: \code{{prefix}_{gene_or_tag}_g{rank}}.
#' Rank is determined within each locus_tag group, ordered by composite_score (descending).
#'
#' @param gRNA_df A data frame containing gRNA data. Must have \code{locus_tag} column.
#'   Optionally includes \code{gene} and \code{composite_score} columns.
#' @param prefix A character string to prepend to each gRNA name.
#'   Default: NULL (auto-detected from nuclease abbreviation if available).
#'   Common prefixes: "Geo" (GeoCas9), "Fn" (FnCas12a), "Fis" (FisCasI_B), "Sp" (SpCas9).
#' @param nuclease A character string specifying the nuclease name for auto-prefix.
#'   Used only if \code{prefix} is NULL. Default: NULL.
#' @param rank_by A character string specifying the column to rank gRNAs by within each gene.
#'   Default: "composite_score".
#' @return The input data frame with added \code{gRNA_name} and \code{gRNA_rank} columns,
#'   sorted by locus_tag and rank.
#' @examples
#' \dontrun{
#' gRNA_df <- data.frame(
#'   locus_tag = c("QT234_RS00005", "QT234_RS00005", "QT234_RS00010"),
#'   gene = c("dnaA", "dnaA", "dnaN"),
#'   protospacer = c("ATCG...", "GCTA...", "TTAA..."),
#'   composite_score = c(92.5, 87.3, 95.1)
#' )
#' named_df <- generate_gRNA_names(gRNA_df, prefix = "Geo")
#' # gRNA_name: Geo_dnaA_g1, Geo_dnaA_g2, Geo_dnaN_g1
#' }
#' @export
generate_gRNA_names <- function(
    gRNA_df,
    prefix = NULL,
    nuclease = NULL,
    rank_by = "composite_score"
) {
  # Validate inputs
  if (!base::is.data.frame(gRNA_df) || base::nrow(gRNA_df) == 0) {
    base::warning("Empty or invalid gRNA data frame. Returning as-is.")
    gRNA_df$gRNA_name <- base::character(0)
    gRNA_df$gRNA_rank <- base::integer(0)
    return(gRNA_df)
  }

  if (!"locus_tag" %in% base::colnames(gRNA_df)) {
    base::stop("gRNA data frame must contain 'locus_tag' column.")
  }

  # Auto-detect prefix from nuclease name
  if (base::is.null(prefix)) {
    nuclease_prefix_map <- base::c(
      "GeoCas9"   = "Geo",
      "FnCas12a"  = "Fn",
      "FisCasI_B" = "Fis",
      "SpCas9"    = "Sp",
      "GtCas9"    = "Gt"
    )
    if (!base::is.null(nuclease) && nuclease %in% base::names(nuclease_prefix_map)) {
      prefix <- nuclease_prefix_map[nuclease]
    } else {
      prefix <- "gRNA"
    }
  }

  # Determine gene identifier for naming
  # Priority: gene name > locus_tag (cleaned)
  if ("gene" %in% base::colnames(gRNA_df)) {
    gRNA_df$`.name_id` <- base::ifelse(
      !base::is.na(gRNA_df$gene) & gRNA_df$gene != "",
      gRNA_df$gene,
      base::gsub("[^a-zA-Z0-9]", "", gRNA_df$locus_tag)
    )
  } else {
    gRNA_df$`.name_id` <- base::gsub("[^a-zA-Z0-9]", "", gRNA_df$locus_tag)
  }

  # Determine ranking column
  if (!rank_by %in% base::colnames(gRNA_df)) {
    base::warning("rank_by column '", rank_by, "' not found. Using row order for ranking.")
    gRNA_df$`.rank_value` <- base::seq_len(base::nrow(gRNA_df)) * -1  # preserve original order
  } else {
    gRNA_df$`.rank_value` <- gRNA_df[[rank_by]]
  }

  # Assign rank within each locus_tag group (descending by score)
  gRNA_df <- gRNA_df %>%
    dplyr::group_by(locus_tag) %>%
    dplyr::arrange(dplyr::desc(`.rank_value`), .by_group = TRUE) %>%
    dplyr::mutate(
      gRNA_rank = dplyr::row_number(),
      gRNA_name = base::paste0(prefix, "_", `.name_id`, "_g", gRNA_rank)
    ) %>%
    dplyr::ungroup()

  # Clean up temp columns
  gRNA_df$`.name_id` <- NULL
  gRNA_df$`.rank_value` <- NULL

  # Reorder: gRNA_name and gRNA_rank first
  col_order <- base::c("gRNA_name", "gRNA_rank",
                  base::setdiff(base::colnames(gRNA_df), base::c("gRNA_name", "gRNA_rank")))
  gRNA_df <- gRNA_df[, col_order, drop = FALSE]

  return(gRNA_df)
}


#' Format Primer Name from User-Defined Pattern
#'
#' Replaces placeholder tokens in a name pattern string with actual values.
#' Used internally by \code{design_cloning_primers} to generate primer names.
#'
#' @param pattern Character string with placeholders: \code{{gene}}, \code{{locus_tag}},
#'   \code{{rank}}, \code{{dir}} (F/R), \code{{nuclease}}, \code{{prefix}},
#'   \code{{method}} (OA/IA), \code{{spacer_num}}.
#' @param gene Gene name.
#' @param locus_tag Locus tag identifier.
#' @param rank gRNA rank within gene.
#' @param direction Primer direction ("F" or "R").
#' @param nuclease Nuclease abbreviation.
#' @param prefix User-defined prefix.
#' @param method Cloning method abbreviation ("OA" for Oligo Annealing/Golden Gate,
#'   "IA" for Inverse PCR Amplification/Gibson).
#' @param spacer_num Global spacer number (zero-padded string).
#' @return A formatted primer name string.
#' @examples
#' \dontrun{
#' format_primer_name("sgRNA_{locus_tag}_g{rank}_{method}_{dir}",
#'   gene = "dnaA", locus_tag = "QT234_RS00005", rank = "1",
#'   direction = "F", nuclease = "Geo", prefix = "Geo",
#'   method = "OA", spacer_num = "001")
#' # Returns: "sgRNA_QT234_RS00005_g1_OA_F"
#' }
#' @export
format_primer_name <- function(pattern, gene = "", locus_tag = "", rank = "",
                                direction = "", nuclease = "", prefix = "",
                                method = "", spacer_num = "") {
  # Use locus_tag as fallback for gene if gene is empty
  gene_or_tag <- base::ifelse(base::nchar(gene) > 0, gene, base::gsub("[^a-zA-Z0-9]", "", locus_tag))

  name <- pattern
  name <- base::gsub("\\{gene\\}", gene_or_tag, name)
  name <- base::gsub("\\{locus_tag\\}", locus_tag, name)
  name <- base::gsub("\\{rank\\}", rank, name)
  name <- base::gsub("\\{dir\\}", direction, name)
  name <- base::gsub("\\{nuclease\\}", nuclease, name)
  name <- base::gsub("\\{prefix\\}", prefix, name)
  name <- base::gsub("\\{method\\}", method, name)
  name <- base::gsub("\\{spacer_num\\}", spacer_num, name)

  return(name)
}
