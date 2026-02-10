#' Design Cloning Primers for gRNA Library
#'
#' Designs cloning primers (oligos) for inserting gRNA spacers into an expression vector.
#' Supports Golden Gate Assembly (BbsI, BsaI, BsmBI), Gibson Assembly, and
#' Deletion (4-primer Gibson arm) methods.
#'
#' For Golden Gate and Gibson (simple 2-primer spacer insertion), each gRNA in the
#' data frame receives a forward and reverse primer pair.
#'
#' For Deletion mode (\code{cloning_method = "deletion"}), the function designs
#' 4 primers per locus_tag (upstream_F/R + downstream_F/R) with Tm-optimized arms
#' and overlap regions. This mode requires additional parameters: \code{genome_seq},
#' \code{genbank_table}, \code{upstream_bp}, \code{downstream_bp}.
#'
#' @param gRNA_df A data frame containing gRNA data (output from \code{run_gRNA_list_generator}
#'   or \code{generate_gRNA_for_locus_tags}). Must contain \code{protospacer} column.
#'   For deletion mode, must contain \code{locus_tag} column.
#' @param vector_file Path to the vector file (.dna, .gb, .fasta). Required for Gibson
#'   and Deletion modes. Optional for Golden Gate (not used).
#' @param start Integer. Vector insertion site start position (1-based). Required for Gibson/Deletion.
#' @param end Integer. Vector insertion site end position (1-based). Required for Gibson/Deletion.
#' @param cloning_method Character string: \code{"golden_gate"} (default), \code{"gibson"},
#'   or \code{"deletion"} (4-primer arm-based Gibson Assembly).
#' @param enzyme Character string for Golden Gate enzyme: "BbsI" (default), "BsaI", or "BsmBI".
#' @param custom_overhangs Named list with custom overhang sequences:
#'   \code{list(F_5prime = "CACC", R_5prime = "AAAC")}. Overrides enzyme defaults.
#' @param tm_target Numeric. Target melting temperature for Gibson/Deletion primers (default: 60).
#' @param overlap_length Integer. Gibson/Deletion initial overlap length in bp (default: 20).
#' @param add_g Logical. Whether to add a 5' G for U6 promoter compatibility (default: TRUE).
#'   Adds G to forward oligo if spacer doesn't start with G. Golden Gate only.
#' @param name_pattern Character string. Template for primer naming with placeholders.
#'   Available tokens: \code{{gene}}, \code{{locus_tag}}, \code{{rank}}, \code{{dir}} (F/R),
#'   \code{{nuclease}}, \code{{prefix}}, \code{{method}} (OA/IA), \code{{spacer_num}}.
#'   Default: \code{"sgRNA_{locus_tag}_g{rank}_{method}_{dir}"} where
#'   OA = Oligo Annealing (Golden Gate), IA = Inverse PCR Amplification (Gibson).
#' @param output_file Character string. Path to save Excel output (default: NULL, no file saved).
#' @param output_dir Character string. Directory to save gRNA vector GenBank files.
#'   For golden_gate/gibson: generates one \code{{gRNA_name}_vector.gbk} per gRNA,
#'   showing the spacer inserted into the vector. Requires \code{vector_file}, \code{start},
#'   and \code{end}. For deletion: passed through to \code{batch_deletion_primers()}.
#'   Default: NULL (no GenBank output).
#' @param genome_seq Character string. Full genome sequence. Required for deletion mode
#'   (unless \code{genbank_file} is provided).
#' @param genbank_table Data frame with columns: locus_tag, start, end. Required for deletion mode
#'   (unless \code{genbank_file} is provided).
#' @param genbank_file Character string. Path to a genome GenBank file (.gb, .gbk, .gbff).
#'   When provided for deletion mode, \code{genome_seq} and \code{genbank_table} are
#'   automatically extracted from this file.
#' @param upstream_bp Integer. Upstream arm length for deletion mode (default: 500).
#' @param downstream_bp Integer. Downstream arm length for deletion mode (default: 500).
#' @param kill_snapgene Logical. Kill running SnapGene before reading .dna files (default: FALSE).
#' @return For golden_gate/gibson: the input data frame with added primer columns.
#'   For deletion: a data frame with 4 primers per unique locus_tag.
#' @examples
#' \dontrun{
#' # Golden Gate with BbsI (OA = Oligo Annealing)
#' # Primer names: sgRNA_[locus_tag]_g[rank]_OA_F / _OA_R
#' primers <- design_cloning_primers(
#'   gRNA_df = gRNA_lib,
#'   cloning_method = "golden_gate",
#'   enzyme = "BbsI"
#' )
#'
#' # Gibson Assembly (IA = Inverse PCR Amplification)
#' # Primer names: sgRNA_[locus_tag]_g[rank]_IA_F / _IA_R
#' primers <- design_cloning_primers(
#'   gRNA_df = gRNA_lib,
#'   vector_file = "my_vector.dna",
#'   start = 1234, end = 1256,
#'   cloning_method = "gibson"
#' )
#'
#' # Deletion (4-primer Gibson arms) — with separate genome_seq + genbank_table
#' del_primers <- design_cloning_primers(
#'   gRNA_df = gRNA_lib,
#'   vector_file = "my_vector.dna",
#'   start = 1234, end = 1256,
#'   cloning_method = "deletion",
#'   genome_seq = genome_sequence,
#'   genbank_table = genbank_df,
#'   output_file = "deletion_primers.xlsx"
#' )
#'
#' # Deletion — with genome GenBank file (auto-extract genome_seq + genbank_table)
#' del_primers <- design_cloning_primers(
#'   gRNA_df = gRNA_lib,
#'   vector_file = "my_vector.dna",
#'   start = 1234, end = 1256,
#'   cloning_method = "deletion",
#'   genbank_file = "my_genome.gbk",
#'   output_file = "deletion_primers.xlsx"
#' )
#' }
#' @export
design_cloning_primers <- function(
    gRNA_df,
    vector_file = NULL,
    start = NULL,
    end = NULL,
    cloning_method = "golden_gate",
    enzyme = "BbsI",
    custom_overhangs = NULL,
    tm_target = 60,
    overlap_length = 20,
    add_g = TRUE,
    name_pattern = "sgRNA_{locus_tag}_{method}_{dir}",
    output_file = NULL,
    output_dir = NULL,
    genome_seq = NULL,
    genbank_table = NULL,
    genbank_file = NULL,
    upstream_bp = 500,
    downstream_bp = 500,
    kill_snapgene = TRUE
) {
  # --- Input validation ---
  if (!base::is.data.frame(gRNA_df) || base::nrow(gRNA_df) == 0) {
    base::stop("gRNA_df must be a non-empty data frame.")
  }
  if (!cloning_method %in% base::c("golden_gate", "gibson", "deletion")) {
    base::stop("cloning_method must be 'golden_gate', 'gibson', or 'deletion'.")
  }

  # --- Deletion mode: dispatch to batch_deletion_primers ---
  if (cloning_method == "deletion") {
    if (base::is.null(vector_file)) base::stop("vector_file is required for deletion mode.")
    if (base::is.null(start) || base::is.null(end)) {
      base::stop("start and end positions are required for deletion mode.")
    }

    # Auto-extract from genbank_file if genome_seq/genbank_table not provided
    if (base::is.null(genome_seq) && base::is.null(genbank_table) &&
        base::is.null(genbank_file)) {
      base::stop("For deletion mode, provide genome_seq + genbank_table, or genbank_file.")
    }

    if (!"locus_tag" %in% base::colnames(gRNA_df)) {
      base::stop("gRNA_df must contain 'locus_tag' column for deletion mode.")
    }

    locus_tags <- base::unique(gRNA_df$locus_tag[!base::is.na(gRNA_df$locus_tag)])
    base::cat("Deletion mode: designing 4-primer Gibson arms for",
              base::length(locus_tags), "locus tags\n")

    return(batch_deletion_primers(
      vector_file = vector_file,
      genome_seq = genome_seq,
      genbank_table = genbank_table,
      locus_tags = locus_tags,
      start = start,
      end = end,
      upstream_bp = upstream_bp,
      downstream_bp = downstream_bp,
      tm_target = tm_target,
      overlap_length = overlap_length,
      output_file = output_file,
      output_dir = output_dir,
      genbank_file = genbank_file,
      kill_snapgene = kill_snapgene
    ))
  }

  # --- Golden Gate / Gibson modes ---
  if (!"protospacer" %in% base::colnames(gRNA_df)) {
    base::stop("gRNA_df must contain a 'protospacer' column.")
  }

  base::cat("Designing", cloning_method, "cloning primers for", base::nrow(gRNA_df), "gRNAs\n")

  # --- Dispatch to method-specific function ---
  if (cloning_method == "golden_gate") {
    gRNA_df <- .design_golden_gate_oligos(
      gRNA_df = gRNA_df,
      enzyme = enzyme,
      custom_overhangs = custom_overhangs,
      add_g = add_g,
      vector_file = vector_file,
      start = start,
      end = end,
      kill_snapgene = kill_snapgene
    )
  } else if (cloning_method == "gibson") {
    # Gibson requires vector info
    if (base::is.null(vector_file)) base::stop("vector_file is required for Gibson Assembly.")
    if (base::is.null(start) || base::is.null(end)) base::stop("start and end positions are required for Gibson Assembly.")

    gRNA_df <- .design_gibson_primers(
      gRNA_df = gRNA_df,
      vector_file = vector_file,
      start = base::as.integer(start),
      end = base::as.integer(end),
      tm_target = tm_target,
      overlap_length = overlap_length,
      kill_snapgene = kill_snapgene
    )
  }

  # --- Apply user-defined naming ---
  # OA = Oligo Annealing (Golden Gate), IA = Inverse PCR Amplification (Gibson)
  method_abbrev <- base::ifelse(cloning_method == "golden_gate", "OA", "IA")

  # Determine nuclease prefix from gRNA_name if available
  nuclease_prefix <- ""
  if ("gRNA_name" %in% base::colnames(gRNA_df) && base::nrow(gRNA_df) > 0) {
    first_name <- gRNA_df$gRNA_name[1]
    nuclease_prefix <- base::strsplit(first_name, "_")[[1]][1]
  }

  gRNA_df$primer_F_name <- base::sapply(base::seq_len(base::nrow(gRNA_df)), function(i) {
    format_primer_name(
      pattern   = name_pattern,
      gene      = .safe_get(gRNA_df, "gene", i),
      locus_tag = .safe_get(gRNA_df, "locus_tag", i),
      rank      = .safe_get(gRNA_df, "gRNA_rank", i),
      direction = "F",
      nuclease  = nuclease_prefix,
      prefix    = nuclease_prefix,
      method    = method_abbrev,
      spacer_num = base::sprintf("%03d", i)
    )
  })

  gRNA_df$primer_R_name <- base::sapply(base::seq_len(base::nrow(gRNA_df)), function(i) {
    format_primer_name(
      pattern   = name_pattern,
      gene      = .safe_get(gRNA_df, "gene", i),
      locus_tag = .safe_get(gRNA_df, "locus_tag", i),
      rank      = .safe_get(gRNA_df, "gRNA_rank", i),
      direction = "R",
      nuclease  = nuclease_prefix,
      prefix    = nuclease_prefix,
      method    = method_abbrev,
      spacer_num = base::sprintf("%03d", i)
    )
  })

  # --- Add oligo length columns (useful for ordering) ---
  if ("primer_F" %in% base::colnames(gRNA_df)) {
    gRNA_df$oligo_F_length <- base::nchar(gRNA_df$primer_F)
  }
  if ("primer_R" %in% base::colnames(gRNA_df)) {
    gRNA_df$oligo_R_length <- base::nchar(gRNA_df$primer_R)
  }

  # --- Reorder columns: comprehensive clean layout ---
  # Priority order for Excel readability:
  #   Gene info → gRNA identity → Primers → gRNA reference → Quality → Off-target → Method → Rest
  col_priority <- base::c(
    # Gene identification
    "locus_tag", "gene", "product",
    # gRNA identity
    "gRNA_name", "gRNA_rank",
    # Primer info (what you order)
    "primer_F_name", "primer_R_name",
    "primer_F", "primer_R",
    "primer_F_Tm", "primer_R_Tm",
    "oligo_F_length", "oligo_R_length",
    # gRNA reference (for verification)
    "protospacer", "pam", "spacer_context", "spacer_length",
    "strand", "pam_site",
    # Gene coordinates
    "gene_start", "gene_end", "direction",
    # Quality metrics
    "composite_score", "percentGC", "Tm", "position_percent",
    # Off-target
    "n0", "n1",
    # Methylation
    "methylation_conflict", "methylation_motifs",
    # Cloning method
    "cloning_method", "enzyme", "vector_name"
  )

  existing_priority <- base::intersect(col_priority, base::colnames(gRNA_df))
  remaining_cols <- base::setdiff(base::colnames(gRNA_df), existing_priority)
  # Remove non-serializable columns
  remove_cols <- base::c("alignments", "score")
  remaining_cols <- base::setdiff(remaining_cols, remove_cols)
  gRNA_df <- gRNA_df[, base::c(existing_priority, remaining_cols), drop = FALSE]

  base::cat("Primer design complete.", base::nrow(gRNA_df), "primer pairs generated.\n")

  # --- Save to Excel if output_file specified ---
  if (!base::is.null(output_file)) {
    base::tryCatch({
      excel_dir <- base::dirname(output_file)
      if (!base::dir.exists(excel_dir) && excel_dir != ".") {
        base::dir.create(excel_dir, recursive = TRUE)
      }
      if (base::requireNamespace("openxlsx", quietly = TRUE)) {
        openxlsx::write.xlsx(gRNA_df, file = output_file, overwrite = TRUE)
      } else if (base::requireNamespace("writexl", quietly = TRUE)) {
        writexl::write_xlsx(gRNA_df, path = output_file)
      } else if (base::requireNamespace("xlsx", quietly = TRUE)) {
        xlsx::write.xlsx(base::as.data.frame(gRNA_df), file = output_file, row.names = FALSE)
      } else {
        base::warning("No Excel writing package found. Install openxlsx, writexl, or xlsx.")
      }
      base::cat("Primers saved to:", output_file, "\n")
    }, error = function(e) {
      base::warning("Failed to save Excel file: ", base::conditionMessage(e))
    })
  }

  # --- Generate GenBank files (gRNA expression vector) ---
  # Only top-1 gRNA per gene gets a GenBank file (best candidate only)
  if (!base::is.null(output_dir) && !base::is.null(vector_file) &&
      !base::is.null(start) && !base::is.null(end)) {

    # Filter to rank 1 only for GenBank output
    genbank_idx <- base::seq_len(base::nrow(gRNA_df))
    if ("gRNA_rank" %in% base::colnames(gRNA_df)) {
      genbank_idx <- base::which(gRNA_df$gRNA_rank == 1)
      if (base::length(genbank_idx) == 0) {
        genbank_idx <- base::seq_len(base::min(1, base::nrow(gRNA_df)))
      }
    }

    base::cat("Generating gRNA expression vector GenBank files for",
              base::length(genbank_idx), "top-1 gRNA(s) in:", output_dir, "\n")

    for (i in genbank_idx) {
      base::tryCatch({
        spacer_i <- base::toupper(gRNA_df$protospacer[i])

        # Determine gRNA name for file naming
        grna_name_i <- ""
        if ("gRNA_name" %in% base::colnames(gRNA_df)) {
          grna_name_i <- base::as.character(gRNA_df$gRNA_name[i])
        }
        if (base::is.na(grna_name_i) || base::nchar(grna_name_i) == 0) {
          grna_name_i <- base::paste0("gRNA_", base::sprintf("%03d", i))
        }

        # Get gene/locus_tag info
        gene_i <- .safe_get(gRNA_df, "gene", i)
        locus_i <- .safe_get(gRNA_df, "locus_tag", i)

        # Build primer_info list
        pinfo <- base::list(
          primer_F = .safe_get(gRNA_df, "primer_F", i),
          primer_R = .safe_get(gRNA_df, "primer_R", i),
          primer_F_Tm = base::ifelse("primer_F_Tm" %in% base::colnames(gRNA_df),
                                      gRNA_df$primer_F_Tm[i], NA_real_),
          primer_R_Tm = base::ifelse("primer_R_Tm" %in% base::colnames(gRNA_df),
                                      gRNA_df$primer_R_Tm[i], NA_real_),
          primer_F_name = .safe_get(gRNA_df, "primer_F_name", i),
          primer_R_name = .safe_get(gRNA_df, "primer_R_name", i),
          cloning_method = base::ifelse(cloning_method == "golden_gate",
                                         "Golden_Gate", "Gibson"),
          enzyme = base::ifelse(cloning_method == "golden_gate", enzyme, ""),
          gibson_overlap = base::ifelse("gibson_overlap" %in% base::colnames(gRNA_df),
                                         gRNA_df$gibson_overlap[i], NA_integer_)
        )

        # Sanitize filename - use locus_tag as primary label
        gbk_label <- ""
        if (!base::is.null(locus_i) && !base::is.na(locus_i) && base::nchar(locus_i) > 0) {
          gbk_label <- locus_i
        } else {
          gbk_label <- grna_name_i
        }
        safe_name <- base::gsub("[^A-Za-z0-9_.-]", "_", gbk_label)
        gbk_path <- base::file.path(output_dir,
                                      base::paste0(safe_name, "_vector.gbk"))

        write_grna_vector_genbank(
          vector_file  = vector_file,
          spacer_seq   = spacer_i,
          start        = base::as.integer(start),
          end          = base::as.integer(end),
          gRNA_name    = grna_name_i,
          locus_tag    = locus_i,
          gene         = gene_i,
          nuclease     = nuclease_prefix,
          primer_info  = pinfo,
          output_path  = gbk_path,
          kill_snapgene = kill_snapgene
        )
      }, error = function(e) {
        base::warning("Failed to write GenBank for gRNA #", i, ": ",
                       base::conditionMessage(e))
      })
    }
    base::cat("GenBank files generated:", base::length(genbank_idx), "files (top-1 per gene)\n")
  } else if (!base::is.null(output_dir) && base::is.null(vector_file)) {
    base::warning("output_dir specified but vector_file is NULL. ",
                   "GenBank output requires vector_file + start + end.")
  }

  return(gRNA_df)
}


# ============================================================
# Golden Gate Assembly
# ============================================================

#' @keywords internal
.design_golden_gate_oligos <- function(gRNA_df, enzyme, custom_overhangs, add_g,
                                        vector_file = NULL, start = NULL,
                                        end = NULL, kill_snapgene = TRUE) {
  # --- Type IIS enzyme cut pattern database ---
  # Verified from NEB (neb.com) and Thermo Fisher product pages.
  # Notation: recognition(cut_fwd/cut_rev) where cut_fwd = distance from
  # recognition 3'-end to cut on the recognition-bearing strand,
  # cut_rev = distance to cut on the complementary strand.
  # Overhang length = cut_rev - cut_fwd (dynamically computed, NOT hardcoded).
  #
  # Source references:
  #   BbsI:  NEB R0539  — GAAGAC(2/6)   → 4nt 5' overhang
  #   BpiI:  Thermo ER1011 (isoschizomer of BbsI)
  #   BsaI:  NEB R0535  — GGTCTC(1/5)   → 4nt 5' overhang
  #   Eco31I: Thermo ER0291 (isoschizomer of BsaI)
  #   BsmBI: NEB R0739  — CGTCTC(1/5)   → 4nt 5' overhang
  #   Esp3I: NEB R0734  (isoschizomer of BsmBI)
  #   SapI:  NEB R0569  — GCTCTTC(1/4)  → 3nt 5' overhang (7bp recognition)
  #   BspQI: NEB R0712  (isoschizomer of SapI)
  #   PaqCI: NEB R0745  — CACCTGC(4/8)  → 4nt 5' overhang (7bp recognition)
  type_IIS_db <- base::list(
    BbsI   = base::list(recognition = "GAAGAC",  cut_fwd = 2L, cut_rev = 6L),
    BpiI   = base::list(recognition = "GAAGAC",  cut_fwd = 2L, cut_rev = 6L),
    BsaI   = base::list(recognition = "GGTCTC",  cut_fwd = 1L, cut_rev = 5L),
    Eco31I = base::list(recognition = "GGTCTC",  cut_fwd = 1L, cut_rev = 5L),
    BsmBI  = base::list(recognition = "CGTCTC",  cut_fwd = 1L, cut_rev = 5L),
    Esp3I  = base::list(recognition = "CGTCTC",  cut_fwd = 1L, cut_rev = 5L),
    SapI   = base::list(recognition = "GCTCTTC", cut_fwd = 1L, cut_rev = 4L),
    BspQI  = base::list(recognition = "GCTCTTC", cut_fwd = 1L, cut_rev = 4L),
    PaqCI  = base::list(recognition = "CACCTGC", cut_fwd = 4L, cut_rev = 8L)
  )

  # Default hardcoded overhangs (fallback when vector_file not provided).
  # These are VECTOR-SPECIFIC conventions for common CRISPR gRNA expression
  # vectors (e.g. pSpCas9, px330). For non-standard vectors, always provide
  # vector_file + start + end to compute real overhangs from the sequence.
  # Enzymes without a default require vector_file or custom_overhangs.
  enzyme_overhangs <- base::list(
    BbsI   = base::list(F_5prime = "CACC", R_5prime = "AAAC"),
    BpiI   = base::list(F_5prime = "CACC", R_5prime = "AAAC"),
    BsaI   = base::list(F_5prime = "CTTC", R_5prime = "CAAA"),
    Eco31I = base::list(F_5prime = "CTTC", R_5prime = "CAAA"),
    BsmBI  = base::list(F_5prime = "CACC", R_5prime = "AAAC"),
    Esp3I  = base::list(F_5prime = "CACC", R_5prime = "AAAC")
    # SapI, BspQI, PaqCI: no universal default — vector_file required
  )

  # --- Determine overhangs + fill ---
  use_vector_overhangs <- FALSE
  fill_left  <- ""
  fill_right <- ""

  if (!base::is.null(custom_overhangs)) {
    # Mode 1: User-specified custom overhangs
    if (!base::all(base::c("F_5prime", "R_5prime") %in% base::names(custom_overhangs))) {
      base::stop("custom_overhangs must contain 'F_5prime' and 'R_5prime' elements.")
    }
    overhangs <- custom_overhangs
    base::cat("Using custom overhangs: F=", overhangs$F_5prime, ", R=", overhangs$R_5prime, "\n")

  } else if (!base::is.null(vector_file) && !base::is.null(start) && !base::is.null(end)) {
    # Mode 2: Vector-based overhang detection
    if (!enzyme %in% base::names(type_IIS_db)) {
      base::stop("Unsupported enzyme for vector-based OA: ", enzyme,
                 "\nSupported: ", base::paste(base::names(type_IIS_db), collapse = ", "))
    }

    base::cat("Reading vector file for OA primer design...\n")
    vec_info <- read_vector_file(vector_file, kill_snapgene = kill_snapgene)
    vec_seq  <- base::toupper(base::as.character(vec_info$sequence))

    enzyme_info <- type_IIS_db[[enzyme]]
    site_result <- .find_enzyme_sites_and_compute_overhangs(
      vec_seq     = vec_seq,
      start       = base::as.integer(start),
      end         = base::as.integer(end),
      enzyme_info = enzyme_info,
      enzyme_name = enzyme
    )

    if (site_result$found) {
      overhangs <- base::list(
        F_5prime = site_result$left_oh,
        R_5prime = site_result$right_oh
      )
      fill_left  <- site_result$fill_left
      fill_right <- site_result$fill_right
      use_vector_overhangs <- TRUE
      base::cat("Using vector-derived ", site_result$overhang_len, "nt ",
                site_result$overhang_type, " overhangs: F=", overhangs$F_5prime,
                ", R=", overhangs$R_5prime,
                base::ifelse(base::nchar(fill_left) > 0, base::paste0(", fill_L=", fill_left), ""),
                base::ifelse(base::nchar(fill_right) > 0, base::paste0(", fill_R=", fill_right), ""),
                "\n")
    } else {
      # Fallback to default hardcoded overhangs
      if (enzyme %in% base::names(enzyme_overhangs)) {
        overhangs <- enzyme_overhangs[[enzyme]]
        base::cat("Fallback to default ", enzyme, " overhangs: F=", overhangs$F_5prime,
                  ", R=", overhangs$R_5prime, "\n")
      } else {
        base::stop("No default overhangs for enzyme: ", enzyme)
      }
    }

  } else if (enzyme %in% base::names(enzyme_overhangs)) {
    # Mode 3: Default hardcoded overhangs (no vector provided)
    overhangs <- enzyme_overhangs[[enzyme]]
    base::cat("Using", enzyme, "default overhangs: F=", overhangs$F_5prime,
              ", R=", overhangs$R_5prime, "\n")

  } else if (enzyme %in% base::names(type_IIS_db)) {
    # Enzyme is in DB but has no default overhangs → vector_file required
    base::stop("Enzyme ", enzyme, " has no default overhangs (overhang = ",
               type_IIS_db[[enzyme]]$cut_rev - type_IIS_db[[enzyme]]$cut_fwd,
               "nt). Please provide vector_file + start + end, ",
               "or use custom_overhangs = list(F_5prime = '...', R_5prime = '...').")
  } else {
    base::stop("Unsupported enzyme: ", enzyme,
               "\nSupported (with defaults): ", base::paste(base::names(enzyme_overhangs), collapse = ", "),
               "\nSupported (vector required): SapI, BspQI, PaqCI",
               "\nOr provide custom_overhangs = list(F_5prime = '...', R_5prime = '...')")
  }

  # Enzyme recognition sites for internal site checking
  recognition_site <- NULL
  if (enzyme %in% base::names(type_IIS_db)) {
    recognition_site <- type_IIS_db[[enzyme]]$recognition
  }

  # Design oligos for each spacer
  spacers <- base::toupper(gRNA_df$protospacer)
  n <- base::length(spacers)

  primer_F <- base::character(n)
  primer_R <- base::character(n)
  primer_F_Tm <- base::numeric(n)
  primer_R_Tm <- base::numeric(n)

  for (i in base::seq_len(n)) {
    sp <- spacers[i]

    # Check for internal enzyme recognition site
    if (!base::is.null(recognition_site)) {
      rc_site <- base::as.character(Biostrings::reverseComplement(Biostrings::DNAString(recognition_site)))
      if (base::grepl(recognition_site, sp, ignore.case = TRUE) ||
          base::grepl(rc_site, sp, ignore.case = TRUE)) {
        base::warning("gRNA #", i, " (", sp, ") contains internal ", enzyme,
                " recognition site. Cloning may fail.")
      }
    }

    # G/C addition for U6 promoter compatibility
    g_add <- ""
    c_add <- ""
    if (add_g && base::substring(sp, 1, 1) != "G") {
      g_add <- "G"
      c_add <- "C"
    }

    # Forward oligo: LEFT_OH + fill_left + G + spacer + fill_right
    primer_F[i] <- base::paste0(overhangs$F_5prime, fill_left, g_add, sp, fill_right)

    # Reverse oligo: RIGHT_OH + RC(fill_right) + RC(spacer) + C + RC(fill_left)
    sp_rc <- base::as.character(Biostrings::reverseComplement(Biostrings::DNAString(sp)))

    rc_fill_left  <- ""
    rc_fill_right <- ""
    if (base::nchar(fill_left) > 0) {
      rc_fill_left <- base::as.character(
        Biostrings::reverseComplement(Biostrings::DNAString(fill_left))
      )
    }
    if (base::nchar(fill_right) > 0) {
      rc_fill_right <- base::as.character(
        Biostrings::reverseComplement(Biostrings::DNAString(fill_right))
      )
    }

    primer_R[i] <- base::paste0(overhangs$R_5prime, rc_fill_right, sp_rc, c_add, rc_fill_left)

    # Tm for annealing region (fill + G + spacer + fill, excluding overhang)
    anneal_F <- base::paste0(fill_left, g_add, sp, fill_right)
    anneal_R <- base::paste0(rc_fill_right, sp_rc, c_add, rc_fill_left)

    primer_F_Tm[i] <- base::tryCatch({
      TmCalculator::Tm_NN(anneal_F, Na = 50)$Tm
    }, error = function(e) NA_real_)

    primer_R_Tm[i] <- base::tryCatch({
      TmCalculator::Tm_NN(anneal_R, Na = 50)$Tm
    }, error = function(e) NA_real_)
  }

  gRNA_df$primer_F <- primer_F
  gRNA_df$primer_R <- primer_R
  gRNA_df$primer_F_Tm <- base::round(primer_F_Tm, 1)
  gRNA_df$primer_R_Tm <- base::round(primer_R_Tm, 1)
  gRNA_df$cloning_method <- "Golden_Gate"
  gRNA_df$enzyme <- enzyme

  # Add fill info columns if vector-based
  if (use_vector_overhangs) {
    gRNA_df$overhang_F <- overhangs$F_5prime
    gRNA_df$overhang_R <- overhangs$R_5prime
    if (base::nchar(fill_left) > 0)  gRNA_df$fill_left  <- fill_left
    if (base::nchar(fill_right) > 0) gRNA_df$fill_right <- fill_right
  }

  return(gRNA_df)
}


# ============================================================
# Gibson Assembly
# ============================================================

#' @keywords internal
.design_gibson_primers <- function(gRNA_df, vector_file, start, end,
                                    tm_target, overlap_length,
                                    kill_snapgene = TRUE) {
  # Read vector
  vector <- read_vector_file(vector_file, kill_snapgene = kill_snapgene)
  vec_seq <- base::toupper(vector$sequence)
  vec_len <- base::nchar(vec_seq)

  if (start < 1 || end > vec_len || start > end) {
    base::stop("Invalid start/end positions. Vector length: ", vec_len, " bp")
  }

  # === IA (Inverse PCR Amplification) Primer Design ===
  # Primers amplify the vector backbone (everything except the stuffer).
  # The gRNA spacer is encoded as a 5' tail on BOTH primers (= Gibson overlap).
  #
  # F primer: 5'-[spacer]-------[backbone_head_F]-3'
  #               overhang (tail)   annealing (head), Tm 57-65
  # R primer: 5'-[RC(spacer)]---[RC(backbone_head_R)]-3'
  #               overhang (tail)   annealing (head), Tm 57-65
  #
  # Gibson overlap = full spacer sequence (both primers share it).
  # After inverse PCR, the linear product self-circularizes via Gibson
  # at the complementary spacer overhangs, inserting the gRNA.

  # Tm range for backbone-annealing HEAD

  tm_head_min <- 57
  tm_head_max <- 65
  min_head_len <- 15L
  max_head_len <- 30L

  # --- F primer HEAD: vector sequence starting at end+1, extending rightward ---
  # Extends 1bp at a time. Stops as soon as Tm >= 57 (cost-effective: shortest head)
  f_head <- ""
  f_head_tm <- NA_real_

  for (len in base::seq.int(min_head_len, max_head_len)) {
    f_start <- end + 1L
    f_end_pos <- end + len
    if (f_end_pos <= vec_len) {
      candidate <- base::substring(vec_seq, f_start, f_end_pos)
    } else {
      candidate <- base::paste0(
        base::substring(vec_seq, f_start, vec_len),
        base::substring(vec_seq, 1L, f_end_pos - vec_len)
      )
    }
    f_head <- candidate
    f_head_tm <- .safe_tm(candidate)
    if (!base::is.na(f_head_tm) && f_head_tm >= tm_head_min) break
  }

  # --- R primer HEAD: vector sequence ending at start-1, extending leftward ---
  # Store in sense (top-strand) orientation; will RC for primer later
  r_head_sense <- ""
  r_head_tm <- NA_real_

  for (len in base::seq.int(min_head_len, max_head_len)) {
    r_start <- start - len
    r_end_pos <- start - 1L
    if (r_start >= 1L) {
      candidate <- base::substring(vec_seq, r_start, r_end_pos)
    } else {
      candidate <- base::paste0(
        base::substring(vec_seq, vec_len + r_start, vec_len),
        base::substring(vec_seq, 1L, r_end_pos)
      )
    }
    r_head_sense <- candidate
    r_head_tm <- .safe_tm(candidate)
    if (!base::is.na(r_head_tm) && r_head_tm >= tm_head_min) break
  }

  # RC the R head for the actual primer sequence
  r_head_rc <- .rc(r_head_sense)

  base::cat("IA backbone annealing: F_head=", base::nchar(f_head), "bp (Tm=",
            base::round(f_head_tm, 1), "\u00b0C), R_head=", base::nchar(r_head_sense),
            "bp (Tm=", base::round(r_head_tm, 1), "\u00b0C)\n")

  # --- Assemble primers for each spacer ---
  # Asymmetric Gibson overlap design:
  #   F primer: 5'-[FULL spacer]-[backbone head F]-3'
  #     -> Carries full spacer (template for polymerase fill-in)
  #   R primer: 5'-[RC(spacer first N bp)]-[backbone head R RC]-3'
  #     -> Carries shortened spacer tail (Gibson overlap only, >= 15bp)
  #
  # After inverse PCR + Gibson:
  #   1. T5 exonuclease creates 3' single-stranded overhangs
  #   2. The N bp overlap (from R) anneals with N bp of F's spacer tail
  #   3. Polymerase fills in the remaining (spacer_len - N) bases using F's
  #      full spacer as template
  #   4. Result: full spacer correctly inserted
  #
  # This saves oligo cost by shortening R primer while preserving function.

  spacers <- base::toupper(gRNA_df$protospacer)
  n <- base::length(spacers)

  primer_F <- base::character(n)
  primer_R <- base::character(n)
  primer_F_Tm <- base::numeric(n)
  primer_R_Tm <- base::numeric(n)
  gibson_overlap <- base::integer(n)

  gibson_overlap_tm_min <- 48
  gibson_overlap_min_len <- 15L

  for (i in base::seq_len(n)) {
    sp <- spacers[i]
    sp_len <- base::nchar(sp)

    # F primer: full spacer + backbone head (spacer = insert template, cannot shorten)
    primer_F[i] <- base::paste0(sp, f_head)
    primer_F_Tm[i] <- base::round(f_head_tm, 1)

    # R primer: optimize Gibson overlap length
    # Start from 15bp (Gibson minimum), stop when overlap Tm >= 48C
    r_tail_len <- sp_len  # default: full spacer
    for (tl in base::seq.int(gibson_overlap_min_len, sp_len)) {
      tail_candidate <- base::substring(sp, 1L, tl)
      tail_tm <- .safe_tm(tail_candidate)
      if (!base::is.na(tail_tm) && tail_tm >= gibson_overlap_tm_min) {
        r_tail_len <- tl
        break
      }
    }

    r_sp_tail <- base::substring(sp, 1L, r_tail_len)
    r_sp_tail_rc <- .rc(r_sp_tail)

    primer_R[i] <- base::paste0(r_sp_tail_rc, r_head_rc)
    primer_R_Tm[i] <- base::round(r_head_tm, 1)
    gibson_overlap[i] <- r_tail_len
  }

  base::cat("IA primer: F = spacer(", sp_len, "bp) + head(",
            base::nchar(f_head), "bp) = ", base::nchar(primer_F[1]), "bp\n")
  base::cat("IA primer: R = overlap(", r_tail_len, "bp) + head(",
            base::nchar(r_head_sense), "bp) = ", base::nchar(primer_R[1]), "bp\n")
  base::cat("Gibson overlap = ", r_tail_len, "bp (saved ",
            sp_len - r_tail_len, "bp on R primer)\n")

  gRNA_df$primer_F <- primer_F
  gRNA_df$primer_R <- primer_R
  gRNA_df$primer_F_Tm <- primer_F_Tm
  gRNA_df$primer_R_Tm <- primer_R_Tm
  gRNA_df$gibson_overlap <- gibson_overlap
  gRNA_df$cloning_method <- "Gibson"
  gRNA_df$vector_name <- vector$name

  return(gRNA_df)
}


# ============================================================
# Deletion Gibson Assembly (4-primer arm design)
# ============================================================

#' Design Deletion Primers (Gibson Assembly, 4-Primer)
#'
#' Designs 4 Gibson Assembly primers for gene deletion vector construction.
#' The upstream and downstream arms flanking the target gene are joined with
#' an optimized overlap, and outer primers include vector homology arms.
#'
#' Two usage modes:
#' \enumerate{
#'   \item \strong{Direct mode}: Provide \code{insert_seq} (pre-extracted upstream + downstream arm).
#'   \item \strong{Locus tag mode}: Provide \code{locus_tag} + \code{genome_seq} + \code{genbank_table},
#'     and the function auto-extracts arms from the genome (supports circular genomes).
#' }
#'
#' This function faithfully ports the primer design logic from
#' \code{PrimerDesigner_for_Deletion.py}, including:
#' \itemize{
#'   \item Outer primer Tm optimization (extend to 25bp, break if >60°C)
#'   \item Vector overlap extension when overlap Tm < 55°C (up to 20bp)
#'   \item Inner primer Tm optimization (55-65°C range, 60°C optimal, ±1°C tolerance)
#'   \item Overlap region optimization (≥20bp combined, Tm ≥55°C, alternate extension by 2bp)
#'   \item Max primer length 70bp constraint
#'   \item Self-dimerization checking
#' }
#'
#' @param vector_file Path to the vector file (.dna, .gb, .fasta).
#' @param insert_seq Character string. Combined arm sequence (upstream_arm + downstream_arm).
#'   Not required if \code{locus_tag}, \code{genome_seq}, and \code{genbank_table} are provided.
#' @param start Integer. Vector insertion site start position (1-based, inclusive).
#' @param end Integer. Vector insertion site end position (1-based, inclusive).
#' @param upstream_bp Integer. Length of the upstream arm in bp (default: 500).
#' @param downstream_bp Integer. Length of the downstream arm in bp (default: 500).
#' @param tm_target Numeric. Target melting temperature for outer primers (default: 60).
#' @param overlap_length Integer. Initial vector overlap length in bp (default: 16).
#' @param min_target_length Integer. Minimum target binding region length (default: 18).
#' @param max_target_length Integer. Maximum target binding region length (default: 50).
#' @param max_primer_length Integer. Maximum total primer length (default: 70).
#' @param max_iterations Integer. Maximum optimization iterations (default: 100).
#' @param locus_tag Character. A single locus tag (e.g., \code{"QT234_RS00005"}) or
#'   a range (\code{"TAG1-TAG2"}). When provided with \code{genome_seq} and
#'   \code{genbank_table} (or \code{genbank_file}), arms are auto-extracted from the genome.
#' @param genome_seq Character string. Full genome sequence. Required for locus_tag mode
#'   (unless \code{genbank_file} is provided).
#' @param genbank_table Data frame with columns: \code{locus_tag}, \code{start}, \code{end}.
#'   Required for locus_tag mode (unless \code{genbank_file} is provided).
#'   Optionally: \code{gene}, \code{product}.
#' @param genbank_file Character string. Path to a genome GenBank file (.gb, .gbk, .gbff).
#'   When provided, \code{genome_seq} and \code{genbank_table} are automatically extracted
#'   from this file using \code{read_genome_genbank()}. This overrides any separately
#'   provided \code{genome_seq} and \code{genbank_table}.
#' @param output_dir Character. Directory to save output GenBank files. If provided,
#'   a \code{{locus_tag}_deletion_vector.gbk} file is generated for SnapGene viewing.
#'   Default: NULL (no GenBank output).
#' @param kill_snapgene Logical. Kill running SnapGene process before reading .dna
#'   files to avoid file lock conflicts (default: FALSE).
#' @return A named list with 4 primer sequences, Tm values (target and full),
#'   primer positions in the modified vector, overlap Tm, and overlap length.
#' @examples
#' \dontrun{
#' # --- Mode 1: Single locus_tag with genome_seq + genbank_table ---
#' result <- design_deletion_primers(
#'   vector_file = "my_vector.dna",
#'   locus_tag = "QT234_RS00005",
#'   genome_seq = genome_sequence,
#'   genbank_table = genbank_df,
#'   start = 1234, end = 1256,
#'   upstream_bp = 500,
#'   downstream_bp = 500,
#'   output_dir = "deletion_results",
#'   kill_snapgene = TRUE
#' )
#'
#' # --- Mode 1b: Single locus_tag with genbank_file (auto-extract) ---
#' result <- design_deletion_primers(
#'   vector_file = "my_vector.dna",
#'   locus_tag = "QT234_RS00005",
#'   genbank_file = "my_genome.gbk",
#'   start = 1234, end = 1256,
#'   output_dir = "deletion_results"
#' )
#'
#' # --- Mode 2: Direct insert_seq ---
#' insert_seq <- paste0(upstream_arm, downstream_arm)
#' result <- design_deletion_primers(
#'   vector_file = "my_vector.dna",
#'   insert_seq = insert_seq,
#'   start = 1234, end = 1256,
#'   upstream_bp = 500,
#'   downstream_bp = 500
#' )
#'
#' cat("Upstream F:", result$upstream_forward_primer, "\n")
#' cat("Upstream R:", result$upstream_reverse_primer, "\n")
#' cat("Downstream F:", result$downstream_forward_primer, "\n")
#' cat("Downstream R:", result$downstream_reverse_primer, "\n")
#' }
#' @export
design_deletion_primers <- function(
    vector_file,
    insert_seq = NULL,
    start,
    end,
    upstream_bp = 500,
    downstream_bp = 500,
    tm_target = 60,
    overlap_length = 16,
    min_target_length = 18,
    max_target_length = 50,
    max_primer_length = 70,
    max_iterations = 100,
    locus_tag = NULL,
    genome_seq = NULL,
    genbank_table = NULL,
    genbank_file = NULL,
    output_dir = NULL,
    kill_snapgene = TRUE
) {
  # --- Input validation ---
  if (base::missing(vector_file) || !base::is.character(vector_file)) {
    base::stop("vector_file must be a valid file path.")
  }

  # --- Auto-extract genome_seq and genbank_table from genbank_file ---
  if (!base::is.null(genbank_file)) {
    base::cat("Parsing genome GenBank file for genome_seq and genbank_table...\n")
    parsed <- read_genome_genbank(genbank_file)
    genome_seq <- parsed$genome_seq
    genbank_table <- parsed$genbank_table
    base::cat("  Genome:", parsed$length, "bp,", base::nrow(genbank_table),
              "features extracted\n")
  }

  # --- Locus tag mode: auto-extract arms from genome ---
  if (!base::is.null(locus_tag) && !base::is.null(genome_seq) &&
      !base::is.null(genbank_table)) {
    required_cols <- base::c("locus_tag", "start", "end")
    missing_cols <- base::setdiff(required_cols, base::colnames(genbank_table))
    if (base::length(missing_cols) > 0) {
      base::stop("genbank_table missing required columns: ", base::paste(missing_cols, collapse = ", "))
    }

    genome_seq_upper <- base::toupper(genome_seq)
    seq_length <- base::nchar(genome_seq_upper)

    # Handle range notation (e.g., "TAG1-TAG2")
    if (base::grepl("-", locus_tag)) {
      parts <- base::strsplit(locus_tag, "-", fixed = TRUE)[[1]]
      start_row <- genbank_table[genbank_table$locus_tag == parts[1], ]
      end_row <- genbank_table[genbank_table$locus_tag == parts[2], ]
      if (base::nrow(start_row) == 0 || base::nrow(end_row) == 0) {
        base::stop("One or both locus tags not found: ", locus_tag)
      }
      gene_start <- base::as.integer(start_row$start[1])
      gene_end <- base::as.integer(end_row$end[1])
    } else {
      target_row <- genbank_table[genbank_table$locus_tag == locus_tag, ]
      if (base::nrow(target_row) == 0) {
        base::stop("Locus tag not found in genbank_table: ", locus_tag)
      }
      gene_start <- base::as.integer(target_row$start[1])
      gene_end <- base::as.integer(target_row$end[1])
    }

    base::cat("Extracting deletion arms for:", locus_tag,
        "(gene pos:", gene_start, "-", gene_end, ")\n")

    # Extract upstream arm (handle circular genome)
    up_end <- gene_start - 1
    up_start <- gene_start - upstream_bp
    if (up_start <= 0) {
      up_start_circ <- seq_length + up_start
      upstream_arm <- base::paste0(
        base::substring(genome_seq_upper, up_start_circ, seq_length),
        base::substring(genome_seq_upper, 1, up_end)
      )
    } else {
      upstream_arm <- base::substring(genome_seq_upper, up_start, up_end)
    }

    # Extract downstream arm (handle circular genome)
    dn_start <- gene_end + 1
    dn_end <- gene_end + downstream_bp
    if (dn_end > seq_length) {
      dn_end_circ <- dn_end %% seq_length
      downstream_arm <- base::paste0(
        base::substring(genome_seq_upper, dn_start, seq_length),
        base::substring(genome_seq_upper, 1, dn_end_circ)
      )
    } else {
      downstream_arm <- base::substring(genome_seq_upper, dn_start, dn_end)
    }

    insert_seq <- base::paste0(upstream_arm, downstream_arm)
    base::cat("  Upstream arm:", base::nchar(upstream_arm), "bp, Downstream arm:",
        base::nchar(downstream_arm), "bp\n")
  }

  # --- Validate insert_seq ---
  if (base::is.null(insert_seq) || !base::is.character(insert_seq) || base::nchar(insert_seq) == 0) {
    base::stop("insert_seq is required. Provide it directly, or provide locus_tag + genome_seq + genbank_table.")
  }
  if (base::nchar(insert_seq) != (upstream_bp + downstream_bp)) {
    base::stop("Length of insert_seq (", base::nchar(insert_seq),
         ") must equal upstream_bp + downstream_bp (", upstream_bp + downstream_bp, ").")
  }

  # Read vector
  vector <- read_vector_file(vector_file)
  vec_seq <- base::tolower(vector$sequence)

  base::cat("Designing deletion Gibson Assembly primers (4-primer)\n")
  base::cat("  Vector:", vector$name, "(", vector$length, "bp)\n")
  base::cat("  Insert:", base::nchar(insert_seq), "bp (upstream:", upstream_bp, "bp + downstream:", downstream_bp, "bp)\n")
  base::cat("  Insertion site:", start, "-", end, "(1-based)\n")

  # Convert to 0-based internally (matching Python logic)
  start_0 <- base::as.integer(start) - 1L
  end_0 <- base::as.integer(end)

  insert_seq <- base::toupper(insert_seq)

  # Split insert into upstream and downstream arms
  upstream_seq <- base::substring(insert_seq, 1, upstream_bp)
  downstream_seq <- base::substring(insert_seq, upstream_bp + 1, base::nchar(insert_seq))

  # Extract vector overlaps (0-based slicing → R 1-based)
  # Python: vector_seq[start_0 - overlap_length : start_0]
  upstream_vector_overlap <- base::toupper(base::substring(vec_seq,
                                                start_0 - overlap_length + 1,
                                                start_0))
  # Python: vector_seq[end_0 : end_0 + overlap_length]
  downstream_vector_overlap <- base::toupper(base::substring(vec_seq,
                                                  end_0 + 1,
                                                  end_0 + overlap_length))

  # ===========================================================
  # 1. UPSTREAM FORWARD PRIMER (outer — vector overlap + arm start)
  # ===========================================================
  upstream_left <- base::substring(upstream_seq, 1, min_target_length)
  upstream_fwd_primer <- base::paste0(upstream_vector_overlap, upstream_left)
  upstream_fwd_tm_target <- .safe_tm(upstream_left)
  upstream_fwd_overlap_tm <- .safe_tm(upstream_vector_overlap)
  upstream_fwd_tm_full <- .safe_tm(upstream_fwd_primer)

  best_upstream_fwd <- base::list(
    primer = upstream_fwd_primer,
    tm_target = upstream_fwd_tm_target,
    tm_full = upstream_fwd_tm_full,
    diff = base::abs(.na0(upstream_fwd_tm_target) - tm_target)
  )

  # Extend target region until Tm ≥ 55°C (max 25bp)
  iteration <- 0
  while (.na0(upstream_fwd_tm_target) < 55 && iteration < max_iterations) {
    if (base::nchar(upstream_left) < 25) {
      next_len <- base::nchar(upstream_left) + 1
      if (next_len > base::nchar(upstream_seq)) break
      upstream_left <- base::substring(upstream_seq, 1, next_len)
      upstream_fwd_primer <- base::paste0(upstream_vector_overlap, upstream_left)
      upstream_fwd_tm_target <- .safe_tm(upstream_left)
      upstream_fwd_tm_full <- .safe_tm(upstream_fwd_primer)

      # If overshooting > 60°C, trim back 1bp and stop
      if (.na0(upstream_fwd_tm_target) > 60) {
        upstream_left <- base::substring(upstream_left, 1, base::nchar(upstream_left) - 1)
        upstream_fwd_primer <- base::paste0(upstream_vector_overlap, upstream_left)
        upstream_fwd_tm_target <- .safe_tm(upstream_left)
        upstream_fwd_tm_full <- .safe_tm(upstream_fwd_primer)
        break
      }

      if (!.check_self_dimerization(upstream_fwd_primer)) {
        best_upstream_fwd <- base::list(
          primer = upstream_fwd_primer,
          tm_target = upstream_fwd_tm_target,
          tm_full = upstream_fwd_tm_full,
          diff = base::abs(.na0(upstream_fwd_tm_target) - tm_target)
        )
      }
    }
    iteration <- iteration + 1
  }

  # Extend vector overlap if overlap Tm < 55°C (up to 20bp)
  if (.na0(upstream_fwd_overlap_tm) < 55 && base::nchar(upstream_vector_overlap) < 20) {
    while (.na0(upstream_fwd_overlap_tm) < 55 && base::nchar(upstream_vector_overlap) < 20) {
      new_ol_len <- base::nchar(upstream_vector_overlap) + 1
      upstream_vector_overlap <- base::toupper(base::substring(vec_seq,
                                                    start_0 - new_ol_len + 1,
                                                    start_0))
      upstream_fwd_primer <- base::paste0(upstream_vector_overlap, upstream_left)
      upstream_fwd_overlap_tm <- .safe_tm(upstream_vector_overlap)
      upstream_fwd_tm_full <- .safe_tm(upstream_fwd_primer)
      if (!.check_self_dimerization(upstream_fwd_primer)) {
        best_upstream_fwd <- base::list(
          primer = upstream_fwd_primer,
          tm_target = upstream_fwd_tm_target,
          tm_full = upstream_fwd_tm_full,
          diff = base::abs(.na0(upstream_fwd_tm_target) - tm_target)
        )
      }
    }
  }

  # ===========================================================
  # 2. DOWNSTREAM REVERSE PRIMER (outer — arm end + vector overlap)
  # ===========================================================
  downstream_right <- base::substring(downstream_seq,
                                 base::nchar(downstream_seq) - min_target_length + 1,
                                 base::nchar(downstream_seq))
  downstream_rev_primer <- .rc(base::paste0(downstream_right, downstream_vector_overlap))
  downstream_rev_target <- .rc(downstream_right)
  downstream_rev_tm_target <- .safe_tm(downstream_rev_target)
  downstream_rev_overlap_tm <- .safe_tm(downstream_vector_overlap)
  downstream_rev_tm_full <- .safe_tm(downstream_rev_primer)

  best_downstream_rev <- base::list(
    primer = downstream_rev_primer,
    tm_target = downstream_rev_tm_target,
    tm_full = downstream_rev_tm_full,
    diff = base::abs(.na0(downstream_rev_tm_target) - tm_target)
  )

  # Extend target region until Tm ≥ 55°C (max 25bp)
  iteration <- 0
  while (.na0(downstream_rev_tm_target) < 55 && iteration < max_iterations) {
    if (base::nchar(downstream_right) < 25) {
      new_start <- base::nchar(downstream_seq) - base::nchar(downstream_right)
      if (new_start < 1) break
      downstream_right <- base::substring(downstream_seq, new_start, base::nchar(downstream_seq))
      downstream_rev_primer <- .rc(base::paste0(downstream_right, downstream_vector_overlap))
      downstream_rev_target <- .rc(downstream_right)
      downstream_rev_tm_target <- .safe_tm(downstream_rev_target)
      downstream_rev_tm_full <- .safe_tm(downstream_rev_primer)

      # If overshooting > 60°C, trim 1bp from front and stop
      if (.na0(downstream_rev_tm_target) > 60) {
        downstream_right <- base::substring(downstream_right, 2, base::nchar(downstream_right))
        downstream_rev_primer <- .rc(base::paste0(downstream_right, downstream_vector_overlap))
        downstream_rev_target <- .rc(downstream_right)
        downstream_rev_tm_target <- .safe_tm(downstream_rev_target)
        downstream_rev_tm_full <- .safe_tm(downstream_rev_primer)
        break
      }

      if (!.check_self_dimerization(downstream_rev_primer)) {
        best_downstream_rev <- base::list(
          primer = downstream_rev_primer,
          tm_target = downstream_rev_tm_target,
          tm_full = downstream_rev_tm_full,
          diff = base::abs(.na0(downstream_rev_tm_target) - tm_target)
        )
      }
    }
    iteration <- iteration + 1
  }

  # Extend vector overlap if overlap Tm < 55°C (up to 20bp)
  if (.na0(downstream_rev_overlap_tm) < 55 && base::nchar(downstream_vector_overlap) < 20) {
    while (.na0(downstream_rev_overlap_tm) < 55 && base::nchar(downstream_vector_overlap) < 20) {
      new_ol_len <- base::nchar(downstream_vector_overlap) + 1
      downstream_vector_overlap <- base::toupper(base::substring(vec_seq,
                                                      end_0 + 1,
                                                      end_0 + new_ol_len))
      downstream_rev_primer <- .rc(base::paste0(downstream_right, downstream_vector_overlap))
      downstream_rev_overlap_tm <- .safe_tm(downstream_vector_overlap)
      downstream_rev_tm_full <- .safe_tm(downstream_rev_primer)
      if (!.check_self_dimerization(downstream_rev_primer)) {
        best_downstream_rev <- base::list(
          primer = downstream_rev_primer,
          tm_target = downstream_rev_tm_target,
          tm_full = downstream_rev_tm_full,
          diff = base::abs(.na0(downstream_rev_tm_target) - tm_target)
        )
      }
    }
  }

  # ===========================================================
  # 3. UPSTREAM REVERSE PRIMER (inner — Tm 55-65°C optimization)
  # ===========================================================
  target_tm_min <- 55
  target_tm_max <- 65
  target_tm_optimal <- 60

  upstream_right_inner <- base::substring(upstream_seq,
                                     base::nchar(upstream_seq) - min_target_length + 1,
                                     base::nchar(upstream_seq))
  ur_tm <- .safe_tm(.rc(upstream_right_inner))

  best_ur_target <- upstream_right_inner
  best_ur_tm <- ur_tm

  iteration <- 0
  while (iteration < max_iterations) {
    if (base::is.na(ur_tm)) break
    if (ur_tm >= target_tm_min && base::abs(ur_tm - target_tm_optimal) <= 1) break

    if (base::nchar(upstream_right_inner) >= max_target_length) break

    if (ur_tm < target_tm_min) {
      # Extend left (add 1bp from upstream_seq)
      new_start <- base::nchar(upstream_seq) - base::nchar(upstream_right_inner)
      if (new_start < 1) break
      upstream_right_inner <- base::substring(upstream_seq, new_start, base::nchar(upstream_seq))
    } else if (ur_tm > target_tm_max) {
      # Trim from left
      upstream_right_inner <- base::substring(upstream_right_inner, 2, base::nchar(upstream_right_inner))
    } else {
      # Within range but not optimal — fine-tune
      if (ur_tm < target_tm_optimal && base::nchar(upstream_right_inner) < max_target_length) {
        new_start <- base::nchar(upstream_seq) - base::nchar(upstream_right_inner)
        if (new_start < 1) break
        upstream_right_inner <- base::substring(upstream_seq, new_start, base::nchar(upstream_seq))
      } else if (ur_tm > target_tm_optimal && base::nchar(upstream_right_inner) > min_target_length) {
        upstream_right_inner <- base::substring(upstream_right_inner, 2, base::nchar(upstream_right_inner))
      } else {
        break
      }
    }

    ur_tm <- .safe_tm(.rc(upstream_right_inner))
    if (!base::is.na(ur_tm) && ur_tm >= target_tm_min && ur_tm <= target_tm_max) {
      best_ur_target <- upstream_right_inner
      best_ur_tm <- ur_tm
    }
    iteration <- iteration + 1
  }

  # ===========================================================
  # 4. DOWNSTREAM FORWARD PRIMER (inner — Tm 55-65°C optimization)
  # ===========================================================
  downstream_left_inner <- base::substring(downstream_seq, 1, min_target_length)
  df_tm <- .safe_tm(downstream_left_inner)

  best_df_target <- downstream_left_inner
  best_df_tm <- df_tm

  iteration <- 0
  while (iteration < max_iterations) {
    if (base::is.na(df_tm)) break
    if (df_tm >= target_tm_min && base::abs(df_tm - target_tm_optimal) <= 1) break

    if (base::nchar(downstream_left_inner) >= max_target_length) break

    if (df_tm < target_tm_min) {
      # Extend right
      next_len <- base::nchar(downstream_left_inner) + 1
      if (next_len > base::nchar(downstream_seq)) break
      downstream_left_inner <- base::substring(downstream_seq, 1, next_len)
    } else if (df_tm > target_tm_max) {
      # Trim from right
      downstream_left_inner <- base::substring(downstream_left_inner, 1, base::nchar(downstream_left_inner) - 1)
    } else {
      # Within range but not optimal — fine-tune
      if (df_tm < target_tm_optimal && base::nchar(downstream_left_inner) < max_target_length) {
        next_len <- base::nchar(downstream_left_inner) + 1
        if (next_len > base::nchar(downstream_seq)) break
        downstream_left_inner <- base::substring(downstream_seq, 1, next_len)
      } else if (df_tm > target_tm_optimal && base::nchar(downstream_left_inner) > min_target_length) {
        downstream_left_inner <- base::substring(downstream_left_inner, 1, base::nchar(downstream_left_inner) - 1)
      } else {
        break
      }
    }

    df_tm <- .safe_tm(downstream_left_inner)
    if (!base::is.na(df_tm) && df_tm >= target_tm_min && df_tm <= target_tm_max) {
      best_df_target <- downstream_left_inner
      best_df_tm <- df_tm
    }
    iteration <- iteration + 1
  }

  # ===========================================================
  # 5. OVERLAP OPTIMIZATION between inner primers
  # ===========================================================
  max_overlap <- base::min(20,
                      max_primer_length - base::nchar(best_ur_target),
                      max_primer_length - base::nchar(best_df_target))
  tail_ur_len <- base::min(10, max_overlap)
  tail_df_len <- base::min(10, max_overlap)

  # tail_upstream_reverse = first N bp of downstream_seq (extends into downstream from junction)
  tail_ur <- base::substring(downstream_seq, 1, tail_ur_len)
  # tail_downstream_forward = last N bp of upstream_seq (extends into upstream from junction)
  tail_df <- base::substring(upstream_seq, base::nchar(upstream_seq) - tail_df_len + 1, base::nchar(upstream_seq))

  # Assemble inner primers
  ur_primer <- .rc(base::paste0(best_ur_target, tail_ur))
  df_primer <- base::paste0(tail_df, best_df_target)

  ur_tm_full <- .safe_tm(ur_primer)
  df_tm_full <- .safe_tm(df_primer)
  overlap_combined <- base::paste0(tail_ur, tail_df)
  ol_tm <- .safe_tm(overlap_combined)
  ol_len_sum <- base::nchar(tail_ur) + base::nchar(tail_df)

  best_ur <- base::list(primer = ur_primer, tm_target = best_ur_tm,
                   tm_full = ur_tm_full, overlap_tm = ol_tm)
  best_df <- base::list(primer = df_primer, tm_target = best_df_tm,
                   tm_full = df_tm_full, overlap_tm = ol_tm)

  # Iteratively extend overlap (alternate by 2bp) until ≥20bp AND Tm ≥55°C
  iteration <- 0
  while ((ol_len_sum < 20 || .na0(ol_tm) < 55) && iteration < max_iterations) {
    if (base::nchar(ur_primer) < max_primer_length && base::nchar(df_primer) < max_primer_length) {
      if (iteration %% 2 == 0 && tail_ur_len < max_overlap) {
        tail_ur_len <- tail_ur_len + 2
        tail_ur <- base::substring(downstream_seq, 1, base::min(tail_ur_len, base::nchar(downstream_seq)))
      } else if (tail_df_len < max_overlap) {
        tail_df_len <- tail_df_len + 2
        tail_df <- base::substring(upstream_seq,
                              base::max(1, base::nchar(upstream_seq) - tail_df_len + 1),
                              base::nchar(upstream_seq))
      }

      ur_primer <- .rc(base::paste0(best_ur_target, tail_ur))
      df_primer <- base::paste0(tail_df, best_df_target)

      ur_tm_full <- .safe_tm(ur_primer)
      df_tm_full <- .safe_tm(df_primer)
      overlap_combined <- base::paste0(tail_ur, tail_df)
      ol_tm <- .safe_tm(overlap_combined)
      ol_len_sum <- base::nchar(tail_ur) + base::nchar(tail_df)

      if (base::nchar(ur_primer) <= max_primer_length &&
          base::nchar(df_primer) <= max_primer_length &&
          !.check_self_dimerization(ur_primer) &&
          !.check_self_dimerization(df_primer)) {
        best_ur <- base::list(primer = ur_primer, tm_target = best_ur_tm,
                         tm_full = ur_tm_full, overlap_tm = ol_tm)
        best_df <- base::list(primer = df_primer, tm_target = best_df_tm,
                         tm_full = df_tm_full, overlap_tm = ol_tm)
        if (!base::is.na(ol_tm) && ol_tm >= 55 && ol_len_sum >= 20) break
      }
    }
    iteration <- iteration + 1
  }

  # ===========================================================
  # Build result list
  # ===========================================================
  # Generate primer names: [locus_tag]_UF, _UR, _DF, _DR
  tag_label <- base::ifelse(!base::is.null(locus_tag), locus_tag, "unknown")

  result <- base::list(
    # Primer names
    upstream_forward_name   = base::paste0(tag_label, "_UF"),
    upstream_reverse_name   = base::paste0(tag_label, "_UR"),
    downstream_forward_name = base::paste0(tag_label, "_DF"),
    downstream_reverse_name = base::paste0(tag_label, "_DR"),
    # Primer sequences
    upstream_forward_primer   = best_upstream_fwd$primer,
    upstream_forward_tm_target = .safe_round(best_upstream_fwd$tm_target),
    upstream_forward_tm_full   = .safe_round(best_upstream_fwd$tm_full),
    upstream_reverse_primer   = best_ur$primer,
    upstream_reverse_tm_target = .safe_round(best_ur$tm_target),
    upstream_reverse_tm_full   = .safe_round(best_ur$tm_full),
    downstream_forward_primer = best_df$primer,
    downstream_forward_tm_target = .safe_round(best_df$tm_target),
    downstream_forward_tm_full   = .safe_round(best_df$tm_full),
    downstream_reverse_primer   = best_downstream_rev$primer,
    downstream_reverse_tm_target = .safe_round(best_downstream_rev$tm_target),
    downstream_reverse_tm_full   = .safe_round(best_downstream_rev$tm_full),
    # Positions in modified vector (1-based)
    upstream_forward_start = start_0 - base::nchar(upstream_vector_overlap) + 1,
    upstream_forward_end   = start_0 + base::nchar(upstream_left),
    upstream_reverse_start = start_0 + upstream_bp - base::nchar(best_ur_target) + 1,
    upstream_reverse_end   = start_0 + upstream_bp + base::nchar(tail_ur),
    downstream_forward_start = start_0 + upstream_bp - base::nchar(tail_df) + 1,
    downstream_forward_end   = start_0 + upstream_bp + base::nchar(best_df_target),
    downstream_reverse_start = start_0 + base::nchar(insert_seq) - base::nchar(downstream_right) + 1,
    downstream_reverse_end   = start_0 + base::nchar(insert_seq) + base::nchar(downstream_vector_overlap),
    # Overlap info
    overlap_tm     = .safe_round(ol_tm),
    overlap_length = ol_len_sum,
    # Metadata
    vector_name = vector$name,
    locus_tag   = locus_tag
  )

  base::cat("\n=== Deletion Primer Design Results ===\n")
  if (!base::is.null(locus_tag)) base::cat("Locus tag:", locus_tag, "\n")
  base::cat("  ", result$upstream_forward_name, ":   Tm(target)=", result$upstream_forward_tm_target,
      "°C, Tm(full)=", result$upstream_forward_tm_full, "°C,",
      base::nchar(best_upstream_fwd$primer), "bp\n")
  base::cat("  ", result$upstream_reverse_name, ":   Tm(target)=", result$upstream_reverse_tm_target,
      "°C, Tm(full)=", result$upstream_reverse_tm_full, "°C,",
      base::nchar(best_ur$primer), "bp\n")
  base::cat("  ", result$downstream_forward_name, ": Tm(target)=", result$downstream_forward_tm_target,
      "°C, Tm(full)=", result$downstream_forward_tm_full, "°C,",
      base::nchar(best_df$primer), "bp\n")
  base::cat("  ", result$downstream_reverse_name, ": Tm(target)=", result$downstream_reverse_tm_target,
      "°C, Tm(full)=", result$downstream_reverse_tm_full, "°C,",
      base::nchar(best_downstream_rev$primer), "bp\n")
  base::cat("  Overlap:", result$overlap_length, "bp, Tm=", result$overlap_tm, "°C\n")

  # --- GenBank output ---
  if (!base::is.null(output_dir)) {
    tag_label <- base::ifelse(!base::is.null(locus_tag), locus_tag, "unknown")
    gbk_path <- base::file.path(output_dir, base::paste0(tag_label, "_deletion_vector.gbk"))
    base::tryCatch({
      write_deletion_genbank(
        vector_file    = vector_file,
        insert_seq     = insert_seq,
        start          = start,
        end            = end,
        locus_tag      = tag_label,
        primers        = result,
        upstream_bp    = upstream_bp,
        downstream_bp  = downstream_bp,
        output_path    = gbk_path,
        kill_snapgene  = kill_snapgene
      )
    }, error = function(e) {
      base::warning("Failed to write GenBank file: ", base::conditionMessage(e))
    })
  }

  return(result)
}


#' Batch Design Deletion Primers for Multiple Locus Tags
#'
#' Wrapper that calls \code{design_deletion_primers} for each locus tag and
#' returns a combined data frame. Extracts upstream/downstream arms from the
#' genome sequence based on gene coordinates in a genbank table.
#'
#' @param vector_file Path to the vector file (.dna, .gb, .fasta).
#' @param genome_seq Character string of the full genome sequence.
#'   Not required if \code{genbank_file} is provided.
#' @param genbank_table Data frame with columns: \code{locus_tag}, \code{start}, \code{end}.
#'   Optionally: \code{gene}, \code{product}. Not required if \code{genbank_file} is provided.
#' @param locus_tags Character vector of locus tags to process. Use "all" for all tags.
#' @param start Integer. Vector insertion site start (1-based).
#' @param end Integer. Vector insertion site end (1-based).
#' @param upstream_bp Integer. Upstream arm length in bp (default: 500).
#' @param downstream_bp Integer. Downstream arm length in bp (default: 500).
#' @param tm_target Numeric. Target Tm (default: 60).
#' @param overlap_length Integer. Initial vector overlap (default: 16).
#' @param output_file Character. Path to save Excel output (default: NULL).
#' @param output_dir Character. Directory to save GenBank files per locus tag (default: NULL).
#'   Each locus tag gets a \code{{tag}_deletion_vector.gbk} file with primer annotations.
#' @param genbank_file Character string. Path to a genome GenBank file (.gb, .gbk, .gbff).
#'   When provided, \code{genome_seq} and \code{genbank_table} are automatically extracted
#'   from this file. This overrides any separately provided \code{genome_seq} and \code{genbank_table}.
#' @param kill_snapgene Logical. Kill SnapGene before reading .dna files (default: FALSE).
#' @return A data frame with one row per locus tag, containing all 4 primer
#'   sequences, Tm values, and metadata.
#' @examples
#' \dontrun{
#' # With separate genome_seq + genbank_table
#' result_df <- batch_deletion_primers(
#'   vector_file = "my_vector.dna",
#'   genome_seq = genome_sequence,
#'   genbank_table = genbank_df,
#'   locus_tags = c("QT234_RS00005", "QT234_RS00010"),
#'   start = 1234, end = 1256,
#'   output_file = "deletion_primers.xlsx",
#'   output_dir = "deletion_results",
#'   kill_snapgene = TRUE
#' )
#'
#' # With genome GenBank file (auto-extract genome_seq + genbank_table)
#' result_df <- batch_deletion_primers(
#'   vector_file = "my_vector.dna",
#'   genbank_file = "my_genome.gbk",
#'   locus_tags = c("QT234_RS00005", "QT234_RS00010"),
#'   start = 1234, end = 1256,
#'   output_file = "deletion_primers.xlsx",
#'   output_dir = "deletion_results"
#' )
#' }
#' @export
batch_deletion_primers <- function(
    vector_file,
    genome_seq = NULL,
    genbank_table = NULL,
    locus_tags,
    start,
    end,
    upstream_bp = 500,
    downstream_bp = 500,
    tm_target = 60,
    overlap_length = 16,
    output_file = NULL,
    output_dir = NULL,
    genbank_file = NULL,
    kill_snapgene = TRUE
) {
  # --- Auto-extract from genbank_file if provided ---
  if (!base::is.null(genbank_file)) {
    base::cat("Parsing genome GenBank file for genome_seq and genbank_table...\n")
    parsed <- read_genome_genbank(genbank_file)
    genome_seq <- parsed$genome_seq
    genbank_table <- parsed$genbank_table
    base::cat("  Genome:", parsed$length, "bp,", base::nrow(genbank_table),
              "features extracted\n")
  }

  # Validate inputs
  if (base::is.null(genome_seq) || !base::is.character(genome_seq) ||
      base::nchar(genome_seq) == 0) {
    base::stop("genome_seq is required. Provide it directly or via genbank_file parameter.")
  }
  if (base::is.null(genbank_table) || !base::is.data.frame(genbank_table)) {
    base::stop("genbank_table is required. Provide it directly or via genbank_file parameter.")
  }

  required_cols <- base::c("locus_tag", "start", "end")
  missing_cols <- base::setdiff(required_cols, base::colnames(genbank_table))
  if (base::length(missing_cols) > 0) {
    base::stop("genbank_table missing required columns: ", base::paste(missing_cols, collapse = ", "))
  }

  if (base::length(locus_tags) == 1 && base::tolower(locus_tags) == "all") {
    locus_tags <- base::unique(genbank_table$locus_tag[!base::is.na(genbank_table$locus_tag)])
    base::cat("Processing all", base::length(locus_tags), "locus tags\n")
  }

  genome_seq <- base::toupper(genome_seq)
  seq_length <- base::nchar(genome_seq)
  results <- base::list()

  for (tag in locus_tags) {
    base::tryCatch({
      base::cat("\n--- Processing:", tag, "---\n")

      # Handle range notation (e.g., "TAG1-TAG2")
      if (base::grepl("-", tag)) {
        parts <- base::strsplit(tag, "-", fixed = TRUE)[[1]]
        start_row <- genbank_table[genbank_table$locus_tag == parts[1], ]
        end_row <- genbank_table[genbank_table$locus_tag == parts[2], ]
        if (base::nrow(start_row) == 0 || base::nrow(end_row) == 0) {
          base::warning("Locus tag(s) not found: ", tag)
          next
        }
        gene_start <- base::as.integer(start_row$start[1])
        gene_end <- base::as.integer(end_row$end[1])
      } else {
        target_row <- genbank_table[genbank_table$locus_tag == tag, ]
        if (base::nrow(target_row) == 0) {
          base::warning("Locus tag not found: ", tag)
          next
        }
        gene_start <- base::as.integer(target_row$start[1])
        gene_end <- base::as.integer(target_row$end[1])
      }

      # Extract upstream arm (handle circular genome)
      up_end <- gene_start - 1
      up_start <- gene_start - upstream_bp
      if (up_start <= 0) {
        up_start_circ <- seq_length + up_start
        upstream_arm <- base::paste0(
          base::substring(genome_seq, up_start_circ, seq_length),
          base::substring(genome_seq, 1, up_end)
        )
      } else {
        upstream_arm <- base::substring(genome_seq, up_start, up_end)
      }

      # Extract downstream arm (handle circular genome)
      dn_start <- gene_end + 1
      dn_end <- gene_end + downstream_bp
      if (dn_end > seq_length) {
        dn_end_circ <- dn_end %% seq_length
        downstream_arm <- base::paste0(
          base::substring(genome_seq, dn_start, seq_length),
          base::substring(genome_seq, 1, dn_end_circ)
        )
      } else {
        downstream_arm <- base::substring(genome_seq, dn_start, dn_end)
      }

      insert_seq <- base::paste0(upstream_arm, downstream_arm)

      # Design primers (+ GenBank output if output_dir specified)
      primers <- design_deletion_primers(
        vector_file = vector_file,
        insert_seq = insert_seq,
        start = start,
        end = end,
        upstream_bp = upstream_bp,
        downstream_bp = downstream_bp,
        tm_target = tm_target,
        overlap_length = overlap_length,
        locus_tag = tag,
        output_dir = output_dir,
        kill_snapgene = kill_snapgene
      )

      # Build result row
      gene_name <- ""
      product_name <- ""
      if ("gene" %in% base::colnames(genbank_table) && !base::grepl("-", tag)) {
        gene_name <- base::as.character(target_row$gene[1])
        if (base::is.na(gene_name)) gene_name <- ""
      }
      if ("product" %in% base::colnames(genbank_table) && !base::grepl("-", tag)) {
        product_name <- base::as.character(target_row$product[1])
        if (base::is.na(product_name)) product_name <- ""
      }

      row <- base::data.frame(
        locus_tag = tag,
        gene = gene_name,
        product = product_name,
        upstream_forward_name = primers$upstream_forward_name,
        upstream_forward_primer = primers$upstream_forward_primer,
        upstream_forward_tm_target = primers$upstream_forward_tm_target,
        upstream_forward_tm_full = primers$upstream_forward_tm_full,
        upstream_reverse_name = primers$upstream_reverse_name,
        upstream_reverse_primer = primers$upstream_reverse_primer,
        upstream_reverse_tm_target = primers$upstream_reverse_tm_target,
        upstream_reverse_tm_full = primers$upstream_reverse_tm_full,
        downstream_forward_name = primers$downstream_forward_name,
        downstream_forward_primer = primers$downstream_forward_primer,
        downstream_forward_tm_target = primers$downstream_forward_tm_target,
        downstream_forward_tm_full = primers$downstream_forward_tm_full,
        downstream_reverse_name = primers$downstream_reverse_name,
        downstream_reverse_primer = primers$downstream_reverse_primer,
        downstream_reverse_tm_target = primers$downstream_reverse_tm_target,
        downstream_reverse_tm_full = primers$downstream_reverse_tm_full,
        overlap_tm = primers$overlap_tm,
        overlap_length = primers$overlap_length,
        stringsAsFactors = FALSE
      )
      results[[tag]] <- row

    }, error = function(e) {
      base::warning("Error processing locus_tag ", tag, ": ", base::conditionMessage(e))
    })
  }

  if (base::length(results) == 0) {
    base::warning("No deletion primers were successfully designed.")
    return(base::data.frame())
  }

  result_df <- base::do.call(base::rbind, results)
  base::rownames(result_df) <- NULL

  base::cat("\n=== Batch complete:", base::nrow(result_df), "of", base::length(locus_tags), "locus tags processed ===\n")

  # Save to Excel
  if (!base::is.null(output_file)) {
    base::tryCatch({
      excel_dir <- base::dirname(output_file)
      if (!base::dir.exists(excel_dir) && excel_dir != ".") {
        base::dir.create(excel_dir, recursive = TRUE)
      }
      if (base::requireNamespace("openxlsx", quietly = TRUE)) {
        openxlsx::write.xlsx(result_df, file = output_file, overwrite = TRUE)
      } else if (base::requireNamespace("writexl", quietly = TRUE)) {
        writexl::write_xlsx(result_df, path = output_file)
      }
      base::cat("Results saved to:", output_file, "\n")
    }, error = function(e) {
      base::warning("Failed to save Excel file: ", base::conditionMessage(e))
    })
  }

  return(result_df)
}


# ============================================================
# Helpers
# ============================================================

#' @keywords internal
.safe_tm <- function(seq) {
  base::tryCatch({
    if (base::is.na(seq) || base::nchar(seq) < 6) return(NA_real_)
    TmCalculator::Tm_NN(seq, Na = 50)$Tm
  }, error = function(e) NA_real_)
}

#' @keywords internal
.safe_get <- function(df, col, i) {
  if (col %in% base::colnames(df)) {
    val <- df[[col]][i]
    if (base::is.na(val) || val == "") return("")
    return(base::as.character(val))
  }
  return("")
}

#' Find Type IIS enzyme recognition sites and compute overhangs + fill sequences
#'
#' Scans +/-20bp around the stuffer region in a vector for enzyme recognition sites
#' on BOTH strands. Determines which orientation has the enzyme cutting TOWARD the
#' stuffer, then computes overhangs and fill sequences.
#'
#' Cut pattern handling:
#'   (a/b) notation: a = recognition-strand cut offset, b = complementary-strand offset.
#'   - a < b (e.g. BbsI 2/6): 5' overhang on complementary strand, length = b - a
#'   - a > b (e.g. hypothetical 6/2): 3' overhang on recognition strand, length = a - b
#'   Overhang length = abs(cut_rev - cut_fwd). Type depends on which is larger.
#'
#' For the LEFT site (upstream): sense recognition → cuts rightward (toward stuffer).
#' For the RIGHT site (downstream): antisense recognition → cuts leftward (toward stuffer).
#' Both strands are searched at both positions; only the correct orientation is used.
#'
#' @param vec_seq Character string. Full vector sequence (uppercase).
#' @param start Integer. Stuffer start position (1-based).
#' @param end Integer. Stuffer end position (1-based).
#' @param enzyme_info List with elements: recognition, cut_fwd, cut_rev.
#' @param enzyme_name Character. Enzyme name for messages.
#' @return List with: left_oh, right_oh, fill_left, fill_right, found (logical),
#'   left_site_pos, right_site_pos, overhang_type ("5prime" or "3prime"),
#'   overhang_len (integer).
#' @keywords internal
.find_enzyme_sites_and_compute_overhangs <- function(vec_seq, start, end, enzyme_info, enzyme_name) {
  recognition <- enzyme_info$recognition
  recog_len   <- base::nchar(recognition)
  cut_fwd     <- enzyme_info$cut_fwd   # recognition-strand cut offset from recog end
  cut_rev     <- enzyme_info$cut_rev   # complementary-strand cut offset from recog end
  vec_len     <- base::nchar(vec_seq)
  oh_len      <- base::abs(cut_rev - cut_fwd)

  # Determine overhang type from cut pattern
  if (cut_fwd < cut_rev) {
    # e.g. BbsI (2/6): recognition strand cut closer → 5' overhang on complementary strand
    overhang_type <- "5prime"
    cut_near <- cut_fwd   # closer cut (on recognition strand)
    cut_far  <- cut_rev   # farther cut (on complementary strand)
  } else if (cut_fwd > cut_rev) {
    # Hypothetical (6/2): recognition strand cut farther → 3' overhang
    overhang_type <- "3prime"
    cut_near <- cut_rev   # closer cut (on complementary strand)
    cut_far  <- cut_fwd   # farther cut (on recognition strand)
  } else {
    # cut_fwd == cut_rev → blunt end, no overhang
    base::warning(enzyme_name, " has cut_fwd == cut_rev (", cut_fwd,
                  ") → blunt end. Cannot use for Golden Gate oligo annealing.")
    return(base::list(found = FALSE))
  }

  base::cat("  [OA] ", enzyme_name, ": ", recognition, " (", cut_fwd, "/", cut_rev,
            ") → ", oh_len, "nt ", overhang_type, " overhang\n", sep = "")

  rc_recognition <- base::as.character(
    Biostrings::reverseComplement(Biostrings::DNAString(recognition))
  )

  # ============================================================
  # LEFT SITE: search BOTH strands in upstream region
  # ============================================================
  left_search_start <- base::max(1L, start - 20L)
  left_search_end   <- start - 1L
  left_region <- base::substr(vec_seq, left_search_start, left_search_end)

  # Search sense strand: recognition found → enzyme cuts RIGHTWARD (toward stuffer) ✓
  left_sense_hits <- base::gregexpr(recognition, left_region, fixed = TRUE)[[1]]
  # Search antisense strand: RC found → enzyme cuts LEFTWARD (away from stuffer) ✗
  left_anti_hits  <- base::gregexpr(rc_recognition, left_region, fixed = TRUE)[[1]]

  left_found <- FALSE
  p <- NA_integer_

  if (left_sense_hits[1] > 0) {
    # Sense recognition found upstream → cuts toward stuffer ✓
    best_hit <- left_sense_hits[base::length(left_sense_hits)]
    p <- left_search_start + best_hit - 1L
    left_found <- TRUE
    base::cat("    Left: sense ", recognition, " at pos ", p,
              " → cuts rightward toward stuffer\n", sep = "")
  }

  if (left_anti_hits[1] > 0 && !left_found) {
    # Only RC found upstream → enzyme cuts AWAY from stuffer
    anti_pos <- left_search_start + left_anti_hits[base::length(left_anti_hits)] - 1L
    base::warning(enzyme_name, " RC site (", rc_recognition, ") found at pos ", anti_pos,
                  " UPSTREAM of stuffer — enzyme cuts AWAY from stuffer. ",
                  "Check vector design. Site ignored.")
  }

  # ============================================================
  # RIGHT SITE: search BOTH strands in downstream region
  # ============================================================
  right_search_start <- end + 1L
  right_search_end   <- base::min(vec_len, end + 20L)
  right_region <- base::substr(vec_seq, right_search_start, right_search_end)

  # Search antisense: RC on sense → enzyme cuts LEFTWARD (toward stuffer) ✓
  right_anti_hits  <- base::gregexpr(rc_recognition, right_region, fixed = TRUE)[[1]]
  # Search sense: recognition found → enzyme cuts RIGHTWARD (away from stuffer) ✗
  right_sense_hits <- base::gregexpr(recognition, right_region, fixed = TRUE)[[1]]

  right_found <- FALSE
  q <- NA_integer_

  if (right_anti_hits[1] > 0) {
    # RC found downstream → cuts toward stuffer ✓
    best_hit <- right_anti_hits[1]
    q <- right_search_start + best_hit - 1L
    right_found <- TRUE
    base::cat("    Right: antisense (", rc_recognition, " on sense) at pos ", q,
              " → cuts leftward toward stuffer\n", sep = "")
  }

  if (right_sense_hits[1] > 0 && !right_found) {
    # Only sense recognition found downstream → enzyme cuts AWAY from stuffer
    sense_pos <- right_search_start + right_sense_hits[1] - 1L
    base::warning(enzyme_name, " sense site (", recognition, ") found at pos ", sense_pos,
                  " DOWNSTREAM of stuffer — enzyme cuts AWAY from stuffer. ",
                  "Check vector design. Site ignored.")
  }

  # --- Check if both sites found ---
  if (!left_found || !right_found) {
    if (!left_found) {
      base::warning(enzyme_name, " recognition site (", recognition,
                    ") not found on sense strand within ",
                    left_search_end - left_search_start + 1,
                    "bp upstream of stuffer (pos ", left_search_start, "-",
                    left_search_end, "). Falling back to default overhangs.")
    }
    if (!right_found) {
      base::warning(enzyme_name, " RC recognition (", rc_recognition,
                    ") not found on sense strand within ",
                    right_search_end - right_search_start + 1,
                    "bp downstream of stuffer (pos ", right_search_start, "-",
                    right_search_end, "). Falling back to default overhangs.")
    }
    return(base::list(found = FALSE))
  }

  # ============================================================
  # Compute overhangs and fill based on overhang type
  # ============================================================
  #
  # LEFT site: sense recognition at position p (1-based)
  #   Sense strand cut at:       p + recog_len + cut_fwd - 1 (last retained base on sense)
  #   Antisense strand cut at:   p + recog_len + cut_rev - 1 (last retained base on anti)
  #
  # RIGHT site: RC on sense at position q (1-based), recognition on antisense
  #   Sense strand cut at:       q - cut_rev (first retained base on sense, backbone side)
  #   Antisense strand cut at:   q - cut_fwd (first retained base on antisense)
  #
  # For (a/b) with a < b (5' overhang, e.g. BbsI 2/6):
  #   LEFT:  antisense extends further → 5' bottom overhang
  #          overhang on sense: positions [p+recog_len+cut_fwd] to [p+recog_len+cut_rev-1]
  #   RIGHT: sense extends further → 5' top overhang
  #          overhang on sense: positions [q-cut_rev] to [q-cut_fwd-1]
  #
  # For (a/b) with a > b (3' overhang, hypothetical 6/2):
  #   LEFT:  sense extends further → 3' sense extension
  #          overhang on sense: positions [p+recog_len+cut_rev] to [p+recog_len+cut_fwd-1]
  #   RIGHT: antisense extends further → 3' antisense extension
  #          overhang on sense: positions [q-cut_fwd] to [q-cut_rev-1]

  if (overhang_type == "5prime") {
    # --- 5' overhang case (standard: BbsI 2/6, BsaI 1/5, etc.) ---

    # LEFT overhang: sense positions between the two cuts
    left_oh_start <- p + recog_len + cut_fwd
    left_oh_end   <- p + recog_len + cut_rev - 1L
    left_oh <- base::substr(vec_seq, left_oh_start, left_oh_end)

    # LEFT fill: bases removed between overhang end and stuffer start
    fill_left_start <- p + recog_len + cut_rev
    fill_left_end   <- start - 1L
    fill_left <- ""
    if (fill_left_start <= fill_left_end) {
      fill_left <- base::substr(vec_seq, fill_left_start, fill_left_end)
    }

    # RIGHT overhang: sense positions between the two cuts (then RC for R oligo)
    right_oh_sense_start <- q - cut_rev
    right_oh_sense_end   <- q - cut_fwd - 1L
    right_oh_sense <- base::substr(vec_seq, right_oh_sense_start, right_oh_sense_end)
    right_oh <- base::as.character(
      Biostrings::reverseComplement(Biostrings::DNAString(right_oh_sense))
    )

    # RIGHT fill: bases removed between stuffer end and right overhang start
    fill_right_start <- end + 1L
    fill_right_end   <- q - cut_rev - 1L
    fill_right <- ""
    if (fill_right_start <= fill_right_end) {
      fill_right <- base::substr(vec_seq, fill_right_start, fill_right_end)
    }

    # F oligo 5' = left_oh (pairs with backbone bottom 5' overhang)
    # R oligo 5' = right_oh = RC(right_oh_sense) (pairs with backbone top 5' overhang)

  } else {
    # --- 3' overhang case (a > b, e.g. hypothetical 6/2) ---
    # Here the overhang polarity is reversed:
    #   LEFT:  sense strand extends further (3' overhang on top)
    #          → insert needs bottom 5' extension at left = R oligo provides this
    #   RIGHT: antisense extends further (3' overhang on bottom)
    #          → insert needs top 5' extension at right = F oligo provides this
    #
    # For oligo annealing: we swap which oligo gets which overhang.
    # F oligo 5' = right overhang,  R oligo 5' = left overhang (swapped vs 5' case)

    # LEFT overhang on sense strand (sense extends further)
    left_oh_start <- p + recog_len + cut_rev
    left_oh_end   <- p + recog_len + cut_fwd - 1L
    left_oh_sense <- base::substr(vec_seq, left_oh_start, left_oh_end)
    # For 3' overhang at left: the insert's bottom strand needs the complement
    left_oh <- base::as.character(
      Biostrings::reverseComplement(Biostrings::DNAString(left_oh_sense))
    )

    # LEFT fill: bases between the far cut and stuffer start
    fill_left_start <- p + recog_len + cut_fwd
    fill_left_end   <- start - 1L
    fill_left <- ""
    if (fill_left_start <= fill_left_end) {
      fill_left <- base::substr(vec_seq, fill_left_start, fill_left_end)
    }

    # RIGHT overhang on sense strand (antisense extends further)
    right_oh_sense_start <- q - cut_fwd
    right_oh_sense_end   <- q - cut_rev - 1L
    right_oh_sense <- base::substr(vec_seq, right_oh_sense_start, right_oh_sense_end)
    # For 3' overhang at right: sense positions directly → F oligo
    right_oh <- right_oh_sense

    # RIGHT fill: bases between stuffer end and right overhang start
    fill_right_start <- end + 1L
    fill_right_end   <- q - cut_fwd - 1L
    fill_right <- ""
    if (fill_right_start <= fill_right_end) {
      fill_right <- base::substr(vec_seq, fill_right_start, fill_right_end)
    }

    # NOTE: For 3' overhang, the oligo construction in .design_golden_gate_oligos()
    # uses the same formula F = left_oh + fill + G + spacer + fill,
    # R = right_oh + RC(fill_right) + RC(spacer) + C + RC(fill_left)
    # but left_oh and right_oh are already adjusted here to go on the correct oligo.
    base::cat("    NOTE: 3' overhang pattern (", cut_fwd, "/", cut_rev,
              ") — overhangs adjusted for oligo annealing compatibility.\n", sep = "")
  }

  # Warn if fill is unusually long
  if (base::nchar(fill_left) > 5L) {
    base::warning("Left fill sequence is ", base::nchar(fill_left), "bp — unusually long. ",
                  "Check stuffer boundaries (start=", start, ").")
  }
  if (base::nchar(fill_right) > 5L) {
    base::warning("Right fill sequence is ", base::nchar(fill_right), "bp — unusually long. ",
                  "Check stuffer boundaries (end=", end, ").")
  }

  # Summary log
  base::cat("    → F overhang: ", left_oh, " (", base::nchar(left_oh), "nt)",
            base::ifelse(base::nchar(fill_left) > 0,
                         base::paste0(", fill_L=", fill_left), ""), "\n",
            "    → R overhang: ", right_oh, " (", base::nchar(right_oh), "nt)",
            base::ifelse(base::nchar(fill_right) > 0,
                         base::paste0(", fill_R=", fill_right), ""), "\n",
            sep = "")

  return(base::list(
    left_oh        = left_oh,
    right_oh       = right_oh,
    fill_left      = fill_left,
    fill_right     = fill_right,
    found          = TRUE,
    left_site_pos  = p,
    right_site_pos = q,
    overhang_type  = overhang_type,
    overhang_len   = oh_len
  ))
}


#' Reverse complement shorthand
#' @keywords internal
.rc <- function(seq_str) {
  base::as.character(Biostrings::reverseComplement(Biostrings::DNAString(base::toupper(seq_str))))
}

#' Safe NA-to-0 for comparisons
#' @keywords internal
.na0 <- function(x) {
  base::ifelse(base::is.na(x), 0, x)
}

#' Safe round that handles NA
#' @keywords internal
.safe_round <- function(x, digits = 2) {
  if (base::is.na(x)) return(NA_real_)
  base::round(x, digits)
}

#' Check self-dimerization of a primer
#'
#' Checks whether a primer sequence can form a self-dimer by aligning
#' with its own reverse complement. Returns TRUE if a dimerization
#' region of at least \code{threshold} bp is found.
#'
#' @param primer_seq Character string of the primer sequence.
#' @param threshold Integer. Minimum match length to report (default: 8).
#' @return Logical. TRUE if self-dimerization detected.
#' @keywords internal
.check_self_dimerization <- function(primer_seq, threshold = 8) {
  base::tryCatch({
    primer_seq <- base::toupper(primer_seq)
    rc <- .rc(primer_seq)
    # Sliding window check for complementary stretches
    p_len <- base::nchar(primer_seq)
    for (i in base::seq_len(p_len - threshold + 1)) {
      window <- base::substring(primer_seq, i, i + threshold - 1)
      if (base::grepl(window, rc, fixed = TRUE)) {
        return(TRUE)
      }
    }
    return(FALSE)
  }, error = function(e) FALSE)
}
