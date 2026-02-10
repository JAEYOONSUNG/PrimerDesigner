#' Design gRNA + Deletion Arm Primers (Unified Pipeline)
#'
#' All-in-one pipeline for CRISPR gene knockout primer design.
#' Combines gRNA expression construct primer design and deletion arm primer
#' design into a single function call.
#'
#' For each target gene, this function:
#' \enumerate{
#'   \item Generates and selects the best gRNA(s) (via \code{design_grna_construct})
#'   \item Designs gRNA cloning primers (Gibson or Golden Gate, 2 per gRNA)
#'   \item Designs deletion arm primers (4-primer Gibson for homologous recombination)
#'   \item Produces a single combined GenBank showing the final construct
#' }
#'
#' @section Combined vs Separate mode:
#' \describe{
#'   \item{Combined (all-in-one vector)}{When \code{deletion_vector_file} is NULL
#'     (same vector for both), \code{deletion_start/end} must be different from
#'     \code{grna_start/end}. The gRNA spacer is inserted first, coordinates
#'     are auto-adjusted, and deletion primers are designed against the
#'     intermediate (gRNA-modified) vector. One combined GenBank per gRNA.}
#'   \item{Separate vectors}{When \code{deletion_vector_file} is a different file,
#'     each construct gets its own GenBank file.}
#' }
#'
#' @param gRNA_df Optional. Pre-generated gRNA library data frame.
#' @param genbank_file Character. Path to local genome GenBank file (.gbff, .gb).
#' @param genbank_accession Character. GenBank accession for auto-generation.
#' @param locus_tags Character vector. Target locus tag(s).
#' @param nuclease Character. Nuclease name (default: "GeoCas9").
#' @param top_n Integer. Number of best gRNAs to select (default: 1).
#' @param position_range Numeric vector. gRNA position range (default: c(5, 95)).
#' @param tm_range Numeric vector. gRNA Tm filter range (default: c(50, 70)).
#' @param strand Character. Strand filter (default: "both").
#' @param methylation_patterns Character vector. IUPAC methylation motifs.
#' @param dir Character. Working directory for auto-generation.
#'
#' @param grna_vector_file Character. Path to gRNA expression vector.
#' @param grna_start Integer. gRNA stuffer start position (1-based).
#' @param grna_end Integer. gRNA stuffer end position (1-based).
#' @param grna_cloning_method Character. "gibson" or "golden_gate".
#' @param grna_enzyme Character. Golden Gate enzyme.
#' @param grna_custom_overhangs Named list. Custom overhangs.
#'
#' @param deletion_vector_file Character or NULL. NULL = same as gRNA vector.
#' @param deletion_start Integer. Deletion arm stuffer start (1-based).
#'   Required when using all-in-one vector. Must differ from grna_start.
#' @param deletion_end Integer. Deletion arm stuffer end (1-based).
#'   Required when using all-in-one vector. Must differ from grna_end.
#'
#' @param upstream_bp Integer. Upstream arm length (default: 500).
#' @param downstream_bp Integer. Downstream arm length (default: 500).
#' @param tm_target Numeric. Target Tm (default: 60).
#' @param overlap_length Integer. Gibson overlap length (default: 20).
#' @param add_g Logical. Add 5' G for transcription (default: TRUE).
#' @param grna_name_pattern Character. Primer naming pattern.
#' @param output_file Character or NULL. Excel output (.xlsx).
#' @param output_dir Character or NULL. GenBank output directory.
#' @param kill_snapgene Logical. Kill SnapGene before reading .dna.
#' @param resume Logical. If TRUE, skip Phase 1 by loading cached gRNA results
#'   from \code{output_dir}, and skip locus tags whose GenBank files already exist.
#'   Useful when a previous run was interrupted. (default: FALSE)
#'
#' @return A list with three components:
#'   \describe{
#'     \item{grna_primers}{Wide-format gRNA primer data (2 per gRNA).}
#'     \item{deletion_primers}{Wide-format deletion arm primer data (4 per locus_tag).}
#'     \item{unified_primers}{Long-format unified primer table (all primers):
#'       locus_tag, type (sgRNA/homologous_arm), direction (F/R),
#'       primer_name, sequence, Tm_target, Tm_full, length, ...}
#'   }
#'
#' @examples
#' \dontrun{
#' # All-in-one vector: two different stuffer sites
#' result <- design_grna_and_deletion(
#'   genbank_file = "~/genomes/GCF_030376765.1.gbff",
#'   locus_tags = "QT235_RS00005",
#'   nuclease = "GeoCas9",
#'   grna_vector_file = "~/vectors/pG1Kt-GeoCas9EF-OA-sfGFP-ACrec.dna",
#'   grna_start = 8811, grna_end = 8840,          # gRNA stuffer
#'   deletion_start = 5001, deletion_end = 5030,   # arm stuffer
#'   grna_cloning_method = "gibson",
#'   upstream_bp = 500, downstream_bp = 500,
#'   output_file = "all_primers.xlsx",
#'   output_dir  = "constructs/"
#' )
#' }
#' @export
design_grna_and_deletion <- function(
    # ========== Genome & gRNA Selection ==========
    gRNA_df = NULL,
    genbank_file = NULL,
    genbank_accession = NULL,
    locus_tags,
    nuclease = "GeoCas9",
    top_n = 1,
    position_range = c(5, 95),
    tm_range = c(50, 70),
    strand = "both",
    methylation_patterns = NULL,
    dir = base::getwd(),

    # ========== gRNA Expression Vector ==========
    grna_vector_file,
    grna_start,
    grna_end,
    grna_cloning_method = "gibson",
    grna_enzyme = "BbsI",
    grna_custom_overhangs = NULL,

    # ========== Deletion Vector / Stuffer Site ==========
    deletion_vector_file = NULL,
    deletion_start = NULL,
    deletion_end = NULL,

    # ========== Deletion Arm Parameters ==========
    upstream_bp = 500,
    downstream_bp = 500,

    # ========== Shared Primer Parameters ==========
    tm_target = 60,
    overlap_length = 20,
    add_g = TRUE,
    grna_name_pattern = "sgRNA_{locus_tag}_{method}_{dir}",

    # ========== Output ==========
    output_file = NULL,
    output_dir = NULL,
    kill_snapgene = TRUE,
    resume = FALSE
) {

  # ============================================================
  # Phase 0: Validation & Mode Detection
  # ============================================================
  base::cat("============================================================\n")
  base::cat(" Unified gRNA + Deletion Arm Primer Design\n")
  base::cat("============================================================\n\n")

  # --- Smart detection: genbank_accession may be a file path ---
  if (!base::is.null(genbank_accession) && base::is.null(genbank_file)) {
    if (base::grepl("\\.(gbff|gb|gbk)$", genbank_accession, ignore.case = TRUE) ||
        (base::file.exists(genbank_accession) &&
         !base::grepl("^(GCF|GCA)_", genbank_accession))) {
      base::cat("  Auto-detected: genbank_accession is a file path, treating as genbank_file\n")
      genbank_file <- genbank_accession
      genbank_accession <- NULL
    }
  }

  if (base::missing(grna_vector_file) || !base::is.character(grna_vector_file)) {
    base::stop("grna_vector_file is required.")
  }
  if (!base::file.exists(grna_vector_file)) {
    base::stop("grna_vector_file not found: ", grna_vector_file)
  }
  if (base::missing(grna_start) || base::missing(grna_end)) {
    base::stop("grna_start and grna_end are required.")
  }
  if (base::missing(locus_tags)) {
    base::stop("locus_tags is required.")
  }

  # --- Determine mode: combined (same vector) or separate ---
  same_vector <- base::is.null(deletion_vector_file)
  if (!same_vector) {
    same_vector <- base::tryCatch({
      base::normalizePath(grna_vector_file, mustWork = FALSE) ==
        base::normalizePath(deletion_vector_file, mustWork = FALSE)
    }, error = function(e) FALSE)
  }

  if (same_vector) {
    # All-in-one vector: require separate deletion stuffer site
    if (base::is.null(deletion_start) || base::is.null(deletion_end)) {
      base::stop(
        "All-in-one vector mode: deletion_start and deletion_end are required.\n",
        "  These specify the arm insertion stuffer site (must differ from gRNA site)."
      )
    }
    if (deletion_start == grna_start && deletion_end == grna_end) {
      base::stop(
        "gRNA and deletion stuffer sites cannot be identical on one vector.\n",
        "  gRNA site:     ", grna_start, "-", grna_end, "\n",
        "  Deletion site: ", deletion_start, "-", deletion_end
      )
    }
    # Check non-overlap
    if (!(deletion_end < grna_start || deletion_start > grna_end)) {
      base::stop(
        "gRNA and deletion stuffer sites must not overlap.\n",
        "  gRNA site:     ", grna_start, "-", grna_end, "\n",
        "  Deletion site: ", deletion_start, "-", deletion_end
      )
    }
    combined_mode <- TRUE
    deletion_vector_file <- grna_vector_file
    base::cat("  Mode: COMBINED (all-in-one vector, two stuffer sites)\n")
  } else {
    combined_mode <- FALSE
    if (base::is.null(deletion_start)) deletion_start <- grna_start
    if (base::is.null(deletion_end))   deletion_end   <- grna_end
    base::cat("  Mode: SEPARATE (gRNA + deletion on different vectors)\n")
  }

  # --- Create output_dir ---
  if (!base::is.null(output_dir)) {
    base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # --- Summary ---
  grna_start_i  <- base::as.integer(grna_start)
  grna_end_i    <- base::as.integer(grna_end)
  del_start_i   <- base::as.integer(deletion_start)
  del_end_i     <- base::as.integer(deletion_end)

  base::cat("  Nuclease:         ", nuclease, "\n")
  base::cat("  Target:           ", base::paste(locus_tags, collapse = ", "), "\n")
  base::cat("  gRNA vector:      ", base::basename(grna_vector_file), "\n")
  base::cat("  gRNA stuffer:     ", grna_start_i, "-", grna_end_i, "\n")
  base::cat("  gRNA cloning:     ", grna_cloning_method, "\n")
  base::cat("  Deletion vector:  ", base::basename(deletion_vector_file), "\n")
  base::cat("  Deletion stuffer: ", del_start_i, "-", del_end_i, "\n")
  base::cat("  Arm lengths:      ", upstream_bp, "bp (up) / ", downstream_bp, "bp (down)\n")
  base::cat("\n")

  # ============================================================
  # Phase 1: gRNA Selection + Cloning Primers
  # ============================================================

  # --- Resume: try to load cached gRNA results ---
  cache_path <- NULL
  if (!base::is.null(output_dir)) {
    cache_path <- base::file.path(output_dir, ".grna_primers_cache.rds")
  }

  if (resume && !base::is.null(cache_path) && base::file.exists(cache_path)) {
    base::cat("=== Phase 1: Loading cached gRNA results (resume mode) ===\n\n")
    grna_primers <- base::readRDS(cache_path)
    base::cat("  Loaded ", base::nrow(grna_primers), " cached gRNA construct(s)\n\n", sep = "")
  } else {
    base::cat("=== Phase 1: gRNA Selection + Cloning Primers ===\n\n")

    grna_primers <- design_grna_construct(
      gRNA_df              = gRNA_df,
      genbank_file         = genbank_file,
      genbank_accession    = genbank_accession,
      vector_file          = grna_vector_file,
      start                = grna_start_i,
      end                  = grna_end_i,
      locus_tags           = locus_tags,
      nuclease             = nuclease,
      top_n                = top_n,
      cloning_method       = grna_cloning_method,
      enzyme               = grna_enzyme,
      custom_overhangs     = grna_custom_overhangs,
      tm_target            = tm_target,
      overlap_length       = overlap_length,
      add_g                = add_g,
      name_pattern         = grna_name_pattern,
      output_file          = NULL,
      output_dir           = if (!combined_mode && !base::is.null(output_dir))
                               base::file.path(output_dir, "gRNA_constructs") else NULL,
      kill_snapgene        = kill_snapgene,
      methylation_patterns = methylation_patterns,
      dir                  = dir,
      position_range       = position_range,
      tm_range             = tm_range,
      strand               = strand
    )

    base::cat("\n  gRNA primers: ", base::nrow(grna_primers), " construct(s)\n\n", sep = "")

    # Cache gRNA results for resume
    if (!base::is.null(cache_path)) {
      base::tryCatch(base::saveRDS(grna_primers, cache_path),
                      error = function(e) NULL)
    }
  }

  # ============================================================
  # Phase 2: Genome Data Extraction
  # ============================================================
  base::cat("=== Phase 2: Genome Data Extraction ===\n\n")

  resolved_gbff <- NULL
  if (!base::is.null(genbank_file) && base::file.exists(genbank_file)) {
    resolved_gbff <- genbank_file
  } else {
    # Look in seqs_srcdir (placed by Phase 1)
    accession_for_search <- genbank_accession
    if (base::is.null(accession_for_search) && !base::is.null(genbank_file)) {
      accession_for_search <- base::sub("\\.(gbff|gb|gbk)$", "",
                                          base::basename(genbank_file),
                                          ignore.case = TRUE)
    }
    if (!base::is.null(accession_for_search)) {
      candidate <- base::file.path(dir, "seqs_srcdir",
                                     base::paste0(accession_for_search, ".gbff"))
      if (base::file.exists(candidate)) resolved_gbff <- candidate
    }
  }
  if (base::is.null(resolved_gbff)) {
    base::stop("Cannot locate genome GenBank file for deletion arm extraction.\n",
               "  Provide genbank_file or genbank_accession.")
  }

  base::cat("  Parsing: ", base::basename(resolved_gbff), "\n")

  # Try to reuse cached genbank_table from Phase 1 (.rds)
  cached_gb_table <- NULL
  accession_id <- genbank_accession
  if (base::is.null(accession_id) && !base::is.null(genbank_file)) {
    accession_id <- base::sub("\\.(gbff|gb|gbk)$", "",
                                base::basename(genbank_file), ignore.case = TRUE)
  }
  if (!base::is.null(accession_id)) {
    rds_path <- base::file.path(dir, base::paste0(accession_id, "_genbank_table.rds"))
    if (base::file.exists(rds_path)) {
      cached_gb_table <- base::tryCatch(base::readRDS(rds_path), error = function(e) NULL)
      if (!base::is.null(cached_gb_table)) {
        base::cat("  Loaded cached genbank_table (", base::nrow(cached_gb_table),
                  " features) from Phase 1\n", sep = "")
      }
    }
  }

  if (!base::is.null(cached_gb_table)) {
    # Extract only genome sequence (skip feature parsing)
    genbank_table <- cached_gb_table
    genome_seq <- .extract_genome_seq_only(resolved_gbff)
  } else {
    parsed <- read_genome_genbank(resolved_gbff)
    genome_seq    <- parsed$genome_seq
    genbank_table <- parsed$genbank_table
  }
  genome_len <- base::nchar(genome_seq)
  base::cat("  Genome: ", genome_len, " bp, ",
            base::nrow(genbank_table), " features\n\n", sep = "")

  # ============================================================
  # Phase 3: Deletion Arm Primer Design
  # ============================================================
  base::cat("=== Phase 3: Deletion Arm Primer Design ===\n\n")

  # Resolve deletion locus tags
  deletion_tags <- locus_tags
  if (base::length(deletion_tags) == 1 &&
      base::tolower(deletion_tags[1]) == "all" &&
      "locus_tag" %in% base::colnames(grna_primers)) {
    deletion_tags <- base::unique(grna_primers$locus_tag)
  }

  if (combined_mode) {
    # ========== COMBINED MODE ==========
    # Read original vector once
    vec_info <- read_vector_file(grna_vector_file, kill_snapgene = kill_snapgene)
    original_vec_seq <- base::toupper(base::as.character(vec_info$sequence))
    original_vec_len <- base::nchar(original_vec_seq)

    # Determine insertion order: gRNA first, then adjust deletion coords
    grna_stuffer_len <- grna_end_i - grna_start_i + 1
    del_stuffer_len  <- del_end_i  - del_start_i  + 1

    # Track per-locus_tag deletion results (arms are same for all gRNAs of same tag)
    deletion_data <- base::list()

    # Resume: detect already-completed locus_tags from existing GenBank files
    completed_tags <- base::character(0)
    if (resume && !base::is.null(output_dir) && base::dir.exists(output_dir)) {
      existing_gbk <- base::list.files(output_dir, pattern = "_combined_construct\\.gbk$",
                                         full.names = FALSE)
      for (gbk_f in existing_gbk) {
        # Extract locus_tag from filename: {tag}_g{rank}_combined_construct.gbk
        m <- base::regmatches(gbk_f,
               base::regexpr("^(.+)_g\\d+_combined_construct\\.gbk$", gbk_f))
        if (base::length(m) > 0) {
          extracted_tag <- base::sub("_g\\d+_combined_construct\\.gbk$", "", m)
          completed_tags <- base::c(completed_tags, extracted_tag)
        }
      }
      completed_tags <- base::unique(completed_tags)
      if (base::length(completed_tags) > 0) {
        base::cat("  Resume: skipping ", base::length(completed_tags),
                  " already-completed tag(s): ",
                  base::paste(utils::head(completed_tags, 5), collapse = ", "),
                  if (base::length(completed_tags) > 5) ", ..." else "", "\n", sep = "")
      }
    }

    for (tag in deletion_tags) {
      # Skip if already completed (resume mode)
      if (tag %in% completed_tags) {
        base::cat("  --- ", tag, " --- SKIPPED (resume)\n")
        next
      }

      base::tryCatch({
        base::cat("  --- ", tag, " ---\n")

        target_row <- genbank_table[genbank_table$locus_tag == tag, ]
        if (base::nrow(target_row) == 0) {
          base::warning("Locus tag not found: ", tag)
          next
        }
        gene_start <- base::as.integer(target_row$start[1])
        gene_end   <- base::as.integer(target_row$end[1])

        # Extract upstream arm (circular genome aware)
        up_start <- gene_start - upstream_bp
        up_end   <- gene_start - 1
        if (up_start <= 0) {
          upstream_arm <- base::paste0(
            base::substring(genome_seq, genome_len + up_start, genome_len),
            base::substring(genome_seq, 1, up_end)
          )
        } else {
          upstream_arm <- base::substring(genome_seq, up_start, up_end)
        }

        # Extract downstream arm (circular genome aware)
        dn_start <- gene_end + 1
        dn_end   <- gene_end + downstream_bp
        if (dn_end > genome_len) {
          downstream_arm <- base::paste0(
            base::substring(genome_seq, dn_start, genome_len),
            base::substring(genome_seq, 1, dn_end %% genome_len)
          )
        } else {
          downstream_arm <- base::substring(genome_seq, dn_start, dn_end)
        }

        insert_seq <- base::paste0(upstream_arm, downstream_arm)

        # Use first gRNA for this tag to build intermediate vector
        # (all spacers of same nuclease have same length → same delta)
        tag_rows <- base::which(grna_primers$locus_tag == tag)
        spacer_example <- base::as.character(grna_primers$protospacer[tag_rows[1]])
        delta_grna <- base::nchar(spacer_example) - grna_stuffer_len

        # Build intermediate vector (gRNA spacer inserted)
        intermediate_seq <- base::paste0(
          base::substring(original_vec_seq, 1, grna_start_i - 1),
          spacer_example,
          base::substring(original_vec_seq, grna_end_i + 1, original_vec_len)
        )

        # Adjust deletion coordinates on intermediate vector
        if (del_start_i > grna_end_i) {
          adj_del_start <- del_start_i + delta_grna
          adj_del_end   <- del_end_i   + delta_grna
        } else {
          # Deletion site is upstream of gRNA site → no shift
          adj_del_start <- del_start_i
          adj_del_end   <- del_end_i
        }

        base::cat("    gRNA delta: ", delta_grna, " bp -> deletion site adjusted: ",
                  adj_del_start, "-", adj_del_end, "\n", sep = "")

        # Write intermediate to temp file for deletion primer design
        temp_fasta <- base::tempfile(fileext = ".fasta")
        base::on.exit(
          if (base::file.exists(temp_fasta)) base::file.remove(temp_fasta),
          add = TRUE
        )
        base::writeLines(
          base::c(base::paste0(">", vec_info$name, "_intermediate"),
                   intermediate_seq),
          temp_fasta
        )

        # Design deletion primers against INTERMEDIATE vector
        # Note: returned primer coordinates are in modified-intermediate space,
        # which is identical to final construct space (both insertions applied).
        del_result <- design_deletion_primers(
          vector_file     = temp_fasta,
          insert_seq      = insert_seq,
          start           = adj_del_start,
          end             = adj_del_end,
          upstream_bp     = upstream_bp,
          downstream_bp   = downstream_bp,
          tm_target       = tm_target,
          overlap_length  = base::min(overlap_length, 16L),
          locus_tag       = tag,
          output_dir      = NULL,
          kill_snapgene   = FALSE
        )

        if (base::file.exists(temp_fasta)) base::file.remove(temp_fasta)

        # Store for later GenBank generation
        gene_name    <- ""
        product_name <- ""
        if ("gene"    %in% base::colnames(genbank_table)) gene_name    <- base::as.character(target_row$gene[1])
        if ("product" %in% base::colnames(genbank_table)) product_name <- base::as.character(target_row$product[1])
        if (base::is.na(gene_name))    gene_name    <- ""
        if (base::is.na(product_name)) product_name <- ""

        deletion_data[[tag]] <- base::list(
          primers       = del_result,
          insert_seq    = insert_seq,
          adj_del_start = adj_del_start,
          adj_del_end   = adj_del_end,
          delta_grna    = delta_grna,
          gene          = gene_name,
          product       = product_name
        )

      }, error = function(e) {
        base::warning("Error processing ", tag, ": ", base::conditionMessage(e))
      })
    }

    # Build deletion_primers data frame
    del_rows <- base::list()
    for (tag in base::names(deletion_data)) {
      dd <- deletion_data[[tag]]
      p  <- dd$primers
      del_rows[[tag]] <- base::data.frame(
        locus_tag = tag, gene = dd$gene, product = dd$product,
        upstream_forward_name       = p$upstream_forward_name,
        upstream_forward_primer     = p$upstream_forward_primer,
        upstream_forward_tm_target  = p$upstream_forward_tm_target,
        upstream_forward_tm_full    = p$upstream_forward_tm_full,
        upstream_reverse_name       = p$upstream_reverse_name,
        upstream_reverse_primer     = p$upstream_reverse_primer,
        upstream_reverse_tm_target  = p$upstream_reverse_tm_target,
        upstream_reverse_tm_full    = p$upstream_reverse_tm_full,
        downstream_forward_name     = p$downstream_forward_name,
        downstream_forward_primer   = p$downstream_forward_primer,
        downstream_forward_tm_target = p$downstream_forward_tm_target,
        downstream_forward_tm_full  = p$downstream_forward_tm_full,
        downstream_reverse_name     = p$downstream_reverse_name,
        downstream_reverse_primer   = p$downstream_reverse_primer,
        downstream_reverse_tm_target = p$downstream_reverse_tm_target,
        downstream_reverse_tm_full  = p$downstream_reverse_tm_full,
        overlap_tm = p$overlap_tm, overlap_length = p$overlap_length,
        stringsAsFactors = FALSE
      )
    }
    deletion_primers <- if (base::length(del_rows) > 0) {
      base::do.call(base::rbind, del_rows)
    } else {
      base::data.frame()
    }
    base::rownames(deletion_primers) <- NULL

    # ===== Write COMBINED GenBank files (one per gRNA) =====
    if (!base::is.null(output_dir)) {
      base::cat("\n  Writing combined GenBank files...\n")

      # Pre-cache vector GenBank lines + sequence (avoid SnapGene CLI per gRNA)
      cache_start <- base::Sys.time()
      ext_vec <- base::tolower(tools::file_ext(grna_vector_file))
      vec_norm <- base::normalizePath(grna_vector_file, mustWork = TRUE)
      if (ext_vec == "dna") {
        cached_gb_lines_vec <- .get_genbank_from_dna(vec_norm, kill_snapgene)
      } else if (ext_vec %in% base::c("gb", "gbk", "genbank")) {
        cached_gb_lines_vec <- base::readLines(vec_norm, warn = FALSE)
      } else {
        vec_tmp <- read_vector_file(vec_norm)
        cached_gb_lines_vec <- .build_minimal_genbank(vec_tmp$name, vec_tmp$sequence)
      }
      base::cat("  Vector cached in ",
                base::round(base::difftime(base::Sys.time(), cache_start, units = "secs"), 1),
                " sec\n", sep = "")

      for (i in base::seq_len(base::nrow(grna_primers))) {
        tag_i     <- base::as.character(grna_primers$locus_tag[i])
        spacer_i  <- base::as.character(grna_primers$protospacer[i])
        dd        <- deletion_data[[tag_i]]
        if (base::is.null(dd)) next

        grna_name_i <- ""
        if ("gRNA_name" %in% base::colnames(grna_primers)) {
          grna_name_i <- base::as.character(grna_primers$gRNA_name[i])
        }
        rank_i <- if ("gRNA_rank" %in% base::colnames(grna_primers)) {
          grna_primers$gRNA_rank[i]
        } else i

        # Build gRNA primer_info list
        grna_pinfo <- base::list(
          primer_F      = if ("primer_F" %in% base::colnames(grna_primers))
                            grna_primers$primer_F[i] else NULL,
          primer_R      = if ("primer_R" %in% base::colnames(grna_primers))
                            grna_primers$primer_R[i] else NULL,
          primer_F_Tm   = if ("primer_F_Tm" %in% base::colnames(grna_primers))
                            grna_primers$primer_F_Tm[i] else NA,
          primer_R_Tm   = if ("primer_R_Tm" %in% base::colnames(grna_primers))
                            grna_primers$primer_R_Tm[i] else NA,
          primer_F_name = if ("primer_F_name" %in% base::colnames(grna_primers))
                            grna_primers$primer_F_name[i] else "gRNA_F",
          primer_R_name = if ("primer_R_name" %in% base::colnames(grna_primers))
                            grna_primers$primer_R_name[i] else "gRNA_R",
          cloning_method = grna_cloning_method,
          enzyme         = grna_enzyme,
          gibson_overlap = if ("gibson_overlap" %in% base::colnames(grna_primers))
                             grna_primers$gibson_overlap[i] else NA
        )

        gbk_label <- tag_i
        if (base::nchar(gbk_label) == 0) gbk_label <- grna_name_i
        out_path <- base::file.path(
          output_dir,
          base::paste0(gbk_label, "_g", rank_i, "_combined_construct.gbk")
        )

        # Skip if GenBank already exists (resume mode)
        if (resume && base::file.exists(out_path)) next

        .write_combined_construct_genbank(
          vector_file     = grna_vector_file,
          spacer_seq      = spacer_i,
          grna_start      = grna_start_i,
          grna_end        = grna_end_i,
          insert_seq      = dd$insert_seq,
          deletion_start  = del_start_i,
          deletion_end    = del_end_i,
          locus_tag       = tag_i,
          gene            = dd$gene,
          gRNA_name       = grna_name_i,
          nuclease        = nuclease,
          grna_primer_info = grna_pinfo,
          deletion_primers = dd$primers,
          upstream_bp     = upstream_bp,
          downstream_bp   = downstream_bp,
          output_path     = out_path,
          kill_snapgene   = FALSE,
          cached_gb_lines = cached_gb_lines_vec,
          cached_vec_seq  = original_vec_seq
        )
      }
    }

  } else {
    # ========== SEPARATE MODE ==========
    deletion_primers <- batch_deletion_primers(
      vector_file    = deletion_vector_file,
      genome_seq     = genome_seq,
      genbank_table  = genbank_table,
      locus_tags     = deletion_tags,
      start          = del_start_i,
      end            = del_end_i,
      upstream_bp    = upstream_bp,
      downstream_bp  = downstream_bp,
      tm_target      = tm_target,
      overlap_length = base::min(overlap_length, 16L),
      output_file    = NULL,
      output_dir     = if (!base::is.null(output_dir))
                         base::file.path(output_dir, "deletion_constructs") else NULL,
      kill_snapgene  = kill_snapgene
    )
  }

  base::cat("\n  Deletion primers: ", base::nrow(deletion_primers),
            " construct(s)\n\n", sep = "")

  # ============================================================
  # Phase 4: Combined Excel Output
  # ============================================================
  # Build unified primer table:
  #   locus_tag | type | direction | primer_name | sequence | Tm_target | Tm_full | ...
  unified_rows <- base::list()
  row_idx <- 0L

  for (i in base::seq_len(base::nrow(grna_primers))) {
    tag_i <- base::as.character(grna_primers$locus_tag[i])
    gname <- if ("gRNA_name" %in% base::colnames(grna_primers))
               base::as.character(grna_primers$gRNA_name[i]) else ""
    spacer_i <- if ("protospacer" %in% base::colnames(grna_primers))
                  base::as.character(grna_primers$protospacer[i]) else ""
    score_i  <- if ("composite_score" %in% base::colnames(grna_primers))
                  grna_primers$composite_score[i] else NA

    # gRNA Forward primer
    if ("primer_F" %in% base::colnames(grna_primers)) {
      row_idx <- row_idx + 1L
      unified_rows[[row_idx]] <- base::data.frame(
        locus_tag  = tag_i,
        type       = "sgRNA",
        direction  = "F",
        primer_name = if ("primer_F_name" %in% base::colnames(grna_primers))
                        grna_primers$primer_F_name[i] else base::paste0(tag_i, "_gRNA_F"),
        sequence   = grna_primers$primer_F[i],
        Tm_target  = if ("primer_F_Tm" %in% base::colnames(grna_primers))
                       grna_primers$primer_F_Tm[i] else NA,
        Tm_full    = if ("primer_F_Tm_full" %in% base::colnames(grna_primers))
                       grna_primers$primer_F_Tm_full[i] else NA,
        length     = base::nchar(grna_primers$primer_F[i]),
        gRNA_name  = gname,
        protospacer = spacer_i,
        composite_score = score_i,
        stringsAsFactors = FALSE
      )
    }

    # gRNA Reverse primer
    if ("primer_R" %in% base::colnames(grna_primers)) {
      row_idx <- row_idx + 1L
      unified_rows[[row_idx]] <- base::data.frame(
        locus_tag  = tag_i,
        type       = "sgRNA",
        direction  = "R",
        primer_name = if ("primer_R_name" %in% base::colnames(grna_primers))
                        grna_primers$primer_R_name[i] else base::paste0(tag_i, "_gRNA_R"),
        sequence   = grna_primers$primer_R[i],
        Tm_target  = if ("primer_R_Tm" %in% base::colnames(grna_primers))
                       grna_primers$primer_R_Tm[i] else NA,
        Tm_full    = if ("primer_R_Tm_full" %in% base::colnames(grna_primers))
                       grna_primers$primer_R_Tm_full[i] else NA,
        length     = base::nchar(grna_primers$primer_R[i]),
        gRNA_name  = gname,
        protospacer = spacer_i,
        composite_score = score_i,
        stringsAsFactors = FALSE
      )
    }
  }

  # Deletion arm primers (4 per locus_tag)
  if (base::nrow(deletion_primers) > 0) {
    for (i in base::seq_len(base::nrow(deletion_primers))) {
      dp <- deletion_primers[i, ]
      tag_i <- base::as.character(dp$locus_tag)

      # UF (upstream forward)
      row_idx <- row_idx + 1L
      unified_rows[[row_idx]] <- base::data.frame(
        locus_tag  = tag_i,
        type       = "homologous_arm",
        direction  = "F",
        primer_name = dp$upstream_forward_name,
        sequence   = dp$upstream_forward_primer,
        Tm_target  = dp$upstream_forward_tm_target,
        Tm_full    = dp$upstream_forward_tm_full,
        length     = base::nchar(dp$upstream_forward_primer),
        gRNA_name  = NA,
        protospacer = NA,
        composite_score = NA,
        stringsAsFactors = FALSE
      )
      # UR (upstream reverse)
      row_idx <- row_idx + 1L
      unified_rows[[row_idx]] <- base::data.frame(
        locus_tag  = tag_i,
        type       = "homologous_arm",
        direction  = "R",
        primer_name = dp$upstream_reverse_name,
        sequence   = dp$upstream_reverse_primer,
        Tm_target  = dp$upstream_reverse_tm_target,
        Tm_full    = dp$upstream_reverse_tm_full,
        length     = base::nchar(dp$upstream_reverse_primer),
        gRNA_name  = NA,
        protospacer = NA,
        composite_score = NA,
        stringsAsFactors = FALSE
      )
      # DF (downstream forward)
      row_idx <- row_idx + 1L
      unified_rows[[row_idx]] <- base::data.frame(
        locus_tag  = tag_i,
        type       = "homologous_arm",
        direction  = "F",
        primer_name = dp$downstream_forward_name,
        sequence   = dp$downstream_forward_primer,
        Tm_target  = dp$downstream_forward_tm_target,
        Tm_full    = dp$downstream_forward_tm_full,
        length     = base::nchar(dp$downstream_forward_primer),
        gRNA_name  = NA,
        protospacer = NA,
        composite_score = NA,
        stringsAsFactors = FALSE
      )
      # DR (downstream reverse)
      row_idx <- row_idx + 1L
      unified_rows[[row_idx]] <- base::data.frame(
        locus_tag  = tag_i,
        type       = "homologous_arm",
        direction  = "R",
        primer_name = dp$downstream_reverse_name,
        sequence   = dp$downstream_reverse_primer,
        Tm_target  = dp$downstream_reverse_tm_target,
        Tm_full    = dp$downstream_reverse_tm_full,
        length     = base::nchar(dp$downstream_reverse_primer),
        gRNA_name  = NA,
        protospacer = NA,
        composite_score = NA,
        stringsAsFactors = FALSE
      )
    }
  }

  unified_primer_table <- if (base::length(unified_rows) > 0) {
    base::do.call(base::rbind, unified_rows)
  } else {
    base::data.frame()
  }
  base::rownames(unified_primer_table) <- NULL

  if (!base::is.null(output_file)) {
    base::cat("=== Phase 4: Saving Excel ===\n\n")
    base::tryCatch({
      excel_dir <- base::dirname(output_file)
      if (excel_dir != "." && !base::dir.exists(excel_dir)) {
        base::dir.create(excel_dir, recursive = TRUE)
      }
      if (base::requireNamespace("openxlsx", quietly = TRUE)) {
        wb <- openxlsx::createWorkbook()

        # Sheet 1: Unified primer table (main deliverable)
        openxlsx::addWorksheet(wb, "all_primers")
        openxlsx::writeData(wb, "all_primers", unified_primer_table)

        # Sheet 2: gRNA library (raw data with scores)
        openxlsx::addWorksheet(wb, "gRNA_library")
        openxlsx::writeData(wb, "gRNA_library", grna_primers)

        # Sheet 3: deletion primer detail (wide format)
        openxlsx::addWorksheet(wb, "deletion_detail")
        openxlsx::writeData(wb, "deletion_detail", deletion_primers)

        # Styling
        hs <- openxlsx::createStyle(textDecoration = "bold",
                                     fgFill = "#DCE6F1", border = "Bottom")
        for (sn in base::c("all_primers", "gRNA_library", "deletion_detail")) {
          nc <- base::tryCatch(
            base::ncol(openxlsx::readWorkbook(wb, sheet = sn)),
            error = function(e) 10L)
          # Use conditional ncol for each sheet
          sn_data <- base::switch(sn,
            "all_primers"     = unified_primer_table,
            "gRNA_library"    = grna_primers,
            "deletion_detail" = deletion_primers)
          if (base::nrow(sn_data) > 0) {
            openxlsx::addStyle(wb, sn, hs, rows = 1,
                                cols = base::seq_len(base::ncol(sn_data)))
            openxlsx::freezePane(wb, sn, firstRow = TRUE)
            base::tryCatch(
              openxlsx::setColWidths(wb, sn,
                                      cols = base::seq_len(base::ncol(sn_data)),
                                      widths = "auto"),
              error = function(e) NULL)
          }
        }

        # Conditional formatting: color-code type column in all_primers
        if (base::nrow(unified_primer_table) > 0) {
          sgRNA_style <- openxlsx::createStyle(fgFill = "#E2EFDA")  # light green
          arm_style   <- openxlsx::createStyle(fgFill = "#FCE4D6")  # light orange
          type_col <- base::which(base::colnames(unified_primer_table) == "type")
          sgRNA_rows <- base::which(unified_primer_table$type == "sgRNA") + 1L
          arm_rows   <- base::which(unified_primer_table$type != "sgRNA") + 1L
          if (base::length(sgRNA_rows) > 0) {
            openxlsx::addStyle(wb, "all_primers", sgRNA_style,
                                rows = sgRNA_rows, cols = type_col, stack = TRUE)
          }
          if (base::length(arm_rows) > 0) {
            openxlsx::addStyle(wb, "all_primers", arm_style,
                                rows = arm_rows, cols = type_col, stack = TRUE)
          }
        }

        openxlsx::saveWorkbook(wb, file = output_file, overwrite = TRUE)
        base::cat("  Saved: ", output_file, "\n", sep = "")
      } else if (base::requireNamespace("writexl", quietly = TRUE)) {
        writexl::write_xlsx(
          base::list(all_primers = unified_primer_table,
                     gRNA_library = grna_primers,
                     deletion_detail = deletion_primers),
          path = output_file
        )
        base::cat("  Saved: ", output_file, "\n", sep = "")
      }
    }, error = function(e) {
      base::warning("Excel save failed: ", base::conditionMessage(e))
    })
  }

  # ============================================================
  # Phase 5: Summary
  # ============================================================
  base::cat("\n============================================================\n")
  base::cat(" Design Complete\n")
  base::cat("============================================================\n")
  base::cat("  gRNA constructs:     ", base::nrow(grna_primers), "\n", sep = "")
  base::cat("  Deletion constructs: ", base::nrow(deletion_primers), "\n", sep = "")
  if (combined_mode) {
    base::cat("  Mode: COMBINED (single GenBank per gRNA with both spacer + arms)\n")
  }
  if (!base::is.null(output_file)) base::cat("  Excel: ", output_file, "\n", sep = "")
  if (!base::is.null(output_dir))  base::cat("  GenBank: ", output_dir, "/\n", sep = "")

  if ("gRNA_name" %in% base::colnames(grna_primers)) {
    base::cat("\n  Selected gRNAs:\n")
    for (i in base::seq_len(base::nrow(grna_primers))) {
      tag_i <- base::as.character(grna_primers$locus_tag[i])
      nm_i  <- base::as.character(grna_primers$gRNA_name[i])
      sp_i  <- base::as.character(grna_primers$protospacer[i])
      sc_i  <- if ("composite_score" %in% base::colnames(grna_primers))
                 base::paste0(", score=", base::round(grna_primers$composite_score[i], 1)) else ""
      base::cat("    ", tag_i, " -> ", nm_i, ": ", sp_i,
                " (", base::nchar(sp_i), "bp", sc_i, ")\n", sep = "")
    }
  }
  base::cat("\n")

  base::invisible(base::list(
    grna_primers       = grna_primers,
    deletion_primers   = deletion_primers,
    unified_primers    = unified_primer_table
  ))
}


# ==============================================================================
# Internal: Write combined construct GenBank (gRNA spacer + deletion arms)
# ==============================================================================

#' @keywords internal
.write_combined_construct_genbank <- function(
    vector_file,
    spacer_seq,
    grna_start,
    grna_end,
    insert_seq,
    deletion_start,
    deletion_end,
    locus_tag,
    gene,
    gRNA_name,
    nuclease,
    grna_primer_info,
    deletion_primers,
    upstream_bp,
    downstream_bp,
    output_path,
    kill_snapgene = FALSE,
    cached_gb_lines = NULL,
    cached_vec_seq = NULL
) {
  # --- Step 1: Get GenBank text and vector sequence ---
  ext <- base::tolower(tools::file_ext(vector_file))
  vector_file <- base::normalizePath(vector_file, mustWork = TRUE)

  if (!base::is.null(cached_gb_lines)) {
    gb_lines <- cached_gb_lines
  } else if (ext == "dna") {
    gb_lines <- .get_genbank_from_dna(vector_file, kill_snapgene)
  } else if (ext %in% base::c("gb", "gbk", "genbank")) {
    gb_lines <- base::readLines(vector_file, warn = FALSE)
  } else {
    vec <- read_vector_file(vector_file)
    gb_lines <- .build_minimal_genbank(vec$name, vec$sequence)
  }

  if (!base::is.null(cached_vec_seq)) {
    original_seq <- base::toupper(cached_vec_seq)
  } else {
    vec_info <- read_vector_file(vector_file, kill_snapgene = FALSE)
    original_seq <- base::toupper(base::as.character(vec_info$sequence))
  }
  spacer_seq   <- base::toupper(spacer_seq)
  insert_seq   <- base::toupper(insert_seq)

  grna_start_i <- base::as.integer(grna_start)
  grna_end_i   <- base::as.integer(grna_end)
  del_start_i  <- base::as.integer(deletion_start)
  del_end_i    <- base::as.integer(deletion_end)

  grna_stuffer_len <- grna_end_i - grna_start_i + 1
  del_stuffer_len  <- del_end_i  - del_start_i  + 1
  spacer_len       <- base::nchar(spacer_seq)
  insert_len       <- base::nchar(insert_seq)

  delta_grna <- spacer_len - grna_stuffer_len
  delta_del  <- insert_len - del_stuffer_len

  # --- Step 2: Build final sequence with both modifications ---
  # Apply from right-to-left (downstream first) to preserve upstream coordinates
  if (grna_start_i < del_start_i) {
    # gRNA site is upstream: apply deletion first (at original coords)
    step1_seq <- base::paste0(
      base::substring(original_seq, 1, del_start_i - 1),
      insert_seq,
      base::substring(original_seq, del_end_i + 1, base::nchar(original_seq))
    )
    # Then apply gRNA (original coords still valid since it's upstream)
    final_seq <- base::paste0(
      base::substring(step1_seq, 1, grna_start_i - 1),
      spacer_seq,
      base::substring(step1_seq, grna_end_i + 1, base::nchar(step1_seq))
    )
  } else {
    # Deletion site is upstream: apply gRNA first (at original coords)
    step1_seq <- base::paste0(
      base::substring(original_seq, 1, grna_start_i - 1),
      spacer_seq,
      base::substring(original_seq, grna_end_i + 1, base::nchar(original_seq))
    )
    # Then apply deletion (original coords still valid since it's upstream)
    final_seq <- base::paste0(
      base::substring(step1_seq, 1, del_start_i - 1),
      insert_seq,
      base::substring(step1_seq, del_end_i + 1, base::nchar(step1_seq))
    )
  }
  new_length <- base::nchar(final_seq)

  # --- Step 3: Parse GenBank header and features ---
  header_lines  <- base::character(0)
  feature_lines <- base::character(0)
  in_features   <- FALSE
  in_origin     <- FALSE

  for (line in gb_lines) {
    if (base::startsWith(line, "FEATURES")) {
      in_features <- TRUE
      header_lines <- base::c(header_lines, line)
      next
    } else if (base::startsWith(line, "ORIGIN")) {
      in_features <- FALSE
      in_origin   <- TRUE
      next
    } else if (base::startsWith(line, "//")) {
      next
    }
    if (in_features)  feature_lines <- base::c(feature_lines, line)
    else if (!in_origin) header_lines <- base::c(header_lines, line)
  }

  # --- Step 4: Update LOCUS line ---
  construct_label <- base::ifelse(base::nchar(locus_tag) > 0, locus_tag,
                     base::ifelse(base::nchar(gRNA_name) > 0, gRNA_name, "construct"))

  for (i in base::seq_along(header_lines)) {
    line <- header_lines[i]
    if (base::startsWith(line, "LOCUS")) {
      parts <- base::strsplit(base::trimws(line), "\\s+")[[1]]
      parts[2] <- base::paste0(parts[2], "_", construct_label)
      parts[3] <- base::paste0(new_length, " bp")
      header_lines[i] <- base::paste(parts, collapse = "  ")
    } else if (base::grepl("^\\s*REFERENCE.*\\(bases", line)) {
      header_lines[i] <- base::sub("\\(bases.*\\)",
                                     base::paste0("(bases 1 to ", new_length, ")"),
                                     line)
    }
  }

  # --- Step 5: Adjust existing features for BOTH insertions ---
  # We need to apply two rounds of adjustment.
  # Determine which is upstream/downstream to apply adjustments in order.
  if (grna_start_i < del_start_i) {
    # gRNA upstream, deletion downstream
    first_start <- grna_start_i; first_end <- grna_end_i; first_diff <- delta_grna
    second_start_orig <- del_start_i; second_end_orig <- del_end_i; second_diff <- delta_del
    # After first adjustment, second insertion site shifts by first_diff
    second_start_adj <- second_start_orig + first_diff
    second_end_adj   <- second_end_orig   + first_diff
  } else {
    # deletion upstream, gRNA downstream
    first_start <- del_start_i; first_end <- del_end_i; first_diff <- delta_del
    second_start_orig <- grna_start_i; second_end_orig <- grna_end_i; second_diff <- delta_grna
    second_start_adj <- second_start_orig + first_diff
    second_end_adj   <- second_end_orig   + first_diff
  }

  # Parse feature blocks
  adjusted_features <- base::character(0)
  current_block <- base::character(0)

  .flush_block <- function(block) {
    if (base::length(block) == 0) return(NULL)
    # First adjustment (upstream insertion)
    adj1 <- .adjust_feature_block(block, first_start, first_end, first_diff)
    if (base::is.null(adj1)) return(NULL)
    # Second adjustment (downstream insertion, at adjusted coordinates)
    adj2 <- .adjust_feature_block(adj1, second_start_adj, second_end_adj, second_diff)
    return(adj2)
  }

  for (line in feature_lines) {
    trimmed <- base::trimws(line)
    is_new_feature <- (base::nchar(trimmed) > 0 &&
                        !base::startsWith(trimmed, "/") &&
                        (base::grepl("\\.\\.|\\.\\^\\.", line) ||
                         (base::grepl("^     \\S", line) && !base::grepl("^     /", line))))

    if (is_new_feature && base::length(current_block) > 0) {
      result <- .flush_block(current_block)
      if (!base::is.null(result)) adjusted_features <- base::c(adjusted_features, result)
      current_block <- line
    } else {
      current_block <- base::c(current_block, line)
    }
  }
  if (base::length(current_block) > 0) {
    result <- .flush_block(current_block)
    if (!base::is.null(result)) adjusted_features <- base::c(adjusted_features, result)
  }

  # --- Step 6: Compute final positions for new features ---
  # In the final construct, both insertions have been applied.
  # Spacer position in final construct:
  if (grna_start_i < del_start_i) {
    # gRNA is upstream: spacer at original grna_start_i
    spacer_final_start <- grna_start_i
    spacer_final_end   <- grna_start_i + spacer_len - 1
    # Arms position: after gRNA delta adjustment
    arms_final_start   <- del_start_i + delta_grna
    arms_final_end     <- del_start_i + delta_grna + insert_len - 1
  } else {
    # Deletion is upstream: arms at original del_start_i
    arms_final_start   <- del_start_i
    arms_final_end     <- del_start_i + insert_len - 1
    # Spacer position: after deletion delta adjustment
    spacer_final_start <- grna_start_i + delta_del
    spacer_final_end   <- grna_start_i + delta_del + spacer_len - 1
  }

  up_arm_end   <- arms_final_start + upstream_bp - 1
  dn_arm_start <- up_arm_end + 1

  # --- Step 7: gRNA spacer annotation ---
  spacer_label <- base::ifelse(base::nchar(locus_tag) > 0, locus_tag,
                  base::ifelse(base::nchar(gRNA_name) > 0, gRNA_name, "gRNA_spacer"))
  new_features <- base::c(
    "",
    base::sprintf("     misc_feature    %d..%d", spacer_final_start, spacer_final_end),
    base::sprintf("                     /label=%s", spacer_label)
  )
  if (base::nchar(locus_tag) > 0) {
    new_features <- base::c(new_features,
      base::sprintf("                     /locus_tag=%s", locus_tag))
  }
  if (base::nchar(gene) > 0) {
    new_features <- base::c(new_features,
      base::sprintf("                     /gene=%s", gene))
  }
  if (base::nchar(nuclease) > 0) {
    new_features <- base::c(new_features,
      base::sprintf("                     /note=nuclease: %s", nuclease))
  }
  new_features <- base::c(new_features,
    base::sprintf("                     /note=spacer: %s (%d bp)", spacer_seq, spacer_len))
  if (base::nchar(gRNA_name) > 0) {
    new_features <- base::c(new_features,
      base::sprintf("                     /note=gRNA_name: %s", gRNA_name))
  }

  # --- Step 8: Upstream/downstream arm annotations ---
  new_features <- base::c(new_features,
    "",
    base::sprintf("     misc_feature    %d..%d", arms_final_start, up_arm_end),
    base::sprintf("                     /label=%s_upstream_arm", locus_tag),
    base::sprintf("                     /note=upstream homology arm (%d bp)", upstream_bp),
    "",
    base::sprintf("     misc_feature    %d..%d", dn_arm_start, arms_final_end),
    base::sprintf("                     /label=%s_downstream_arm", locus_tag),
    base::sprintf("                     /note=downstream homology arm (%d bp)", downstream_bp)
  )

  # --- Step 9: gRNA primer annotations ---
  if (!base::is.null(grna_primer_info$primer_F)) {
    is_gibson <- base::grepl("gibson", grna_primer_info$cloning_method, ignore.case = TRUE)
    if (is_gibson) {
      f_head_len <- base::nchar(grna_primer_info$primer_F) - spacer_len
      r_overlap  <- spacer_len
      if (!base::is.na(grna_primer_info$gibson_overlap)) {
        r_overlap <- base::as.integer(grna_primer_info$gibson_overlap)
      }
      r_head_len <- base::nchar(grna_primer_info$primer_R) - r_overlap

      f_start <- spacer_final_start
      f_end   <- spacer_final_end + base::max(0L, f_head_len)
      r_start <- spacer_final_start - base::max(0L, r_head_len)
      r_end   <- spacer_final_start + r_overlap - 1L
    } else {
      f_start <- spacer_final_start; f_end <- spacer_final_end
      r_start <- spacer_final_start; r_end <- spacer_final_end
    }

    f_label <- base::ifelse(base::nchar(grna_primer_info$primer_F_name) > 0,
                              grna_primer_info$primer_F_name, "gRNA_F")
    r_label <- base::ifelse(base::nchar(grna_primer_info$primer_R_name) > 0,
                              grna_primer_info$primer_R_name, "gRNA_R")

    new_features <- base::c(new_features,
      "",
      base::sprintf("     primer          %d..%d", f_start, f_end),
      base::sprintf("                     /label=%s", f_label),
      base::sprintf("                     /note=seq: %s (%d bp)",
                     grna_primer_info$primer_F, base::nchar(grna_primer_info$primer_F)),
      "",
      base::sprintf("     primer          complement(%d..%d)", r_start, r_end),
      base::sprintf("                     /label=%s", r_label),
      base::sprintf("                     /note=seq: %s (%d bp)",
                     grna_primer_info$primer_R, base::nchar(grna_primer_info$primer_R))
    )
  }

  # --- Step 10: Deletion arm primer annotations ---
  # Primer coordinates from design_deletion_primers() are in modified-intermediate
  # vector space, which equals the final construct coordinate space.
  dp <- deletion_primers
  if (!base::is.null(dp$upstream_forward_primer)) {
    new_features <- base::c(new_features,
      "",
      base::sprintf("     primer          %d..%d",
                     dp$upstream_forward_start, dp$upstream_forward_end),
      base::sprintf("                     /label=%s", dp$upstream_forward_name),
      base::sprintf("                     /note=Tm(target): %.2f C, Tm(full): %.2f C",
                     dp$upstream_forward_tm_target, dp$upstream_forward_tm_full),
      base::sprintf("                     /note=seq: %s", dp$upstream_forward_primer),
      "",
      base::sprintf("     primer          complement(%d..%d)",
                     dp$upstream_reverse_start, dp$upstream_reverse_end),
      base::sprintf("                     /label=%s", dp$upstream_reverse_name),
      base::sprintf("                     /note=Tm(target): %.2f C, Tm(full): %.2f C",
                     dp$upstream_reverse_tm_target, dp$upstream_reverse_tm_full),
      base::sprintf("                     /note=seq: %s", dp$upstream_reverse_primer),
      "",
      base::sprintf("     primer          %d..%d",
                     dp$downstream_forward_start, dp$downstream_forward_end),
      base::sprintf("                     /label=%s", dp$downstream_forward_name),
      base::sprintf("                     /note=Tm(target): %.2f C, Tm(full): %.2f C",
                     dp$downstream_forward_tm_target, dp$downstream_forward_tm_full),
      base::sprintf("                     /note=seq: %s", dp$downstream_forward_primer),
      "",
      base::sprintf("     primer          complement(%d..%d)",
                     dp$downstream_reverse_start, dp$downstream_reverse_end),
      base::sprintf("                     /label=%s", dp$downstream_reverse_name),
      base::sprintf("                     /note=Tm(target): %.2f C, Tm(full): %.2f C",
                     dp$downstream_reverse_tm_target, dp$downstream_reverse_tm_full),
      base::sprintf("                     /note=seq: %s", dp$downstream_reverse_primer)
    )
  }

  # --- Step 11: Build ORIGIN ---
  origin_lines <- "ORIGIN"
  seq_lower <- base::tolower(final_seq)
  for (i in base::seq(1, base::nchar(seq_lower), by = 60)) {
    chunk <- base::substring(seq_lower, i,
                              base::min(i + 59, base::nchar(seq_lower)))
    blocks <- base::character(0)
    for (j in base::seq(1, base::nchar(chunk), by = 10)) {
      blocks <- base::c(blocks,
                          base::substring(chunk, j,
                                           base::min(j + 9, base::nchar(chunk))))
    }
    formatted <- base::paste(blocks, collapse = " ")
    pos_str <- base::as.character(i)
    padding <- base::paste(base::rep(" ", 9 - base::nchar(pos_str)), collapse = "")
    origin_lines <- base::c(origin_lines, base::paste0(padding, pos_str, " ", formatted))
  }
  origin_lines <- base::c(origin_lines, "//")

  # --- Step 12: Assemble and write ---
  all_lines <- base::c(header_lines, adjusted_features, new_features, origin_lines)

  out_dir <- base::dirname(output_path)
  if (!base::dir.exists(out_dir) && out_dir != ".") {
    base::dir.create(out_dir, recursive = TRUE)
  }
  base::writeLines(all_lines, output_path)
  base::cat("  Combined GenBank: ", output_path, " (",
            new_length, " bp)\n", sep = "")

  base::invisible(output_path)
}


# ==============================================================================
# Internal: Fast genome sequence extractor (skips feature parsing)
# ==============================================================================

#' Extract only the genome sequence from a GenBank file (ORIGIN section).
#' Much faster than full read_genome_genbank() when genbank_table is cached.
#' @keywords internal
.extract_genome_seq_only <- function(gbff_path) {
  lines <- base::readLines(gbff_path, warn = FALSE)
  in_origin <- FALSE
  seq_parts <- base::character(0)
  for (line in lines) {
    if (base::startsWith(line, "ORIGIN")) {
      in_origin <- TRUE
      next
    }
    if (in_origin) {
      if (base::startsWith(line, "//")) break
      cleaned <- base::gsub("[^a-zA-Z]", "", line)
      if (base::nchar(cleaned) > 0) {
        seq_parts <- base::c(seq_parts, cleaned)
      }
    }
  }
  genome_seq <- base::paste0(seq_parts, collapse = "")
  base::toupper(genome_seq)
}
