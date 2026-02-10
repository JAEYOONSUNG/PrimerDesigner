#' Design gRNA Expression Construct (End-to-End Pipeline)
#'
#' All-in-one pipeline for CRISPR gRNA construct design. Supports three modes:
#'
#' \strong{Mode 1a — From GenBank file (recommended):} Provide \code{genbank_file}
#' (path to a local .gbff file) and \code{locus_tags}. The accession is auto-extracted
#' from the filename. BSgenome/Bowtie are built in \code{dir}.
#'
#' \strong{Mode 1b — From accession:} Provide \code{genbank_accession}
#' and \code{locus_tags}. Files are downloaded if not present.
#'
#' \strong{Mode 2 — Pre-generated library:} Provide a pre-generated
#' \code{gRNA_df} from \code{run_gRNA_list_generator()}.
#'
#' For multi-gene deletion targets (consecutive locus_tags, e.g.,
#' \code{c("RS00005", "RS00010")}), gRNAs from all specified genes are pooled
#' and the best \code{top_n} gRNA(s) are selected globally. This is the
#' standard workflow for region deletions: the gRNA cuts somewhere within
#' the deletion target to provide selection pressure.
#'
#' @param gRNA_df Optional. A pre-generated gRNA library data frame
#'   (from \code{run_gRNA_list_generator}). If NULL (default), gRNAs are
#'   auto-generated using \code{genbank_file} or \code{genbank_accession}.
#' @param genbank_file Character. Path to a local GenBank file (.gbff, .gb, .gbk).
#'   The accession is extracted from the filename (e.g., \code{GCF_030376765.1.gbff}
#'   → \code{"GCF_030376765.1"}). The file is copied to \code{dir/seqs_srcdir/}
#'   if not already there. Preferred over \code{genbank_accession}.
#' @param genbank_accession Character. GenBank accession number for auto-generation
#'   (e.g., "GCF_030376765.1"). Files are downloaded if not present.
#'   Used when \code{genbank_file} is not provided.
#' @param vector_file Path to the gRNA expression vector (.dna, .gb, .fasta).
#' @param start Integer. Stuffer/insertion site start position in the vector
#'   (1-based, inclusive).
#' @param end Integer. Stuffer/insertion site end position in the vector
#'   (1-based, inclusive).
#' @param locus_tags Character vector of target locus tags.
#'   \itemize{
#'     \item \code{"all"} — genome-wide, all genes, top_n per gene.
#'     \item Single gene: \code{"QT234_RS00005"} — best gRNA within that gene.
#'     \item Multi-gene region: \code{c("QT234_RS00005", "QT234_RS00010")} —
#'       all gRNAs pooled and ranked globally (for region deletion).
#'   }
#' @param nuclease Character. Nuclease name (default: "GeoCas9"). Also used for
#'   methylation PAM-side determination.
#' @param top_n Integer. Number of best gRNAs to select (default: 1).
#'   For multi-gene regions, top_n is global (not per gene).
#' @param cloning_method Character: "gibson" (IA, default) or "golden_gate" (OA).
#' @param enzyme Character. Golden Gate enzyme: "BbsI" (default), "BsaI", "BsmBI".
#' @param custom_overhangs Named list for custom overhangs (optional).
#' @param tm_target Numeric. Target Tm for Gibson primers (default: 60).
#' @param overlap_length Integer. Gibson overlap length (default: 20).
#' @param add_g Logical. Add 5' G for transcription (default: TRUE).
#' @param name_pattern Character. Primer naming template.
#'   Default: \code{"sgRNA_{locus_tag}_g{rank}_{method}_{dir}"}.
#' @param output_file Character. Path to save primer Excel output (default: NULL).
#' @param output_dir Character. Directory to save GenBank files (default: NULL).
#' @param kill_snapgene Logical. Kill SnapGene before reading .dna (default: FALSE).
#' @param methylation_patterns Character vector of methylation recognition
#'   sequences using IUPAC codes (default: NULL). Applied before gRNA selection.
#' @param dir Character. Working directory for auto-generation mode (default: getwd()).
#'   BSgenome, bowtie index, and intermediate library are saved here.
#' @param position_range Numeric vector of length 2. gRNA position range as
#'   percentage within target gene (default: c(5, 95)).
#' @param tm_range Numeric vector of length 2. gRNA Tm range (default: c(50, 70)).
#' @param strand Character. Strand filter ("both", "5", "3"; default: "both").
#' @return A data frame with the selected gRNAs and their cloning primers.
#' @examples
#' \dontrun{
#' # ====== From local GenBank file (recommended) ======
#' result <- design_grna_construct(
#'   genbank_file = "~/genomes/GCF_030376765.1.gbff",
#'   locus_tags  = "QT235_RS00005",
#'   nuclease    = "GeoCas9",
#'   vector_file = "~/vectors/pG1Kt-GeoCas9EF-OA-sfGFP-ACrec.dna",
#'   start = 8811, end = 8840,
#'   methylation_patterns = c("GCCAT", "CCANNNNNTTG"),
#'   dir = "~/Desktop/gRNA_project",
#'   output_file = "primers.xlsx",
#'   output_dir  = "constructs/"
#' )
#'
#' # ====== Genome-wide from file ======
#' result <- design_grna_construct(
#'   genbank_file = "~/genomes/GCF_030376765.1.gbff",
#'   locus_tags  = "all",
#'   vector_file = "~/vectors/pG1Kt-GeoCas9EF-OA-sfGFP-ACrec.dna",
#'   start = 8811, end = 8840,
#'   top_n = 3,
#'   dir = "~/Desktop/gRNA_project",
#'   output_file = "genome_wide_primers.xlsx"
#' )
#'
#' # ====== From accession (downloads if needed) ======
#' result <- design_grna_construct(
#'   genbank_accession = "GCF_030376765.1",
#'   locus_tags  = c("QT235_RS00005", "QT235_RS00010"),
#'   nuclease    = "GeoCas9",
#'   vector_file = "~/vectors/pG1Kt-GeoCas9EF-OA-sfGFP-ACrec.dna",
#'   start = 8811, end = 8840,
#'   top_n = 2,
#'   dir = "~/Desktop/gRNA_project",
#'   output_file = "region_primers.xlsx",
#'   output_dir  = "region_constructs/"
#' )
#'
#' # ====== With pre-generated library (backward compatible) ======
#' result <- design_grna_construct(
#'   gRNA_df     = my_grna_lib,
#'   vector_file = "my_vector.dna",
#'   start = 8811, end = 8840,
#'   locus_tags  = "QT235_RS00005",
#'   output_file = "primers.xlsx"
#' )
#' }
#' @export
design_grna_construct <- function(
    gRNA_df = NULL,
    genbank_file = NULL,
    genbank_accession = NULL,
    vector_file,
    start,
    end,
    locus_tags = NULL,
    nuclease = "GeoCas9",
    top_n = 1,
    cloning_method = "gibson",
    enzyme = "BbsI",
    custom_overhangs = NULL,
    tm_target = 60,
    overlap_length = 20,
    add_g = TRUE,
    name_pattern = "sgRNA_{locus_tag}_{method}_{dir}",
    output_file = NULL,
    output_dir = NULL,
    kill_snapgene = TRUE,
    methylation_patterns = NULL,
    dir = base::getwd(),
    position_range = c(5, 95),
    tm_range = c(50, 70),
    strand = "both"
) {
  # ============================================================
  # Phase 0: Determine mode & validate common inputs
  # ============================================================
  auto_mode <- base::is.null(gRNA_df)

  # --- Smart detection: genbank_accession may be a file path ---
  # If user accidentally passes a file path as genbank_accession, auto-correct
  if (!base::is.null(genbank_accession) && base::is.null(genbank_file)) {
    if (base::grepl("\\.(gbff|gb|gbk)$", genbank_accession, ignore.case = TRUE) ||
        (base::file.exists(genbank_accession) &&
         !base::grepl("^(GCF|GCA)_", genbank_accession))) {
      base::cat("  Auto-detected: genbank_accession is a file path, treating as genbank_file\n")
      genbank_file <- genbank_accession
      genbank_accession <- NULL
    }
  }

  # --- Resolve genbank_file → genbank_accession ---
  # If user provides a local .gbff file, extract accession from filename
  # and ensure files are in the expected dir/seqs_srcdir/ location.

  if (!base::is.null(genbank_file) && base::is.null(genbank_accession)) {
    if (!base::file.exists(genbank_file)) {
      base::stop("genbank_file not found: ", genbank_file)
    }

    # Extract accession from filename
    # e.g., "GCF_030376765.1.gbff" → "GCF_030376765.1"
    gb_basename <- base::basename(genbank_file)
    genbank_accession <- base::sub("\\.(gbff|gb|gbk)$", "", gb_basename,
                                    ignore.case = TRUE)
    base::cat("  GenBank file: ", genbank_file, "\n", sep = "")
    base::cat("  Extracted accession: ", genbank_accession, "\n", sep = "")

    # Ensure .gbff is in dir/seqs_srcdir/
    seqs_dir <- base::file.path(dir, "seqs_srcdir")
    expected_gbff <- base::file.path(seqs_dir,
                                      base::paste0(genbank_accession, ".gbff"))
    if (!base::file.exists(expected_gbff)) {
      base::dir.create(seqs_dir, recursive = TRUE, showWarnings = FALSE)
      base::file.copy(genbank_file, expected_gbff, overwrite = FALSE)
      base::cat("  Copied .gbff → ", expected_gbff, "\n", sep = "")
    }

    # Also check for .fna file alongside the .gbff
    gb_dir <- base::dirname(genbank_file)
    fna_candidates <- base::list.files(
      gb_dir,
      pattern = base::paste0("^", base::gsub("\\.", "\\\\.", genbank_accession),
                              "\\.fna$"),
      full.names = TRUE
    )
    if (base::length(fna_candidates) > 0) {
      expected_fna <- base::file.path(seqs_dir,
                                       base::paste0(genbank_accession, ".fna"))
      if (!base::file.exists(expected_fna)) {
        base::file.copy(fna_candidates[1], expected_fna, overwrite = FALSE)
        base::cat("  Copied .fna  → ", expected_fna, "\n", sep = "")
      }
    }
  }

  if (auto_mode && base::is.null(genbank_accession)) {
    base::stop("Either gRNA_df, genbank_file, or genbank_accession must be provided.")
  }
  if (auto_mode && base::is.null(locus_tags)) {
    base::stop("locus_tags is required when using auto-generation mode.")
  }
  if (!cloning_method %in% base::c("golden_gate", "gibson")) {
    base::stop("cloning_method must be 'golden_gate' or 'gibson'. ",
               "For deletion mode, use design_cloning_primers() directly.")
  }
  if (base::missing(vector_file) || !base::is.character(vector_file)) {
    base::stop("vector_file is required.")
  }
  if (base::missing(start) || base::missing(end)) {
    base::stop("start and end positions (stuffer region) are required.")
  }

  start_i <- base::as.integer(start)
  end_i <- base::as.integer(end)
  stuffer_len <- end_i - start_i + 1
  method_label <- base::ifelse(cloning_method == "golden_gate",
                                "Golden Gate (OA)", "Gibson (IA)")

  # Track whether methylation was already applied
  methylation_applied <- FALSE
  lib_file <- NULL

  # ============================================================
  # Phase 1: Get gRNA candidates
  # ============================================================
  if (auto_mode) {
    # ----- Auto-generation mode -----
    base::cat("============================================================\n")
    base::cat(" gRNA Construct Design — Auto-Generation Mode\n")
    base::cat("============================================================\n")
    base::cat("  Accession:  ", genbank_accession, "\n")

    if (base::length(locus_tags) > 1) {
      base::cat("  Target:      ", base::paste(locus_tags, collapse = " → "),
                " (multi-gene region)\n")
    } else {
      base::cat("  Target:      ", locus_tags, "\n")
    }
    base::cat("  Nuclease:    ", nuclease, "\n")
    base::cat("  Stuffer:     ", start_i, "-", end_i, " (", stuffer_len, " bp)\n")
    base::cat("  Method:      ", method_label, "\n")
    if (!base::is.null(methylation_patterns)) {
      base::cat("  Methylation: ", base::paste(methylation_patterns, collapse = ", "), "\n")
    }
    base::cat("\n")

    # Generate gRNA library
    # (includes methylation filtering if patterns provided)
    lib_tag_str <- base::substr(
      base::paste(locus_tags, collapse = "_"), 1, 80
    )
    lib_file <- base::paste0("gRNA_lib_", lib_tag_str, ".xlsx")

    gRNA_df <- run_gRNA_list_generator(
      genbank_accession    = genbank_accession,
      target_sequence      = locus_tags,
      nuclease             = nuclease,
      dir                  = dir,
      position_range       = position_range,
      tm_range             = tm_range,
      strand               = strand,
      methylation_patterns = methylation_patterns,
      output_file          = lib_file
    )

    if (!base::is.data.frame(gRNA_df) || base::nrow(gRNA_df) == 0) {
      base::stop("No gRNAs found for the specified locus_tags. ",
                 "Check genbank_accession and locus_tag names.")
    }

    base::cat("\ngRNA library: ", base::nrow(gRNA_df), " candidates generated\n", sep = "")
    methylation_applied <- TRUE

  } else {
    # ----- Pre-generated library mode -----
    if (!base::is.data.frame(gRNA_df) || base::nrow(gRNA_df) == 0) {
      base::stop("gRNA_df must be a non-empty data frame.")
    }
    if (!"protospacer" %in% base::colnames(gRNA_df)) {
      base::stop("gRNA_df must contain a 'protospacer' column.")
    }

    base::cat("=== gRNA Expression Construct Design ===\n")
    base::cat("  Vector stuffer region: ", start_i, " - ", end_i,
              " (", stuffer_len, " bp)\n", sep = "")
    base::cat("  Cloning method: ", method_label, "\n", sep = "")
  }

  # ============================================================
  # Phase 2: Filter to specified locus_tags (if provided)
  # ============================================================
  # Detect "all" mode — no filtering needed
  is_all <- base::length(locus_tags) == 1 &&
            base::tolower(locus_tags[1]) == "all"

  if (!is_all && !base::is.null(locus_tags) &&
      "locus_tag" %in% base::colnames(gRNA_df)) {
    available_tags <- base::unique(gRNA_df$locus_tag[!base::is.na(gRNA_df$locus_tag)])
    missing_tags <- base::setdiff(locus_tags, available_tags)
    if (base::length(missing_tags) > 0) {
      base::warning("Locus tag(s) not found in gRNA_df: ",
                     base::paste(missing_tags, collapse = ", "))
    }
    gRNA_df <- gRNA_df[gRNA_df$locus_tag %in% locus_tags, , drop = FALSE]
    if (base::nrow(gRNA_df) == 0) {
      base::stop("No gRNAs found for the specified locus_tags.")
    }
    base::cat("  Filtered to ", base::length(locus_tags), " locus tag(s): ",
              base::nrow(gRNA_df), " gRNAs\n", sep = "")
  } else if (is_all) {
    base::cat("  Genome-wide mode: ", base::nrow(gRNA_df), " gRNAs from all genes\n", sep = "")
  }

  # ============================================================
  # Phase 3: Methylation filtering (Mode 2 only — Mode 1 already done)
  # ============================================================
  if (!methylation_applied &&
      !base::is.null(methylation_patterns) &&
      base::length(methylation_patterns) > 0 &&
      base::nrow(gRNA_df) > 0) {
    base::cat("\nFiltering methylation site conflicts...\n")
    n_before <- base::nrow(gRNA_df)
    gRNA_df <- filter_methylation_sites(
      gRNA_df = gRNA_df,
      methylation_patterns = methylation_patterns,
      nuclease = nuclease,
      exclude = TRUE
    )
    n_after <- base::nrow(gRNA_df)
    base::cat("  Methylation filter: ", n_before, " -> ", n_after,
              " gRNAs (removed ", n_before - n_after, ")\n", sep = "")
    if (base::nrow(gRNA_df) == 0) {
      base::stop("No gRNAs remaining after methylation filtering. ",
                 "Consider relaxing methylation_patterns or expanding the library.")
    }
  }

  # ============================================================
  # Phase 4: Select best gRNA(s)
  # ============================================================
  multi_gene_region <- !base::is.null(locus_tags) && base::length(locus_tags) > 1

  if (multi_gene_region) {
    # ----- Multi-gene region: pool all gRNAs, rank globally -----
    region_label <- base::paste0(locus_tags[1], "--",
                                  locus_tags[base::length(locus_tags)])
    base::cat("\n  Multi-gene region (", region_label, "): pooling ",
              base::nrow(gRNA_df), " gRNAs for global ranking\n", sep = "")

    # Rank globally by composite_score (best first)
    if ("composite_score" %in% base::colnames(gRNA_df)) {
      gRNA_df <- gRNA_df[base::order(gRNA_df$composite_score, decreasing = TRUE), ]
    }

    # Assign new global rank
    gRNA_df$gRNA_rank <- base::seq_len(base::nrow(gRNA_df))

    # Select top_n from the combined pool
    n_select <- base::min(top_n, base::nrow(gRNA_df))
    selected_df <- gRNA_df[base::seq_len(n_select), , drop = FALSE]

  } else {
    # ----- Single gene or per-gene selection -----
    all_tags <- base::unique(gRNA_df$locus_tag[!base::is.na(gRNA_df$locus_tag)])
    if (base::length(all_tags) == 0) all_tags <- "unknown"

    selected_rows <- base::list()
    for (tag in all_tags) {
      tag_df <- gRNA_df[gRNA_df$locus_tag == tag, , drop = FALSE]

      if ("gRNA_rank" %in% base::colnames(tag_df)) {
        tag_df <- tag_df[base::order(tag_df$gRNA_rank, decreasing = FALSE), ]
      } else if ("composite_score" %in% base::colnames(tag_df)) {
        tag_df <- tag_df[base::order(tag_df$composite_score, decreasing = TRUE), ]
      }

      n_select <- base::min(top_n, base::nrow(tag_df))
      selected_rows[[tag]] <- tag_df[base::seq_len(n_select), , drop = FALSE]
    }

    selected_df <- base::do.call(base::rbind, selected_rows)
    base::rownames(selected_df) <- NULL
  }

  base::cat("  Selected ", base::nrow(selected_df), " gRNA(s)",
            base::ifelse(multi_gene_region, " (global)", " (per gene)"), "\n", sep = "")

  # Log selected spacers
  for (i in base::seq_len(base::nrow(selected_df))) {
    spacer_i <- selected_df$protospacer[i]
    tag_i <- base::ifelse("locus_tag" %in% base::colnames(selected_df),
                           base::as.character(selected_df$locus_tag[i]), "")
    rank_i <- base::ifelse("gRNA_rank" %in% base::colnames(selected_df),
                            selected_df$gRNA_rank[i], i)
    score_i <- NA
    if ("composite_score" %in% base::colnames(selected_df)) {
      score_i <- base::round(selected_df$composite_score[i], 1)
    }
    gene_i <- ""
    if ("gene" %in% base::colnames(selected_df)) {
      gene_i <- base::as.character(selected_df$gene[i])
      if (base::is.na(gene_i)) gene_i <- ""
    }
    gene_str <- base::ifelse(base::nchar(gene_i) > 0,
                              base::paste0(" (", gene_i, ")"), "")

    base::cat("    ", tag_i, gene_str, " g", rank_i,
              ": ", spacer_i, " (", base::nchar(spacer_i), "bp",
              base::ifelse(!base::is.na(score_i),
                            base::paste0(", score=", score_i), ""),
              ")\n", sep = "")
  }

  # Report spacer vs stuffer size
  spacer_lengths <- base::nchar(selected_df$protospacer)
  unique_lengths <- base::unique(spacer_lengths)
  if (base::length(unique_lengths) == 1) {
    base::cat("  Spacer length: ", unique_lengths[1], " bp -> replaces ",
              stuffer_len, " bp stuffer (delta: ",
              unique_lengths[1] - stuffer_len, " bp)\n", sep = "")
  } else {
    base::cat("  Spacer lengths: ", base::paste(unique_lengths, collapse = "/"),
              " bp -> replaces ", stuffer_len, " bp stuffer\n", sep = "")
  }

  # ============================================================
  # Phase 4.5: Enrich with gRNA reference columns
  # ============================================================
  base::cat("\nAdding gRNA reference information...\n")

  # spacer_length: quick reference for oligo ordering
  selected_df$spacer_length <- base::nchar(selected_df$protospacer)

  # spacer_context: visual representation of the targeting site
  # 3prime PAM (GeoCas9, SpCas9): spacer|PAM  (e.g., ATCGATCG...|NNNNCAAA)
  # 5prime PAM (FnCas12a):        PAM|spacer  (e.g., TTV|ATCGATCG...)
  pam_side <- .get_pam_side(nuclease)
  if ("pam" %in% base::colnames(selected_df)) {
    selected_df$spacer_context <- base::sapply(
      base::seq_len(base::nrow(selected_df)),
      function(i) {
        sp <- base::as.character(selected_df$protospacer[i])
        pm <- base::as.character(selected_df$pam[i])
        if (base::is.na(sp) || base::is.na(pm)) return(NA_character_)
        if (pam_side == "3prime") {
          base::paste0(sp, "|", pm)
        } else {
          base::paste0(pm, "|", sp)
        }
      }
    )
  }

  # target_gene_start / target_gene_end: gene coordinates from GenBank
  # (already present from merge in run_gRNA_list_generator: start, end columns)
  # Rename to avoid confusion with vector start/end
  if ("start" %in% base::colnames(selected_df) &&
      !"gene_start" %in% base::colnames(selected_df)) {
    base::colnames(selected_df)[base::colnames(selected_df) == "start"] <- "gene_start"
  }
  if ("end" %in% base::colnames(selected_df) &&
      !"gene_end" %in% base::colnames(selected_df)) {
    base::colnames(selected_df)[base::colnames(selected_df) == "end"] <- "gene_end"
  }

  # pam_site_position: genomic position of the PAM (from crisprDesign)
  # Already present as pam_site if available

  # direction: gene strand from GenBank (already present from merge)

  base::cat("  Reference columns added: spacer_length, spacer_context",
            base::ifelse("gene_start" %in% base::colnames(selected_df),
                          ", gene_start, gene_end", ""), "\n")

  # ============================================================
  # Phase 5: Design cloning primers
  # ============================================================
  base::cat("\nDesigning cloning primers...\n")

  result <- design_cloning_primers(
    gRNA_df          = selected_df,
    vector_file      = vector_file,
    start            = start_i,
    end              = end_i,
    cloning_method   = cloning_method,
    enzyme           = enzyme,
    custom_overhangs = custom_overhangs,
    tm_target        = tm_target,
    overlap_length   = overlap_length,
    add_g            = add_g,
    name_pattern     = name_pattern,
    output_file      = output_file,
    output_dir       = output_dir,
    kill_snapgene    = kill_snapgene
  )

  # ============================================================
  # Phase 6: Summary
  # ============================================================
  base::cat("\n=== Construct Design Complete ===\n")
  base::cat("  Total constructs: ", base::nrow(result), "\n", sep = "")
  if (!base::is.null(output_file)) {
    base::cat("  Primers saved to: ", output_file, "\n", sep = "")
  }
  if (!base::is.null(output_dir)) {
    base::cat("  GenBank files in: ", output_dir, "\n", sep = "")
    base::cat("  (Open .gbk files in SnapGene to visualize constructs)\n")
  }
  if (auto_mode) {
    lib_path <- base::file.path(dir, lib_file)
    base::cat("  gRNA library: ", lib_path, "\n", sep = "")
  }

  return(result)
}
