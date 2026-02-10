#' Write gRNA Expression Vector GenBank File
#'
#' Generates a modified GenBank file with a gRNA spacer sequence inserted into
#' the expression vector at the specified insertion site. The output file can be
#' opened in SnapGene to visualize the final gRNA expression construct with
#' spacer and primer annotations.
#'
#' This function replaces the stuffer/placeholder sequence at the insertion site
#' (defined by \code{start} and \code{end}) with the gRNA spacer, adjusts all
#' existing feature locations accordingly, and adds annotations for:
#' \itemize{
#'   \item The gRNA spacer (misc_feature with label, nuclease, target info)
#'   \item Forward and reverse primer binding sites (if primer_info provided)
#' }
#'
#' @param vector_file Path to the original vector file (.dna, .gb, .fasta).
#' @param spacer_seq Character string. The gRNA spacer sequence to insert.
#' @param start Integer. Vector insertion site start position (1-based, inclusive).
#' @param end Integer. Vector insertion site end position (1-based, inclusive).
#' @param gRNA_name Character. Label for the gRNA (e.g., "Geo_dnaA_g1"). Default: NULL.
#' @param locus_tag Character. Target locus tag (optional, for annotation).
#' @param gene Character. Target gene name (optional, for annotation).
#' @param nuclease Character. Nuclease name (optional, for annotation).
#' @param primer_info Named list with primer details (optional):
#'   \code{primer_F}, \code{primer_R}, \code{primer_F_Tm}, \code{primer_R_Tm},
#'   \code{primer_F_name}, \code{primer_R_name}, \code{cloning_method}, \code{enzyme}.
#' @param output_path Character. Path for the output GenBank file.
#' @param kill_snapgene Logical. Kill SnapGene before reading .dna (default: FALSE).
#' @return Invisible. The output file path.
#' @examples
#' \dontrun{
#' # Simple usage: insert spacer into vector
#' write_grna_vector_genbank(
#'   vector_file = "my_vector.dna",
#'   spacer_seq = "ATCGATCGATCGATCGATCGA",
#'   start = 3977, end = 3987,
#'   gRNA_name = "Geo_dnaA_g1",
#'   gene = "dnaA",
#'   nuclease = "GeoCas9",
#'   output_path = "Geo_dnaA_g1_vector.gbk"
#' )
#'
#' # With primer info (from design_cloning_primers output)
#' write_grna_vector_genbank(
#'   vector_file = "my_vector.dna",
#'   spacer_seq = "ATCGATCGATCGATCGATCGA",
#'   start = 3977, end = 3987,
#'   gRNA_name = "Geo_dnaA_g1",
#'   nuclease = "GeoCas9",
#'   primer_info = list(
#'     primer_F = "CACCGATCGATCGATCGATCGATCGA",
#'     primer_R = "AAACTCGATCGATCGATCGATCGATC",
#'     primer_F_Tm = 62.1,
#'     primer_R_Tm = 61.8,
#'     primer_F_name = "dnaA_g1_F",
#'     primer_R_name = "dnaA_g1_R",
#'     cloning_method = "Golden_Gate",
#'     enzyme = "BbsI"
#'   ),
#'   output_path = "Geo_dnaA_g1_vector.gbk"
#' )
#' }
#' @export
write_grna_vector_genbank <- function(
    vector_file,
    spacer_seq,
    start,
    end,
    gRNA_name = NULL,
    locus_tag = NULL,
    gene = NULL,
    nuclease = NULL,
    primer_info = NULL,
    output_path,
    kill_snapgene = FALSE
) {
  # --- Input validation ---
  if (base::missing(vector_file) || !base::is.character(vector_file)) {
    base::stop("vector_file must be a valid file path.")
  }
  if (base::missing(spacer_seq) || !base::is.character(spacer_seq) ||
      base::nchar(spacer_seq) == 0) {
    base::stop("spacer_seq must be a non-empty DNA sequence string.")
  }
  if (base::missing(start) || base::missing(end)) {
    base::stop("start and end positions are required.")
  }

  # --- Step 1: Get vector GenBank content ---
  ext <- base::tolower(tools::file_ext(vector_file))
  vector_file <- base::normalizePath(vector_file, mustWork = TRUE)

  if (ext == "dna") {
    gb_lines <- .get_genbank_from_dna(vector_file, kill_snapgene)
  } else if (ext %in% base::c("gb", "gbk", "genbank")) {
    gb_lines <- base::readLines(vector_file, warn = FALSE)
  } else {
    # For FASTA or other formats, build a minimal GenBank template
    vec <- read_vector_file(vector_file)
    gb_lines <- .build_minimal_genbank_grna(vec$name, vec$sequence)
  }

  # --- Step 2: Parse GenBank sections ---
  header_lines <- base::character(0)
  feature_lines <- base::character(0)
  in_features <- FALSE
  in_origin <- FALSE

  for (line in gb_lines) {
    if (base::startsWith(line, "FEATURES")) {
      in_features <- TRUE
      header_lines <- base::c(header_lines, line)
      next
    } else if (base::startsWith(line, "ORIGIN")) {
      in_features <- FALSE
      in_origin <- TRUE
      next
    } else if (base::startsWith(line, "//")) {
      next
    }

    if (in_features) {
      feature_lines <- base::c(feature_lines, line)
    } else if (!in_origin) {
      header_lines <- base::c(header_lines, line)
    }
    # skip ORIGIN sequence lines — we rebuild them
  }

  # --- Step 3: Replace insertion site with spacer ---
  vec <- read_vector_file(vector_file, kill_snapgene = FALSE)
  vec_seq <- vec$sequence
  spacer_seq <- base::toupper(spacer_seq)

  start_i <- base::as.integer(start)
  end_i <- base::as.integer(end)

  if (start_i < 1 || end_i > base::nchar(vec_seq) || start_i > end_i) {
    base::stop("Invalid start/end positions. Vector length: ", base::nchar(vec_seq),
               " bp, given: ", start_i, "-", end_i)
  }

  new_seq <- base::paste0(
    base::substring(vec_seq, 1, start_i - 1),
    spacer_seq,
    base::substring(vec_seq, end_i + 1, base::nchar(vec_seq))
  )
  new_length <- base::nchar(new_seq)

  # --- Step 4: Update header (LOCUS line length, construct name) ---
  # Prioritize locus_tag for construct label (user preference)
  construct_label <- ""
  if (!base::is.null(locus_tag) && base::nchar(locus_tag) > 0) {
    construct_label <- locus_tag
  } else if (!base::is.null(gRNA_name) && base::nchar(gRNA_name) > 0) {
    construct_label <- gRNA_name
  } else if (!base::is.null(gene) && base::nchar(gene) > 0) {
    construct_label <- gene
  } else {
    construct_label <- "gRNA"
  }

  for (i in base::seq_along(header_lines)) {
    line <- header_lines[i]
    if (base::startsWith(line, "LOCUS")) {
      parts <- base::strsplit(base::trimws(line), "\\s+")[[1]]
      parts[2] <- base::paste0(parts[2], "_", construct_label)
      parts[3] <- base::paste0(new_length, " bp")
      header_lines[i] <- base::paste(parts, collapse = "  ")
    } else if (base::grepl("^\\s*REFERENCE.*\\(bases", line)) {
      header_lines[i] <- base::sub(
        "\\(bases.*\\)",
        base::paste0("(bases 1 to ", new_length, ")"),
        line
      )
    }
  }

  # --- Step 5: Adjust existing feature locations ---
  length_diff <- base::nchar(spacer_seq) - (end_i - start_i + 1)
  adjusted_features <- base::character(0)
  current_block <- base::character(0)

  for (line in feature_lines) {
    trimmed <- base::trimws(line)
    # Detect new feature line (not a qualifier)
    if (base::nchar(trimmed) > 0 && !base::startsWith(trimmed, "/") &&
        base::grepl("\\.\\.|\\.\\^\\.", line) && !base::startsWith(trimmed, "/")) {
      # Save previous block
      if (base::length(current_block) > 0) {
        adjusted <- .adjust_feature_block(current_block, start_i, end_i, length_diff)
        if (!base::is.null(adjusted)) {
          adjusted_features <- base::c(adjusted_features, adjusted)
        }
      }
      current_block <- line
    } else if (base::nchar(trimmed) > 0 && !base::startsWith(trimmed, "/") &&
               base::grepl("^     \\S", line) && !base::grepl("^     /", line)) {
      # Another new feature line (no ..)
      if (base::length(current_block) > 0) {
        adjusted <- .adjust_feature_block(current_block, start_i, end_i, length_diff)
        if (!base::is.null(adjusted)) {
          adjusted_features <- base::c(adjusted_features, adjusted)
        }
      }
      current_block <- line
    } else {
      current_block <- base::c(current_block, line)
    }
  }
  # Process last block
  if (base::length(current_block) > 0) {
    adjusted <- .adjust_feature_block(current_block, start_i, end_i, length_diff)
    if (!base::is.null(adjusted)) {
      adjusted_features <- base::c(adjusted_features, adjusted)
    }
  }

  # --- Step 6: Add gRNA spacer annotation ---
  spacer_start <- start_i
  spacer_end <- start_i + base::nchar(spacer_seq) - 1

  # Spacer label: prioritize locus_tag
  spacer_label <- "gRNA_spacer"
  if (!base::is.null(locus_tag) && base::nchar(locus_tag) > 0) {
    spacer_label <- locus_tag
  } else if (!base::is.null(gRNA_name) && base::nchar(gRNA_name) > 0) {
    spacer_label <- gRNA_name
  }

  # Build spacer feature with all target info
  spacer_features <- base::c(
    base::sprintf("     misc_feature    %d..%d", spacer_start, spacer_end),
    base::sprintf("                     /label=%s", spacer_label)
  )

  # /locus_tag qualifier
  if (!base::is.null(locus_tag) && base::nchar(locus_tag) > 0) {
    spacer_features <- base::c(spacer_features,
      base::sprintf("                     /locus_tag=%s", locus_tag))
  }

  # /gene qualifier
  if (!base::is.null(gene) && base::nchar(gene) > 0) {
    spacer_features <- base::c(spacer_features,
      base::sprintf("                     /gene=%s", gene))
  }

  # /note with nuclease, spacer length, target info
  note_parts <- base::character(0)
  if (!base::is.null(nuclease) && base::nchar(nuclease) > 0) {
    spacer_features <- base::c(spacer_features,
      base::sprintf("                     /note=nuclease: %s", nuclease))
  }
  spacer_features <- base::c(spacer_features,
    base::sprintf("                     /note=spacer: %s (%d bp)",
                   spacer_seq, base::nchar(spacer_seq))
  )
  if (!base::is.null(gRNA_name) && base::nchar(gRNA_name) > 0) {
    spacer_features <- base::c(spacer_features,
      base::sprintf("                     /note=gRNA_name: %s", gRNA_name))
  }

  # --- Step 7: Add primer annotations (if provided) ---
  primer_features <- base::character(0)
  if (!base::is.null(primer_info)) {
    method_str <- ""
    if (!base::is.null(primer_info$cloning_method)) {
      method_str <- primer_info$cloning_method
    }
    enzyme_str <- ""
    if (!base::is.null(primer_info$enzyme)) {
      enzyme_str <- base::paste0(" (", primer_info$enzyme, ")")
    }

    # Determine primer binding positions based on cloning method
    spacer_len <- base::nchar(spacer_seq)
    is_gibson <- (!base::is.null(primer_info$cloning_method) &&
                  base::grepl("Gibson", primer_info$cloning_method, ignore.case = TRUE))

    if (is_gibson) {
      # IA primers (asymmetric Gibson overlap):
      #   F primer: 5'-[FULL spacer]-[backbone head]-3' (sense strand)
      #   R primer: 5'-[RC(spacer first N bp)]-[backbone head RC]-3' (complement)
      #
      # In the final construct, the spacer replaced the stuffer. On the construct:
      #   F covers: spacer + downstream backbone head (sense strand)
      #   R covers: upstream backbone head + partial spacer (complement strand)
      #
      # Annotation spans the FULL primer extent so SnapGene shows correct region.
      f_head_len <- base::nchar(primer_info$primer_F) - spacer_len

      # R primer: get Gibson overlap length (may be < spacer_len)
      r_overlap <- spacer_len  # default: full spacer
      if (!base::is.null(primer_info$gibson_overlap) &&
          !base::is.na(primer_info$gibson_overlap)) {
        r_overlap <- base::as.integer(primer_info$gibson_overlap)
      }
      r_head_len <- base::nchar(primer_info$primer_R) - r_overlap

      # F primer: spacer_start to spacer_end + head (full extent, sense)
      f_primer_start <- spacer_start
      f_primer_end <- spacer_end + base::max(0L, f_head_len)
      # R primer: (spacer_start - head) to (spacer_start + overlap - 1) (complement)
      r_primer_start <- spacer_start - base::max(0L, r_head_len)
      r_primer_end <- spacer_start + r_overlap - 1L
    } else {
      # OA primers: overhang + spacer (bind within spacer region only)
      f_primer_start <- spacer_start
      f_primer_end <- spacer_end
      r_primer_start <- spacer_start
      r_primer_end <- spacer_end
      r_overlap <- spacer_len
      f_head_len <- 0L
      r_head_len <- 0L
    }

    # Forward primer annotation
    if (!base::is.null(primer_info$primer_F) && base::nchar(primer_info$primer_F) > 0) {
      f_label <- base::ifelse(
        !base::is.null(primer_info$primer_F_name) && base::nchar(primer_info$primer_F_name) > 0,
        primer_info$primer_F_name,
        "Forward_primer"
      )
      f_tm_str <- base::ifelse(
        !base::is.null(primer_info$primer_F_Tm) && !base::is.na(primer_info$primer_F_Tm),
        base::sprintf("%.1f", primer_info$primer_F_Tm),
        "NA"
      )

      primer_features <- base::c(
        primer_features,
        "",
        base::sprintf("     primer          %d..%d", f_primer_start, f_primer_end),
        base::sprintf("                     /label=%s", f_label),
        base::sprintf("                     /note=Tm(head): %s C", f_tm_str),
        base::sprintf("                     /note=seq: %s (%d bp)",
                       primer_info$primer_F, base::nchar(primer_info$primer_F)),
        base::sprintf("                     /note=method: %s%s", method_str, enzyme_str)
      )
    }

    # Reverse primer annotation
    if (!base::is.null(primer_info$primer_R) && base::nchar(primer_info$primer_R) > 0) {
      r_label <- base::ifelse(
        !base::is.null(primer_info$primer_R_name) && base::nchar(primer_info$primer_R_name) > 0,
        primer_info$primer_R_name,
        "Reverse_primer"
      )
      r_tm_str <- base::ifelse(
        !base::is.null(primer_info$primer_R_Tm) && !base::is.na(primer_info$primer_R_Tm),
        base::sprintf("%.1f", primer_info$primer_R_Tm),
        "NA"
      )

      primer_features <- base::c(
        primer_features,
        "",
        base::sprintf("     primer          complement(%d..%d)", r_primer_start, r_primer_end),
        base::sprintf("                     /label=%s", r_label),
        base::sprintf("                     /note=Tm(head): %s C", r_tm_str),
        base::sprintf("                     /note=seq: %s (%d bp)",
                       primer_info$primer_R, base::nchar(primer_info$primer_R)),
        base::sprintf("                     /note=method: %s%s", method_str, enzyme_str)
      )
    }
  }

  # --- Step 8: Build ORIGIN section ---
  origin_lines <- "ORIGIN"
  seq_lower <- base::tolower(new_seq)
  for (i in base::seq(1, base::nchar(seq_lower), by = 60)) {
    chunk <- base::substring(seq_lower, i,
                              base::min(i + 59, base::nchar(seq_lower)))
    # Split into 10-char blocks
    blocks <- base::character(0)
    for (j in base::seq(1, base::nchar(chunk), by = 10)) {
      blocks <- base::c(blocks,
                          base::substring(chunk, j,
                                           base::min(j + 9, base::nchar(chunk))))
    }
    formatted <- base::paste(blocks, collapse = " ")
    pos_str <- base::as.character(i)
    padding <- base::paste(base::rep(" ", 9 - base::nchar(pos_str)), collapse = "")
    origin_lines <- base::c(origin_lines,
                              base::paste0(padding, pos_str, " ", formatted))
  }
  origin_lines <- base::c(origin_lines, "//")

  # --- Step 9: Assemble and write ---
  all_lines <- base::c(
    header_lines,
    adjusted_features,
    spacer_features,
    primer_features,
    origin_lines
  )

  output_dir_path <- base::dirname(output_path)
  if (!base::dir.exists(output_dir_path) && output_dir_path != ".") {
    base::dir.create(output_dir_path, recursive = TRUE)
  }
  base::writeLines(all_lines, output_path)
  base::cat("gRNA vector GenBank saved:", output_path, "(", new_length, "bp)\n")

  invisible(output_path)
}


# ============================================================
# Internal helper
# ============================================================

#' Build minimal GenBank template for gRNA expression vector
#' @keywords internal
.build_minimal_genbank_grna <- function(name, sequence) {
  len <- base::nchar(sequence)
  base::c(
    base::sprintf("LOCUS       %s  %d bp  DNA  circular  UNK", name, len),
    base::sprintf("DEFINITION  %s gRNA expression vector.", name),
    base::sprintf("ACCESSION   %s", name),
    base::sprintf("VERSION     %s", name),
    "FEATURES             Location/Qualifiers",
    base::sprintf("     source          1..%d", len),
    base::sprintf("                     /organism=%s", name),
    "ORIGIN"
  )
}
