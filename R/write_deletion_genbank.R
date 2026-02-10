#' Write Deletion Vector GenBank File
#'
#' Generates a modified GenBank file for a deletion vector construct.
#' Takes the original vector GenBank (converted from .dna if needed),
#' replaces the insertion site with upstream + downstream arms,
#' adjusts all existing feature locations, and adds annotations for:
#' arm regions (misc_feature) and 4 primers (primer feature with Tm info).
#'
#' The output .gbk file can be opened directly in SnapGene to visualize
#' the deletion construct with all primer binding sites.
#'
#' @param vector_file Path to the original vector file (.dna, .gb, .fasta).
#' @param insert_seq Character string. Combined arm sequence (upstream + downstream).
#' @param start Integer. Vector insertion site start position (1-based, inclusive).
#' @param end Integer. Vector insertion site end position (1-based, inclusive).
#' @param locus_tag Character. Target locus tag for annotation labels.
#' @param primers Named list. Output from \code{design_deletion_primers()}.
#' @param upstream_bp Integer. Length of upstream arm.
#' @param downstream_bp Integer. Length of downstream arm.
#' @param output_path Character. Path for the output GenBank file.
#' @param kill_snapgene Logical. Kill SnapGene before reading .dna (default: FALSE).
#' @return Invisible. The output file path.
#' @examples
#' \dontrun{
#' write_deletion_genbank(
#'   vector_file = "my_vector.dna",
#'   insert_seq = base::paste0(upstream_arm, downstream_arm),
#'   start = 1234, end = 1256,
#'   locus_tag = "QT234_RS00005",
#'   primers = primer_result,
#'   upstream_bp = 500, downstream_bp = 500,
#'   output_path = "QT234_RS00005_deletion_vector.gbk"
#' )
#' }
#' @export
write_deletion_genbank <- function(
    vector_file,
    insert_seq,
    start,
    end,
    locus_tag,
    primers,
    upstream_bp,
    downstream_bp,
    output_path,
    kill_snapgene = FALSE
) {
  # --- Step 1: Get vector GenBank content ---
  # We need the original GenBank text (with features) as a template.
  # If .dna file, convert to GenBank first via SnapGene CLI or read sequence.
  ext <- base::tolower(tools::file_ext(vector_file))
  vector_file <- base::normalizePath(vector_file, mustWork = TRUE)

  if (ext == "dna") {
    gb_lines <- .get_genbank_from_dna(vector_file, kill_snapgene)
  } else if (ext %in% base::c("gb", "gbk", "genbank")) {
    gb_lines <- base::readLines(vector_file, warn = FALSE)
  } else {
    # For FASTA, build a minimal GenBank template
    vec <- read_vector_file(vector_file)
    gb_lines <- .build_minimal_genbank(vec$name, vec$sequence)
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

  # --- Step 3: Replace sequence ---
  # Read original vector sequence
  vec <- read_vector_file(vector_file, kill_snapgene = FALSE)
  vec_seq <- vec$sequence
  insert_seq <- base::toupper(insert_seq)

  start_i <- base::as.integer(start)
  end_i <- base::as.integer(end)

  new_seq <- base::paste0(
    base::substring(vec_seq, 1, start_i - 1),
    insert_seq,
    base::substring(vec_seq, end_i + 1, base::nchar(vec_seq))
  )
  new_length <- base::nchar(new_seq)

  # --- Step 4: Update header (LOCUS line length, name) ---
  for (i in base::seq_along(header_lines)) {
    line <- header_lines[i]
    if (base::startsWith(line, "LOCUS")) {
      parts <- base::strsplit(base::trimws(line), "\\s+")[[1]]
      parts[2] <- base::paste0(parts[2], "_", locus_tag, "_deletion")
      parts[3] <- base::paste0(new_length, " bp")
      header_lines[i] <- base::paste(parts, collapse = "  ")
    } else if (base::grepl("^\\s*REFERENCE.*\\(bases", line)) {
      header_lines[i] <- base::sub("\\(bases.*\\)", base::paste0("(bases 1 to ", new_length, ")"), line)
    }
  }

  # --- Step 5: Adjust existing feature locations ---
  length_diff <- base::nchar(insert_seq) - (end_i - start_i + 1)
  adjusted_features <- base::character(0)
  current_block <- base::character(0)

  for (line in feature_lines) {
    trimmed <- base::trimws(line)
    # New feature starts (not a qualifier, not empty)
    if (base::nchar(trimmed) > 0 && !base::startsWith(trimmed, "/") && base::grepl("\\.\\.|\\.\\^\\.", line) &&
        !base::startsWith(trimmed, "/")) {
      # Process previous block
      if (base::length(current_block) > 0) {
        adjusted <- .adjust_feature_block(current_block, start_i, end_i, length_diff)
        if (!base::is.null(adjusted)) adjusted_features <- base::c(adjusted_features, adjusted)
      }
      current_block <- line
    } else if (base::nchar(trimmed) > 0 && !base::startsWith(trimmed, "/") &&
               base::grepl("^     \\S", line) && !base::grepl("^     /", line)) {
      # Also a new feature line (no ..)
      if (base::length(current_block) > 0) {
        adjusted <- .adjust_feature_block(current_block, start_i, end_i, length_diff)
        if (!base::is.null(adjusted)) adjusted_features <- base::c(adjusted_features, adjusted)
      }
      current_block <- line
    } else {
      current_block <- base::c(current_block, line)
    }
  }
  # Process last block
  if (base::length(current_block) > 0) {
    adjusted <- .adjust_feature_block(current_block, start_i, end_i, length_diff)
    if (!base::is.null(adjusted)) adjusted_features <- base::c(adjusted_features, adjusted)
  }

  # --- Step 6: Add arm annotations (misc_feature) ---
  # In the new vector, insert starts at position start_i (1-based)
  arm_up_start <- start_i
  arm_up_end <- start_i + upstream_bp - 1
  arm_dn_start <- start_i + upstream_bp
  arm_dn_end <- start_i + upstream_bp + downstream_bp - 1

  arm_features <- base::c(
    base::sprintf("     misc_feature    %d..%d", arm_up_start, arm_up_end),
    base::sprintf("                     /label=%s_upstream", locus_tag),
    base::sprintf("                     /note=upstream arm for %s deletion (%d bp)", locus_tag, upstream_bp),
    "",
    base::sprintf("     misc_feature    %d..%d", arm_dn_start, arm_dn_end),
    base::sprintf("                     /label=%s_downstream", locus_tag),
    base::sprintf("                     /note=downstream arm for %s deletion (%d bp)", locus_tag, downstream_bp)
  )

  # --- Step 7: Add primer annotations ---
  # Use primer names from result list (e.g., [locus_tag]_UF, _UR, _DF, _DR)
  uf_label <- base::ifelse(!base::is.null(primers$upstream_forward_name),
                            primers$upstream_forward_name, "Upstream_Forward")
  ur_label <- base::ifelse(!base::is.null(primers$upstream_reverse_name),
                            primers$upstream_reverse_name, "Upstream_Reverse")
  df_label <- base::ifelse(!base::is.null(primers$downstream_forward_name),
                            primers$downstream_forward_name, "Downstream_Forward")
  dr_label <- base::ifelse(!base::is.null(primers$downstream_reverse_name),
                            primers$downstream_reverse_name, "Downstream_Reverse")

  primer_features <- base::c(
    # Upstream Forward (sense strand)
    base::sprintf("     primer          %d..%d", primers$upstream_forward_start, primers$upstream_forward_end),
    base::sprintf("                     /label=%s", uf_label),
    base::sprintf("                     /note=Tm(target): %.2f C, Tm(full): %.2f C",
            primers$upstream_forward_tm_target, primers$upstream_forward_tm_full),
    base::sprintf("                     /note=seq: %s", primers$upstream_forward_primer),
    "",
    # Upstream Reverse (complement)
    base::sprintf("     primer          complement(%d..%d)", primers$upstream_reverse_start, primers$upstream_reverse_end),
    base::sprintf("                     /label=%s", ur_label),
    base::sprintf("                     /note=Tm(target): %.2f C, Tm(full): %.2f C",
            primers$upstream_reverse_tm_target, primers$upstream_reverse_tm_full),
    base::sprintf("                     /note=seq: %s", primers$upstream_reverse_primer),
    "",
    # Downstream Forward (sense strand)
    base::sprintf("     primer          %d..%d", primers$downstream_forward_start, primers$downstream_forward_end),
    base::sprintf("                     /label=%s", df_label),
    base::sprintf("                     /note=Tm(target): %.2f C, Tm(full): %.2f C",
            primers$downstream_forward_tm_target, primers$downstream_forward_tm_full),
    base::sprintf("                     /note=seq: %s", primers$downstream_forward_primer),
    "",
    # Downstream Reverse (complement)
    base::sprintf("     primer          complement(%d..%d)", primers$downstream_reverse_start, primers$downstream_reverse_end),
    base::sprintf("                     /label=%s", dr_label),
    base::sprintf("                     /note=Tm(target): %.2f C, Tm(full): %.2f C",
            primers$downstream_reverse_tm_target, primers$downstream_reverse_tm_full),
    base::sprintf("                     /note=seq: %s", primers$downstream_reverse_primer)
  )

  # --- Step 8: Build ORIGIN section ---
  origin_lines <- "ORIGIN"
  seq_lower <- base::tolower(new_seq)
  for (i in base::seq(1, base::nchar(seq_lower), by = 60)) {
    chunk <- base::substring(seq_lower, i, base::min(i + 59, base::nchar(seq_lower)))
    # Split into 10-char blocks
    blocks <- base::character(0)
    for (j in base::seq(1, base::nchar(chunk), by = 10)) {
      blocks <- base::c(blocks, base::substring(chunk, j, base::min(j + 9, base::nchar(chunk))))
    }
    formatted <- base::paste(blocks, collapse = " ")
    pos_str <- base::as.character(i)
    padding <- base::paste(base::rep(" ", 9 - base::nchar(pos_str)), collapse = "")
    origin_lines <- base::c(origin_lines, base::paste0(padding, pos_str, " ", formatted))
  }
  origin_lines <- base::c(origin_lines, "//")

  # --- Step 9: Assemble and write ---
  all_lines <- base::c(
    header_lines,
    adjusted_features,
    arm_features,
    primer_features,
    origin_lines
  )

  output_dir <- base::dirname(output_path)
  if (!base::dir.exists(output_dir) && output_dir != ".") {
    base::dir.create(output_dir, recursive = TRUE)
  }
  base::writeLines(all_lines, output_path)
  base::cat("GenBank file saved:", output_path, "(", new_length, "bp)\n")

  invisible(output_path)
}


# ============================================================
# Internal helpers for GenBank writing
# ============================================================

#' Get GenBank lines from .dna file (via SnapGene CLI)
#' @keywords internal
.get_genbank_from_dna <- function(dna_file, kill_snapgene = FALSE) {
  if (kill_snapgene) .kill_snapgene()

  # Try SnapGene CLI conversion first

  snapgene_paths <- base::c(
    "/Applications/SnapGene.app/Contents/MacOS/SnapGene",
    "C:/Program Files/SnapGene/snapgene.exe",
    base::Sys.which("snapgene")
  )
  snapgene_bin <- NULL
  for (p in snapgene_paths) {
    if (base::nchar(p) > 0 && base::file.exists(p)) {
      snapgene_bin <- p
      break
    }
  }

  if (!base::is.null(snapgene_bin)) {
    temp_gb <- tempfile(fileext = ".gbk")
    cmd <- base::sprintf('%s --convert "GenBank - SnapGene" --input "%s" --output "%s"',
                   snapgene_bin, dna_file, temp_gb)
    exit_code <- base::system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

    if (exit_code == 0 && base::file.exists(temp_gb) && base::file.info(temp_gb)$size > 0) {
      lines <- base::readLines(temp_gb, warn = FALSE)
      base::file.remove(temp_gb)
      return(lines)
    }
  }

  # Fallback: build minimal GenBank from sequence
  vec <- read_vector_file(dna_file, kill_snapgene = FALSE)
  return(.build_minimal_genbank(vec$name, vec$sequence))
}


#' Build minimal GenBank from name + sequence
#' @keywords internal
.build_minimal_genbank <- function(name, sequence) {
  len <- base::nchar(sequence)
  base::c(
    base::sprintf("LOCUS       %s  %d bp  DNA  circular  UNK", name, len),
    base::sprintf("DEFINITION  %s deletion vector.", name),
    base::sprintf("ACCESSION   %s", name),
    base::sprintf("VERSION     %s", name),
    "FEATURES             Location/Qualifiers",
    base::sprintf("     source          1..%d", len),
    base::sprintf("                     /organism=%s", name),
    "ORIGIN"
  )
}


#' Adjust a single coordinate based on insertion site
#' @keywords internal
.adjust_one_coord <- function(coord, ins_start, ins_end, length_diff) {
  # coord is AFTER insertion site → shift
  if (coord > ins_end) return(coord + length_diff)
  # coord is BEFORE insertion site → keep
  if (coord < ins_start) return(coord)
  # coord is INSIDE insertion site → mark as removed (-1)
  return(-1L)
}

#' Adjust a feature block's location based on insertion length difference
#'
#' Handles simple ranges, complement(), join(), and circular wrap-around
#' features (where start > end, e.g., 12000..27). Each coordinate is
#' adjusted independently, so wrap-around features are preserved correctly.
#'
#' @keywords internal
.adjust_feature_block <- function(block, start, end, length_diff) {
  first_line <- block[1]
  qualifiers <- if (base::length(block) > 1) block[2:base::length(block)] else base::character(0)

  # Extract feature type
  feat_type <- base::trimws(base::sub("\\s+.*", "", base::trimws(first_line)))

  is_complement <- base::grepl("complement", first_line)
  is_join <- base::grepl("join", first_line)

  # ---- Helper: adjust a single "start..end" segment ----
  .adjust_segment <- function(seg_str) {
    nums <- base::regmatches(seg_str, base::gregexpr("[0-9]+", seg_str))[[1]]
    if (base::length(nums) < 2) return(seg_str)

    s <- base::as.integer(nums[1])
    e <- base::as.integer(nums[2])
    is_wraparound <- (s > e)

    new_s <- .adjust_one_coord(s, start, end, length_diff)
    new_e <- .adjust_one_coord(e, start, end, length_diff)

    # If either coordinate hit the insertion site → remove segment
    if (new_s < 0 || new_e < 0) return(NULL)

    # Preserve partial-end markers (< >)
    s_prefix <- base::ifelse(base::grepl("^<", seg_str), "<", "")
    e_prefix <- base::ifelse(base::grepl("\\.\\.<", seg_str) ||
                              base::grepl("\\.\\.>", seg_str),
                              base::sub(".*\\.\\.", "", base::sub("[0-9]+$", "", seg_str)),
                              "")

    base::paste0(s_prefix, new_s, "..", e_prefix, new_e)
  }

  # ---- Join features: parse each segment ----
  if (is_join) {
    # Extract full location string (may span multiple lines)
    full_loc <- first_line
    # Find all "number..number" segments
    segments <- base::regmatches(full_loc, base::gregexpr("[0-9<>]+\\.\\.[0-9<>]+", full_loc))[[1]]
    if (base::length(segments) == 0) return(base::c(block))

    new_segments <- base::character(0)
    for (seg in segments) {
      adj <- .adjust_segment(seg)
      if (!base::is.null(adj)) new_segments <- base::c(new_segments, adj)
    }

    if (base::length(new_segments) == 0) return(NULL)  # all segments removed

    if (base::length(new_segments) == 1) {
      new_loc <- new_segments[1]
    } else {
      new_loc <- base::paste0("join(", base::paste(new_segments, collapse = ","), ")")
    }
    if (is_complement) new_loc <- base::paste0("complement(", new_loc, ")")
    new_first <- base::sprintf("     %-16s%s", feat_type, new_loc)
    return(base::c(new_first, qualifiers))
  }

  # ---- Simple location (possibly wrap-around) ----
  loc_match <- base::regmatches(first_line, base::regexpr("[0-9<>]+\\.\\.[0-9<>]+", first_line))
  if (base::length(loc_match) == 0) return(base::c(block))

  nums <- base::as.integer(base::regmatches(loc_match, base::gregexpr("[0-9]+", loc_match))[[1]])
  if (base::length(nums) < 2) return(base::c(block))

  loc_start <- nums[1]
  loc_end   <- nums[2]
  is_wraparound <- (loc_start > loc_end)

  if (is_wraparound) {
    # Circular wrap-around: adjust each coordinate independently
    new_s <- .adjust_one_coord(loc_start, start, end, length_diff)
    new_e <- .adjust_one_coord(loc_end,   start, end, length_diff)
    if (new_s < 0 || new_e < 0) return(NULL)
    new_loc <- base::paste0(new_s, "..", new_e)
    if (is_complement) new_loc <- base::paste0("complement(", new_loc, ")")
    new_first <- base::sprintf("     %-16s%s", feat_type, new_loc)
    return(base::c(new_first, qualifiers))
  }

  # Normal linear range
  # Feature entirely before insertion → keep as-is
  if (loc_end < start) return(block)

  # Feature entirely after insertion → shift both
  if (loc_start > end) {
    new_start <- loc_start + length_diff
    new_end   <- loc_end   + length_diff
    new_loc <- base::paste0(new_start, "..", new_end)
    if (is_complement) new_loc <- base::paste0("complement(", new_loc, ")")
    new_first <- base::sprintf("     %-16s%s", feat_type, new_loc)
    return(base::c(new_first, qualifiers))
  }

  # Feature overlaps insertion site → remove
  return(NULL)
}
