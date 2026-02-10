#' Read Genome GenBank File and Extract Sequence + Feature Table
#'
#' Parses a genome GenBank file (.gb, .gbk, .gbff) to extract both the
#' nucleotide sequence and a feature annotation table. This allows users
#' to provide a single GenBank file instead of separate \code{genome_seq}
#' and \code{genbank_table} arguments.
#'
#' The function extracts CDS features to build a data frame with columns:
#' \code{locus_tag}, \code{start}, \code{end}, \code{strand}, \code{gene},
#' \code{product}, and \code{nt_seq} (nucleotide sequence of each CDS).
#'
#' @param genbank_file Character string. Path to a genome GenBank file
#'   (.gb, .gbk, .gbff).
#' @return A list with components:
#'   \describe{
#'     \item{genome_seq}{Character string. The full genome nucleotide sequence (uppercase).}
#'     \item{genbank_table}{Data frame with columns: locus_tag, start, end, strand, gene, product, nt_seq.}
#'     \item{organism}{Character string. Organism name from the GenBank header.}
#'     \item{accession}{Character string. Accession from the GenBank header.}
#'     \item{length}{Integer. Genome length in bp.}
#'   }
#' @examples
#' \dontrun{
#' # Read a genome GenBank file
#' genome <- read_genome_genbank("my_genome.gbk")
#'
#' # Use directly with deletion primer design
#' result <- batch_deletion_primers(
#'   vector_file = "my_vector.dna",
#'   genome_seq = genome$genome_seq,
#'   genbank_table = genome$genbank_table,
#'   locus_tags = c("QT234_RS00005", "QT234_RS00010"),
#'   start = 1234, end = 1256
#' )
#'
#' # Or use via the genbank_file parameter
#' result <- batch_deletion_primers(
#'   vector_file = "my_vector.dna",
#'   genbank_file = "my_genome.gbk",
#'   locus_tags = c("QT234_RS00005", "QT234_RS00010"),
#'   start = 1234, end = 1256
#' )
#' }
#' @export
read_genome_genbank <- function(genbank_file) {
  # --- Input validation ---
  if (base::missing(genbank_file) || !base::is.character(genbank_file) ||
      base::nchar(genbank_file) == 0) {
    base::stop("genbank_file must be a non-empty file path.")
  }
  if (!base::file.exists(genbank_file)) {
    base::stop("GenBank file not found: ", genbank_file)
  }

  ext <- base::tolower(tools::file_ext(genbank_file))
  if (!ext %in% c("gb", "gbk", "genbank", "gbff")) {
    base::warning("Unexpected file extension '.", ext,
                  "'. Expected: .gb, .gbk, .gbff, .genbank. Attempting to parse anyway.")
  }

  base::cat("Reading genome GenBank file:", genbank_file, "\n")
  lines <- base::readLines(genbank_file, warn = FALSE)

  if (base::length(lines) == 0) {
    base::stop("GenBank file is empty: ", genbank_file)
  }

  # --- Parse LOCUS line ---
  locus_line <- base::grep("^LOCUS", lines, value = TRUE)
  accession_val <- ""
  organism_val <- ""

  if (base::length(locus_line) > 0) {
    locus_parts <- base::strsplit(base::trimws(locus_line[1]), "\\s+")[[1]]
    if (base::length(locus_parts) >= 2) {
      accession_val <- locus_parts[2]
    }
  }

  # --- Parse ACCESSION line ---
  accession_line <- base::grep("^ACCESSION", lines, value = TRUE)
  if (base::length(accession_line) > 0) {
    acc_parts <- base::strsplit(base::trimws(accession_line[1]), "\\s+")[[1]]
    if (base::length(acc_parts) >= 2) {
      accession_val <- acc_parts[2]
    }
  }

  # --- Parse organism name ---
  organism_idx <- base::grep("^\\s+/organism=", lines)
  if (base::length(organism_idx) > 0) {
    org_line <- lines[organism_idx[1]]
    organism_val <- base::sub("^\\s+/organism=\"?", "", org_line)
    organism_val <- base::sub("\"?\\s*$", "", organism_val)
  } else {
    # Try from DEFINITION line
    def_line <- base::grep("^DEFINITION", lines, value = TRUE)
    if (base::length(def_line) > 0) {
      organism_val <- base::sub("^DEFINITION\\s+", "", def_line[1])
      organism_val <- base::sub(",.*$", "", organism_val)
    }
  }

  # ======================================================
  # Parse nucleotide sequence from ORIGIN section
  # ======================================================
  origin_idx <- base::grep("^ORIGIN", lines)
  end_idx <- base::grep("^//", lines)

  if (base::length(origin_idx) == 0 || base::length(end_idx) == 0) {
    base::stop(
      "This GenBank file does not contain a nucleotide sequence (no ORIGIN section).\n",
      "A full GenBank file is required.\n",
      "Please re-download in 'GenBank (full)' format from NCBI,\n",
      "or use download_genbank_gbff() to fetch the complete file.\n",
      "File: ", genbank_file
    )
  }

  # Find the end marker that comes after ORIGIN
  valid_ends <- end_idx[end_idx > origin_idx[1]]
  if (base::length(valid_ends) == 0) {
    base::stop("GenBank file has ORIGIN but no terminator (//) found after it.")
  }

  seq_lines <- lines[(origin_idx[1] + 1):(valid_ends[1] - 1)]
  # Remove line numbers and spaces, keep only nucleotide characters
  genome_seq <- base::toupper(base::gsub("[^a-zA-Z]", "",
                                          base::paste(seq_lines, collapse = "")))

  if (base::nchar(genome_seq) == 0) {
    base::stop(
      "The ORIGIN section in this GenBank file contains no nucleotide sequence.\n",
      "This is not a full GenBank file. Please provide a file with sequence data.\n",
      "File: ", genbank_file
    )
  }

  genome_length <- base::nchar(genome_seq)
  base::cat("  Genome sequence:", genome_length, "bp\n")

  # ======================================================
  # Parse FEATURES section - extract CDS features
  # ======================================================
  features_idx <- base::grep("^FEATURES", lines)
  if (base::length(features_idx) == 0) {
    base::warning("No FEATURES section found. Returning empty genbank_table.")
    genbank_table <- base::data.frame(
      locus_tag = base::character(0),
      start = base::integer(0),
      end = base::integer(0),
      strand = base::character(0),
      gene = base::character(0),
      product = base::character(0),
      nt_seq = base::character(0),
      stringsAsFactors = FALSE
    )
    return(base::list(
      genome_seq    = genome_seq,
      genbank_table = genbank_table,
      organism      = organism_val,
      accession     = accession_val,
      length        = genome_length
    ))
  }

  # Feature section: from FEATURES line to ORIGIN (or next section)
  feat_start <- features_idx[1] + 1
  feat_end <- base::ifelse(base::length(origin_idx) > 0,
                            origin_idx[1] - 1, base::length(lines))

  # Also check for other sections that might come between FEATURES and ORIGIN
  other_sections <- base::grep("^[A-Z]{2,}", lines)
  other_sections <- other_sections[other_sections > feat_start &
                                     other_sections < origin_idx[1]]
  if (base::length(other_sections) > 0) {
    feat_end <- base::min(other_sections) - 1
  }

  feat_lines <- lines[feat_start:feat_end]

  # Parse features into blocks (each feature starts at column 6 with no leading /)
  feature_blocks <- base::list()
  current_block <- base::character(0)
  current_type <- ""

  for (fl in feat_lines) {
    # A new feature line starts with 5 spaces followed by a non-space character
    # and is NOT a qualifier (doesn't start with /)
    if (base::grepl("^     [A-Za-z]", fl) && !base::grepl("^     /", fl)) {
      # Save previous block
      if (base::length(current_block) > 0 && base::nchar(current_type) > 0) {
        feature_blocks[[base::length(feature_blocks) + 1]] <- base::list(
          type = current_type,
          lines = current_block
        )
      }
      # Parse feature type and location from this line
      trimmed <- base::trimws(fl)
      parts <- base::strsplit(trimmed, "\\s+")[[1]]
      current_type <- parts[1]
      current_block <- fl
    } else {
      current_block <- base::c(current_block, fl)
    }
  }
  # Save last block
  if (base::length(current_block) > 0 && base::nchar(current_type) > 0) {
    feature_blocks[[base::length(feature_blocks) + 1]] <- base::list(
      type = current_type,
      lines = current_block
    )
  }

  # Filter for CDS features and parse qualifiers
  cds_blocks <- Filter(function(x) x$type == "CDS", feature_blocks)
  # Also collect gene features as fallback for locus_tag if CDS count is 0
  gene_blocks <- Filter(function(x) x$type == "gene", feature_blocks)

  # If no CDS features found, use gene features
  blocks_to_parse <- if (base::length(cds_blocks) > 0) cds_blocks else gene_blocks

  if (base::length(blocks_to_parse) == 0) {
    base::warning("No CDS or gene features found in GenBank file.")
    genbank_table <- base::data.frame(
      locus_tag = base::character(0),
      start = base::integer(0),
      end = base::integer(0),
      strand = base::character(0),
      gene = base::character(0),
      product = base::character(0),
      nt_seq = base::character(0),
      stringsAsFactors = FALSE
    )
    return(base::list(
      genome_seq    = genome_seq,
      genbank_table = genbank_table,
      organism      = organism_val,
      accession     = accession_val,
      length        = genome_length
    ))
  }

  base::cat("  Parsing", base::length(blocks_to_parse), "CDS/gene features...\n")

  # Parse each CDS/gene block
  parsed_features <- base::lapply(blocks_to_parse, function(block) {
    block_lines <- block$lines
    all_text <- base::paste(block_lines, collapse = "\n")

    # --- Parse location from first line ---
    first_line <- base::trimws(block_lines[1])
    loc_str <- base::sub("^\\S+\\s+", "", first_line)

    # Handle continuation lines for location (multi-line join)
    i <- 2
    while (i <= base::length(block_lines) &&
           !base::grepl("^\\s+/", block_lines[i])) {
      loc_str <- base::paste0(loc_str, base::trimws(block_lines[i]))
      i <- i + 1
    }

    # Determine strand
    strand <- "+"
    if (base::grepl("complement", loc_str)) {
      strand <- "-"
    }

    # Extract numeric positions
    nums <- base::as.integer(
      base::regmatches(loc_str, base::gregexpr("[0-9]+", loc_str))[[1]]
    )
    if (base::length(nums) < 2) return(NULL)

    feat_start <- base::min(nums)
    feat_end <- base::max(nums)

    # --- Parse qualifiers ---
    qualifiers <- base::list()
    current_qual <- NULL
    current_val <- ""

    for (j in base::seq_along(block_lines)) {
      line <- block_lines[j]
      if (base::grepl("^\\s+/([^=]+)=(.*)", line)) {
        # Save previous qualifier
        if (!base::is.null(current_qual)) {
          qualifiers[[current_qual]] <- base::gsub("^\"|\"$", "", current_val)
        }
        # New qualifier
        m <- base::regmatches(line, base::regexec("^\\s+/([^=]+)=(.*)", line))[[1]]
        current_qual <- m[2]
        current_val <- m[3]
      } else if (base::grepl("^\\s+/([^=]+)$", line)) {
        # Flag qualifier (no value)
        if (!base::is.null(current_qual)) {
          qualifiers[[current_qual]] <- base::gsub("^\"|\"$", "", current_val)
        }
        m <- base::regmatches(line, base::regexec("^\\s+/([^=]+)$", line))[[1]]
        current_qual <- m[2]
        current_val <- "TRUE"
      } else if (!base::is.null(current_qual) && base::grepl("^\\s{21}", line)) {
        # Continuation of qualifier value
        current_val <- base::paste0(current_val, base::trimws(line))
      }
    }
    # Save last qualifier
    if (!base::is.null(current_qual)) {
      qualifiers[[current_qual]] <- base::gsub("^\"|\"$", "", current_val)
    }

    locus_tag <- qualifiers[["locus_tag"]]
    gene_name <- qualifiers[["gene"]]
    product_name <- qualifiers[["product"]]

    if (base::is.null(locus_tag)) locus_tag <- NA_character_
    if (base::is.null(gene_name)) gene_name <- NA_character_
    if (base::is.null(product_name)) product_name <- NA_character_

    # Remove quotes
    locus_tag <- base::gsub("\"", "", locus_tag)
    gene_name <- base::gsub("\"", "", gene_name)
    product_name <- base::gsub("\"", "", product_name)

    # --- Extract nucleotide sequence for this feature ---
    nt_seq <- ""
    if (feat_start >= 1 && feat_end <= genome_length) {
      raw_seq <- base::substring(genome_seq, feat_start, feat_end)
      if (strand == "-") {
        # Reverse complement for minus strand
        nt_seq <- base::tryCatch({
          base::as.character(
            Biostrings::reverseComplement(Biostrings::DNAString(raw_seq))
          )
        }, error = function(e) raw_seq)
      } else {
        nt_seq <- raw_seq
      }
    }

    base::data.frame(
      locus_tag = locus_tag,
      start     = feat_start,
      end       = feat_end,
      strand    = strand,
      gene      = gene_name,
      product   = product_name,
      nt_seq    = nt_seq,
      stringsAsFactors = FALSE
    )
  })

  # Remove NULLs and combine
  parsed_features <- Filter(base::Negate(base::is.null), parsed_features)
  genbank_table <- base::do.call(base::rbind, parsed_features)
  base::rownames(genbank_table) <- NULL

  # Remove rows with no locus_tag (e.g., some hypothetical features)
  has_tag <- !base::is.na(genbank_table$locus_tag) &
    genbank_table$locus_tag != ""
  n_removed <- base::sum(!has_tag)
  if (n_removed > 0) {
    base::cat("  Removed", n_removed, "features without locus_tag\n")
    genbank_table <- genbank_table[has_tag, , drop = FALSE]
  }

  # Remove duplicates (same locus_tag from CDS)
  genbank_table <- genbank_table[!base::duplicated(genbank_table$locus_tag),
                                  , drop = FALSE]
  base::rownames(genbank_table) <- NULL

  base::cat("  Extracted", base::nrow(genbank_table), "annotated features\n")
  base::cat("  Organism:", organism_val, "\n")
  base::cat("  Accession:", accession_val, "\n")

  return(base::list(
    genome_seq    = genome_seq,
    genbank_table = genbank_table,
    organism      = organism_val,
    accession     = accession_val,
    length        = genome_length
  ))
}
