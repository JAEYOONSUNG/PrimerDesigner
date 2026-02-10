#' Generate gRNA for Specific Locus Tags
#'
#' Generates gRNA lists for one or more locus tags using the sequences from genbank_table and optionally saves the result to an Excel file.
#'
#' @param locus_tags A character vector of locus tags (e.g., "QT234_RS00005").
#' @param genbank_table A data frame containing GenBank table data with columns including locus_tag and nt_seq.
#' @param genbank_accession A character string specifying the GenBank accession number (e.g., "GCF_030376745.1").
#' @param nuclease A character string specifying the nuclease to use ("GeoCas9", "FnCas12a", or "FisCasI_B"; default: "GeoCas9").
#' @param bowtie_index A character string specifying the path to the Bowtie index.
#' @param genome A BSgenome object specifying the reference genome (default: NULL, auto-loaded).
#' @param position_range A numeric vector of length 2 specifying the position range (as percentage) for gRNA filtering (default: c(10, 90)).
#' @param tm_range A numeric vector of length 2 specifying the Tm range for gRNA filtering (default: c(50, 70)).
#' @param strand A character string specifying the strand to filter ("both", "5", or "3"; default: "both").
#' @param dir A character string specifying the working directory (default: current working directory).
#' @param output_file A character string specifying the output Excel file path (default: NULL, no file saved).
#' @param name_prefix A character string to prepend to gRNA names (default: NULL, auto-detected from nuclease).
#' @return A data frame containing the gRNA list for the specified locus tags with gRNA_name and gRNA_rank columns.
#' @examples
#' \dontrun{
#' genbank_table <- data.frame(locus_tag = "QT234_RS00005", nt_seq = "atcg")
#' gRNA_list <- generate_gRNA_for_locus_tags(
#'   locus_tags = "QT234_RS00005",
#'   genbank_table = genbank_table,
#'   genbank_accession = "GCF_030376745.1",
#'   nuclease = "GeoCas9",
#'   bowtie_index = "path/to/bowtie_index",
#'   genome = NULL,
#'   position_range = c(10, 90),
#'   tm_range = c(50, 70),
#'   output_file = "gRNA_list.xlsx"
#' )
#' print(gRNA_list)
#' }
#' @export
generate_gRNA_for_locus_tags <- function(
    locus_tags,
    genbank_table,
    genbank_accession,
    nuclease = "GeoCas9",
    bowtie_index,
    genome = NULL,
    position_range = c(10, 90),
    tm_range = c(50, 70),
    strand = "both",
    dir = base::getwd(),
    output_file = NULL,
    name_prefix = NULL
) {
  # Step 1: Validate genbank_accession
  if (base::missing(genbank_accession)) {
    base::stop("Error: 'genbank_accession' is missing. Please provide a valid GenBank accession number.")
  }

  # Step 2: Validate locus_tags
  if (!base::all(locus_tags %in% genbank_table$locus_tag)) {
    missing_tags <- locus_tags[!(locus_tags %in% genbank_table$locus_tag)]
    base::stop("Error: The following locus_tags were not found in genbank_table: ", base::paste(missing_tags, collapse = ", "))
  }

  # Step 3: Set default genome if NULL
  if (base::is.null(genome)) {
    clean_accession <- base::gsub("[^a-zA-Z0-9]", "", genbank_accession)
    pkg_name <- base::paste0("BSgenome.", clean_accession)
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      base::cat("Building BSgenome package for accession:", genbank_accession, "\n")
      build_bsgenome_from_accession(genbank_accession = genbank_accession, dir = dir)
    }
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      base::stop("BSgenome package not found:", pkg_name, "\nPlease ensure it is installed.")
    }
    genome <- base::get(base::gsub(".*\\.", "", pkg_name), envir = base::asNamespace(pkg_name))
    base::cat("Using default genome package:", pkg_name, "\n")
  }

  # Step 4: Filter genbank_table for the given locus_tags and extract nt_seq + gene info
  target_info <- genbank_table %>%
    dplyr::filter(locus_tag %in% locus_tags)
  target_seqs <- target_info %>% dplyr::pull(nt_seq)
  # Extract gene names if available (for naming)
  target_genes <- if ("gene" %in% base::colnames(target_info)) {
    target_info %>% dplyr::pull(gene)
  } else {
    base::rep(NA_character_, base::length(locus_tags))
  }

  # Step 5: Clean nt_seq data
  target_seqs <- base::sapply(base::seq_along(target_seqs), function(i) {
    seq <- target_seqs[i]
    locus_tag <- locus_tags[i]
    # Check for invalid characters and log them
    if (!is_valid_dna(seq)) {
      log_invalid_sequence(seq, locus_tag)
    }
    # Replace invalid characters with 'N'
    seq <- base::gsub("[^ATGCNatgcnRYMKSWHBVDrymkswbdvh-]", "N", seq)
    return(seq)
  })

  # Step 6: Generate gRNA for each locus_tag
  gRNA_list <- base::lapply(base::seq_along(locus_tags), function(i) {
    locus_tag <- locus_tags[i]
    target_seq <- target_seqs[i]

    base::cat("Generating gRNA list for locus_tag:", locus_tag, "\n")
    base::tryCatch({
      gRNA_df <- nuclease_with_parameter(
        target_sequence = target_seq,
        strand = strand,
        bowtie_index = bowtie_index,
        genome = genome,
        nuclease = nuclease,
        position_range = position_range,
        tm_range = tm_range
      )
      gRNA_df$locus_tag <- locus_tag
      # Add gene name for downstream naming
      gene_name <- target_genes[i]
      gRNA_df$gene <- base::ifelse(base::is.na(gene_name) || gene_name == "", NA_character_, gene_name)
      return(gRNA_df)
    }, error = function(e) {
      base::cat("Error for locus_tag", locus_tag, ":", base::conditionMessage(e), "\n")
      return(NULL)
    })
  })

  # Step 7: Remove NULL entries and combine
  gRNA_list <- base::Filter(base::Negate(base::is.null), gRNA_list)
  if (base::length(gRNA_list) == 0) {
    base::stop("No gRNAs generated for any locus_tag")
  }
  gRNA_list <- base::do.call(rbind, gRNA_list)

  # Step 8: Generate systematic gRNA names and ranks
  base::cat("Step 8: Generating systematic gRNA names\n")
  gRNA_list <- generate_gRNA_names(
    gRNA_df = gRNA_list,
    prefix = name_prefix,
    nuclease = nuclease,
    rank_by = "composite_score"
  )
  base::cat("gRNA naming complete. Sample names:", utils::head(gRNA_list$gRNA_name, 3), "\n")

  # Step 9: Save to Excel if output_file is specified
  if (!base::is.null(output_file)) {
    base::cat("Step 9: Saving gRNA list to Excel\n")
    base::tryCatch({
      # Ensure the directory exists
      output_dir <- base::dirname(output_file)
      if (!base::dir.exists(output_dir)) {
        base::dir.create(output_dir, recursive = TRUE)
        base::cat("Created output directory:", output_dir, "\n")
      }
      # Save the file using xlsx::write.xlsx
      xlsx::write.xlsx(gRNA_list, file = output_file, row.names = FALSE)
      base::cat("Successfully saved gRNA list to:", output_file, "\n")
    }, error = function(e) {
      base::cat("Error saving gRNA list to Excel:", base::conditionMessage(e), "\n")
      base::stop("Failed to save gRNA list to ", output_file, ". Please check file path, permissions, or xlsx package installation.")
    })
  }

  return(gRNA_list)
}
