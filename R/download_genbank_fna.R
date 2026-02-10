#' Download GenBank .fna File
#'
#' Downloads a GenBank .fna file for a given accession number from NCBI.
#' Creates a directory to store the downloaded file and retrieves metadata such as organism name and release date.
#'
#' @param dir A character string specifying the directory to save the .fna file (default: current working directory).
#' @param accession A character string specifying the GenBank accession number (e.g., "GCF_030376745.1").
#' @return A list containing the file path of the downloaded .fna file, organism name, release date, and sequence directory.
#' @examples
#' \dontrun{
#' metadata <- download_genbank_fna(dir = getwd(), accession = "GCF_030376745.1")
#' print(metadata)
#' }

download_genbank_fna <- function(dir = getwd(), accession) {
  low_dir <- "seqs_srcdir"
  seqs_dir <- base::file.path(dir, low_dir)
  base::dir.create(seqs_dir, showWarnings = FALSE, recursive = TRUE)
  base::cat("Ensured seqs_srcdir exists:", seqs_dir, "\n")

  # Validate accession format (GCF_xxxxxx.x or GCA_xxxxxx.x)
  if (!base::grepl("^(GCF|GCA)_[0-9]+\\.[0-9]$", accession)) {
    base::stop("Invalid GenBank accession format. Expected: GCF_xxxxxx.x or GCA_xxxxxx.x")
  }

  # Search NCBI assembly database
  base::cat("Searching for accession:", accession, "\n")
  search_result <- rentrez::entrez_search(
    db = "assembly",
    term = accession,
    retmode = "xml",
    retmax = 1,
    use_history = TRUE
  )

  # Debug: Log search results
  base::cat("Search result count:", base::length(search_result$ids), "\n")
  if (base::length(search_result$ids) > 0) {
    base::cat("Found IDs:", search_result$ids, "\n")
  }

  if (base::length(search_result$ids) == 0) {
    # Fallback query
    base::cat("Retrying with alternative query...\n")
    search_result <- rentrez::entrez_search(
      db = "assembly",
      term = base::paste0(accession, " AND (complete_genome[filter] OR scaffold[filter] OR contig[filter])"),
      retmode = "xml",
      retmax = 1,
      use_history = TRUE
    )
    base::cat("Retry result count:", base::length(search_result$ids), "\n")
    if (base::length(search_result$ids) == 0) {
      base::stop("No results found for accession:", accession, "\nQueries attempted:", accession, " and ", base::paste0(accession, " AND (complete_genome[filter] OR scaffold[filter] OR contig[filter])"))
    }
  }

  # Fetch metadata
  summary <- rentrez::entrez_summary(db = "assembly", id = search_result$ids[1])
  organism_name <- base::ifelse(
    base::is.null(summary$speciesname),
    "Unknown",
    summary$speciesname
  )
  release_date <- base::ifelse(
    base::is.null(summary$lastupdatedate),
    base::format(base::Sys.Date(), "%Y-%m-%d"),
    summary$lastupdatedate
  )

  # Fetch .fna file
  recs <- rentrez::entrez_link(dbfrom = "assembly", id = search_result$ids[1], db = "nuccore")
  nuccore_ids <- NULL
  if (base::grepl("^GCF_", accession)) {
    # Prefer RefSeq for GCF_
    nuccore_ids <- recs$links$assembly_nuccore_refseq
    if (base::is.null(nuccore_ids) || base::length(nuccore_ids) == 0) {
      base::cat("No RefSeq sequences found, falling back to GenBank/INSDC\n")
      nuccore_ids <- recs$links$assembly_nuccore_insdc
    }
  } else {
    # Prefer GenBank/INSDC for GCA_
    nuccore_ids <- recs$links$assembly_nuccore_insdc
    if (base::is.null(nuccore_ids) || base::length(nuccore_ids) == 0) {
      base::cat("No GenBank/INSDC sequences found, falling back to RefSeq\n")
      nuccore_ids <- recs$links$assembly_nuccore_refseq
    }
  }

  if (base::is.null(nuccore_ids) || base::length(nuccore_ids) == 0) {
    base::stop("No nucleotide sequences found for accession:", accession)
  }

  base::cat("Found nuccore IDs:", nuccore_ids, "\n")
  fna_file <- base::file.path(seqs_dir, base::paste0(accession, ".fna"))
  base::cat("Downloading .fna to:", fna_file, "\n")
  base::sink(fna_file)
  for (id in nuccore_ids) {
    fna <- rentrez::entrez_fetch(db = "nuccore", id = id, rettype = "fasta")
    base::cat(fna)
  }
  base::sink()

  if (!base::file.exists(fna_file)) {
    base::stop("Failed to download .fna file for accession:", accession)
  }

  # Return metadata for seed file
  base::list(
    fna_file = fna_file,
    organism_name = organism_name,
    release_date = release_date,
    seqs_dir = seqs_dir
  )
}
