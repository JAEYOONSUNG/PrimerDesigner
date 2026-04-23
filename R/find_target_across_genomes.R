#' Find Target Across Multiple GenBank Files
#'
#' Scans one or more directories of GenBank files and returns a
#' \code{target_table} suitable for \code{design_shared_grna_and_deletion()}.
#' The user supplies a single query (locus_tag, gene name, or product keyword)
#' and the function discovers matching CDSes per genome.
#'
#' When a genome has more than one match (e.g. a gene family) the function
#' prompts the user to pick which locus_tag to use, unless \code{interactive}
#' is \code{FALSE}, in which case a data frame of candidates is returned so
#' the caller can resolve ambiguity programmatically.
#'
#' @param genbank_dirs Character vector of directory paths holding genome
#'   files. Accepts GenBank (\code{.gbk}, \code{.gb}, \code{.gbff},
#'   \code{.genbank}) and SnapGene binary (\code{.dna}); \code{.dna} files are
#'   converted to GenBank on the fly via the SnapGene CLI. Each file is
#'   treated as one genome; \code{genome_id} defaults to the file basename
#'   (no extension).
#' @param query Character. Locus tag, gene symbol, or product keyword to search.
#' @param query_type One of \code{"auto"} (default), \code{"locus_tag"},
#'   \code{"gene"}, or \code{"product"}. In \code{"auto"} mode the query is
#'   matched against all three fields and the field with the strongest hit is
#'   used.
#' @param ignore_case Logical. Case-insensitive match for \code{gene} and
#'   \code{product} searches. Default \code{TRUE}. \code{locus_tag} is always
#'   exact-match.
#' @param regex Logical. If \code{TRUE}, treat \code{query} as a regex for
#'   \code{gene}/\code{product} searches. Default \code{FALSE} (fixed substring).
#' @param interactive Logical. If \code{TRUE} (default when R is interactive)
#'   and a genome has more than one matching CDS, prompt the user to choose
#'   which one to use. If \code{FALSE}, all candidates are returned and
#'   \code{one_per_genome} is ignored.
#' @param one_per_genome Logical. If \code{TRUE}, require exactly one locus_tag
#'   per genome after resolution. If any genome ends up with zero matches the
#'   call stops. Default \code{TRUE}.
#' @param genome_id_from One of \code{"file"} (default, strip extension) or
#'   \code{"accession"} (use LOCUS accession from the file header).
#' @return A data frame with columns \code{genome_id}, \code{genbank_file},
#'   \code{locus_tag}, \code{gene}, \code{product}, ready to pass to
#'   \code{design_shared_grna_and_deletion()}. If \code{interactive = FALSE}
#'   and matches are ambiguous, additional rows may share a \code{genome_id};
#'   a \code{match_field} column then records which field matched.
#' @examples
#' \dontrun{
#' # One folder, product keyword
#' tt <- find_target_across_genomes(
#'   genbank_dirs = "/data/geobacillus_gbk",
#'   query = "jetD", query_type = "gene")
#'
#' # Multiple folders, product keyword with interactive disambiguation
#' tt <- find_target_across_genomes(
#'   genbank_dirs = c("/data/typeA", "/data/typeB"),
#'   query = "DUF2625 domain-containing",
#'   query_type = "product")
#'
#' # Non-interactive; caller filters manually
#' cand <- find_target_across_genomes(
#'   genbank_dirs = "/data/gbk",
#'   query = "recA", interactive = FALSE, one_per_genome = FALSE)
#' }
#' @export
find_target_across_genomes <- function(
    genbank_dirs,
    query,
    query_type = c("auto", "locus_tag", "gene", "product"),
    ignore_case = TRUE,
    regex = FALSE,
    interactive = base::interactive(),
    one_per_genome = TRUE,
    genome_id_from = c("file", "accession"),
    kill_snapgene = FALSE,
    sequence_fallback = TRUE,
    sequence_fallback_min_identity = 0.35,
    sequence_fallback_min_qcov = 0.60,
    sequence_fallback_evalue = 1e-5
) {
  query_type <- base::match.arg(query_type)
  genome_id_from <- base::match.arg(genome_id_from)

  if (base::missing(query) || !base::is.character(query) || base::length(query) != 1L ||
      base::nchar(query) == 0) {
    base::stop("query must be a non-empty single string.")
  }

  # --- Collect candidate files across all dirs ---
  # Accepts GenBank (.gb/.gbk/.gbff/.genbank) and SnapGene binary (.dna).
  # .dna files are converted to GenBank on the fly via SnapGene CLI.
  exts <- c("gbk", "gb", "gbff", "genbank", "dna")
  gb_files <- base::unlist(base::lapply(genbank_dirs, function(d) {
    if (!base::dir.exists(d)) base::stop("Directory not found: ", d)
    base::list.files(d, pattern = base::paste0("\\.(",
                                                base::paste(exts, collapse = "|"),
                                                ")$"),
                     full.names = TRUE, ignore.case = TRUE, recursive = FALSE)
  }))
  if (base::length(gb_files) == 0) {
    base::stop("No genome files (", base::paste(exts, collapse = ", "),
               ") found in: ", base::paste(genbank_dirs, collapse = ", "))
  }
  base::cat("Scanning", base::length(gb_files), "genome file(s)...\n")

  # --- Per-file search (parse once per file; CDS-only feature extraction) ---
  rows <- base::list()
  parsed_cache <- base::list()   # cache parsed records to reuse on fallback pass
  file_map <- base::list()        # genome_id -> genbank_file
  missed <- base::list()          # genome_id -> list(parsed, gf) with no hit
  for (gf in gb_files) {
    base::cat("  ", base::basename(gf), "\n")
    parsed <- base::tryCatch(.find_target_read_any(gf, kill_snapgene),
                             error = function(e) { base::warning(e$message); NULL })
    if (base::is.null(parsed)) next
    gt <- parsed$genbank_table
    if (base::is.null(gt) || base::nrow(gt) == 0) next

    gid <- base::switch(genome_id_from,
      file = base::sub("\\.[^.]+$", "", base::basename(gf)),
      accession = if (base::nzchar(parsed$accession)) parsed$accession
                   else base::sub("\\.[^.]+$", "", base::basename(gf)))
    parsed_cache[[gid]] <- parsed
    file_map[[gid]] <- gf

    hits <- .find_target_hits(gt, query, query_type, ignore_case, regex)
    if (base::nrow(hits) == 0) {
      missed[[gid]] <- base::list(parsed = parsed, gf = gf)
      next
    }
    hits$genome_id    <- gid
    hits$genbank_file <- gf
    rows[[base::length(rows) + 1]] <- hits
  }

  # --- Sequence-based fallback for genomes with no annotation match ----
  # Seeds are the CDS records of every annotation hit we already have.
  # For each missing genome, run MMseqs2 (or Biostrings) to locate a
  # homolog and append a row with `match_field = "sequence_homolog"`.
  seq_rows <- base::list()
  if (base::isTRUE(sequence_fallback) && base::length(missed) > 0L &&
      base::length(rows) > 0L) {
    all_annotation_hits <- base::do.call(base::rbind, rows)
    seeds <- base::list()
    for (k in base::seq_len(base::nrow(all_annotation_hits))) {
      r <- all_annotation_hits[k, , drop = FALSE]
      nt <- if ("nt_seq" %in% base::colnames(r)) base::as.character(r$nt_seq[1]) else NA_character_
      aa <- if ("protein" %in% base::colnames(r)) base::as.character(r$protein[1]) else NA_character_
      if (base::is.na(aa) && !base::is.na(nt))
        aa <- .shared_translate_cds(nt)
      if (!base::is.na(aa) && base::nchar(aa) > 10L) {
        seed_id <- base::paste0(r$genome_id[1], "__", r$locus_tag[1])
        seeds[[seed_id]] <- base::list(nt = nt, aa = aa)
      }
    }
    if (base::length(seeds) > 0L) {
      base::cat("\nSequence-based homolog fallback for ",
                 base::length(missed), " genome(s) without annotation match:\n",
                 sep = "")
      for (gid in base::names(missed)) {
        parsed <- missed[[gid]]$parsed
        gf <- missed[[gid]]$gf
        row <- base::tryCatch(
          .shared_find_homolog_by_sequence(
            seeds = seeds, genome_id = gid,
            genome_record = parsed, genbank_table = parsed$genbank_table,
            genbank_file = gf,
            evalue = sequence_fallback_evalue,
            min_aa_identity = sequence_fallback_min_identity,
            min_q_cov = sequence_fallback_min_qcov, verbose = TRUE),
          error = function(e) {
            base::warning(sprintf("homolog fallback failed for %s: %s",
                                   gid, base::conditionMessage(e)), call. = FALSE)
            NULL
          })
        if (!base::is.null(row)) seq_rows[[base::length(seq_rows) + 1L]] <- row
      }
    }
  }

  if (base::length(rows) == 0 && base::length(seq_rows) == 0) {
    base::stop("No GenBank file contained a match for query '", query, "'.")
  }
  all_hits <- if (base::length(rows) > 0)
    base::do.call(base::rbind, rows) else base::data.frame()
  if (base::length(seq_rows) > 0) {
    seq_df <- dplyr::bind_rows(seq_rows)
    # Reconcile columns: annotation rows have nt_seq/protein etc.; homolog
    # rows carry the extra diagnostic fields. Use dplyr::bind_rows to align.
    all_hits <- dplyr::bind_rows(all_hits, seq_df)
  }
  base::rownames(all_hits) <- NULL

  # --- Resolve ambiguity per genome ---
  desired_cols <- base::c("genome_id", "genbank_file", "locus_tag",
                           "gene", "product", "match_field",
                           "sequence_hit", "pident", "qcov", "bits",
                           "homology_seed", "pseudo_locus")
  return_cols <- base::intersect(desired_cols, base::colnames(all_hits))
  if (!interactive) {
    return(all_hits[, return_cols, drop = FALSE])
  }

  resolved <- base::do.call(base::rbind, base::lapply(
    base::split(all_hits, all_hits$genome_id), function(df) {
      if (base::nrow(df) == 1L) return(df)
      .prompt_pick_locus(df)
    }))
  base::rownames(resolved) <- NULL

  if (one_per_genome) {
    dup_gids <- base::names(base::which(base::table(resolved$genome_id) != 1L))
    if (base::length(dup_gids) > 0) {
      base::stop("one_per_genome = TRUE but ambiguous after resolution for: ",
                 base::paste(dup_gids, collapse = ", "))
    }
  }

  # --- Summary print ---
  base::cat("\nSelected target across", base::nrow(resolved), "genome(s):\n")
  base::print(resolved[, base::intersect(base::c("genome_id", "locus_tag",
                                                    "gene", "product",
                                                    "match_field"),
                                            base::colnames(resolved))],
              row.names = FALSE)

  final_cols <- base::intersect(return_cols, base::colnames(resolved))
  resolved[, final_cols, drop = FALSE]
}

# ---- Internal helpers -------------------------------------------------------

# Read a genome file; transparently handle .dna via SnapGene CLI conversion.
.find_target_read_any <- function(path, kill_snapgene = FALSE) {
  ext <- base::tolower(tools::file_ext(path))
  if (ext == "dna") {
    gb_lines <- .get_genbank_from_dna(path, kill_snapgene = kill_snapgene)
    tmp_gbk <- base::tempfile(fileext = ".gbk")
    base::on.exit(if (base::file.exists(tmp_gbk)) base::unlink(tmp_gbk),
                  add = TRUE)
    base::writeLines(gb_lines, tmp_gbk)
    read_genome_genbank(tmp_gbk)
  } else {
    read_genome_genbank(path)
  }
}

.find_target_hits <- function(gt, query, query_type, ignore_case, regex) {
  if (!"locus_tag" %in% base::names(gt)) gt$locus_tag <- NA_character_
  if (!"gene" %in% base::names(gt))      gt$gene      <- NA_character_
  if (!"product" %in% base::names(gt))   gt$product   <- NA_character_

  match_field <- function(vals, qry, fixed) {
    if (fixed) {
      # grepl() rejects fixed + ignore.case; apply case-folding manually.
      if (ignore_case) {
        base::grepl(base::tolower(qry), base::tolower(vals), fixed = TRUE)
      } else {
        base::grepl(qry, vals, fixed = TRUE)
      }
    } else {
      base::grepl(qry, vals, perl = TRUE, ignore.case = ignore_case)
    }
  }

  fixed <- !regex

  idx_lt <- base::logical(base::nrow(gt))
  idx_gn <- base::logical(base::nrow(gt))
  idx_pr <- base::logical(base::nrow(gt))

  # locus_tag: always exact (case-sensitive) match when that field is queried
  if (query_type %in% c("auto", "locus_tag")) {
    idx_lt <- !base::is.na(gt$locus_tag) & gt$locus_tag == query
  }
  if (query_type %in% c("auto", "gene")) {
    idx_gn <- match_field(base::ifelse(base::is.na(gt$gene), "", gt$gene),
                          query, fixed)
  }
  if (query_type %in% c("auto", "product")) {
    idx_pr <- match_field(base::ifelse(base::is.na(gt$product), "", gt$product),
                          query, fixed)
  }

  # "auto": prefer the most specific field (locus_tag > gene > product)
  if (query_type == "auto") {
    if (base::any(idx_lt))      { keep <- idx_lt; mf <- "locus_tag" }
    else if (base::any(idx_gn)) { keep <- idx_gn; mf <- "gene"      }
    else                        { keep <- idx_pr; mf <- "product"   }
  } else {
    keep <- base::switch(query_type,
      locus_tag = idx_lt, gene = idx_gn, product = idx_pr)
    mf <- query_type
  }

  out <- gt[keep, base::intersect(c("locus_tag", "start", "end", "strand",
                                    "gene", "product", "nt_seq", "protein",
                                    "contig"),
                                   base::names(gt)),
            drop = FALSE]
  if (base::nrow(out) > 0) out$match_field <- mf
  out
}

.prompt_pick_locus <- function(df) {
  gid <- df$genome_id[1]
  base::cat(base::sprintf("\n%d matches in '%s':\n", base::nrow(df), gid))
  display <- base::data.frame(
    n = base::seq_len(base::nrow(df)),
    locus_tag = df$locus_tag,
    gene = df$gene,
    product = df$product,
    pos = base::ifelse(base::is.na(df$start), "",
                        base::sprintf("%s:%d-%d", df$strand, df$start, df$end)),
    stringsAsFactors = FALSE
  )
  base::print(display, row.names = FALSE)
  repeat {
    ans <- base::readline(base::sprintf(
      "  Pick 1-%d (or 0 to skip, a = accept all, q = abort): ",
      base::nrow(df)))
    ans <- base::trimws(ans)
    if (base::identical(ans, "q")) base::stop("Aborted by user.")
    if (base::identical(ans, "a")) return(df)
    if (base::identical(ans, "0") || base::identical(ans, "")) {
      base::message("  skipped ", gid)
      return(df[0, , drop = FALSE])
    }
    sel <- base::suppressWarnings(base::as.integer(ans))
    if (!base::is.na(sel) && sel >= 1L && sel <= base::nrow(df)) {
      return(df[sel, , drop = FALSE])
    }
    base::cat("  invalid input\n")
  }
}
