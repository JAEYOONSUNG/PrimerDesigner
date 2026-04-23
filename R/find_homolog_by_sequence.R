# Sequence-based homolog search for find_target_across_genomes
#
# When annotation-based search misses a genome (because the ortholog is
# annotated under a different product / hypothetical protein name), locate
# the homolog by aligning the existing annotation hits' CDS sequences
# against the missing genome and map the hit back to a CDS row. Primary
# engine is MMseqs2 (fast translated AA->NT search); pure-R Biostrings
# CDS-level pairwiseAlignment is the fallback.

# Bacterial genetic code table 11, used when GenBank /transl_table is absent.
.shared_bacterial_codon_table <- function() {
  c(
    TTT = "F", TTC = "F", TTA = "L", TTG = "L",
    TCT = "S", TCC = "S", TCA = "S", TCG = "S",
    TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
    TGT = "C", TGC = "C", TGA = "*", TGG = "W",
    CTT = "L", CTC = "L", CTA = "L", CTG = "L",
    CCT = "P", CCC = "P", CCA = "P", CCG = "P",
    CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
    CGT = "R", CGC = "R", CGA = "R", CGG = "R",
    ATT = "I", ATC = "I", ATA = "I", ATG = "M",
    ACT = "T", ACC = "T", ACA = "T", ACG = "T",
    AAT = "N", AAC = "N", AAA = "K", AAG = "K",
    AGT = "S", AGC = "S", AGA = "R", AGG = "R",
    GTT = "V", GTC = "V", GTA = "V", GTG = "V",
    GCT = "A", GCC = "A", GCA = "A", GCG = "A",
    GAT = "D", GAC = "D", GAA = "E", GAG = "E",
    GGT = "G", GGC = "G", GGA = "G", GGG = "G"
  )
}

#' @keywords internal
.shared_translate_cds <- function(nt_seq) {
  if (base::is.na(nt_seq) || base::nchar(nt_seq) < 3L) return(NA_character_)
  if (base::requireNamespace("Biostrings", quietly = TRUE)) {
    aa <- base::tryCatch({
      dna <- Biostrings::DNAString(base::toupper(nt_seq))
      base::as.character(Biostrings::translate(dna, if.fuzzy.codon = "X",
                                                no.init.codon = TRUE))
    }, error = function(e) NA_character_)
    if (!base::is.na(aa)) return(base::sub("\\*+$", "", aa))
  }
  # Pure-R fallback
  s <- base::toupper(base::as.character(nt_seq))
  n <- base::nchar(s)
  n3 <- n %/% 3L
  if (n3 == 0L) return(NA_character_)
  codons <- base::substring(s, base::seq.int(1L, n3 * 3L, by = 3L),
                             base::seq.int(3L, n3 * 3L, by = 3L))
  tbl <- .shared_bacterial_codon_table()
  aa <- tbl[codons]
  aa[base::is.na(aa)] <- "X"
  base::sub("\\*+$", "", base::paste(aa, collapse = ""))
}

#' @keywords internal
.shared_mmseqs_bin <- function() {
  # Respects the path cached on load (zzz.R :: .shared_cache_mmseqs_path)
  # and the public `find_mmseqs_binary()` helper; both write into
  # `options("PrimerDesigner.mmseqs_path")`.
  bin <- base::getOption("PrimerDesigner.mmseqs_path", NA_character_)
  if (!base::is.na(bin) && base::nzchar(bin) && base::file.exists(bin))
    return(bin)
  # Cache miss (e.g. devtools::load_all without .onLoad) -- re-scan.
  if (base::exists("find_mmseqs_binary", mode = "function")) {
    return(find_mmseqs_binary(refresh = TRUE))
  }
  NA_character_
}

.shared_mmseqs_available <- function() {
  !base::is.na(.shared_mmseqs_bin())
}

# Write a FASTA file from a named character vector of sequences.
.shared_write_fasta <- function(seqs, path) {
  nm <- base::names(seqs)
  if (base::is.null(nm)) nm <- base::paste0("seq_", base::seq_along(seqs))
  lines <- base::character(2L * base::length(seqs))
  lines[base::seq.int(1L, 2L * base::length(seqs), by = 2L)] <-
    base::paste0(">", nm)
  lines[base::seq.int(2L, 2L * base::length(seqs), by = 2L)] <- seqs
  base::writeLines(lines, path)
  path
}

# ---- MMseqs2 primary -----------------------------------------------------
#
# Runs `mmseqs easy-search` with AA seeds against a genome nucleotide FASTA
# (translated target search, bacterial table 11). Returns a data frame of
# hits with contig / start / end / strand / identity / qcov / bitscore.
.shared_mmseqs_homolog_search <- function(seeds_aa, genome_fna, tmp_dir,
                                           sensitivity = 7.0, evalue = 1e-5,
                                           max_seqs = 20L) {
  if (base::length(seeds_aa) == 0L) return(base::data.frame())
  seed_faa <- base::file.path(tmp_dir, "seeds.faa")
  .shared_write_fasta(seeds_aa, seed_faa)
  out_tsv <- base::file.path(tmp_dir, "hits.m8")
  mm_tmp <- base::file.path(tmp_dir, "mm_tmp"); base::dir.create(mm_tmp, showWarnings = FALSE)
  fmt_cols <- "query,target,pident,alnlen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen"
  args <- base::c("easy-search", seed_faa, genome_fna, out_tsv, mm_tmp,
                   "--translation-table", "11",
                   "-s", base::as.character(sensitivity),
                   "-e", base::as.character(evalue),
                   "--max-seqs", base::as.character(max_seqs),
                   "--format-output", fmt_cols,
                   "--search-type", "2",   # 2 = translated
                   "--min-seq-id", "0.0")
  mm_bin <- .shared_mmseqs_bin()
  status <- base::suppressWarnings(
    base::system2(mm_bin, args, stdout = FALSE, stderr = FALSE))
  if (status != 0 || !base::file.exists(out_tsv) ||
      base::file.info(out_tsv)$size == 0L) {
    return(base::data.frame())
  }
  tab <- base::tryCatch(
    utils::read.table(out_tsv, sep = "\t", header = FALSE,
                       stringsAsFactors = FALSE,
                       col.names = base::strsplit(fmt_cols, ",")[[1]]),
    error = function(e) base::data.frame())
  if (base::nrow(tab) == 0L) return(base::data.frame())
  tab$pident <- base::as.numeric(tab$pident)
  tab$bits   <- base::as.numeric(tab$bits)
  tab$qcov   <- (base::pmin(tab$qend, tab$qlen) -
                  base::pmax(tab$qstart, 1L) + 1L) / tab$qlen
  tab$strand <- base::ifelse(tab$tstart <= tab$tend, "+", "-")
  tab$contig <- tab$target
  tab$start  <- base::pmin(tab$tstart, tab$tend)
  tab$end    <- base::pmax(tab$tstart, tab$tend)
  tab[, base::c("query", "contig", "start", "end", "strand",
                  "pident", "qcov", "bits", "evalue")]
}

# ---- Biostrings fallback -------------------------------------------------
#
# CDS-level local alignment: we compare each seed AA against every CDS AA
# translation, rank by alignment score, then keep hits meeting identity /
# coverage thresholds. Avoids whole-genome DP.
.shared_biostrings_homolog_search <- function(seeds_aa, genbank_table,
                                                genome_seq = NULL,
                                                min_aa_identity = 0.35,
                                                min_q_cov = 0.60,
                                                top_n_rank = 30L) {
  if (!base::requireNamespace("Biostrings", quietly = TRUE)) {
    base::warning("Biostrings not installed; skipping sequence fallback.",
                  call. = FALSE)
    return(base::data.frame())
  }
  if (!("nt_seq" %in% base::colnames(genbank_table)) ||
      !("start" %in% base::colnames(genbank_table))) {
    return(base::data.frame())
  }
  cds_nt <- base::as.character(genbank_table$nt_seq)
  cds_aa <- base::vapply(cds_nt, .shared_translate_cds, character(1))
  keep <- !base::is.na(cds_aa) & base::nchar(cds_aa) > 10L
  if (!base::any(keep)) return(base::data.frame())
  cds_aa_vec <- cds_aa[keep]
  cds_idx <- base::which(keep)
  cds_aa_set <- Biostrings::AAStringSet(cds_aa_vec)

  hits <- base::list()
  for (s_name in base::names(seeds_aa)) {
    seed <- Biostrings::AAString(seeds_aa[[s_name]])
    scores <- base::suppressWarnings(
      Biostrings::pairwiseAlignment(cds_aa_set, seed,
                                     type = "local", scoreOnly = TRUE))
    top <- base::order(scores, decreasing = TRUE)[base::seq_len(
      base::min(top_n_rank, base::length(scores)))]
    for (j in top) {
      aln <- base::suppressWarnings(
        Biostrings::pairwiseAlignment(cds_aa_set[j], seed, type = "local"))
      pid <- base::as.numeric(Biostrings::pid(aln, type = "PID3")) / 100
      patt_str <- base::as.character(Biostrings::pattern(aln))
      qcov <- base::nchar(patt_str) / base::nchar(base::as.character(seed))
      score_val <- base::as.numeric(Biostrings::score(aln))
      if (!base::is.na(pid) && pid >= min_aa_identity && qcov >= min_q_cov) {
        g_row <- genbank_table[cds_idx[j], , drop = FALSE]
        hits[[base::length(hits) + 1L]] <- base::data.frame(
          query  = s_name,
          contig = base::as.character(
            if ("contig" %in% base::colnames(g_row)) g_row$contig[1] else NA),
          start  = base::as.integer(g_row$start[1]),
          end    = base::as.integer(g_row$end[1]),
          strand = base::as.character(
            if ("strand" %in% base::colnames(g_row)) g_row$strand[1] else "+"),
          pident = pid,
          qcov   = qcov,
          bits   = score_val,
          evalue = NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (base::length(hits) == 0L) return(base::data.frame())
  dplyr::bind_rows(hits)
}

# ---- Map a hit back to a CDS row or synthesize a pseudo locus_tag --------
.shared_map_hit_to_cds <- function(hit, genbank_table, genome_id) {
  if ("contig" %in% base::colnames(genbank_table) &&
      !base::is.null(hit$contig) && !base::is.na(hit$contig) &&
      base::nzchar(hit$contig)) {
    gt <- genbank_table[!base::is.na(genbank_table$contig) &
                          genbank_table$contig == hit$contig, , drop = FALSE]
    if (base::nrow(gt) == 0L) gt <- genbank_table
  } else {
    gt <- genbank_table
  }
  if (base::nrow(gt) == 0L) {
    return(.shared_pseudo_locus_row(hit, genome_id))
  }
  # Overlap calculation
  ov_lo <- base::pmax(gt$start, hit$start)
  ov_hi <- base::pmin(gt$end, hit$end)
  ov_len <- base::pmax(0L, ov_hi - ov_lo + 1L)
  hit_len <- hit$end - hit$start + 1L
  gt_len  <- gt$end - gt$start + 1L
  reciprocal_ov <- ov_len / base::pmin(hit_len, gt_len)
  midpoint <- base::as.integer((hit$start + hit$end) / 2L)
  contains_mid <- gt$start <= midpoint & gt$end >= midpoint
  candidate <- contains_mid | reciprocal_ov >= 0.6
  if (!base::any(candidate)) {
    return(.shared_pseudo_locus_row(hit, genome_id))
  }
  best <- base::which.max(base::ifelse(candidate, ov_len, -1L))
  row <- gt[best, , drop = FALSE]
  base::data.frame(
    genome_id    = genome_id,
    locus_tag    = base::as.character(row$locus_tag[1]),
    gene         = if ("gene" %in% base::colnames(row)) base::as.character(row$gene[1]) else NA_character_,
    product      = if ("product" %in% base::colnames(row)) base::as.character(row$product[1]) else NA_character_,
    match_field  = "sequence_homolog",
    sequence_hit = TRUE,
    homology_seed = hit$query,
    pident       = hit$pident,
    qcov         = hit$qcov,
    bits         = hit$bits,
    evalue       = hit$evalue,
    pseudo_locus = FALSE,
    stringsAsFactors = FALSE
  )
}

.shared_pseudo_locus_row <- function(hit, genome_id) {
  contig_str <- if (!base::is.null(hit$contig) && !base::is.na(hit$contig))
    base::as.character(hit$contig) else "contig1"
  pseudo <- base::sprintf("SEQHOM_%s_%s_%d_%d_%s",
                           genome_id, contig_str,
                           hit$start, hit$end, hit$strand)
  base::data.frame(
    genome_id    = genome_id,
    locus_tag    = pseudo,
    gene         = NA_character_,
    product      = "sequence homolog (no annotated CDS)",
    match_field  = "sequence_homolog",
    sequence_hit = TRUE,
    homology_seed = hit$query,
    pident       = hit$pident,
    qcov         = hit$qcov,
    bits         = hit$bits,
    evalue       = hit$evalue,
    pseudo_locus = TRUE,
    stringsAsFactors = FALSE
  )
}

# ---- Top-level dispatcher -------------------------------------------------
#' Find a homolog of `seed_cds` in `genome_record` using MMseqs2 (preferred)
#' or Biostrings (fallback). Seeds should be a named list of NT + AA
#' entries: list(seed_id = list(nt = "...", aa = "..."))
#' @keywords internal
.shared_find_homolog_by_sequence <- function(seeds, genome_id, genome_record,
                                               genbank_table, genbank_file,
                                               sensitivity = 7.0, evalue = 1e-5,
                                               min_aa_identity = 0.35,
                                               min_q_cov = 0.60, verbose = TRUE) {
  if (base::length(seeds) == 0L) return(NULL)
  # Prepare AA seed vector keyed by seed_id.
  seeds_aa <- base::vapply(seeds, function(s) {
    if (!base::is.null(s$aa) && base::nchar(s$aa) > 10L) base::as.character(s$aa)
    else if (!base::is.null(s$nt)) .shared_translate_cds(s$nt)
    else NA_character_
  }, character(1))
  seeds_aa <- seeds_aa[!base::is.na(seeds_aa) & base::nchar(seeds_aa) > 10L]
  if (base::length(seeds_aa) == 0L) return(NULL)
  seeds_aa <- seeds_aa[!base::duplicated(seeds_aa)]

  hits <- base::data.frame()
  used_tool <- NA_character_
  if (.shared_mmseqs_available()) {
    tmp_dir <- base::tempfile("mmseqs_homolog_"); base::dir.create(tmp_dir)
    base::on.exit(base::unlink(tmp_dir, recursive = TRUE, force = TRUE),
                   add = TRUE)
    # Write genome FASTA for mmseqs
    contig_name <- if (!base::is.null(genome_record$accession) &&
                        base::nzchar(genome_record$accession))
      genome_record$accession else genome_id
    genome_fna <- base::file.path(tmp_dir, "genome.fna")
    .shared_write_fasta(
      stats::setNames(base::as.character(genome_record$genome_seq),
                       contig_name),
      genome_fna)
    hits <- .shared_mmseqs_homolog_search(seeds_aa, genome_fna, tmp_dir,
                                            sensitivity = sensitivity,
                                            evalue = evalue)
    used_tool <- "mmseqs2"
  }

  if (base::nrow(hits) == 0L) {
    if (verbose && base::is.na(used_tool)) {
      base::message("[homolog] mmseqs2 not found, falling back to Biostrings.")
    }
    hits <- .shared_biostrings_homolog_search(
      seeds_aa, genbank_table, genome_record$genome_seq,
      min_aa_identity = min_aa_identity, min_q_cov = min_q_cov)
    used_tool <- if (base::is.na(used_tool)) "biostrings" else used_tool
  }
  if (base::nrow(hits) == 0L) return(NULL)

  # Keep only hits passing AA id / coverage, pick best by bitscore.
  hits <- hits[!base::is.na(hits$pident) &
                 hits$pident >= min_aa_identity &
                 !base::is.na(hits$qcov) & hits$qcov >= min_q_cov, ,
               drop = FALSE]
  if (base::nrow(hits) == 0L) return(NULL)
  hits <- hits[base::order(-hits$bits), , drop = FALSE]
  best <- hits[1, , drop = FALSE]

  if (verbose) {
    base::message(sprintf(
      "[homolog] %s: %s hit for %s at %s:%d-%d (%s) pid=%.2f qcov=%.2f bits=%.1f",
      genome_id, used_tool, best$query,
      if (!base::is.null(best$contig)) best$contig else "?",
      best$start, best$end, best$strand,
      best$pident, best$qcov, best$bits))
  }

  row <- .shared_map_hit_to_cds(best, genbank_table, genome_id)
  row$genbank_file <- genbank_file
  row
}
