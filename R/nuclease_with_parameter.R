#' Design gRNAs with Specific Nuclease Parameters
#'
#' Designs guide RNAs (gRNAs) for a target sequence using a specified nuclease.
#' Performs pre- and post-alignment filtering, including GC content, Tm, position, and off-target analysis.
#'
#' @param target_sequence A character string specifying the target DNA sequence for gRNA design.
#' @param strand A character string specifying the strand to filter ("both", "5", or "3"; default: "both").
#' @param bowtie_index A character string specifying the path to the Bowtie index (required if aligner = "bowtie" or "bwa").
#' @param genome A BSgenome object specifying the reference genome (required if aligner = "bowtie" or "bwa" or txObject is provided).
#' @param nuclease A character string specifying a predefined nuclease ("GeoCas9", "GtCas9", "FnCas12a", "FisCasI_B", "SpCas9") or a CrisprNuclease object for a custom nuclease (default: "GeoCas9").
#' @param position_range A numeric vector of length 2 specifying the position range (as percentage) for gRNA filtering (default: c(10, 90)).
#' @param tm_range A numeric vector of length 2 specifying the Tm range for gRNA filtering (default: c(50, 70)).
#' @param n_mismatches An integer specifying the maximum number of mismatches for off-target analysis (default: 1).
#' @param custom_seq A character string specifying a custom off-target sequence for biostrings alignment (default: NULL).
#' @param txObject A TxDb object specifying transcript annotations (default: NULL).
#' @param aligner A character string specifying the aligner to use ("bowtie", "bwa", or "biostrings"; default: "bowtie").
#' @param return_alignments A logical indicating whether to return alignment details (default: FALSE).
#' @return A data frame containing the filtered gRNA list with sequence features and off-target annotations. If return_alignments = TRUE, returns a list with the data frame and alignment details.
#' @export
nuclease_with_parameter <- function(
    target_sequence,
    strand = "both",
    bowtie_index,
    genome,
    nuclease = "GeoCas9",
    position_range = c(10, 90),
    tm_range = c(50, 70),
    n_mismatches = 1,
    custom_seq = NULL,
    txObject = NULL,
    aligner = "bowtie",
    return_alignments = FALSE
) {
  # Step 1: Define nuclease list
  nuclease_list <- list(
    GeoCas9 = CrisprNuclease(
      "GeoCas9",
      pams = c("NNNNCAAA"),
      pam_side = "3prime",
      spacer_length = 21
    ),
    FnCas12a = CrisprNuclease(
      "FnCas12a",
      pams = c("TTV"),
      pam_side = "5prime",
      spacer_length = 18
    ),
    FisCasI_B = CrisprNuclease(
      "FisCasI_B",
      pams = c("TTCA"),
      pam_side = "5prime",
      spacer_length = 28
    ),
    SpCas9 = CrisprNuclease(
      "SpCas9",
      pams = c("NGG"),
      pam_side = "3prime",
      spacer_length = 20
    )
  )

  # Step 2: Select nuclease
  if (base::is.character(nuclease)) {
    if (!nuclease %in% base::names(nuclease_list)) {
      base::stop("Nuclease not found: ", nuclease, ". Available options: ", base::paste(base::names(nuclease_list), collapse = ", "), " or provide a custom CrisprNuclease object.")
    }
    crispr_nuclease <- nuclease_list[[nuclease]]
  } else if (base::inherits(nuclease, "CrisprNuclease")) {
    # Validate CrisprNuclease object
    if (!"motifs" %in% base::slotNames(nuclease) || base::length(nuclease@motifs) == 0) {
      base::stop("Invalid CrisprNuclease object: missing or empty 'motifs' slot.")
    }
    if (!"pam_side" %in% base::slotNames(nuclease) || !nuclease@pam_side %in% base::c("3prime", "5prime")) {
      base::stop("Invalid CrisprNuclease object: missing or invalid 'pam_side' slot.")
    }
    if (!"spacer_length" %in% base::slotNames(nuclease) || nuclease@spacer_length <= 0) {
      base::stop("Invalid CrisprNuclease object: missing or invalid 'spacer_length' slot.")
    }
    crispr_nuclease <- nuclease
    nuclease <- crispr_nuclease@nucleaseName
  } else {
    base::stop("nuclease must be a character string or a valid CrisprNuclease object.")
  }
  base::cat("Using nuclease:", nuclease, "\n")
  base::cat("Motifs:", base::paste(crispr_nuclease@motifs, collapse = ", "), "\n")
  base::cat("PAM side:", crispr_nuclease@pam_side, "\n")
  base::cat("Spacer length:", crispr_nuclease@spacer_length, "\n")

  # Step 3: Validate aligner and inputs
  if (!aligner %in% base::c("bowtie", "bwa", "biostrings")) {
    base::stop("Invalid aligner: ", aligner, ". Must be 'bowtie', 'bwa', or 'biostrings'.")
  }
  if (aligner %in% base::c("bowtie", "bwa") && (base::is.null(bowtie_index) || base::is.null(genome))) {
    base::stop("bowtie_index and genome must be provided when aligner = 'bowtie' or 'bwa'.")
  }
  if (aligner == "biostrings" && base::is.null(custom_seq)) {
    base::stop("custom_seq must be provided when aligner = 'biostrings'.")
  }
  if (!base::is.numeric(n_mismatches) || n_mismatches < 0) {
    base::stop("n_mismatches must be a non-negative integer.")
  }

  # Step 3.5: Verify UCSC lookups are suppressed for non-model organisms
  # Primary suppression is in run_gRNA_list_generator (Step 5.1).
  # This is a defense-in-depth fallback for direct nuclease_with_parameter calls.
  # Uses unlockBinding + assign to directly replace functions in GenomeInfoDb namespace.
  .ucsc_orig_nwp <- base::list()
  .ucsc_gi_ns_nwp <- base::asNamespace("GenomeInfoDb")
  if (!base::is.null(genome)) {
    base::tryCatch({
      bsg_genome_tag <- genome@seqinfo@genome
      if (!base::all(base::is.na(bsg_genome_tag))) {
        # Test if UCSC suppression is already active (e.g., from run_gRNA_list_generator)
        test_result <- base::tryCatch({
          res <- GenomeInfoDb::getChromInfoFromUCSC("__dummy_test_genome__")
          if (base::is.data.frame(res) && base::nrow(res) == 0) "SUPPRESSED" else "LIVE"
        }, error = function(e) "LIVE")

        if (test_result == "LIVE") {
          base::cat("  Step 3.5: Applying UCSC suppression (fallback)...\n")

          .safe_fn <- function(genome, ...) {
            base::data.frame(chrom = base::character(0), size = base::integer(0),
                             assembled = base::logical(0), circular = base::logical(0),
                             stringsAsFactors = FALSE)
          }

          # Replace getChromInfoFromUCSC
          base::tryCatch({
            .ucsc_orig_nwp[["getChromInfoFromUCSC"]] <- base::get("getChromInfoFromUCSC", envir = .ucsc_gi_ns_nwp)
            unlockBinding("getChromInfoFromUCSC", .ucsc_gi_ns_nwp)
            base::assign("getChromInfoFromUCSC", .safe_fn, envir = .ucsc_gi_ns_nwp)
            lockBinding("getChromInfoFromUCSC", .ucsc_gi_ns_nwp)
            base::cat("    Replaced getChromInfoFromUCSC\n")
          }, error = function(e) {
            base::cat("    Warning: UCSC suppression failed: ", base::conditionMessage(e), "\n")
          })

          # Replace fetchExtendedChromInfoFromUCSC
          base::tryCatch({
            .ucsc_orig_nwp[["fetchExtendedChromInfoFromUCSC"]] <- base::get("fetchExtendedChromInfoFromUCSC", envir = .ucsc_gi_ns_nwp)
            unlockBinding("fetchExtendedChromInfoFromUCSC", .ucsc_gi_ns_nwp)
            base::assign("fetchExtendedChromInfoFromUCSC", function(genome, ...) {
              base::data.frame(UCSC_seqlevel = base::character(0), NCBI_seqlevel = base::character(0),
                               stringsAsFactors = FALSE)
            }, envir = .ucsc_gi_ns_nwp)
            lockBinding("fetchExtendedChromInfoFromUCSC", .ucsc_gi_ns_nwp)
          }, error = function(e) NULL)

          base::cat("    UCSC suppression applied (",
                    base::length(.ucsc_orig_nwp), " functions replaced)\n")
        } else {
          base::cat("  Step 3.5: UCSC lookups already suppressed -- OK\n")
        }
      } else {
        base::cat("  Step 3.5: BSgenome genome tag is NA -- no UCSC risk\n")
      }
    }, error = function(e) {
      base::cat("  Note: genome tag check: ", base::conditionMessage(e), "\n")
    })
  }
  # Cleanup: only restore functions WE replaced (not ones from run_gRNA_list_generator)
  base::on.exit({
    for (.fn_name in base::names(.ucsc_orig_nwp)) {
      base::tryCatch({
        unlockBinding(.fn_name, .ucsc_gi_ns_nwp)
        base::assign(.fn_name, .ucsc_orig_nwp[[.fn_name]], envir = .ucsc_gi_ns_nwp)
        lockBinding(.fn_name, .ucsc_gi_ns_nwp)
      }, error = function(e) NULL)
    }
  }, add = TRUE)

  # Step 4: Validate Bowtie index
  if (aligner %in% base::c("bowtie", "bwa")) {
    index_files <- base::list.files(base::dirname(bowtie_index), pattern = "\\.ebwt$", full.names = TRUE)
    if (base::length(index_files) < 6) {
      base::stop("Incomplete Bowtie index at: ", bowtie_index, ". Expected 6 .ebwt files, found: ", base::length(index_files))
    }
    base::cat("Bowtie index validated:", index_files, "\n")
    # Check sequence names in index and AUTO-FIX mismatch
    fna_file <- base::file.path(base::dirname(bowtie_index), "../seqs_srcdir", base::paste0(base::basename(bowtie_index), ".fna"))
    if (base::file.exists(fna_file)) {
      fna_seq <- Biostrings::readDNAStringSet(fna_file)
      bsgenome_seqnames <- GenomeInfoDb::seqnames(genome)
      fna_seq_prefixes <- base::sub(" .*", "", base::names(fna_seq))
      base::cat("  .fna seqnames:    ", base::paste(fna_seq_prefixes, collapse = ", "), "\n")
      base::cat("  BSgenome seqnames:", base::paste(bsgenome_seqnames, collapse = ", "), "\n")

      if (!base::all(fna_seq_prefixes %in% bsgenome_seqnames)) {
        base::cat("  [!] Seqname MISMATCH detected -- auto-fixing...\n")

        if (base::length(fna_seq_prefixes) == base::length(bsgenome_seqnames)) {
          # Rewrite .fna headers to match BSgenome
          base::names(fna_seq) <- bsgenome_seqnames
          Biostrings::writeXStringSet(fna_seq, fna_file)
          base::cat("  .fna headers rewritten to: ",
                    base::paste(bsgenome_seqnames, collapse = ", "), "\n")

          # Rebuild Bowtie index with corrected .fna
          bowtie_dir <- base::dirname(bowtie_index)
          accession  <- base::basename(bowtie_index)
          parent_dir <- base::dirname(bowtie_dir)

          base::cat("  Rebuilding Bowtie index...\n")
          base::unlink(bowtie_dir, recursive = TRUE)
          base::dir.create(bowtie_dir, recursive = TRUE, showWarnings = FALSE)
          base::tryCatch({
            run_bowtie_build(directory = parent_dir, accession = accession)
            base::cat("  Bowtie index rebuilt successfully\n")
          }, error = function(e) {
            base::cat("  Bowtie index rebuild failed -- trying Rbowtie directly...\n")
            base::tryCatch({
              Rbowtie::bowtie_build(
                fna_file,
                outdir = bowtie_dir,
                force  = TRUE,
                prefix = accession
              )
              base::cat("  Bowtie index rebuilt (direct) successfully\n")
            }, error = function(e2) {
              base::warning("Could not rebuild Bowtie index: ",
                             base::conditionMessage(e2),
                             "\n  n0/n1 may remain 0.")
            })
          })
        } else {
          base::warning(
            "Cannot auto-fix seqname mismatch: .fna has ",
            base::length(fna_seq_prefixes), " sequences but BSgenome has ",
            base::length(bsgenome_seqnames),
            ". n0/n1 may be 0 for all spacers."
          )
        }
      } else {
        base::cat("  Seqnames match -- OK\n")
      }
    }
  }

  # Step 5: Process target sequence
  target_dna <- target_sequence
  target_dna_clean <- base::toupper(target_dna)
  if (!base::grepl("^[ATGCN]+$", target_dna_clean)) {
    base::warning("Target sequence contains invalid characters. Only ATGCN are allowed. Invalid characters will be replaced with N.")
    target_dna_clean <- base::gsub("[^ATGCN]", "N", target_dna_clean)
  }
  target_dna_string <- Biostrings::DNAString(target_dna_clean)
  target_dna_set <- Biostrings::DNAStringSet(target_dna_clean)
  base::names(target_dna_set) <- base::paste0("coding_region")
  target_length <- base::length(target_dna_string)  # Use length() for DNAString, not nchar()


  # Step 6: Verify target sequence exists in genome
  if (aligner %in% base::c("bowtie", "bwa")) {
    base::cat("Verifying target sequence in BSgenome...\n")
    # Convert target_dna to DNAString and normalize to uppercase
    target_dna_clean <- base::toupper(target_dna)
    if (!base::grepl("^[ATGCN]+$", target_dna_clean)) {
      base::warning("Target sequence contains invalid characters. Only ATGCN are allowed. Invalid characters will be replaced with N.")
      target_dna_clean <- base::gsub("[^ATGCN]", "N", target_dna_clean)
    }
    target_dna_string <- Biostrings::DNAString(target_dna_clean)
    base::cat("Target sequence length:", base::length(target_dna_string), "\n")
    base::cat("Target sequence (first 20 chars):", base::as.character(Biostrings::subseq(target_dna_string, 1, base::min(20, base::length(target_dna_string)))), "...\n")
    found <- FALSE
    for (seqname in GenomeInfoDb::seqnames(genome)) {
      genome_seq <- genome[[seqname]]
      base::cat("Checking BSgenome sequence:", seqname, "type:", base::class(genome_seq), "length:", base::length(genome_seq), "\n")
      base::tryCatch({
        matches <- Biostrings::matchPattern(target_dna_string, genome_seq, max.mismatch = 0)
        if (base::length(matches) > 0) {
          found <- TRUE
          base::cat("Target sequence found in BSgenome sequence:", seqname, "with", base::length(matches), "exact matches\n")
          break
        }
      }, error = function(e) {
        base::cat("Error matching target sequence in BSgenome sequence:", seqname, ":", base::conditionMessage(e), "\n")
      })
    }
    if (!found) {
      base::warning("Target sequence not found in BSgenome. This may cause n0 = 0 for all spacers.")
    }
  }

  # Step 7: Verify PAM sites in target sequence
  # Use matchPattern (not vmatchPattern) for single DNAString subject
  # fixed = FALSE enables IUPAC ambiguity code matching (N, V, etc.)
  # NOTE: crispr_nuclease@motifs is an S4 slot -- must convert to character for for-loop
  pam_motifs <- base::as.character(crispr_nuclease@motifs)
  pam_count <- 0L
  target_dna_obj <- Biostrings::DNAString(target_dna)
  for (motif in pam_motifs) {
    pm <- Biostrings::matchPattern(motif, target_dna_obj, fixed = FALSE)
    pam_count <- pam_count + base::length(pm)
  }
  base::cat("Number of PAM sites (",
            base::paste(pam_motifs, collapse = ", "),
            ") in target sequence:", pam_count, "\n")
  if (pam_count == 0) {
    base::warning("No valid PAM sites found in target sequence for ", nuclease, ". This will cause n0 = 0.")
    return(base::data.frame())  # Return empty data frame if no PAM sites
  }

  # Step 8: Find spacers with timing
  start_time <- base::Sys.time()
  initial_guides <- crisprDesign::findSpacers(target_dna_set, crisprNuclease = crispr_nuclease, bsgenome = genome)
  end_time <- base::Sys.time()
  base::cat("Time for findSpacers:", base::difftime(end_time, start_time, units = "secs"), "seconds\n")
  base::cat("Initial spacers found:", base::length(initial_guides), "\n")
  if (base::length(initial_guides) > 0) {
    base::cat("Sample spacer sequences:", base::as.character(utils::head(initial_guides$protospacer)), "\n")
    base::cat("Sample PAM sites:", base::as.character(utils::head(initial_guides$pam)), "\n")
    base::cat("Sample PAM positions:", utils::head(initial_guides$pam_site), "\n")
    # Validate spacers and PAMs
    if (base::length(initial_guides$pam) == 0) {
      base::warning("No valid PAM sites assigned to spacers. This will cause n0 = 0.")
      return(base::data.frame())  # Return empty data frame if no PAMs
    } else {
      n_debug <- base::min(3L, base::length(initial_guides))
      for (i in base::seq_len(n_debug)) {
        spacer_seq <- base::as.character(initial_guides$protospacer[i])
        pam_seq <- base::as.character(initial_guides$pam[i])
        base::cat("Spacer", i, "protospacer:", spacer_seq, "PAM:", pam_seq, "\n")
      }
      if (base::length(initial_guides) > n_debug) {
        base::cat("... and", base::length(initial_guides) - n_debug, "more spacers\n")
      }
    }
  } else {
    base::cat("No spacers found. Returning empty data frame.\n")
    return(base::data.frame())
  }

  # Step 9: Pre-alignment filtering
  start_time <- base::Sys.time()
  filtered_guides <- crisprDesign::addSequenceFeatures(initial_guides)

  # Calculate Tm and position_percent
  filtered_guides$Tm <- base::sapply(base::as.character(filtered_guides$protospacer), function(seq) {
    tm <- TmCalculator::Tm_NN(seq, Na = 50, nn_table = "DNA_NN4", outlist = FALSE)
    return(tm)
  })
  filtered_guides$position_percent <- (filtered_guides$pam_site / target_length) * 100

  # Log ranges for debugging
  base::cat("Tm range before filtering:", base::range(filtered_guides$Tm, na.rm = TRUE), "\n")
  base::cat("Position range before filtering:", base::range(filtered_guides$position_percent, na.rm = TRUE), "\n")

  # --- Palindrome & secondary structure scoring ---
  # Check for palindromic sequences (self-complementary regions >= 6bp)
  # and estimate secondary structure propensity (consecutive self-complement runs)
  spacer_seqs <- base::as.character(filtered_guides$protospacer)

  .rc_simple <- function(s) {
    base::paste0(base::rev(base::sapply(base::strsplit(s, "")[[1]], function(b) {
      switch(b, A = "T", T = "A", G = "C", C = "G", "N")
    })), collapse = "")
  }

  # Palindrome score: fraction of spacer that is self-complementary (0 = no palindrome, 1 = perfect palindrome)
  palindrome_score <- base::sapply(spacer_seqs, function(seq) {
    rc_seq <- .rc_simple(seq)
    slen <- base::nchar(seq)
    if (slen < 6L) return(0)
    max_palin <- 0L
    # Slide window sizes from 6 to half the spacer length
    for (wsize in base::seq.int(6L, base::max(6L, slen %/% 2L))) {
      for (pos in base::seq_len(slen - wsize + 1L)) {
        subseq <- base::substring(seq, pos, pos + wsize - 1L)
        # Check if this subsequence appears in the reverse complement
        if (base::grepl(subseq, rc_seq, fixed = TRUE)) {
          max_palin <- base::max(max_palin, wsize)
        }
      }
    }
    return(max_palin / slen)
  }, USE.NAMES = FALSE)

  # Secondary structure score: longest self-complementary stem (internal hairpin potential)
  # A spacer forms hairpin if its 5' end base-pairs with its 3' end
  hairpin_score <- base::sapply(spacer_seqs, function(seq) {
    slen <- base::nchar(seq)
    if (slen < 8L) return(0)
    bases <- base::strsplit(seq, "")[[1]]
    comp <- base::sapply(bases, function(b) switch(b, A = "T", T = "A", G = "C", C = "G", "N"))
    max_stem <- 0L
    # Check from both ends inward (hairpin stem)
    stem <- 0L
    for (k in base::seq_len(slen %/% 2L - 1L)) {
      if (bases[k] == comp[slen - k + 1L]) {
        stem <- stem + 1L
      } else {
        if (stem >= 4L) max_stem <- base::max(max_stem, stem)
        stem <- 0L
      }
    }
    if (stem >= 4L) max_stem <- base::max(max_stem, stem)
    return(max_stem / (slen / 2))
  }, USE.NAMES = FALSE)

  # Add scores to GuideSet as metadata columns
  S4Vectors::mcols(filtered_guides)$palindrome_score <- palindrome_score
  S4Vectors::mcols(filtered_guides)$hairpin_score <- hairpin_score

  base::cat("Palindrome score range:", base::round(base::range(palindrome_score), 3), "\n")
  base::cat("Hairpin score range:", base::round(base::range(hairpin_score), 3), "\n")

  # Apply filters
  filtered_guides <- filtered_guides[filtered_guides$percentGC >= 20]
  filtered_guides <- filtered_guides[filtered_guides$percentGC <= 67]
  filtered_guides <- filtered_guides[!filtered_guides$polyA]
  filtered_guides <- filtered_guides[!filtered_guides$polyC]
  filtered_guides <- filtered_guides[!filtered_guides$polyT]
  filtered_guides <- filtered_guides[!filtered_guides$polyG]
  filtered_guides <- filtered_guides[!filtered_guides$startingGGGGG]
  filtered_guides <- filtered_guides[filtered_guides$Tm >= tm_range[1] & filtered_guides$Tm <= tm_range[2]]
  filtered_guides <- filtered_guides[filtered_guides$position_percent >= position_range[1] & filtered_guides$position_percent <= position_range[2]]

  # Hard filter: remove strong palindromes (>= 50% self-complementary)
  filtered_guides <- filtered_guides[filtered_guides$palindrome_score < 0.5]

  # Hard filter: remove strong hairpin formers (stem >= 60% of half-length)
  filtered_guides <- filtered_guides[filtered_guides$hairpin_score < 0.6]

  base::cat("Spacers after palindrome/hairpin filtering:", base::length(filtered_guides), "\n")

  end_time <- base::Sys.time()
  base::cat("Time for pre-alignment filtering (including Tm and position):", base::difftime(end_time, start_time, units = "secs"), "seconds\n")
  base::cat("Spacers after pre-filtering:", base::length(filtered_guides), "\n")

  # Step 10: Placeholder score (will be replaced by composite_score after alignment)
  filtered_guides$score <- base::rep(0, base::length(filtered_guides))

  # Step 10.5: Validate BSgenome before alignment
  base::tryCatch({
    bsg_sn <- GenomeInfoDb::seqnames(genome)
    bsg_sl <- GenomeInfoDb::seqlengths(genome)
    base::cat("  BSgenome seqlengths:", base::paste(bsg_sn, "=", bsg_sl, collapse = ", "), "\n")
    # If seqlengths are NA, crisprBowtie will drop all alignments
    if (base::any(base::is.na(bsg_sl))) {
      base::cat("  [!] seqlengths NA detected -- fixing in-place\n")
      # Recalculate from actual sequences
      for (sn_i in bsg_sn) {
        chr_seq <- genome[[sn_i]]
        if (!base::is.na(base::length(chr_seq))) {
          GenomeInfoDb::seqlengths(genome)[sn_i] <- base::length(chr_seq)
        }
      }
      base::cat("  seqlengths after fix:", base::paste(GenomeInfoDb::seqlengths(genome), collapse = ", "), "\n")
    }
  }, error = function(e) {
    base::cat("  [diag] seqlengths check failed:", base::conditionMessage(e), "\n")
  })
  utils::flush.console()

  # Step 11: Compute n0/n1 via crisprBowtie::runCrisprBowtie directly.
  # We bypass crisprDesign::addSpacerAlignments which has a bug where its internal
  # GRanges post-processing drops all alignments for custom (non-model) BSgenome
  # objects, resulting in n0=0 for all spacers. Using runCrisprBowtie directly
  # still uses the CrisprVerse Bowtie alignment engine but avoids the buggy wrapper.
  # mode="spacer" ensures only spacer sequences go to Bowtie (not PAM-expanded).
  start_time <- base::Sys.time()
  aligned_guides <- filtered_guides
  if (base::length(filtered_guides) > 0) {
    base::tryCatch({
      spacer_seqs <- base::as.character(filtered_guides$protospacer)
      spacer_seqs_unique <- base::unique(spacer_seqs)
      base::names(spacer_seqs_unique) <- spacer_seqs_unique

      # Single Bowtie run via crisprBowtie (CrisprVerse alignment engine)
      bowtie_msg <- utils::capture.output({
        direct_aln <- crisprBowtie::runCrisprBowtie(
          spacers = spacer_seqs_unique,
          bowtie_index = bowtie_index,
          bsgenome = genome,
          crisprNuclease = crispr_nuclease,
          n_mismatches = n_mismatches,
          mode = "spacer",
          canonical = FALSE,
          ignore_pam = TRUE
        )
      }, type = "message")

      end_time <- base::Sys.time()
      base::cat("  Bowtie alignment:", base::round(base::difftime(end_time, start_time, units = "secs"), 1), "sec\n")

      # Compute n0/n1 per spacer from crisprBowtie result
      n0_fixed <- FALSE
      if (base::nrow(direct_aln) > 0) {
        # Auto-detect match column (spacer or protospacer)
        match_col <- NULL
        for (try_col in base::c("spacer", "protospacer")) {
          if (try_col %in% base::colnames(direct_aln) &&
              base::any(spacer_seqs %in% base::as.character(direct_aln[[try_col]]))) {
            match_col <- try_col
            break
          }
        }
        if (!base::is.null(match_col)) {
          aln_match_vals <- base::as.character(direct_aln[[match_col]])
          new_n0 <- base::integer(base::length(aligned_guides))
          new_n1 <- base::integer(base::length(aligned_guides))
          for (si in base::seq_along(aligned_guides)) {
            hit_idx <- base::which(aln_match_vals == spacer_seqs[si])
            if (base::length(hit_idx) > 0) {
              mm_vals <- direct_aln$n_mismatches[hit_idx]
              new_n0[si] <- base::sum(mm_vals == 0L)
              new_n1[si] <- base::sum(mm_vals == 1L)
            }
          }
          S4Vectors::mcols(aligned_guides)$n0 <- new_n0
          S4Vectors::mcols(aligned_guides)$n1 <- new_n1
          base::cat("  n0 [", base::min(new_n0), ",", base::max(new_n0), "]",
                    " n1 [", base::min(new_n1), ",", base::max(new_n1), "]\n")
          n0_fixed <- TRUE
        }
      }

      # Biostrings matchPattern fallback
      if (!n0_fixed) {
        base::cat("  crisprBowtie match failed -- Biostrings fallback\n")
        bsg_sn <- GenomeInfoDb::seqnames(genome)
        manual_n0 <- base::integer(base::length(aligned_guides))
        manual_n1 <- base::integer(base::length(aligned_guides))
        for (si in base::seq_along(aligned_guides)) {
          proto_dna <- Biostrings::DNAString(spacer_seqs[si])
          proto_rc <- Biostrings::reverseComplement(proto_dna)
          total_n0 <- 0L
          total_n1_only <- 0L
          for (sn in bsg_sn) {
            chr_seq <- genome[[sn]]
            fwd_0 <- Biostrings::matchPattern(proto_dna, chr_seq, max.mismatch = 0)
            rev_0 <- Biostrings::matchPattern(proto_rc, chr_seq, max.mismatch = 0)
            ct0 <- base::length(fwd_0) + base::length(rev_0)
            total_n0 <- total_n0 + ct0
            fwd_1 <- Biostrings::matchPattern(proto_dna, chr_seq, max.mismatch = 1)
            rev_1 <- Biostrings::matchPattern(proto_rc, chr_seq, max.mismatch = 1)
            ct1 <- base::length(fwd_1) + base::length(rev_1) - ct0
            total_n1_only <- total_n1_only + ct1
          }
          manual_n0[si] <- total_n0
          manual_n1[si] <- total_n1_only
        }
        S4Vectors::mcols(aligned_guides)$n0 <- manual_n0
        S4Vectors::mcols(aligned_guides)$n1 <- manual_n1
        base::cat("  n0 [", base::min(manual_n0), ",", base::max(manual_n0), "]",
                  " n1 [", base::min(manual_n1), ",", base::max(manual_n1), "] (Biostrings)\n")
      }
    }, error = function(e) {
      end_time <- base::Sys.time()
      base::cat("Alignment error:", base::conditionMessage(e), "\n")
      base::stop("Failed to perform spacer alignments.")
    })
  } else {
    base::cat("No spacers to align after pre-filtering.\n")
  }

  # Step 12: Add sequence features (post-alignment)
  start_time <- base::Sys.time()
  aligned_guides_features <- crisprDesign::addSequenceFeatures(aligned_guides)
  end_time <- base::Sys.time()
  base::cat("Time for post-alignment sequence features:", base::difftime(end_time, start_time, units = "secs"), "seconds\n")

  # Post-alignment filtering for uniqueness (optional)
  # Uncomment to enforce n0 = 1
  # aligned_guides_features <- aligned_guides_features[aligned_guides_features$n0 == 1]  # Exactly one exact match
  # aligned_guides_features <- aligned_guides_features[aligned_guides_features$n1 == 0]  # No 1-mismatch off-targets
  base::cat("Spacers after post-alignment filtering:", base::length(aligned_guides_features), "\n")

  # Step 13: Convert to data frame
  final_df <- base::data.frame(aligned_guides_features)

  # Step 14: Filter by strand
  if (strand == "both") {
    base::cat("No strand filtering applied (strand = 'both')\n")
  } else if (strand == "5") {
    final_df <- final_df %>% dplyr::filter(strand == "-")
    base::cat("Filtering for '-' strand\n")
  } else if (strand == "3") {
    final_df <- final_df %>% dplyr::filter(strand == "+")
    base::cat("Filtering for '+' strand\n")
  } else {
    base::cat("Invalid strand argument '", strand, "', defaulting to no strand filtering\n")
  }

  # Step 15: Calculate composite score and sort by it
  if (base::nrow(final_df) > 0) {
    final_df <- calculate_composite_score(final_df)
    final_df <- final_df %>%
      dplyr::arrange(dplyr::desc(composite_score))

    base::cat("Spacers scored and sorted by composite_score (descending)\n")
  }

  base::cat("Final spacers:", base::nrow(final_df), "\n")

  # Step 16: Return results
  if (return_alignments) {
    return(base::list(
      gRNA_df = final_df,
      alignments = crisprBase::alignments(aligned_guides_features)
    ))
  } else {
    return(final_df)
  }
}
