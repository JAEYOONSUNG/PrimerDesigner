#' Package setup / runtime hooks
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Try to auto-patch snapgene_reader on load (best-effort; silent on
  # failure so the package still loads when Python isn't available).
  base::tryCatch(patch_snapgene_reader(verbose = FALSE),
                 error = function(e) NULL)
  # Probe MMseqs2 across the common install paths (homebrew, miniforge,
  # miniconda, anaconda) that Rscript doesn't inherit from the user's
  # shell rc. Cache the result so repeated homolog searches don't re-scan.
  base::tryCatch(.shared_cache_mmseqs_path(verbose = FALSE),
                 error = function(e) NULL)
  base::invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  bin <- base::tryCatch(.shared_cache_mmseqs_path(verbose = FALSE),
                         error = function(e) NA_character_)
  if (!base::is.na(bin) && base::nzchar(bin)) {
    base::packageStartupMessage(
      "[PrimerDesigner] MMseqs2 detected at ", bin,
      " (sequence-based homolog fallback enabled).")
  } else {
    base::packageStartupMessage(
      "[PrimerDesigner] MMseqs2 not found — sequence-based homolog fallback ",
      "will use the slower Biostrings engine.\n",
      "  Install via: `brew install mmseqs2` or `conda install -c bioconda mmseqs2`.")
  }
}

#' Locate the installed MMseqs2 binary (auto-caches the result)
#'
#' Scans the user's PATH plus common conda / miniforge / homebrew install
#' directories and caches the first hit on `options("PrimerDesigner.mmseqs_path")`.
#' Returns the absolute path, or `NA_character_` when no binary can be
#' found. Called automatically on package load; users rarely need to call
#' this directly.
#'
#' @param refresh Force a re-scan instead of using the cached value.
#' @return Path to `mmseqs` (character), or `NA_character_`.
#' @export
find_mmseqs_binary <- function(refresh = FALSE) {
  if (!base::isTRUE(refresh)) {
    cached <- base::getOption("PrimerDesigner.mmseqs_path", NA_character_)
    if (!base::is.na(cached) && base::nzchar(cached) &&
        base::file.exists(cached)) return(cached)
  }
  .shared_cache_mmseqs_path(verbose = FALSE)
}

#' @keywords internal
.shared_cache_mmseqs_path <- function(verbose = TRUE) {
  bin <- base::Sys.which("mmseqs")
  if (base::nzchar(bin) && base::file.exists(bin))
    bin <- base::unname(bin)
  else
    bin <- NA_character_
  if (base::is.na(bin)) {
    home <- base::Sys.getenv("HOME", "")
    candidates <- base::c(
      "/opt/homebrew/bin/mmseqs",
      "/opt/homebrew/opt/mmseqs2/bin/mmseqs",
      "/usr/local/bin/mmseqs",
      "/usr/local/opt/mmseqs2/bin/mmseqs",
      base::file.path(home, "miniforge3/bin/mmseqs"),
      base::file.path(home, "mambaforge/bin/mmseqs"),
      base::file.path(home, "miniconda3/bin/mmseqs"),
      base::file.path(home, "anaconda3/bin/mmseqs"),
      base::file.path(home, "opt/anaconda3/bin/mmseqs"),
      base::file.path(home, "Library/r-miniconda-arm64/bin/mmseqs"),
      base::file.path(home, "Library/r-miniconda/bin/mmseqs")
    )
    for (p in candidates) {
      if (base::file.exists(p)) { bin <- p; break }
    }
  }
  base::options(PrimerDesigner.mmseqs_path = bin)
  if (verbose && !base::is.na(bin))
    base::message("[PrimerDesigner] mmseqs cached at ", bin)
  bin
}

#' Patch the installed `snapgene_reader` Python module
#'
#' snapgene_reader 0.1.23 has two pain points for our workflow:
#'
#' 1. `snapgene_file_to_dict()` crashes on `.dna` files whose feature
#'    qualifiers carry more than two `V` items
#'    (`ValueError: too many values to unpack (expected 2)`).
#' 2. It IGNORES the `Primers` block (block id 5), so every primer tracked
#'    in SnapGene's Primers panel is lost on import.
#'
#' This helper rewrites the installed `snapgene_reader.py` in place to
#' address both issues. Safe to call repeatedly: each fix is applied only
#' if the unpatched signature is present. No-op when `snapgene_reader` is
#' not installed.
#'
#' @param verbose Whether to emit a message when a patch is applied.
#' @return `TRUE` on successful patch (or already-patched), `FALSE` when
#'   the module can't be located.
#' @export
patch_snapgene_reader <- function(verbose = TRUE) {
  if (!base::requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
  if (!reticulate::py_module_available("snapgene_reader")) return(FALSE)

  mod <- base::tryCatch(reticulate::import("snapgene_reader"),
                         error = function(e) NULL)
  if (base::is.null(mod)) return(FALSE)
  py_file <- base::tryCatch(mod$`__file__`, error = function(e) NULL)
  if (base::is.null(py_file) || !base::file.exists(py_file)) return(FALSE)

  src <- base::readLines(py_file, warn = FALSE)
  orig <- src

  # Fix 1: (fmt1, value1), (_, value2) = e_v.items() -> tolerant of >2 items.
  bad_unpack <- "                                (fmt1, value1), (_, value2) = e_v.items()"
  new_unpack <- base::c(
    "                                items = list(e_v.items())",
    "                                if len(items) < 2:",
    "                                    continue",
    "                                (fmt1, value1), (_, value2) = items[0], items[1]"
  )
  idx <- base::which(src == bad_unpack)
  if (base::length(idx) == 1L) {
    src <- base::c(src[base::seq_len(idx - 1L)],
                    new_unpack,
                    src[base::seq.int(idx + 1L, base::length(src))])
  }

  # Fix 2: add a block-5 Primers handler just before the generic else-branch
  # that skips unknown blocks.
  else_marker <- "        else:"
  skip_line <- "            # WE IGNORE THE WHOLE BLOCK"
  if (!base::any(base::grepl("READ PRIMERS \\(SnapGene Primers panel", src))) {
    # locate the `        else:` line that appears RIGHT BEFORE
    # `            # WE IGNORE THE WHOLE BLOCK`.
    for (i in base::seq_along(src)) {
      if (src[i] == else_marker && i < base::length(src) &&
          src[i + 1L] == skip_line) {
        block5 <- base::c(
          "        elif ord(next_byte) == 5:",
          "            # READ PRIMERS (SnapGene Primers panel) and emit",
          "            # each binding site as a primer_bind feature so the",
          "            # info survives a Biopython round-trip.",
          "            try:",
          "                block_content = fileobject.read(block_size).decode('utf-8')",
          "                primers_parsed = xmltodict.parse(block_content)",
          "                primers_list = primers_parsed.get('Primers', {}).get(",
          "                    'Primer', [])",
          "                if not isinstance(primers_list, list):",
          "                    primers_list = [primers_list]",
          "                for primer in primers_list:",
          "                    name = primer.get('@name', 'primer')",
          "                    seq_val = primer.get('@sequence', '')",
          "                    desc = primer.get('@description', '')",
          "                    sites = primer.get('BindingSite', [])",
          "                    if not isinstance(sites, list):",
          "                        sites = [sites]",
          "                    for site in sites:",
          "                        if site.get('@simplified') == '1':",
          "                            continue",
          "                        loc = site.get('@location', '')",
          "                        bound_strand = site.get('@boundStrand', '0')",
          "                        strand_sym = '-' if bound_strand == '1' else '+'",
          "                        try:",
          "                            parts = [int(x) for x in loc.split('-')]",
          "                            if len(parts) < 2:",
          "                                continue",
          "                            start_i, end_i = sorted(parts[:2])",
          "                        except (ValueError, AttributeError):",
          "                            continue",
          "                        seg = {",
          "                            '@range': '%d-%d' % (start_i, end_i),",
          "                            '@color': '#a020f0',",
          "                            '@type': 'standard',",
          "                        }",
          "                        notes = ['color: #a020f0']",
          "                        if seq_val:",
          "                            notes.append('sequence: ' + seq_val)",
          "                        if desc:",
          "                            notes.append(desc)",
          "                        quals = {'label': name, 'note': notes}",
          "                        data['features'].append(dict(",
          "                            type='primer_bind',",
          "                            strand=strand_sym,",
          "                            start=start_i,",
          "                            end=end_i,",
          "                            name=name,",
          "                            color='#a020f0',",
          "                            textColor='black',",
          "                            segments=[seg],",
          "                            row=0,",
          "                            isOrf=False,",
          "                            qualifiers=quals,",
          "                        ))",
          "            except Exception:",
          "                pass",
          ""
        )
        src <- base::c(src[base::seq_len(i - 1L)],
                        block5,
                        src[base::seq.int(i, base::length(src))])
        break
      }
    }
  }

  if (!base::identical(src, orig)) {
    base::writeLines(src, py_file)
    if (verbose) base::message(
      "[PrimerDesigner] patched snapgene_reader at ", py_file)
  }
  TRUE
}
