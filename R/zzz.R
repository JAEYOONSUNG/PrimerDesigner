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
    conda <- base::tryCatch(.shared_find_conda_bin(),
                             error = function(e) NA_character_)
    if (!base::is.na(conda) && base::nzchar(conda)) {
      base::packageStartupMessage(
        "[PrimerDesigner] MMseqs2 not found (Biostrings fallback will be used).\n",
        "  Conda detected at ", conda,
        " — run `install_mmseqs_via_conda()` to install it there.")
    } else {
      base::packageStartupMessage(
        "[PrimerDesigner] MMseqs2 not found — sequence-based homolog fallback ",
        "will use the slower Biostrings engine.\n",
        "  Install via: `brew install mmseqs2` or `conda install -c bioconda mmseqs2`.")
    }
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

#' Locate a conda-family executable (conda / mamba / micromamba)
#'
#' Probes the user's PATH plus every common miniforge / mambaforge /
#' miniconda / anaconda install prefix. Prefers `mamba` / `micromamba`
#' over `conda` when available (much faster solver).
#' @keywords internal
.shared_find_conda_bin <- function() {
  for (exe in base::c("micromamba", "mamba", "conda")) {
    bin <- base::Sys.which(exe)
    if (base::nzchar(bin) && base::file.exists(bin))
      return(base::unname(bin))
  }
  home <- base::Sys.getenv("HOME", "")
  for (exe in base::c("micromamba", "mamba", "conda")) {
    for (prefix in base::c("miniforge3", "mambaforge", "miniconda3",
                             "anaconda3", "opt/anaconda3",
                             "Library/r-miniconda-arm64",
                             "Library/r-miniconda")) {
      cand <- base::file.path(home, prefix, "bin", exe)
      if (base::file.exists(cand)) return(cand)
    }
  }
  NA_character_
}

#' Install MMseqs2 into the user's conda environment
#'
#' Convenience wrapper that calls the detected conda-family binary
#' (`micromamba` / `mamba` / `conda`, picked in that order) to install
#' MMseqs2 from the bioconda channel. When `env_name` is supplied the
#' package lands in that env; otherwise the base env is used.
#'
#' @param env_name Optional conda environment name. `NULL` (default)
#'   installs into the base env; supply a name to target a dedicated env.
#'   The env will be created if it doesn't exist.
#' @param channel Conda channel to pull from. Default "bioconda".
#' @param refresh Whether to re-probe the `mmseqs` binary cache after
#'   install. Default `TRUE`.
#' @return Path to the freshly installed `mmseqs` binary, or
#'   `NA_character_` on failure.
#' @export
install_mmseqs_via_conda <- function(env_name = NULL,
                                      channel = "bioconda",
                                      refresh = TRUE) {
  conda <- .shared_find_conda_bin()
  if (base::is.na(conda)) {
    base::stop("No conda / mamba / micromamba binary found. ",
               "Install miniforge, miniconda or anaconda first, or use ",
               "`brew install mmseqs2`.")
  }
  # Create env if specified and missing.
  if (!base::is.null(env_name) && base::nzchar(env_name)) {
    env_list <- base::suppressWarnings(base::system2(
      conda, base::c("env", "list"), stdout = TRUE, stderr = TRUE))
    if (!base::any(base::grepl(base::paste0("^", env_name, "\\s"), env_list))) {
      base::message("[PrimerDesigner] creating conda env '", env_name, "' ...")
      base::system2(conda, base::c("create", "-y", "-n", env_name,
                                     "-c", channel, "mmseqs2"))
    } else {
      base::system2(conda, base::c("install", "-y", "-n", env_name,
                                     "-c", channel, "mmseqs2"))
    }
  } else {
    base::message("[PrimerDesigner] installing mmseqs2 into base env via ",
                   base::basename(conda), " ...")
    base::system2(conda, base::c("install", "-y", "-c", channel, "mmseqs2"))
  }
  if (base::isTRUE(refresh)) find_mmseqs_binary(refresh = TRUE) else NA_character_
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
