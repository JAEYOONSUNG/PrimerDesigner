#' Read Vector File (SnapGene .dna, GenBank, or FASTA)
#'
#' Reads a vector sequence from various file formats commonly used in molecular cloning.
#' Supports SnapGene .dna files (via reticulate + Python snapgene_reader or SnapGene CLI),
#' GenBank (.gb, .gbk), and FASTA (.fa, .fasta, .fna) formats.
#'
#' For .dna files, the function tries these methods in order:
#' \enumerate{
#'   \item Python \code{snapgene_reader} via \code{reticulate}
#'   \item SnapGene application CLI conversion to GenBank
#' }
#' If SnapGene is running, it can interfere with CLI conversion. The function will
#' detect this and ask the user for permission to close SnapGene before proceeding.
#'
#' @param vector_file A character string specifying the path to the vector file.
#' @param kill_snapgene Logical. If TRUE, automatically kills any running SnapGene process
#'   without asking. If FALSE (default), the function will interactively ask the user
#'   for permission if SnapGene is detected running and CLI conversion fails.
#' @return A list with components:
#'   \describe{
#'     \item{sequence}{A character string of the vector DNA sequence (uppercase).}
#'     \item{name}{The vector name extracted from the file.}
#'     \item{length}{The length of the vector sequence in bp.}
#'     \item{file_type}{The detected file format ("dna", "genbank", or "fasta").}
#'   }
#' @examples
#' \dontrun{
#' vec <- read_vector_file("my_vector.dna")
#' vec <- read_vector_file("pUC19.gb")
#' vec <- read_vector_file("my_vector.dna", kill_snapgene = TRUE)
#' cat("Vector:", vec$name, "Length:", vec$length, "bp\n")
#' }
#' @export
read_vector_file <- function(vector_file, kill_snapgene = FALSE) {
  # Validate input
  if (base::missing(vector_file) || !base::is.character(vector_file) ||
      base::nchar(vector_file) == 0) {
    base::stop("vector_file must be a non-empty file path.")
  }
  if (!base::file.exists(vector_file)) {
    base::stop("Vector file not found: ", vector_file)
  }

  ext <- base::tolower(tools::file_ext(vector_file))
  vector_file <- base::normalizePath(vector_file, mustWork = TRUE)

  # Kill SnapGene upfront if explicitly requested
  if (kill_snapgene && ext == "dna") {
    .kill_snapgene()
  }

  result <- base::switch(ext,
    "dna" = .read_snapgene_dna(vector_file, kill_snapgene = kill_snapgene),
    "gb"  = , "gbk" = , "genbank" = .read_genbank_vector(vector_file),
    "fa"  = , "fasta" = , "fna" = .read_fasta_vector(vector_file),
    base::stop("Unsupported vector file format: .", ext,
         "\nSupported formats: .dna (SnapGene), .gb/.gbk (GenBank), .fa/.fasta (FASTA)")
  )

  # Validate result
  if (base::is.null(result$sequence) || base::nchar(result$sequence) == 0) {
    base::stop("Failed to extract sequence from vector file: ", vector_file)
  }

  base::cat("Vector loaded:", result$name, "(", result$length, "bp,",
            result$file_type, ")\n")
  return(result)
}


# --- Internal reader functions ---

#' @keywords internal
.read_snapgene_dna <- function(file_path, kill_snapgene = FALSE) {
  # --- Method 1: Python snapgene_reader via reticulate ---
  if (base::requireNamespace("reticulate", quietly = TRUE) &&
      reticulate::py_module_available("snapgene_reader")) {
    base::cat("Reading .dna file via snapgene_reader (Python)\n")
    py_result <- base::tryCatch({
      snapgene <- reticulate::import("snapgene_reader")
      snap_dict <- snapgene$snapgene_file_to_dict(file_path)

      seq_str <- base::toupper(base::as.character(snap_dict$seq))
      vec_name <- base::ifelse(
        !base::is.null(snap_dict$name) && base::nchar(snap_dict$name) > 0,
        snap_dict$name,
        tools::file_path_sans_ext(base::basename(file_path))
      )

      base::list(
        sequence  = seq_str,
        name      = vec_name,
        length    = base::nchar(seq_str),
        file_type = "dna"
      )
    }, error = function(e) {
      base::warning("snapgene_reader failed: ", base::conditionMessage(e),
                    "\nTrying SnapGene CLI fallback...")
      NULL
    })

    if (!base::is.null(py_result)) return(py_result)
  }

  # --- Method 2: SnapGene CLI conversion to GenBank ---
  cli_result <- .try_snapgene_cli(file_path)
  if (!base::is.null(cli_result)) return(cli_result)

  # --- CLI failed: check if SnapGene is running ---
  if (!kill_snapgene && .is_snapgene_running()) {
    base::cat("\n")
    base::cat("============================================================\n")
    base::cat("  SnapGene is currently running.\n")
    base::cat("  This can prevent the CLI from converting .dna files.\n")
    base::cat("  Close SnapGene to continue? (y/n): ")
    base::cat("\n============================================================\n")

    answer <- base::tryCatch({
      base::readline(prompt = "Close SnapGene? [y/n]: ")
    }, error = function(e) {
      # Non-interactive session (e.g., script mode)
      base::cat("Non-interactive session detected. Use kill_snapgene = TRUE.\n")
      "n"
    })

    if (base::tolower(base::trimws(answer)) %in% base::c("y", "yes")) {
      .kill_snapgene()
      # Retry CLI after killing SnapGene
      cli_result <- .try_snapgene_cli(file_path)
      if (!base::is.null(cli_result)) return(cli_result)
    } else {
      base::cat("SnapGene was not closed. CLI conversion may fail.\n")
      base::cat("Tip: use kill_snapgene = TRUE to auto-close SnapGene.\n")
    }
  } else if (!kill_snapgene) {
    # SnapGene not running but CLI still failed
    base::cat("SnapGene CLI conversion failed (SnapGene is not running).\n")
  }

  # --- Neither method available ---
  base::stop(
    "Cannot read SnapGene .dna file. Try one of:\n",
    "  1) Install Python snapgene_reader: reticulate::py_install('snapgene-reader')\n",
    "  2) Close SnapGene and retry, or use kill_snapgene = TRUE\n",
    "  3) Convert your .dna file to GenBank (.gb) format in SnapGene manually.\n",
    "     (File > Export > GenBank)"
  )
}


#' Try SnapGene CLI conversion to GenBank
#' @keywords internal
.try_snapgene_cli <- function(file_path) {
  snapgene_bin <- .find_snapgene_bin()
  if (base::is.null(snapgene_bin)) return(NULL)

  base::cat("Reading .dna file via SnapGene CLI:", snapgene_bin, "\n")
  temp_gb <- base::tempfile(fileext = ".gbk")
  base::on.exit(
    if (base::file.exists(temp_gb)) base::file.remove(temp_gb),
    add = TRUE
  )

  # Try multiple SnapGene CLI conversion command formats
  commands <- base::list(
    base::sprintf('"%s" --convert "GenBank - SnapGene" --input "%s" --output "%s"',
                  snapgene_bin, file_path, temp_gb),
    base::sprintf('"%s" --export-dna "%s" --format genbank --output "%s"',
                  snapgene_bin, file_path, temp_gb),
    base::sprintf('"%s" "%s" --convert "%s"',
                  snapgene_bin, file_path, temp_gb)
  )

  for (cmd in commands) {
    exit_code <- base::tryCatch(
      base::system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE, timeout = 30),
      error = function(e) -1,
      warning = function(w) -1
    )

    if (exit_code == 0 && base::file.exists(temp_gb) &&
        base::file.info(temp_gb)$size > 100) {
      result <- .read_genbank_vector(temp_gb)
      result$file_type <- "dna"
      if (result$name == tools::file_path_sans_ext(base::basename(temp_gb))) {
        result$name <- tools::file_path_sans_ext(base::basename(file_path))
      }
      return(result)
    }
    # Clean up temp file for next attempt
    if (base::file.exists(temp_gb)) base::file.remove(temp_gb)
  }

  return(NULL)
}


#' Find SnapGene binary path
#' @keywords internal
.find_snapgene_bin <- function() {
  snapgene_paths <- base::c(
    "/Applications/SnapGene.app/Contents/MacOS/SnapGene",
    "C:/Program Files/SnapGene/snapgene.exe",
    "C:/Program Files (x86)/SnapGene/snapgene.exe",
    base::Sys.which("snapgene")
  )

  for (p in snapgene_paths) {
    if (base::nchar(p) > 0 && base::file.exists(p)) {
      return(p)
    }
  }
  return(NULL)
}


#' Check if SnapGene is currently running
#' @keywords internal
.is_snapgene_running <- function() {
  os_type <- .Platform$OS.type
  if (os_type == "unix") {
    # macOS / Linux: check with pgrep
    res <- base::tryCatch(
      base::system("pgrep -f SnapGene", intern = FALSE,
                   ignore.stdout = TRUE, ignore.stderr = TRUE),
      error = function(e) 1
    )
    return(res == 0)
  } else {
    # Windows: check with tasklist
    res <- base::tryCatch({
      output <- base::system("tasklist /FI \"IMAGENAME eq snapgene.exe\"",
                             intern = TRUE, ignore.stderr = TRUE)
      base::any(base::grepl("snapgene", output, ignore.case = TRUE))
    }, error = function(e) FALSE)
    return(res)
  }
}


#' Kill running SnapGene processes
#' @keywords internal
.kill_snapgene <- function() {
  os_type <- .Platform$OS.type
  if (os_type == "unix") {
    # macOS / Linux
    res <- base::tryCatch(
      base::system("pkill -f SnapGene 2>/dev/null", intern = FALSE),
      error = function(e) 1
    )
    if (res == 0) {
      base::cat("SnapGene process terminated.\n")
      base::Sys.sleep(2)  # wait for file locks to release
    }
  } else {
    # Windows
    res <- base::tryCatch(
      base::system("taskkill /F /IM snapgene.exe 2>NUL", intern = FALSE),
      error = function(e) 1
    )
    if (res == 0) {
      base::cat("SnapGene process terminated.\n")
      base::Sys.sleep(2)
    }
  }
}


#' @keywords internal
.read_genbank_vector <- function(file_path) {
  # Read GenBank file — manual ORIGIN parsing (most robust)
  lines <- base::readLines(file_path, warn = FALSE)

  # Extract LOCUS name
  locus_line <- base::grep("^LOCUS", lines, value = TRUE)
  vec_name <- if (base::length(locus_line) > 0) {
    base::trimws(base::strsplit(locus_line[1], "\\s+")[[1]][2])
  } else {
    tools::file_path_sans_ext(base::basename(file_path))
  }

  # Extract sequence from ORIGIN section
  origin_start <- base::grep("^ORIGIN", lines)
  end_marker <- base::grep("^//", lines)
  if (base::length(origin_start) == 0 || base::length(end_marker) == 0) {
    base::stop("Cannot parse GenBank file: no ORIGIN/sequence section found in ",
               file_path)
  }

  valid_end <- end_marker[end_marker > origin_start[1]]
  if (base::length(valid_end) == 0) {
    base::stop("Cannot parse GenBank file: no terminator (//) after ORIGIN in ",
               file_path)
  }

  seq_lines <- lines[(origin_start[1] + 1):(valid_end[1] - 1)]
  seq_str <- base::toupper(base::gsub("[^a-zA-Z]", "",
                                       base::paste(seq_lines, collapse = "")))

  if (base::nchar(seq_str) == 0) {
    base::stop("No sequence data found in ORIGIN section of: ", file_path)
  }

  base::list(
    sequence  = seq_str,
    name      = vec_name,
    length    = base::nchar(seq_str),
    file_type = "genbank"
  )
}


#' @keywords internal
.read_fasta_vector <- function(file_path) {
  if (!base::requireNamespace("Biostrings", quietly = TRUE)) {
    base::stop("Package 'Biostrings' is required to read FASTA files.")
  }

  dna <- Biostrings::readDNAStringSet(file_path)
  if (base::length(dna) == 0) {
    base::stop("No sequences found in FASTA file: ", file_path)
  }

  # Use first sequence (vector should be a single sequence)
  if (base::length(dna) > 1) {
    base::warning("Multiple sequences found in FASTA file. Using the first one: ",
                  base::names(dna)[1])
  }

  seq_str <- base::toupper(base::as.character(dna[[1]]))
  vec_name <- base::names(dna)[1]

  base::list(
    sequence  = seq_str,
    name      = vec_name,
    length    = base::nchar(seq_str),
    file_type = "fasta"
  )
}
