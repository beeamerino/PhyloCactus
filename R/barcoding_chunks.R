# ------------------------------------------------------------------------------
# Molecular diagnostic branch (11_barcoding/7_classifier/), step 7 in chunks. PhyloCactus 0.5.0,
# Phase 6B.
#
# The full IdTaxa run is about 34 hours of one core (projection of 23-09) and runs on a cluster as a
# job array: the folds are dealt to N chunks, each chunk is one task, and a last job merges the chunk
# tables into the tables a single run would have written, byte for byte, or refuses (decisions K4
# and K5 of BMM, 2026-09-26).
# ------------------------------------------------------------------------------

#' Arguments of a chunked run, checked before anything is read
#' @noRd
.bc_check_chunk <- function(chunk, n_chunks) {
  if (is.null(chunk) && is.null(n_chunks)) return(invisible(TRUE))
  if (is.null(n_chunks) || length(n_chunks) != 1L || is.na(n_chunks) || n_chunks < 1 ||
      n_chunks != round(n_chunks)) {
    stop("n_chunks must be a whole number of 1 or more, given together with chunk.", call. = FALSE)
  }
  if (is.null(chunk) || length(chunk) != 1L || is.na(chunk) || chunk < 1 || chunk > n_chunks ||
      chunk != round(chunk)) {
    stop("chunk must be a whole number between 1 and n_chunks (", n_chunks, ").", call. = FALSE)
  }
  invisible(TRUE)
}

#' Chunk of each fold of one locus and scheme: round-robin in the order of step 5
#'
#' Every chunk gets the same share of every locus, within one fold, so the slow loci are spread
#' evenly without a model of their cost.
#' @noRd
.bc_chunk_assign <- function(n, n_chunks) {
  as.integer((seq_len(n) - 1L) %% as.integer(n_chunks) + 1L)
}

#' The queries a whole run classifies, in the order it writes them
#'
#' The same walk as `classify_barcoding_folds()`: schemes, loci in radix order, folds of step 5 with
#' the declared sample, queries in the order of the fold. The merge checks the chunks against it and
#' orders the rows by it.
#' @noRd
.bc_planned_queries <- function(library_dir, folds_dir, schemes, loci, max_folds, seed) {
  lib <- .bc_library_from_dir(library_dir)
  loci_all <- sort(unique(lib$locus), method = "radix")
  loci <- if (is.null(loci)) loci_all else intersect(loci_all, loci)
  out <- list()
  for (sc in schemes) {
    f <- file.path(folds_dir, paste0("TABLE_barcoding_folds_", sc, ".csv"))
    if (!file.exists(f)) {
      stop("No ", basename(f), " in ", folds_dir, ". Run build_barcoding_folds() first.", call. = FALSE)
    }
    tab <- utils::read.csv(f, stringsAsFactors = FALSE)
    rows <- list()
    for (l in loci) {
      d <- lib[lib$locus == l, , drop = FALSE]
      folds <- .bc_sample_folds(.bc_folds_from_table(tab[tab$locus == l, , drop = FALSE], d$sid),
                                max_folds, seed)
      for (p in folds) {
        rows[[length(rows) + 1L]] <- data.frame(locus = l, fold = p$fold, sid = p$test_ids,
                                                stringsAsFactors = FALSE)
      }
    }
    out[[sc]] <- if (length(rows) > 0L) do.call(rbind, rows) else
      data.frame(locus = character(0), fold = integer(0), sid = character(0))
  }
  out
}

.bc_read_chunk <- function(f) {
  x <- tryCatch(utils::read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0L || !all(c("locus", "fold", "sid") %in% names(x))) return(NULL)
  x
}

#' Merge the Chunks of a Chunked Run of Step 7
#'
#' Reads the tables that `classify_barcoding_folds(chunk = k, n_chunks = n)` wrote in `chunks/` of
#' `output_dir`, checks them against the folds of step 5 and writes the tables of the whole run,
#' which are the tables a run without chunks would have written, row for row and in the same order.
#'
#' @details
#' The merge refuses, naming what it found, when a chunk is missing, when a fold appears that is not
#' among the folds planned for the run (the folds of step 5, with the declared sample of
#' `max_folds`), when a query appears more than once, and when a planned query is missing. The
#' arguments that choose the folds (`schemes`, `loci`, `max_folds`, `seed`) must be those of the
#' chunks. The timing table sums the seconds of the chunks by locus and scheme; the chunk tables and
#' their sessions are left in `chunks/`.
#'
#' @param library_dir,folds_dir,output_dir,method,alignment,schemes,loci,max_folds,seed As in
#'   [classify_barcoding_folds()], with the values the chunks were run with.
#' @param n_chunks Integer. Number of chunks of the run.
#' @return Invisibly, a list of prediction tables by scheme. Writes
#'   `TABLE_barcoding_predictions_<scheme>_<method>.csv` and `TABLE_barcoding_timing_<method>.csv`
#'   in `output_dir`.
#' @seealso [classify_barcoding_folds()], [generate_barcoding_job_scripts()].
#' @examples
#' \dontrun{
#' merge_barcoding_chunks(
#'   library_dir = "11_barcoding/4_library",
#'   folds_dir = "11_barcoding/5_folds",
#'   output_dir = "11_barcoding/7_classifier",
#'   method = "idtaxa", n_chunks = 60
#' )
#' }
#' @export
merge_barcoding_chunks <- function(library_dir = file.path("11_barcoding", "4_library"),
                                   folds_dir = file.path("11_barcoding", "5_folds"),
                                   output_dir = file.path("11_barcoding", "7_classifier"),
                                   method = c("nn", "idtaxa"),
                                   n_chunks,
                                   alignment = c("library", "add"),
                                   schemes = c("species", "genus"),
                                   loci = NULL,
                                   max_folds = NULL,
                                   seed = 1L) {
  method <- match.arg(method)
  alignment <- match.arg(alignment)
  .bc_check_chunk(1L, n_chunks)
  n_chunks <- as.integer(n_chunks)
  .bc_assert_output_dir(output_dir)
  suffix <- paste0(method, if (alignment == "add") "_add" else "")
  chunk_dir <- file.path(output_dir, "chunks")
  tags <- paste0("_chunk", seq_len(n_chunks), "of", n_chunks)
  planned <- .bc_planned_queries(library_dir, folds_dir, schemes, loci, max_folds, seed)

  out <- list()
  for (sc in schemes) {
    files <- file.path(chunk_dir, paste0("TABLE_barcoding_predictions_", sc, "_", suffix, tags, ".csv"))
    for (k in seq_len(n_chunks)) {
      if (!file.exists(files[k])) {
        stop("Missing chunk ", k, " of ", n_chunks, " for scheme ", sc, ": ", files[k], call. = FALSE)
      }
    }
    parts <- Filter(Negate(is.null), lapply(files, .bc_read_chunk))
    x <- if (length(parts) > 0L) do.call(rbind, parts) else NULL
    pl <- planned[[sc]]
    if (is.null(x)) {
      if (nrow(pl) > 0L) stop("Every chunk of scheme ", sc, " is empty, and ", nrow(pl),
                              " queries were planned.", call. = FALSE)
      next
    }
    fold_key <- paste(x$locus, x$fold, sep = "\r")
    foreign <- !fold_key %in% paste(pl$locus, pl$fold, sep = "\r")
    if (any(foreign)) {
      i <- which(foreign)[1]
      stop("Fold ", x$fold[i], " of locus ", x$locus[i], ", scheme ", sc,
           ", is not among the folds of step 5 planned for this run.", call. = FALSE)
    }
    key <- paste(x$locus, x$fold, x$sid, sep = "\r")
    if (anyDuplicated(key)) {
      i <- which(duplicated(key))[1]
      stop("Query ", x$sid[i], " of fold ", x$fold[i], ", locus ", x$locus[i], ", scheme ", sc,
           ", appears more than once in the chunks.", call. = FALSE)
    }
    pl_key <- paste(pl$locus, pl$fold, pl$sid, sep = "\r")
    absent <- !pl_key %in% key
    if (any(absent)) {
      i <- which(absent)[1]
      stop("Query ", pl$sid[i], " of fold ", pl$fold[i], ", locus ", pl$locus[i], ", scheme ", sc,
           ", is missing from the chunks (", sum(absent), " planned queries missing).", call. = FALSE)
    }
    pred <- x[match(pl_key, key), , drop = FALSE]
    rownames(pred) <- NULL
    out[[sc]] <- pred
    utils::write.csv(pred, file.path(output_dir, paste0("TABLE_barcoding_predictions_", sc, "_", suffix, ".csv")),
                     row.names = FALSE)
  }

  tim_files <- file.path(chunk_dir, paste0("TABLE_barcoding_timing_", suffix, tags, ".csv"))
  tim <- do.call(rbind, Filter(Negate(is.null), lapply(tim_files[file.exists(tim_files)], function(f) {
    t <- tryCatch(utils::read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(t) || nrow(t) == 0L) NULL else t
  })))
  if (!is.null(tim)) {
    order_key <- unlist(lapply(schemes, function(sc) paste(sc, unique(planned[[sc]]$locus), sep = "\r")))
    k <- paste(tim$scheme, tim$locus, sep = "\r")
    agg <- do.call(rbind, lapply(order_key[order_key %in% k], function(ok) {
      d <- tim[k == ok, , drop = FALSE]
      data.frame(locus = d$locus[1], scheme = d$scheme[1], method = d$method[1], alignment = d$alignment[1],
                 folds = sum(d$folds), queries = sum(d$queries), seconds = sum(d$seconds),
                 stringsAsFactors = FALSE)
    }))
    utils::write.csv(agg, file.path(output_dir, paste0("TABLE_barcoding_timing_", suffix, ".csv")),
                     row.names = FALSE)
  }
  message("Merged ", n_chunks, " chunks of method ", suffix, ": every planned query present once.")
  .bc_classifier_banner(out, output_dir, suffix)
  invisible(out)
}

#' Generate the SLURM Job Scripts of a Chunked Run of Step 7
#'
#' Writes the scripts that run [classify_barcoding_folds()] as a SLURM job array, one chunk per
#' task on one core, and [merge_barcoding_chunks()] once every task has finished without error. It
#' follows the conventions of [generate_ml_search_script()] and [generate_bootstrap_script()]:
#' `cluster_*` arguments, `load_module`, one log per array task, and the versions the job ran with
#' printed at its start.
#'
#' @details
#' Five files are written in `job_dir`: `run_chunk.R` and `run_merge.R`, which hold the arguments of
#' the run; `job_array.sh`, the array (`--array=1-n_chunks`, one core per task); `job_merge.sh`, the
#' merge; and `submit.sh`, which submits the array and then the merge with
#' `--dependency=afterok:<array>`. Submit with `bash submit.sh` from any directory.
#'
#' The paths are written as given and must be those of the cluster: pass absolute paths, or paths
#' relative to the directory the job is submitted from. Nothing about a particular cluster is
#' written unless it is passed: with `load_module = NULL` no module is loaded, with
#' `cluster_queue = NULL` no QoS is requested, and with `r_lib = NULL` the default R library is used.
#'
#' @param library_dir,folds_dir Character. As in [classify_barcoding_folds()], on the cluster.
#' @param classifier_dir Character. `output_dir` of the run, on the cluster.
#' @param method,seed,threshold,schemes,loci,max_folds As in [classify_barcoding_folds()].
#' @param n_chunks Integer. Number of chunks, one array task each.
#' @param job_dir Character. Where the scripts are written. Defaults to `job/` inside
#'   `classifier_dir`.
#' @param cluster_job_name Character. SLURM job name; the merge job adds `_merge`. Defaults to
#'   `"cactus_idtaxa"`.
#' @param cluster_partition Character. SLURM partition. Defaults to `"main"`.
#' @param cluster_mem Character. Memory per task (`--mem`). Defaults to `"4G"`.
#' @param cluster_time Character. Time limit per task (`-t`). Defaults to `"04:00:00"`.
#' @param cluster_queue Character. Optional SLURM QoS (`-q`). Defaults to `NULL`.
#' @param cluster_mail_user Character. Address for SLURM mail. Defaults to
#'   `Sys.getenv("MY_EMAIL", "")`; an empty value writes no mail line. The array mails on failure
#'   only; the merge mails at its end.
#' @param load_module Character vector of environment modules to load, or `NULL`.
#' @param r_lib Character. R library to put first (`R_LIBS_USER`), or `NULL`.
#' @param rscript_exec Character. The `Rscript` command. Defaults to `"Rscript"`.
#' @return Invisibly, the paths of the five files.
#' @seealso [classify_barcoding_folds()], [merge_barcoding_chunks()].
#' @examples
#' \dontrun{
#' generate_barcoding_job_scripts(
#'   library_dir = "/home/user/PhyloCactus_Tutorial/11_barcoding/4_library",
#'   folds_dir = "/home/user/PhyloCactus_Tutorial/11_barcoding/5_folds",
#'   classifier_dir = "/home/user/PhyloCactus_Tutorial/11_barcoding/7_classifier",
#'   method = "idtaxa", n_chunks = 60,
#'   load_module = "R/4.4.0", r_lib = "~/R/phylocactus/4.4"
#' )
#' }
#' @export
generate_barcoding_job_scripts <- function(library_dir, folds_dir, classifier_dir,
                                           method = c("idtaxa", "nn"),
                                           n_chunks,
                                           job_dir = file.path(classifier_dir, "job"),
                                           seed = 1L,
                                           threshold = 60,
                                           schemes = c("species", "genus"),
                                           loci = NULL,
                                           max_folds = NULL,
                                           cluster_job_name = "cactus_idtaxa",
                                           cluster_partition = "main",
                                           cluster_mem = "4G",
                                           cluster_time = "04:00:00",
                                           cluster_queue = NULL,
                                           cluster_mail_user = Sys.getenv("MY_EMAIL", ""),
                                           load_module = NULL,
                                           r_lib = NULL,
                                           rscript_exec = "Rscript") {
  method <- match.arg(method)
  if (missing(n_chunks) || length(n_chunks) != 1L || is.na(n_chunks) || n_chunks < 1 ||
      n_chunks != round(n_chunks)) {
    stop("n_chunks must be a whole number of 1 or more.", call. = FALSE)
  }
  n_chunks <- as.integer(n_chunks)
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
  dq <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
  common <- c(
    paste0("  library_dir = ", dq(library_dir), ","),
    paste0("  folds_dir = ", dq(folds_dir), ","),
    paste0("  output_dir = ", dq(classifier_dir), ","),
    paste0("  method = ", dq(method), ","),
    paste0("  schemes = ", dq(schemes), ","),
    paste0("  loci = ", dq(loci), ","),
    paste0("  max_folds = ", dq(max_folds), ","),
    paste0("  seed = ", dq(as.integer(seed)), ","))
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M")
  run_chunk <- c(
    paste0("# One chunk of step 7. Written by generate_barcoding_job_scripts() on ", stamp, "."),
    "# Called by job_array.sh with the array task number as its only argument.",
    "chunk <- as.integer(commandArgs(trailingOnly = TRUE)[1])",
    "if (is.na(chunk)) stop(\"No chunk number given.\", call. = FALSE)",
    "library(PhyloCactus)",
    "classify_barcoding_folds(",
    common,
    paste0("  threshold = ", dq(threshold), ","),
    "  chunk = chunk,",
    paste0("  n_chunks = ", n_chunks, "L"),
    ")")
  run_merge <- c(
    paste0("# Merge of the chunks of step 7. Written by generate_barcoding_job_scripts() on ", stamp, "."),
    "library(PhyloCactus)",
    "merge_barcoding_chunks(",
    common,
    paste0("  n_chunks = ", n_chunks, "L"),
    ")")

  header <- function(name, array, log, mail_type) {
    h <- c("#!/bin/bash",
           paste0("#SBATCH -J ", name),
           paste0("#SBATCH -p ", cluster_partition),
           "#SBATCH --nodes=1",
           "#SBATCH --ntasks-per-node=1",
           "#SBATCH --cpus-per-task=1",
           paste0("#SBATCH --mem=", cluster_mem),
           paste0("#SBATCH -t ", cluster_time))
    if (array) h <- c(h, paste0("#SBATCH --array=1-", n_chunks))
    h <- c(h, paste0("#SBATCH -o ", log, ".out"), paste0("#SBATCH -e ", log, ".err"))
    if (!is.null(cluster_queue) && nzchar(cluster_queue)) h <- c(h, paste0("#SBATCH -q ", cluster_queue))
    if (!is.null(cluster_mail_user) && nzchar(trimws(cluster_mail_user))) {
      h <- c(h, paste0("#SBATCH --mail-type=", mail_type),
             paste0("#SBATCH --mail-user=", trimws(cluster_mail_user)))
    }
    h <- c(h, "#", "")
    if (!is.null(load_module) && length(load_module) > 0L && nzchar(trimws(paste(load_module, collapse = " ")))) {
      h <- c(h, paste0("module load ", paste(load_module, collapse = " ")), "")
    }
    if (!is.null(r_lib) && nzchar(r_lib)) {
      h <- c(h, paste0("export R_LIBS_USER=", shQuote(r_lib)), "")
    }
    c(h,
      "# Move to the submit directory, from which relative paths are resolved",
      "cd \"${SLURM_SUBMIT_DIR:-$PWD}\"",
      "",
      "echo \"=== Versions used by this job ===\"",
      paste0(rscript_exec, " -e 'cat(R.version.string, \"\\n\"); for (p in c(\"DECIPHER\", \"Biostrings\", \"PhyloCactus\")) ",
             "cat(p, tryCatch(as.character(packageVersion(p)), error = function(e) \"not installed\"), \"\\n\")'"),
      "echo \"================================\"",
      "")
  }
  job_array <- c(header(cluster_job_name, TRUE, "slurm_%A_%a", "FAIL"),
                 paste0("echo \"Chunk ${SLURM_ARRAY_TASK_ID} of ", n_chunks, "\""),
                 paste0(rscript_exec, " ", shQuote(file.path(job_dir, "run_chunk.R")), " \"${SLURM_ARRAY_TASK_ID}\""))
  job_merge <- c(header(paste0(cluster_job_name, "_merge"), FALSE, "slurm_merge_%j", "END,FAIL"),
                 paste0(rscript_exec, " ", shQuote(file.path(job_dir, "run_merge.R"))))
  submit <- c(
    "#!/bin/bash",
    paste0("# Submits the array of ", n_chunks, " chunks and then the merge, which waits for every task"),
    "# to end without error. Written by generate_barcoding_job_scripts().",
    "set -euo pipefail",
    paste0("cd ", shQuote(job_dir)),
    "jid=$(sbatch --parsable job_array.sh | cut -d ';' -f 1)",
    "echo \"Array job: ${jid}\"",
    "sbatch --dependency=afterok:${jid} job_merge.sh")

  files <- file.path(job_dir, c("run_chunk.R", "run_merge.R", "job_array.sh", "job_merge.sh", "submit.sh"))
  writeLines(run_chunk, files[1])
  writeLines(run_merge, files[2])
  writeLines(job_array, files[3])
  writeLines(job_merge, files[4])
  writeLines(submit, files[5])
  Sys.chmod(files[3:5], mode = "0755")
  message("Job scripts written to ", job_dir, ". Submit on the cluster with: bash ", file.path(job_dir, "submit.sh"))
  invisible(files)
}
