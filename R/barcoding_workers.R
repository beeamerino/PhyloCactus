#' Step 7 in parallel on one machine (decisions PA1 to PA3 of BMM, 29-09)
#'
#' Starts `workers` R processes, one chunk of the 6B chunks each, writes the progress to
#' `PROGRESS_<suffix>.txt` of `output_dir` every `getOption("PhyloCactus.poll_seconds", 60)`
#' seconds, merges the chunks with [merge_barcoding_chunks()] when every worker has finished, and
#' sends one email at the end. A worker that fails stops the others and the run. On macOS the machine
#' is kept awake with `caffeinate` while the call lasts (a flag file in `chunks/` removed on exit),
#' when `caffeinate` exists.
#' @param args Named list of the arguments of [classify_barcoding_folds()] passed to every worker.
#' @noRd
.bc_classify_workers <- function(args, workers, suffix, notify, notify_to, notify_credentials) {
  started <- Sys.time()
  out_dir <- args$output_dir
  chunk_dir <- file.path(out_dir, "chunks")
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
  progress_file <- file.path(out_dir, paste0("PROGRESS_", suffix, ".txt"))
  rscript <- file.path(R.home("bin"), "Rscript")
  for (k in seq_len(workers)) {
    unlink(file.path(chunk_dir, paste0(c("STATUS_", "PID_"), suffix, "_worker", k, "of", workers, ".txt")))
  }
  saveRDS(args, file.path(chunk_dir, paste0("ARGS_", suffix, ".rds")))

  # One script per worker: the chunk, then a status file (done, or failed with the error)
  for (k in seq_len(workers)) {
    tag <- paste0(suffix, "_worker", k, "of", workers)
    script <- file.path(chunk_dir, paste0("RUN_", tag, ".R"))
    writeLines(c(
      sprintf('writeLines(as.character(Sys.getpid()), "%s")', file.path(chunk_dir, paste0("PID_", tag, ".txt"))),
      sprintf('a <- readRDS("%s")', file.path(chunk_dir, paste0("ARGS_", suffix, ".rds"))),
      sprintf('a$chunk <- %dL; a$n_chunks <- %dL', k, workers),
      sprintf('st <- "%s"', file.path(chunk_dir, paste0("STATUS_", tag, ".txt"))),
      'r <- tryCatch({',
      sprintf('  if (identical(Sys.getenv("PHYLOCACTUS_TEST_FAIL_CHUNK"), "%d")) stop("failure asked by PHYLOCACTUS_TEST_FAIL_CHUNK")', k),
      '  suppressPackageStartupMessages(library(PhyloCactus)); do.call(classify_barcoding_folds, a); "done"',
      '}, error = function(e) c("failed", conditionMessage(e)))',
      'writeLines(r, st)',
      'if (r[1] == "failed") quit(status = 1)'), script)
    system2(rscript, shQuote(script), stdout = file.path(chunk_dir, paste0("LOG_", tag, ".txt")),
            stderr = file.path(chunk_dir, paste0("LOG_", tag, ".txt")), wait = FALSE)
  }
  # Awake while the run lasts: caffeinate holds a shell that waits on a flag removed when this call
  # ends, finished, failed or interrupted
  flag <- file.path(chunk_dir, paste0("RUNNING_", suffix))
  writeLines(format(started), flag)
  on.exit(unlink(flag), add = TRUE)
  if (nzchar(Sys.which("caffeinate"))) {
    system2("caffeinate", c("-i", "sh", "-c", shQuote(sprintf("while [ -f %s ]; do sleep 30; done", shQuote(normalizePath(flag))))),
            wait = FALSE, stdout = FALSE, stderr = FALSE)
  }

  status_of <- function(k) {
    f <- file.path(chunk_dir, paste0("STATUS_", suffix, "_worker", k, "of", workers, ".txt"))
    if (file.exists(f)) readLines(f, warn = FALSE) else "running"
  }
  steps_done <- function(k) {
    f <- file.path(chunk_dir, paste0("LOG_", suffix, "_worker", k, "of", workers, ".txt"))
    if (!file.exists(f)) return(0L)
    sum(grepl("^Locus '.*predictions written", readLines(f, warn = FALSE)))
  }
  report <- function(state) {
    st <- vapply(seq_len(workers), function(k) status_of(k)[1], "")
    steps <- sum(vapply(seq_len(workers), steps_done, 0L))
    lines <- c(sprintf("Step 7, method %s, %d workers: %s", suffix, workers, state),
               sprintf("Started %s; now %s (%s)", format(started, "%Y-%m-%d %H:%M:%S"), format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                       .format_elapsed(difftime(Sys.time(), started, units = "secs"))),
               sprintf("%d of %d workers finished; %d failed; %d locus and scheme steps written",
                       sum(st == "done"), workers, sum(st == "failed"), steps),
               sprintf("  worker %d: %s, %d steps", seq_len(workers), st, vapply(seq_len(workers), steps_done, 0L)))
    writeLines(lines, progress_file)
    message(lines[3])
    st
  }
  wait <- getOption("PhyloCactus.poll_seconds", 60)
  kill_all <- function() {
    for (k in seq_len(workers)) {
      pf <- file.path(chunk_dir, paste0("PID_", suffix, "_worker", k, "of", workers, ".txt"))
      if (file.exists(pf)) try(tools::pskill(as.integer(readLines(pf, warn = FALSE)[1])), silent = TRUE)
    }
  }
  fail <- function(msg) {
    kill_all()
    report("failed")
    if (isTRUE(notify)) {
      note <- .compose_run_notification(analysis = paste0("barcoding classification of the folds, method ", suffix, ", ", workers, " workers"),
                                        status = "failed", started = started,
                                        outputs = c("Progress" = progress_file, "Worker logs" = chunk_dir),
                                        error_message = msg)
      send_run_notification(note$subject, note$body, to = notify_to, credentials = notify_credentials)
    }
    stop(msg, call. = FALSE)
  }

  repeat {
    Sys.sleep(wait)
    st <- report("running")
    if (any(st == "failed")) {
      k <- which(st == "failed")[1]
      log <- file.path(chunk_dir, paste0("LOG_", suffix, "_worker", k, "of", workers, ".txt"))
      tail_log <- if (file.exists(log)) utils::tail(readLines(log, warn = FALSE), 5) else character(0)
      fail(paste0("Step 7 stopped: chunk ", k, " of ", workers, " failed: ", paste(status_of(k)[-1], collapse = " "),
                  if (length(tail_log)) paste0("\nLast lines of its log:\n", paste(tail_log, collapse = "\n")) else ""))
    }
    if (all(st == "done")) break
  }
  res <- tryCatch(suppressMessages(merge_barcoding_chunks(library_dir = args$library_dir, folds_dir = args$folds_dir,
                                                          output_dir = out_dir, method = args$method, n_chunks = workers,
                                                          alignment = args$alignment, schemes = args$schemes, loci = args$loci,
                                                          max_folds = args$max_folds, seed = args$seed)),
                  error = function(e) e)
  if (inherits(res, "error")) fail(paste0("Step 7 stopped: the merge of the chunks failed: ", conditionMessage(res)))
  report("finished")
  if (isTRUE(notify)) {
    note <- .compose_run_notification(analysis = paste0("barcoding classification of the folds, method ", suffix, ", ", workers, " workers"),
                                      status = "finished", started = started,
                                      outputs = c(stats::setNames(file.path(out_dir, paste0("TABLE_barcoding_predictions_", args$schemes, "_", suffix, ".csv")),
                                                                  paste0("Predictions, scheme ", args$schemes)),
                                                  "Progress" = progress_file))
    send_run_notification(note$subject, note$body, to = notify_to, credentials = notify_credentials)
  }
  invisible(res)
}
