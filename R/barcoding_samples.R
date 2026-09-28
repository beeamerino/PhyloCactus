# ------------------------------------------------------------------------------
# Molecular diagnostic branch (11_barcoding/), the sample sheet. PhyloCactus 0.5.0.
#
# A sheet with one row per sample, routed by the matrix of accepted inputs
# (inst/extdata/barcoding_accepted_inputs.csv): an admitted input goes to its identification function
# and gets its report; any other gets the message of its row. One index table for the whole sheet.
#
# Governing documents: 01_minutes_and_phases/04_barcoding_0.5.0/ of the audit repository, Phase 9
# proposal of 2026-09-28 (decisions IR7 and IR8 of BMM).
# ------------------------------------------------------------------------------

#' The matrix of accepted inputs of the package
#' @noRd
.bc_accepted_inputs <- function() {
  f <- system.file("extdata", "barcoding_accepted_inputs.csv", package = "PhyloCactus")
  if (!nzchar(f)) stop("barcoding_accepted_inputs.csv is not installed with the package.", call. = FALSE)
  utils::read.csv(f, stringsAsFactors = FALSE)
}

#' The answer of a sample against the species declared for it
#' @noRd
.bc_declared_comparison <- function(t, declared) {
  if (is.na(declared) || !nzchar(declared)) return("not declared")
  declared <- gsub("[[:space:]]+", "_", trimws(declared))
  sp <- unique(t$predicted_species[t$state %in% 1L & !is.na(t$predicted_species)])
  ge <- unique(t$predicted_genus[t$state %in% 1:2 & !is.na(t$predicted_genus)])
  if (declared %in% sp) return("species agrees")
  if (length(sp)) return("disagrees")
  if (sub("_.*$", "", declared) %in% ge) return("genus agrees")
  if (length(ge)) return("disagrees")
  "not assigned"
}

#' Identify the samples of a sample sheet
#'
#' Reads a sheet with one row per sample and routes each row by the matrix of accepted inputs of the
#' package (`inst/extdata/barcoding_accepted_inputs.csv`): an admitted input goes to
#' [identify_barcoding_query()] (`sanger`, `plastome`, `genbank`; the sequences of a row are one
#' sample), [identify_barcoding_assembly()] (`assembly`) or [identify_barcoding_reads()]
#' (`reads_paired`), and gets its report ([report_barcoding_identification()]) with the declared
#' species and the voucher; any other input is not run and gets the message of its row. An error in
#' one sample is recorded and the others go on. A sample without a voucher is identified and flagged:
#' without it the leakage rule per specimen cannot be applied to it.
#'
#' Columns of the sheet: `sample_id` (unique), `input_type` (a row of the matrix), `file_1`, and
#' optionally `file_2` (the second file of paired reads), `declared_species` and `voucher`.
#'
#' Writes one folder per identified sample in `output_dir` and `TABLE_barcoding_samples_index.csv`:
#' status (`identified`, `not_run`, `error`), message, loci found, answers per state, species and
#' genera named, the answer against the declared species (`species agrees`, `genus agrees`,
#' `disagrees`, `not assigned`, `not declared`), whether the voucher is missing, and the report.
#'
#' @param sheet Character or data frame. The sample sheet, a CSV path or a data frame.
#' @param library_dir,metrics_dir As in [identify_barcoding_query()].
#' @param output_dir Character. Where the samples and the index are written.
#' @param threads Integer. Threads of the external tools (assemblies and reads).
#' @return Invisibly, the index table.
#' @seealso [identify_barcoding_query()], [identify_barcoding_assembly()], [identify_barcoding_reads()].
#' @examples
#' \dontrun{
#' identify_barcoding_samples("samples.csv")
#' }
#' @export
identify_barcoding_samples <- function(sheet,
                                       library_dir = file.path("11_barcoding", "4_library"),
                                       metrics_dir = file.path("11_barcoding", "11_metrics"),
                                       output_dir = file.path("11_barcoding", "10_identify", "samples"),
                                       threads = 4L) {
  s <- if (is.data.frame(sheet)) sheet else utils::read.csv(sheet, stringsAsFactors = FALSE, colClasses = "character")
  need <- c("sample_id", "input_type", "file_1")
  miss <- setdiff(need, names(s))
  if (length(miss)) stop("The sample sheet lacks the column(s): ", paste(miss, collapse = ", "), ".", call. = FALSE)
  for (col in c("file_2", "declared_species", "voucher")) if (!col %in% names(s)) s[[col]] <- NA_character_
  s[] <- lapply(s, function(v) { v <- as.character(v); v[!is.na(v) & !nzchar(trimws(v))] <- NA; v })
  if (anyNA(s$sample_id) || anyDuplicated(s$sample_id)) stop("Every sample_id must be given and unique.", call. = FALSE)
  m <- .bc_accepted_inputs()
  bad <- setdiff(s$input_type, m$input_type)
  if (length(bad)) {
    stop("Unknown input_type: ", paste(bad, collapse = ", "), ". Valid: ", paste(m$input_type, collapse = ", "), ".", call. = FALSE)
  }
  .bc_assert_output_dir(output_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  rows <- list()
  for (i in seq_len(nrow(s))) {
    r <- s[i, ]
    row <- m[m$input_type == r$input_type, ]
    idx <- data.frame(sample_id = r$sample_id, input_type = r$input_type, status = "not_run", message = row$message,
                      loci_found = NA_integer_, state1 = NA_integer_, state2 = NA_integer_, state3 = NA_integer_,
                      species_named = NA_character_, genera_named = NA_character_, declared_species = r$declared_species,
                      declared_comparison = NA_character_, voucher_missing = is.na(r$voucher), report = NA_character_,
                      stringsAsFactors = FALSE)
    if (row$status != "admitted") {
      message("Sample ", r$sample_id, " (", r$input_type, "): not run. ", row$message)
      rows[[i]] <- idx
      next
    }
    res <- tryCatch({
      files <- c(r$file_1, if (identical(r$input_type, "reads_paired")) r$file_2)
      missing <- files[is.na(files) | !file.exists(files)]
      if (length(missing)) stop("File not found: ", paste(missing, collapse = ", "), ".", call. = FALSE)
      sdir <- file.path(output_dir, r$sample_id)
      invisible(utils::capture.output(tab <- switch(row$function_name,
        identify_barcoding_query = identify_barcoding_query(r$file_1, library_dir = library_dir, metrics_dir = metrics_dir,
                                                            output_dir = sdir, run_name = r$sample_id, sample = r$sample_id,
                                                            report = FALSE),
        identify_barcoding_assembly = identify_barcoding_assembly(r$file_1, library_dir = library_dir, metrics_dir = metrics_dir,
                                                                  output_dir = sdir, run_name = r$sample_id, threads = threads,
                                                                  report = FALSE),
        identify_barcoding_reads = identify_barcoding_reads(r$file_1, r$file_2, library_dir = library_dir, metrics_dir = metrics_dir,
                                                            output_dir = output_dir, run_name = r$sample_id, threads = threads,
                                                            report = FALSE))))
      f_rec <- file.path(sdir, paste0("RUN_", r$sample_id, ".csv"))
      rec <- utils::read.csv(f_rec, stringsAsFactors = FALSE, colClasses = "character")
      rec <- rec[!rec$name %in% c("sample", "declared_species", "voucher"), , drop = FALSE]
      rec <- rbind(rec, data.frame(kind = c("setting", "sample", "sample"), name = c("sample", "declared_species", "voucher"),
                                   value = c(r$sample_id, r$declared_species, r$voucher)))
      utils::write.csv(rec, f_rec, row.names = FALSE)
      p <- report_barcoding_identification(r$sample_id, results_dir = sdir, library_dir = library_dir, metrics_dir = metrics_dir)
      found <- !tab$reason %in% c("no_overlap", "assembly_failed", "single_species_library")
      idx$status <- "identified"; idx$message <- ""
      idx$loci_found <- sum(found)
      idx$state1 <- sum(tab$state == 1L); idx$state2 <- sum(tab$state == 2L); idx$state3 <- sum(tab$state == 3L)
      idx$species_named <- paste(unique(stats::na.omit(tab$predicted_species[tab$state == 1L])), collapse = "; ")
      idx$genera_named <- paste(unique(stats::na.omit(tab$predicted_genus[tab$state %in% 1:2])), collapse = "; ")
      idx$declared_comparison <- .bc_declared_comparison(tab, r$declared_species)
      idx$report <- p[1]
      idx
    }, error = function(e) {
      idx$status <- "error"; idx$message <- conditionMessage(e)
      message("Sample ", r$sample_id, ": error. ", conditionMessage(e))
      idx
    })
    rows[[i]] <- res
  }
  index <- do.call(rbind, rows)
  rownames(index) <- NULL
  f <- file.path(output_dir, "TABLE_barcoding_samples_index.csv")
  utils::write.csv(index, f, row.names = FALSE)
  message("Samples: ", sum(index$status == "identified"), " identified, ", sum(index$status == "not_run"), " not run, ",
          sum(index$status == "error"), " with an error. Index: ", f)
  invisible(index)
}
