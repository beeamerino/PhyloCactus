# AP-07 and AP-25 (closure review of 0.5.4): a metadata batch whose retries all failed was dropped
# without a message, and cp_normalize_marker() replaced a hyphen with a hyphen.

test_that("failed metadata batches are counted and stop the call above max_failed_fraction", {
  dir <- withr::local_tempdir()
  cache <- file.path(dir, "meta.csv")
  calls <- 0L
  testthat::local_mocked_bindings(
    read.GenBank = function(access.nb, ...) {
      calls <<- calls + 1L
      if (any(access.nb == "BAD1")) stop("HTTP 500")
      x <- structure(as.list(access.nb), names = access.nb)
      attr(x, "species") <- rep("Opuntia_ficus-indica", length(access.nb))
      attr(x, "description") <- rep("trnL-trnF", length(access.nb))
      x
    },
    .package = "ape"
  )
  logged <- character(0)
  lg <- function(...) logged <<- c(logged, paste0(...))
  expect_error(
    PhyloCactus:::cp_download_all_metadata(c("A1", "A2", "BAD1"), cache, batch_size = 1L, sleep_time = 0,
                                           max_retries = 2L, log_message = lg),
    "1 of 3 accessions"
  )
  expect_true(any(grepl("unrecovered after 2 attempts", logged)))
  # The batches that worked are cached, and below the fraction the call returns with the failures.
  expect_setequal(utils::read.csv(cache)$sid, c("A1", "A2"))
  out <- PhyloCactus:::cp_download_all_metadata(c("A1", "A2", "BAD1"), cache, batch_size = 1L, sleep_time = 0,
                                                max_retries = 1L, max_failed_fraction = 0.5)
  expect_identical(attr(out, "failed_sids"), "BAD1")
  expect_warning(
    suppressMessages(fetch_genbank_metadata("A1", cache, force_download = TRUE)),
    "deprecated"
  )
})

test_that("typographic dashes in marker names become hyphens", {
  expect_identical(PhyloCactus:::cp_normalize_marker(c("trnL–trnF", "psbA‐trnH", " rbcL ")),
                   c("trnL-trnF", "psbA-trnH", "rbcL"))
})
