# Decision T6-6 of BMM (14-10, Phase 12): the files the package distributes (inst/extdata) and the
# identification reports published on the website (pkgdown/assets/reports) carry no path of the
# machine that made them and no host name, and the examples of Tutorial 6 are public records named by
# their accession. The run record of an identification stores the host name (Sys.info() nodename) and
# the report prints it, so the reports of the website are regenerated without it. Written before the
# reports and the example tables exist.

.df_extdata <- function() {
  d <- system.file("extdata", package = "PhyloCactus")
  skip_if(!nzchar(d), "extdata not installed")
  d
}

.df_reports <- function() {
  d <- testthat::test_path("..", "..", "pkgdown", "assets", "reports")
  skip_if(!dir.exists(d), "pkgdown/assets/reports is not in this tree (it is not part of the built package)")
  d
}

.df_text_files <- function(d) {
  f <- list.files(d, recursive = TRUE, full.names = TRUE)
  f[grepl("\\.(csv|tsv|txt|tree|gb|html|Rmd|R|md)$", f, ignore.case = TRUE)]
}

.df_local_paths <- function(files) {
  bad <- character(0)
  for (f in files) {
    x <- readLines(f, warn = FALSE)
    if (any(grepl("/Users/|/home/|[A-Za-z]:\\\\\\\\Users", x))) bad <- c(bad, basename(f))
  }
  bad
}

.df_accession <- "^([A-Z]{1,2}_?[0-9]{5,8}\\.[0-9]+|GC[AF]_[0-9]{9}\\.[0-9]+|[SEDR]RR[0-9]{6,9})$"

test_that("no file of extdata holds a path of a local machine", {
  expect_identical(.df_local_paths(.df_text_files(.df_extdata())), character(0))
})

test_that("the example tables of Tutorial 6 are named by the accession of a public record", {
  ex <- list.files(.df_extdata(), pattern = "^barcoding_example_.*_(answer|sampling)\\.csv$")
  skip_if(length(ex) == 0L, "no example table of Tutorial 6 yet")
  acc <- sub("^barcoding_example_(.*)_(answer|sampling)\\.csv$", "\\1", ex)
  expect_true(all(grepl(.df_accession, acc)), info = paste(acc[!grepl(.df_accession, acc)], collapse = ", "))
  # Each example has both tables
  a <- sub("_answer\\.csv$", "", ex[grepl("_answer\\.csv$", ex)]); s <- sub("_sampling\\.csv$", "", ex[grepl("_sampling\\.csv$", ex)])
  expect_setequal(a, s)
})

test_that("every query of the example tables is a public accession, and no own sample is among them", {
  ex <- list.files(.df_extdata(), pattern = "^barcoding_example_.*_answer\\.csv$", full.names = TRUE)
  skip_if(length(ex) == 0L, "no example table of Tutorial 6 yet")
  q <- unique(unlist(lapply(ex, function(f) utils::read.csv(f, stringsAsFactors = FALSE)$query)))
  # The regions of an assembly carry the run in their name (run__locus); the contigs that GetOrganelle
  # assembles from reads are named by their organelle (plastome_<n>, nrdna_<n>)
  q <- sub("__.*$", "", q)
  q <- q[!grepl("^(plastome|nrdna)_[0-9]+$", q)]
  expect_true(all(grepl(.df_accession, q)), info = paste(q[!grepl(.df_accession, q)], collapse = ", "))
})

test_that("no example table carries the run record or the provenance of the machine", {
  f <- list.files(.df_extdata(), pattern = "^barcoding_example_")
  expect_false(any(grepl("(^|_)(RUN|provenance|checksums)", f)))
})

test_that("the reports of the website hold no local path and no host name", {
  d <- .df_reports()
  html <- list.files(d, pattern = "\\.html$", full.names = TRUE)
  skip_if(length(html) == 0L, "no report on the website yet")
  expect_identical(.df_local_paths(html), character(0))
  for (f in html) {
    x <- paste(readLines(f, warn = FALSE), collapse = "\n")
    expect_true(grepl("machine not recorded", x, fixed = TRUE), info = basename(f))
    expect_false(grepl("\\.local\\b", x), info = basename(f))
  }
})
