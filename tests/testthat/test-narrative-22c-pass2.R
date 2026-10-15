# Second narrative pass of step 22c (BMM, 15-10): the reader-facing text describes the package and
# its results, not its development. Absent: prescriptions to the reader's Methods, development history,
# the state of the package in earlier versions, the jargon of the audit ("folder 0", "frozen"), the
# construction "rather than", and the figures of 0.4.5 left in Tutorials 1 and 2 and in the help of the
# search and bootstrap functions. Present: the figures of the Phase 12 run (log of the ML search,
# 4468 s, 1208 taxa, 13 partitions, 7899 patterns; bootstopping after 550 replicates; trnT-psbD with
# 57 ingroup sequences) and the closing of Tutorial 4 pointing to the identification branch.
# Written before the text.

.p2_root <- function(...) testthat::test_path("..", "..", ...)

.p2_read <- function(...) {
  f <- .p2_root(...)
  skip_if(!file.exists(f), paste(file.path(...), "is not in this tree"))
  paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

.p2_texts <- function() {
  d <- .p2_root("vignettes"); m <- .p2_root("man")
  skip_if(!dir.exists(d) || !dir.exists(m), "vignettes or man is not in this tree")
  files <- c(list.files(d, "\\.Rmd$", full.names = TRUE), list.files(m, "\\.Rd$", full.names = TRUE), .p2_root("README.md"))
  txt <- vapply(files, function(f) paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n"), character(1))
  news <- readLines(.p2_root("NEWS.md"), warn = FALSE, encoding = "UTF-8")
  txt[["NEWS.md (0.5.0)"]] <- paste(news[seq_len(grep("^# PhyloCactus 0\\.4\\.5", news)[1] - 1L)], collapse = "\n")
  # Line breaks of the help pages and of wrapped paragraphs do not hide a phrase
  vapply(txt, function(x) gsub("\\s+", " ", x), character(1))
}

test_that("the reader-facing text holds no development narrative (second pass of N1)", {
  txt <- .p2_texts()
  stale <- c("rather than", "Methods should", "should be stated in Methods", "Recorded here the development history",
             "under construction", "Section Manifest", "conditions for reopening", "folder 0", "frozen",
             "The manifest of this section", "will be computed by the functions of that section",
             "In the real-file test", "currently available for")
  for (s in stale) {
    where <- names(txt)[grepl(s, txt, fixed = TRUE)]
    expect_length(where, 0L)
    if (length(where)) fail(paste0("\"", s, "\" in ", paste(basename(where), collapse = ", ")))
  }
})

test_that("Tutorials 1 and 2 and the help of the search and bootstrap functions state the Phase 12 run", {
  txt <- .p2_texts()
  old <- c("1024 taxa", "1024 terminals", "12809", "5959 patterns", "11 partitions", "2739 s", "650 of 1000",
           "possessed 50 ingroup sequences")
  for (s in old) {
    where <- names(txt)[grepl(s, txt, fixed = TRUE)]
    expect_length(where, 0L)
    if (length(where)) fail(paste0("\"", s, "\" in ", paste(basename(where), collapse = ", ")))
  }
  t2 <- gsub("\\s+", " ", .p2_read("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd"))
  expect_true(grepl("after 550 of 1000 replicates", t2, fixed = TRUE))
  expect_true(grepl("(1208 taxa, 13 700 sites, 7899 patterns, 13 partitions) took 4468 s", t2, fixed = TRUE))
  t1 <- gsub("\\s+", " ", .p2_read("vignettes", "tutorial-1-cactus-phylogeny-prep.Rmd"))
  expect_true(grepl("`trnT-psbD` held 57 ingroup sequences", t1, fixed = TRUE))
  for (rd in c("generate_bootstrap_script.Rd", "generate_ml_search_script.Rd")) {
    x <- gsub("\\s+", " ", .p2_read("man", rd))
    expect_true(grepl("1208", x, fixed = TRUE), info = rd)
  }
})

test_that("Tutorial 4 closes on the identification branch, as built", {
  t4 <- .p2_read("vignettes", "tutorial-4-cactus-phylogeny-validation.Rmd")
  expect_true(grepl("[Continue to Tutorial 5: Molecular Identification of Cacti (DNA Barcoding)]", t4, fixed = TRUE))
})
