# Version 0.5.1 (closing audit of 0.5.0, 16-10): the release describes what it is and what changed.
# DESCRIPTION names both branches and the external programs; NEWS records the changes of 0.5.0 that
# it left out; the cut of long records in the phylogeny branch is documented; Tutorial 2 reads the
# weak support of the backbone; the website deploy removes the pages of earlier builds; the reports
# of Tutorial 6 explain their two versions; the source of R/ carries no internal code of the audit
# and, with the help pages, no dated decision; the help pages hold no "additionally"; every line
# that names the checklist cites it. Written before the text.

.r51_root <- function(...) testthat::test_path("..", "..", ...)

.r51_read <- function(...) {
  f <- .r51_root(...)
  skip_if(!file.exists(f), paste(file.path(...), "is not in this tree"))
  paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

.r51_flat <- function(x) gsub("\\s+", " ", x)

test_that("the package is version 0.5.1 and NEWS opens with it", {
  expect_identical(as.character(utils::packageVersion("PhyloCactus")), "0.5.1")
  f <- .r51_root("NEWS.md")
  skip_if(!file.exists(f), "NEWS.md is not in this tree")
  expect_identical(readLines(f, n = 1L, warn = FALSE), "# PhyloCactus 0.5.1")
})

test_that("DESCRIPTION names both branches and the external programs (M1, m5)", {
  d <- read.dcf(.r51_root("DESCRIPTION"))
  title <- d[1, "Title"]
  desc <- .r51_flat(d[1, "Description"])
  expect_match(title, "Identification", fixed = TRUE)
  expect_match(desc, "identification", fixed = TRUE)
  expect_match(desc, "DNA barcoding", fixed = TRUE)
  expect_false(grepl("across bootstrap replicates", desc, fixed = TRUE))
  urls <- trimws(strsplit(d[1, "URL"], ",")[[1]])
  expect_identical(anyDuplicated(sub("/$", "", urls)), 0L)
  expect_true("SystemRequirements" %in% colnames(d))
  req <- .r51_flat(d[1, "SystemRequirements"])
  for (p in c("MAFFT", "RAxML-NG", "ModelTest-NG", "treePL", "BLAST+", "minimap2", "samtools", "fastp", "GetOrganelle")) {
    expect_match(req, p, fixed = TRUE, info = p)
  }
})

test_that("NEWS of 0.5.0 records the changes it left out (M2)", {
  news <- readLines(.r51_root("NEWS.md"), warn = FALSE, encoding = "UTF-8")
  from <- grep("^# PhyloCactus 0\\.5\\.0", news)
  to <- grep("^# PhyloCactus 0\\.4\\.5", news)
  expect_length(from, 1L)
  s <- .r51_flat(paste(news[from:(to - 1L)], collapse = "\n"))
  for (k in c("1208", "13 700", "13 partitions", "1e-05", "validate_phylogenies()", "theme_phylocactus()",
              "phylocactus_colours()", "phylocactus_panel_tags()", "extract_barcoding_plastome_loci()",
              "barcoding_outgroup_taxids()", "workers", "expected_tips", "cactus_constraints.csv",
              "assemble_ingroup_phylotar()", "2000 bases")) {
    expect_true(grepl(k, s, fixed = TRUE), info = k)
  }
})

test_that("the cut of long records in the phylogeny branch is documented (M3)", {
  for (rd in c("assemble_ingroup_phylotar.Rd", "assemble_outgroup_phylotar.Rd")) {
    x <- .r51_flat(.r51_read("man", rd))
    expect_match(x, "2000 bases", fixed = TRUE, info = rd)
  }
  expect_match(.r51_flat(.r51_read("man", "assemble_ingroup_phylotar.Rd")), "TABLE_LONG_RECORDS_CUT_INGROUP.csv", fixed = TRUE)
  expect_match(.r51_flat(.r51_read("man", "assemble_outgroup_phylotar.Rd")), "TABLE_LONG_RECORDS_CUT_OUTGROUP.csv", fixed = TRUE)
  t1 <- .r51_flat(.r51_read("vignettes", "tutorial-1-cactus-phylogeny-prep.Rmd"))
  expect_match(t1, "2000 bases", fixed = TRUE)
  expect_match(t1, "TABLE_LONG_RECORDS_CUT_INGROUP.csv", fixed = TRUE)
})

test_that("Tutorial 2 reads the weak support of the backbone (M5)", {
  t2 <- .r51_flat(.r51_read("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd"))
  expect_match(t2, "soft polytomies", fixed = TRUE)
  expect_match(t2, "76.9%", fixed = TRUE)
  expect_match(t2, "single locus", fixed = TRUE)
})

test_that("the website deploy removes the pages of earlier builds (M6)", {
  f <- .r51_root(".github", "workflows", "pkgdown.yaml")
  skip_if(!file.exists(f), "the workflows are not in this tree")
  y <- readLines(f, warn = FALSE)
  expect_true(any(grepl("^\\s*clean:\\s*true\\s*$", y)))
  expect_false(any(grepl("^\\s*clean:\\s*false\\s*$", y)))
})

test_that("Tutorial 6 explains the two versions of its reports (m1)", {
  t6 <- .r51_flat(.r51_read("vignettes", "tutorial-6-cactus-barcoding-identification.Rmd"))
  expect_match(t6, "0.4.5.9000", fixed = TRUE)
})

test_that("the source of R/ carries no internal code of the audit (m2)", {
  d <- .r51_root("R")
  skip_if(!dir.exists(d), "R is not in this tree")
  # The codes of test-narrative-22c.R, applied here to the comments and the internal documentation
  codes <- paste0(
    "\\b[Pp]robes?\\b|\\bPhase [0-9]+|\\b[Dd]ecisions? [A-Z0-9]|\\b[Aa]mendments? [A-Z0-9]|\\bBMM\\b|",
    "\\b(TB[0-9]|R-c|R2|K[0-9]|E[0-9]{1,2}|O[0-9]|I[0-9]|IR[0-9]|J[0-9][a-z]?|C[0-9]|DB[0-9]|D-[ab]|T[0-9][a-z]?(-[0-9a-z]+)?|",
    "W[0-9]|IP[0-9]|L[0-9]|X[0-9]|6[A-E]|5[AB])\\b|\\bCN[1-3]\\b|audit repository")
  for (f in list.files(d, "\\.R$", full.names = TRUE)) {
    x <- readLines(f, warn = FALSE)
    cm <- sub("^[^#]*#", "#", x[grepl("#", x)])
    cm <- gsub("\\\\?[\"']CN[1-3]\\\\?[\"']", "", cm)
    cm <- gsub("Apple M2", "Apple", cm, fixed = TRUE)
    cm <- gsub("\\bP(5|95)\\b", "", cm)
    h <- cm[grepl(codes, cm, perl = TRUE)]
    expect_length(h, 0L)
    if (length(h)) fail(paste0(basename(f), ": ", paste(substr(trimws(h), 1, 90), collapse = " | ")))
  }
  expect_false(any(grepl("1024 tips", readLines(file.path(d, "inference.R"), warn = FALSE), fixed = TRUE)))
})

test_that("the source of R/ and the help pages date no decision (m2)", {
  dated <- "[Dd]ecisions? of 20[0-9]{2}|[Dd]ecided on 20[0-9]{2}|\\b[Uu]ntil 20[0-9]{2}-|\\b[Ss]ince 20[0-9]{2}-"
  for (dir in c("R", "man")) {
    d <- .r51_root(dir)
    skip_if(!dir.exists(d), paste(dir, "is not in this tree"))
    for (f in list.files(d, "\\.(R|Rd)$", full.names = TRUE)) {
      x <- readLines(f, warn = FALSE)
      if (dir == "R") x <- sub("^[^#]*#", "#", x[grepl("#", x)])
      h <- x[grepl(dated, x, perl = TRUE)]
      expect_length(h, 0L)
      if (length(h)) fail(paste0(basename(f), ": ", paste(substr(trimws(h), 1, 90), collapse = " | ")))
    }
  }
})

test_that("the help pages hold no 'additionally' (m3)", {
  d <- .r51_root("man")
  skip_if(!dir.exists(d), "man is not in this tree")
  hits <- Filter(function(f) any(grepl("additionally", readLines(f, warn = FALSE), ignore.case = TRUE)),
                 list.files(d, "\\.Rd$", full.names = TRUE))
  expect_length(hits, 0L)
})

test_that("every prose line that names the checklist cites it (m4)", {
  files <- c(.r51_root("README.md"), list.files(.r51_root("vignettes"), "\\.Rmd$", full.names = TRUE))
  skip_if(!all(file.exists(files[1])), "README.md is not in this tree")
  for (f in files) {
    x <- readLines(f, warn = FALSE, encoding = "UTF-8")
    chunk <- cumsum(grepl("^```", x)) %% 2 == 1 | grepl("^```", x)
    prose <- x[!chunk]
    h <- prose[grepl("checklist|Caryophyllales\\.org|CactaceaeFullList", prose, ignore.case = TRUE) &
                 !grepl("Korotkova", prose, fixed = TRUE)]
    expect_length(h, 0L)
    if (length(h)) fail(paste0(basename(f), ": ", paste(substr(trimws(h), 1, 90), collapse = " | ")))
  }
})

test_that("the homology check states the figures of the reference run of 0.5.0", {
  # TABLE_marker_homology_check.csv of the reference run: the markers with verdict ok (rpL16,
  # trnL_trnF, matK, rbcL) share 0.26 to 0.58 of their outgroup 20-mers; trnT-psbD no longer
  # enters the matrix. The help page, the internal documentation and the warning say so.
  for (f in c(.r51_root("R", "screening.R"), .r51_root("man", "integrate_and_clean_markers.Rd"))) {
    skip_if(!file.exists(f), paste(basename(f), "is not in this tree"))
    x <- .r51_flat(paste(readLines(f, warn = FALSE), collapse = "\n"))
    expect_false(grepl("0.26 to 0.59", x, fixed = TRUE), info = basename(f))
    expect_false(grepl("returns 0.040", x, fixed = TRUE), info = basename(f))
    expect_match(x, "0.26 to 0.58", fixed = TRUE, info = basename(f))
  }
})
