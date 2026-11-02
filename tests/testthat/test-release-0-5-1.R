# Version 0.5.1 (closing audit of 0.5.0, 16-10): the website deploy removes the pages of earlier
# builds; the source of R/ carries no internal code of the audit and, with the help pages, no dated
# decision. Written before the text.

.r51_root <- function(...) testthat::test_path("..", "..", ...)

test_that("the website deploy removes the pages of earlier builds (M6)", {
  f <- .r51_root(".github", "workflows", "pkgdown.yaml")
  skip_if(!file.exists(f), "the workflows are not in this tree")
  y <- readLines(f, warn = FALSE)
  expect_true(any(grepl("^\\s*clean:\\s*true\\s*$", y)))
  expect_false(any(grepl("^\\s*clean:\\s*false\\s*$", y)))
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
