# Decisions N1 to N4 of BMM (15-10, Phase 12, step 22c): the text a reader of the package meets (the
# vignettes, the help pages, the README, NEWS of this version, the comments of the tutorial scripts and
# the strings the functions print or write) holds no em-dash and no internal code of the audit
# (probes, phases, decisions, amendments, rule and figure codes): a control, a rule or a figure is
# named by what it does (N2). The values "CN1", "CN2" and "CN3" of the argument `controls` stay, as
# quoted values. Figure captions are drawn above the figure (N3) and sections are spaced, with a thin
# rule above each section and step (N4-b), by pkgdown/extra.css. Written before the text.

.n_root <- function(...) testthat::test_path("..", "..", ...)

.n_files <- function(dir, pattern) {
  d <- .n_root(dir)
  skip_if(!dir.exists(d), paste(dir, "is not in this tree"))
  list.files(d, pattern = pattern, full.names = TRUE)
}

# Internal codes of the audit. Percentiles (P5, P95), the Apple M2 chip, the schemes E and G, the
# calibration schemes S1 to S4 and the comparison sets A1 to A4 of Tutorial 4 are defined where they
# are used and are not codes of the audit.
.n_codes <- paste0(
  "\\b[Pp]robes?\\b|\\bPhase [0-9]+|\\b[Dd]ecisions? [A-Z0-9]|\\b[Aa]mendments? [A-Z0-9]|\\bBMM\\b|",
  "\\b(TB[0-9]|R-c|R2|K[0-9]|E[0-9]{1,2}|O[0-9]|I[0-9]|IR[0-9]|J[0-9][a-z]?|C[0-9]|DB[0-9]|D-[ab]|T[0-9][a-z]?(-[0-9a-z]+)?|",
  "W[0-9]|IP[0-9]|L[0-9]|X[0-9]|6[A-E]|5[AB])\\b|\\bCN[1-3]\\b")

# Lines with a hit, after removing the quoted values of `controls`
.n_hits <- function(x) {
  y <- gsub("\\\\?[\"']CN[1-3]\\\\?[\"']", "\"\"", x)
  y <- gsub("Apple M2", "Apple", y, fixed = TRUE)
  x[grepl(.n_codes, y, perl = TRUE)]
}

.n_strings <- function(f) {
  x <- readLines(f, warn = FALSE)
  x <- x[!grepl("^\\s*#", x)]
  unlist(regmatches(x, gregexpr("\"([^\"\\\\]|\\\\.)*\"", x)))
}

test_that("no em-dash in the text the package distributes (N1)", {
  dash <- intToUtf8(0x2014)
  files <- c(.n_files("vignettes", "\\.Rmd$"), .n_files("man", "\\.Rd$"), .n_files("R", "\\.R$"),
             .n_files(file.path("inst", "scripts"), "\\.R$"), .n_root("README.md"), .n_root("NEWS.md"))
  for (f in files[file.exists(files)]) {
    x <- readLines(f, warn = FALSE, encoding = "UTF-8")
    expect_false(any(grepl(dash, x, fixed = TRUE)), info = basename(f))
  }
})

test_that("the vignettes and the help pages name controls, rules and figures by what they do (N2)", {
  for (f in c(.n_files("vignettes", "\\.Rmd$"), .n_files("man", "\\.Rd$"), .n_root("README.md"))) {
    h <- .n_hits(readLines(f, warn = FALSE))
    expect_length(h, 0L)
    if (length(h)) fail(paste0(basename(f), ": ", paste(substr(trimws(h), 1, 90), collapse = " | ")))
  }
  news <- readLines(.n_root("NEWS.md"), warn = FALSE)
  this <- news[seq_len(grep("^# PhyloCactus 0\\.4\\.5", news)[1] - 1L)]
  expect_length(.n_hits(this), 0L)
})

test_that("the comments of the tutorial scripts hold no internal code (N2)", {
  for (f in .n_files(file.path("inst", "scripts"), "^tutorial-.*\\.R$")) {
    x <- readLines(f, warn = FALSE)
    # The figure codes of Tutorial 5 name the files the script writes (Figure_TB1_..., Figure_TB9_...)
    h <- .n_hits(gsub("\\bTB[0-9]\\b", "", x[grepl("^\\s*#", x)]))
    expect_length(h, 0L)
    if (length(h)) fail(paste0(basename(f), ": ", paste(substr(trimws(h), 1, 90), collapse = " | ")))
  }
})

test_that("the strings the functions print or write hold no internal code (N2)", {
  for (f in .n_files("R", "\\.R$")) {
    s <- .n_strings(f)
    s <- s[!s %in% c("\"CN1\"", "\"CN2\"", "\"CN3\"")]
    h <- .n_hits(gsub("\\bP(5|95)\\b", "", s))
    expect_length(h, 0L)
    if (length(h)) fail(paste0(basename(f), ": ", paste(substr(h, 1, 90), collapse = " | ")))
  }
})

test_that("pkgdown/extra.css draws captions above figures and spaces the sections with a rule (N3, N4-b)", {
  f <- .n_root("pkgdown", "extra.css")
  expect_true(file.exists(f))
  skip_if(!file.exists(f))
  css <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_match(css, "\\.figure[^{]*\\{[^}]*flex-direction:\\s*column", perl = TRUE)
  expect_match(css, "(\\.caption|figcaption)[^{]*\\{[^}]*order:\\s*-1", perl = TRUE)
  for (h in c("h2", "h3")) {
    rule <- regmatches(css, regexpr(paste0("[^}]*\\b", h, "\\b[^{]*\\{[^}]*\\}"), css, perl = TRUE))
    expect_length(rule, 1L)
    expect_match(rule, "margin-top", fixed = TRUE, info = h)
    expect_match(rule, "border-top:\\s*1px solid #E5E5E5", perl = TRUE, info = h)
  }
})
