# Decision N2 of BMM (15-10, Phase 12, step 22c): the strings the functions print or write hold no
# internal code of the audit (probes, phases, decisions, amendments, rule and figure codes): a
# control, a rule or a figure is named by what it does. The values "CN1", "CN2" and "CN3" of the
# argument `controls` stay, as quoted values. Written before the text.

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

test_that("the strings the functions print or write hold no internal code (N2)", {
  for (f in .n_files("R", "\\.R$")) {
    s <- .n_strings(f)
    s <- s[!s %in% c("\"CN1\"", "\"CN2\"", "\"CN3\"")]
    h <- .n_hits(gsub("\\bP(5|95)\\b", "", s))
    expect_length(h, 0L)
    if (length(h)) fail(paste0(basename(f), ": ", paste(substr(h, 1, 90), collapse = " | ")))
  }
})
