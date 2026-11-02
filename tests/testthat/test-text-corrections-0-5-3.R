# Closure of 0.5.2 into 0.5.3, phase 4a (audit of 26-10-2026, findings M10 and n1): the manifest
# names no path outside the two repositories and no forbidden term. Written before the text.

.t53_extdata <- function(name) {
  f <- system.file("extdata", name, package = "PhyloCactus")
  skip_if(!nzchar(f), paste(name, "not installed"))
  f
}

test_that("the manifest names no path outside the two repositories and no forbidden term (M10, n1)", {
  m <- utils::read.csv(.t53_extdata("MANIFEST.csv"), stringsAsFactors = FALSE)
  txt <- unlist(m, use.names = FALSE)
  expect_false(any(grepl("plantlab-core", txt, fixed = TRUE)))
  expect_false(any(grepl("phylogenom", txt, ignore.case = TRUE)))
  expect_false(any(grepl("^(/|~)", txt)))
})
