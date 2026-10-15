# Third narrative pass of step 22c (BMM, 15-10): Tutorials 3 and 7 and the README read in full.
# Tutorial 3 is the third of four stages, not the last; Figure 12 plots the share of loci, not their
# number; the IUCN data are cited with the Red List version of the run (2026-1, extracted 15-10).
# Tutorial 7 describes automate_treePL() as Tutorial 2 does (parameters chosen on the ML tree, applied
# to the temporal replicates) and gives FBP on the 0 to 1 scale used elsewhere. The README gives its
# references in the short style of the tutorials. Written before the text.

.p3_read <- function(...) {
  f <- testthat::test_path("..", "..", ...)
  skip_if(!file.exists(f), paste(file.path(...), "is not in this tree"))
  gsub("\\s+", " ", paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
}

test_that("Tutorial 3 places itself in the workflow, describes Figure 12 and cites the Red List version", {
  t3 <- .p3_read("vignettes", "tutorial-3-cactus-phylogeny-visualization.Rmd")
  for (s in c("the last stage of the core workflow", "Data Setup & Pre-processing",
              "Number of loci per species in each IUCN category"))
    expect_false(grepl(s, t3, fixed = TRUE), info = s)
  expect_true(grepl("the third stage of the phylogenetic workflow", t3, fixed = TRUE))
  expect_true(grepl("Share of the loci sampled for each species in each IUCN category", t3, fixed = TRUE))
  expect_true(grepl("The reference run used Red List version 2026-1", t3, fixed = TRUE))
  refs <- sub("^.*## References", "", t3)
  expect_true(grepl("- IUCN 2026. *The IUCN Red List of Threatened Species*. Version 2026-1.", refs, fixed = TRUE))
})

test_that("Tutorial 7 describes automate_treePL() and the collapse of nodes as the tutorials do", {
  t7 <- .p3_read("vignettes", "tutorial-7-cactus-phylogeny-functions.Rmd")
  for (s in c("on the maximum-likelihood tree and on each bootstrap replicate", "FBP below 70", "Accessed 2026-09-19"))
    expect_false(grepl(s, t7, fixed = TRUE), info = s)
  expect_true(grepl("FBP below 0.70", t7, fixed = TRUE))
  expect_true(grepl("Accessed 2026-10-15", t7, fixed = TRUE))
  expect_true(grepl("applies them to every temporal bootstrap replicate", t7, fixed = TRUE))
})

test_that("the README gives its references in the short style of the tutorials", {
  f <- testthat::test_path("..", "..", "README.md")
  skip_if(!file.exists(f), "README.md is not in this tree")
  x <- readLines(f, warn = FALSE, encoding = "UTF-8")
  refs <- x[seq(grep("^## References", x), grep("^## License", x) - 1L)]
  refs <- refs[startsWith(refs, "- ")]
  expect_gt(length(refs), 5L)
  expect_false(any(grepl("^- [A-Z][^ ]*, [A-Z]\\.", refs)), info = "full author lists")
  expect_false(any(grepl("[0-9]-[0-9]+\\. <", refs)), info = "page ranges with a hyphen")
  expect_true(any(grepl("Cactaceae at Caryophyllales.org - A dynamic online", refs, fixed = TRUE)))
})
