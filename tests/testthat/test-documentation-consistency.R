# Decision DB3 (c), option c1 of BMM (15-10, Phase 12, D-b): the figures the text of the tutorials
# states are read from the files distributed in inst/extdata, so that the text cannot drift from them
# in later versions. The numbers are taken from the sentences of the vignettes, not written here. Red
# until step 23 publishes the outputs of the Phase 12 run; green after it.

.dc_vignette <- function(name) {
  f <- system.file("doc", paste0(name, ".Rmd"), package = "PhyloCactus")
  if (!nzchar(f)) f <- testthat::test_path("..", "..", "vignettes", paste0(name, ".Rmd"))
  skip_if(!file.exists(f), paste(name, "not found"))
  paste(readLines(f, warn = FALSE), collapse = " ")
}

.dc_extdata <- function(name) {
  f <- system.file("extdata", name, package = "PhyloCactus")
  skip_if(!nzchar(f), paste(name, "not installed"))
  f
}

# The integers captured by the groups of a pattern in a text, thousands separators dropped
.dc_numbers <- function(text, pattern) {
  m <- regmatches(text, regexec(pattern, text))[[1]]
  expect_true(length(m) > 1L, info = pattern)
  as.integer(gsub(",", "", m[-1]))
}

test_that("Tutorial 2 states the bipartitions of the distributed ML and constraint trees", {
  skip_if_not_installed("ape")
  v <- .dc_vignette("tutorial-2-cactus-phylogeny-inference")
  n <- .dc_numbers(v, "fixes ([0-9,]+) of the ([0-9,]+) internal bipartitions")
  ml <- ape::unroot(ape::read.tree(.dc_extdata("phylocactus_ml_tree.tree")))
  cn <- ape::unroot(ape::read.tree(.dc_extdata("phylocactus_constraint_tree.tree")))
  expect_identical(ape::Nnode(ml) - 1L, n[2], info = "internal bipartitions of phylocactus_ml_tree.tree")
  expect_identical(ape::Nnode(cn) - 1L, n[1], info = "bipartitions of phylocactus_constraint_tree.tree")
  t <- .dc_numbers(v, "--mem; [^)]*, ([0-9,]+) terminals\\)")
  expect_identical(ape::Ntip(ml), t[1], info = "terminals of phylocactus_ml_tree.tree")
})

test_that("Tutorial 3 states the registry and the tree of the distributed species table", {
  v <- .dc_vignette("tutorial-3-cactus-phylogeny-visualization")
  sp <- utils::read.csv(.dc_extdata("phylocactus_table_species.csv"), stringsAsFactors = FALSE)
  skip_if(!"status_pipeline" %in% names(sp), "the distributed species table has no status_pipeline column")
  st <- table(sp$status_pipeline)
  n_of <- function(x) if (x %in% names(st)) as.integer(st[[x]]) else 0L
  expect_identical(nrow(sp), .dc_numbers(v, "the registry contains ([0-9,]+) species")[1], info = "registry")
  tips <- .dc_numbers(v, "have ([0-9,]+) tips \\(([0-9,]+) ingroup, ([0-9,]+) outgroup\\)")
  expect_identical(n_of("sampled_and_dated") + n_of("outgroup_reference"), tips[1], info = "tips")
  expect_identical(n_of("sampled_and_dated"), tips[2], info = "ingroup tips")
  expect_identical(n_of("outgroup_reference"), tips[3], info = "outgroup tips")
  expect_identical(n_of("sequenced_collapsed_duplicate"), .dc_numbers(v, "removed ([0-9,]+) ingroup species")[1],
                   info = "identical duplicates")
  expect_identical(n_of("unsampled_checklist_gap"), .dc_numbers(v, "The other ([0-9,]+) accepted species")[1],
                   info = "accepted species without sequences")
})
