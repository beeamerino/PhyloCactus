# M5 (Phase D, decision D5): one GenBank plastome distributed as the example query of Tutorial 5.
# MN517611.1, Mammillaria zephyranthoides, is example X1 of the measured examples of 07-10: identified
# to species from matK against the library of this version.

example_plastome <- "barcoding_example_plastome_MN517611.1.gb"

test_that("the example plastome is distributed and declared as an external input", {
  path <- system.file("extdata", example_plastome, package = "PhyloCactus")
  expect_true(nzchar(path))

  m <- utils::read.csv(system.file("extdata", "MANIFEST.csv", package = "PhyloCactus"),
                       stringsAsFactors = FALSE)
  row <- m[m$canonical == example_plastome, ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$role, "input")
  expect_equal(row$produced_by, "external")
  expect_false(nzchar(row$source_path))
  expect_match(row$description, "MN517611.1", fixed = TRUE)
})

test_that("the identification reads the example plastome as one record named by its accession", {
  path <- system.file("extdata", example_plastome, package = "PhyloCactus")
  skip_if(!nzchar(path), "example plastome not installed")

  q <- PhyloCactus:::.bc_identify_read_query(path)
  expect_equal(names(q), "MN517611.1")
  expect_gt(nchar(q[[1]]), 100000L)
  expect_lt(nchar(q[[1]]), 200000L)
})
