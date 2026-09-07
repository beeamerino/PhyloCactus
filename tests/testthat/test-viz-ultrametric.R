test_that(".validate_treepl_output rejects a non-ultrametric chronogram (real function call)", {
  skip_if_not_installed("ape")

  testdata <- testthat::test_path("testdata")
  nonultra_file <- file.path(testdata, "test_tree_nonultra.nwk")
  skip_if_not(file.exists(nonultra_file), "Test fixture missing")

  # This is the guard that stops a failed treePL run from being propagated downstream
  # as a corrupted chronogram, so the error path is the behaviour that matters.
  expect_error(
    .validate_treepl_output(nonultra_file, label = "test_nonultra", tol = 1e-4),
    "not ultrametric"
  )
})

test_that(".validate_treepl_output accepts an ultrametric chronogram and returns it (real function call)", {
  skip_if_not_installed("ape")

  testdata <- testthat::test_path("testdata")
  ultra_file <- file.path(testdata, "test_tree.nwk")
  skip_if_not(file.exists(ultra_file), "Test fixture missing")
  skip_if_not(ape::is.ultrametric(ape::read.tree(ultra_file), tol = 1e-4),
              "Fixture is not ultrametric at the tolerance used here")

  tr <- .validate_treepl_output(ultra_file, label = "test_ultra", tol = 1e-4)

  expect_s3_class(tr, "phylo")
  expect_equal(length(tr$tip.label), length(ape::read.tree(ultra_file)$tip.label))
})

test_that(".validate_treepl_output rejects a missing or empty output file (real function call)", {
  skip_if_not_installed("ape")

  tmp_dir <- withr::local_tempdir()

  missing_file <- file.path(tmp_dir, "does_not_exist.tre")
  expect_error(
    .validate_treepl_output(missing_file, label = "absent"),
    "missing or empty"
  )

  empty_file <- file.path(tmp_dir, "empty.tre")
  file.create(empty_file)
  expect_error(
    .validate_treepl_output(empty_file, label = "empty"),
    "missing or empty"
  )
})

test_that(".validate_treepl_output rejects a file that is not valid Newick (real function call)", {
  skip_if_not_installed("ape")

  tmp_dir <- withr::local_tempdir()
  junk_file <- file.path(tmp_dir, "junk.tre")
  writeLines("this is not a newick tree", junk_file)

  # ape::read.tree() warns before returning NULL on unparseable input; the warning is
  # expected noise here, the error is the behaviour under test.
  expect_error(
    suppressWarnings(.validate_treepl_output(junk_file, label = "junk")),
    "Newick|ultrametric"
  )
})
