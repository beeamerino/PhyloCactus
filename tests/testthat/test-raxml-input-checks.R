# AP-19 (closure review of 0.5.4): only run_local_bootstraps() checked its input files; the other
# wrappers of RAxML-NG ran the program on a missing file, calculate_rf_distances() did not check the
# binary, and map_branch_supports() reported its failures as "RAxML-NG ML".

test_that("the RAxML-NG wrappers name each missing input file before running the program", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  fake <- file.path(dir, "raxml-ng"); writeLines(c("#!/bin/sh", "exit 3"), fake); Sys.chmod(fake, "0755")
  aln <- file.path(dir, "a.phy"); writeLines("x", aln)
  nope <- file.path(dir, "missing.txt")
  expect_error(calculate_ml_tree(fake, aln, nope, nope, output_dir = dir), "part_file = .*constraint_file =")
  expect_error(map_branch_supports(fake, nope, aln), "best_tree =")
  expect_error(calculate_temporal_bootstraps(fake, aln, aln, nope, output_dir = dir), "best_tree_file =")
  expect_error(calculate_rf_distances(fake, nope), "ml_trees_file =")
  expect_error(calculate_rf_distances(file.path(dir, "no-binary"), aln), "not found in your system's PATH")
  # map_branch_supports() reports its own name when RAxML-NG fails.
  expect_error(map_branch_supports(fake, aln, aln, output_dir = dir), "map_branch_supports")
})
