# Step 23 of Phase 12 (runbook of 13-10, publication of 0.5.0): NEWS keeps the section of 0.5.0,
# and the manifest fixes the number of terminals of every tree the run publishes (column
# expected_tips, read by publish_run_outputs.R of the audit repository, which stops when a tree of
# the run has another number). The trees of the Phase 12 run have 1208 terminals (Module 7, reduced
# matrix; recount after step 19). Written before the code; the last test turns green when step 23
# publishes the trees.

test_that("NEWS keeps the section of 0.5.0 (the version of the release is tested in test-release-0-5-1.R)", {
  f <- testthat::test_path("..", "..", "NEWS.md")
  skip_if(!file.exists(f), "NEWS.md is not in this tree")
  expect_length(grep("^# PhyloCactus 0\\.5\\.0$", readLines(f, warn = FALSE)), 1L)
})

test_that("the manifest fixes the terminals of every tree the run publishes (expected_tips)", {
  path <- system.file("extdata", "MANIFEST.csv", package = "PhyloCactus")
  skip_if(!nzchar(path), "MANIFEST.csv not installed")
  m <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_true("expected_tips" %in% names(m))
  run_trees <- m$role == "run_output" & m$validator %in% c("tree_rooted_clade", "tree_unrooted_clade")
  expect_setequal(m$canonical[run_trees], c("phylocactus_ml_tree.tree", "phylocactus_constraint_tree.tree",
                                            "phylocactus_chronogram_ml.tree", "phylocactus_chronogram_hpd.tree"))
  expect_true(all(m$expected_tips[run_trees] %in% 1208L))
  expect_true(all(is.na(m$expected_tips[!run_trees])))
})

test_that("the distributed trees of the run have the terminals the manifest expects", {
  skip_if_not_installed("ape")
  path <- system.file("extdata", "MANIFEST.csv", package = "PhyloCactus")
  skip_if(!nzchar(path), "MANIFEST.csv not installed")
  m <- utils::read.csv(path, stringsAsFactors = FALSE)
  skip_if(!"expected_tips" %in% names(m), "no expected_tips column")
  for (i in which(!is.na(m$expected_tips))) {
    f <- system.file("extdata", m$canonical[i], package = "PhyloCactus")
    first <- readLines(f, n = 1L, warn = FALSE)
    tr <- if (grepl("^#NEXUS", first, ignore.case = TRUE)) ape::read.nexus(f) else ape::read.tree(f)
    expect_identical(ape::Ntip(tr), as.integer(m$expected_tips[i]), info = m$canonical[i])
  }
})
