# These tests pre-seed the ML_tree results directory with a pre-computed, ultrametric
# "treepl_ML_tree.tre" so automate_treePL() skips Step 1 (which would otherwise require a real
# `treePL` installation) and proceeds directly into the real Step 2/3 smoothing-extraction logic
# and Step 4 config-generation logic being tested here.

test_that("automate_treePL() extracts the smoothing value from configure_smooth_ML_tree and propagates it, together with numsites and calibration lines, into each bootstrap treePL config (real function call)", {
  skip_if_not_installed("ape")
  skip_if(.Platform$OS.type != "unix", "requires a POSIX shell for system()")

  tmp_dir <- withr::local_tempdir()
  results_dir <- file.path(tmp_dir, "results")
  ml_dir <- file.path(results_dir, "ML_tree")
  dir.create(ml_dir, recursive = TRUE)

  ultrametric_tree <- ape::rcoal(4, tip.label = c("A", "B", "C", "Portulaca_fulgens"))
  ape::write.tree(ultrametric_tree, file.path(ml_dir, "treepl_ML_tree.tre"))
  writeLines(c(
    "treefile = supportTree_raxml-ng_fixed.tree",
    "numsites = 100",
    "thorough",
    "smoothing = 0.1"
  ), file.path(ml_dir, "configure_smooth_ML_tree"))

  cfg_file <- file.path(tmp_dir, "calibs.cfg")
  writeLines(c(
    "numsites = 8500",
    "mrca = root A B Portulaca_fulgens",
    "min = root 1",
    "max = root 5"
  ), cfg_file)

  ml_tree_file <- file.path(tmp_dir, "ml_input.nwk")
  ape::write.tree(ultrametric_tree, ml_tree_file)

  bs_all <- c(ultrametric_tree, ultrametric_tree)
  class(bs_all) <- "multiPhylo"
  bs_trees_file <- file.path(tmp_dir, "bs.nwk")
  ape::write.tree(bs_all, bs_trees_file)

  # Step 5 (actually invoking `treePL` on the generated bootstrap config) will fail in any
  # environment without a real treePL installation; that failure is expected and irrelevant
  # here, because Step 4 (config generation, the logic under test) has already run and written
  # its file to disk by the time Step 5 is reached.
  # suppressWarnings because .check_treepl_convergence() fires on these fixtures: the trees
  # have four terminals and an objective around 30, so the printed digits rather than the
  # optimiser set the floor on a measurable improvement. The guard is about real trees; what
  # is under test here is configuration propagation.
  try(suppressWarnings(suppressMessages(
    automate_treePL(
      cfg_file = cfg_file,
      ml_tree_file = ml_tree_file,
      bs_trees_file = bs_trees_file,
      results_dir = results_dir,
      treePL_out = file.path(tmp_dir, "out"),
      outgroup = "Portulaca_fulgens",
      num_bs = 1,
      seed = 1L
    )
  )), silent = TRUE)

  bs_cfg_file <- file.path(results_dir, "BS_tree", "BS_1", "cfg_BS_1.cfg")
  expect_true(file.exists(bs_cfg_file))
  cfg_lines <- readLines(bs_cfg_file)

  # numsites falls back to the value in the main calibration cfg (8500) when not passed
  # explicitly to automate_treePL(), and is then divided by rescale_factor (100 by default),
  # because the branch lengths of every replicate have been multiplied by that same factor and
  # treePL reads a branch as edge.length * numsites expected substitutions. Propagating 8500
  # verbatim alongside rescaled branches is the defect this division exists to prevent.
  expect_true(any(grepl("^numsites = 85$", cfg_lines)))
  expect_false(any(grepl("^numsites = 8500$", cfg_lines)))
  # The smoothing value cross-validated once on the ML tree must be reused, unmodified, and
  # written under the keyword treePL actually reads. This assertion used to require
  # "smoothing = 0.1", which is the spelling treePL discards in silence before falling back to
  # its default of 10: the test passed for months while every chronogram was dated at 10.
  expect_true(any(grepl("^smooth = 0\\.1$", cfg_lines)))
  expect_false(any(grepl("^smoothing *=", cfg_lines)))
  # Calibration constraints from the main cfg must be copied through to each replicate.
  expect_true(any(grepl("^mrca = root A B Portulaca_fulgens$", cfg_lines)))
  expect_true(any(grepl("^min = root 1$", cfg_lines)))
  expect_true(any(grepl("^max = root 5$", cfg_lines)))
})

test_that("numsites and the branch-length rescaling stay coupled whatever the factor (real function call)", {
  skip_if_not_installed("ape")
  skip_if(.Platform$OS.type != "unix", "requires a POSIX shell for system()")

  # treePL reads a branch as edge.length * numsites. The pair therefore has one degree of
  # freedom, not two, and this test holds the product invariant across factors so the two can
  # never drift apart again: a run where the branches were scaled and numsites was not is
  # indistinguishable, in its output, from a correctly scaled run on a matrix that never existed.
  run_with <- function(tmp_dir, factor) {
    results_dir <- file.path(tmp_dir, "results")
    ml_dir <- file.path(results_dir, "ML_tree")
    dir.create(ml_dir, recursive = TRUE)

    # Same seed in both runs: the invariance being tested is over the factor, so the input tree
    # has to be identical and rcoal() is random.
    set.seed(42L)
    ultrametric_tree <- ape::rcoal(4, tip.label = c("A", "B", "C", "Portulaca_fulgens"))
    ape::write.tree(ultrametric_tree, file.path(ml_dir, "treepl_ML_tree.tre"))
    writeLines(c("treefile = x", "numsites = 100", "thorough", "smoothing = 0.1"),
               file.path(ml_dir, "configure_smooth_ML_tree"))

    cfg_file <- file.path(tmp_dir, "calibs.cfg")
    writeLines(c("numsites = 8500", "mrca = root A B Portulaca_fulgens",
                 "min = root 1", "max = root 5"), cfg_file)

    ml_tree_file <- file.path(tmp_dir, "ml_input.nwk")
    ape::write.tree(ultrametric_tree, ml_tree_file)
    bs_all <- c(ultrametric_tree, ultrametric_tree)
    class(bs_all) <- "multiPhylo"
    bs_trees_file <- file.path(tmp_dir, "bs.nwk")
    ape::write.tree(bs_all, bs_trees_file)

    # suppressWarnings because .check_treepl_convergence() fires on these fixtures: the trees
    # have four terminals and an objective around 30, so the printed digits rather than the
    # optimiser set the floor on a measurable improvement. The guard is about real trees; what
    # is under test here is configuration propagation.
    try(suppressWarnings(suppressMessages(automate_treePL(
      cfg_file = cfg_file, ml_tree_file = ml_tree_file,
      bs_trees_file = bs_trees_file, results_dir = results_dir,
      treePL_out = file.path(tmp_dir, "out"), outgroup = "Portulaca_fulgens",
      num_bs = 1, seed = 1L, rescale_factor = factor
    ))), silent = TRUE)

    list(
      cfg = readLines(file.path(results_dir, "BS_tree", "BS_1", "cfg_BS_1.cfg")),
      ml_cfg = readLines(file.path(ml_dir, "calibrations_scaled.cfg")),
      tree = ape::read.tree(file.path(ml_dir, "supportTree_raxml-ng_fixed.tree"))
    )
  }

  unscaled <- run_with(withr::local_tempdir(), 1)
  scaled   <- run_with(withr::local_tempdir(), 100)

  expect_true(any(grepl("^numsites = 8500$", unscaled$cfg)))
  expect_true(any(grepl("^numsites = 85$", scaled$cfg)))

  # The configuration actually handed to the maximum-likelihood run carries the scaled value too,
  # and is written beside the run rather than overwriting the user's declarative cfg.
  expect_true(any(grepl("^numsites = 85$", scaled$ml_cfg)))
  expect_true(any(grepl("^numsites = 8500$", unscaled$ml_cfg)))

  # Product invariance: branches up by the factor, sites down by it.
  expect_equal(
    sum(scaled$tree$edge.length) * 85,
    sum(unscaled$tree$edge.length) * 8500,
    tolerance = 0.01
  )
})

test_that("automate_treePL() falls back to parsing cv_ML_tree and selects the smoothing value with the lowest chi-square score when no smoothing line is present (real function call)", {
  skip_if_not_installed("ape")

  tmp_dir <- withr::local_tempdir()
  results_dir <- file.path(tmp_dir, "results")
  ml_dir <- file.path(results_dir, "ML_tree")
  dir.create(ml_dir, recursive = TRUE)

  ultrametric_tree <- ape::rcoal(4, tip.label = c("A", "B", "C", "Portulaca_fulgens"))
  ape::write.tree(ultrametric_tree, file.path(ml_dir, "treepl_ML_tree.tre"))

  # No "smoothing = " line; forces the fallback branch to parse cv_ML_tree instead.
  writeLines(c("treefile = supportTree_raxml-ng_fixed.tree", "numsites = 100", "thorough"),
             file.path(ml_dir, "configure_smooth_ML_tree"))

  # Lowest chi-square (50.000) is associated with smoothing = 1.0.
  writeLines(c(
    "1\t(0.1)\t123.456",
    "2\t(1.0)\t50.000",
    "3\t(10.0)\t200.000"
  ), file.path(ml_dir, "cv_ML_tree"))

  cfg_file <- file.path(tmp_dir, "calibs.cfg")
  writeLines("numsites = 100", cfg_file)

  ml_tree_file <- file.path(tmp_dir, "ml_input.nwk")
  ape::write.tree(ultrametric_tree, ml_tree_file)

  bs_trees_file <- file.path(tmp_dir, "bs.nwk")
  ape::write.tree(ultrametric_tree, bs_trees_file)

  # num_bs = 0 keeps the run entirely within Steps 0-3 (no bootstrap loop, no treePL execution
  # needed), so the extracted smoothing value can be observed directly via the printed message.
  expect_output(
    automate_treePL(
      cfg_file = cfg_file,
      ml_tree_file = ml_tree_file,
      bs_trees_file = bs_trees_file,
      results_dir = results_dir,
      treePL_out = file.path(tmp_dir, "out"),
      outgroup = "Portulaca_fulgens",
      num_bs = 0
    ),
    "Best smoothing strategy chosen: 1\\b"
  )
})

test_that("automate_treePL() stops with an informative error when no smoothing value can be extracted from either source (real function call)", {
  skip_if_not_installed("ape")

  tmp_dir <- withr::local_tempdir()
  results_dir <- file.path(tmp_dir, "results")
  ml_dir <- file.path(results_dir, "ML_tree")
  dir.create(ml_dir, recursive = TRUE)

  ultrametric_tree <- ape::rcoal(4, tip.label = c("A", "B", "C", "Portulaca_fulgens"))
  ape::write.tree(ultrametric_tree, file.path(ml_dir, "treepl_ML_tree.tre"))

  # Neither a "smoothing = " line nor a cv_ML_tree file is present.
  writeLines(c("treefile = supportTree_raxml-ng_fixed.tree", "numsites = 100", "thorough"),
             file.path(ml_dir, "configure_smooth_ML_tree"))

  cfg_file <- file.path(tmp_dir, "calibs.cfg")
  writeLines("numsites = 100", cfg_file)

  ml_tree_file <- file.path(tmp_dir, "ml_input.nwk")
  ape::write.tree(ultrametric_tree, ml_tree_file)

  bs_trees_file <- file.path(tmp_dir, "bs.nwk")
  ape::write.tree(ultrametric_tree, bs_trees_file)

  expect_error(
    automate_treePL(
      cfg_file = cfg_file,
      ml_tree_file = ml_tree_file,
      bs_trees_file = bs_trees_file,
      results_dir = results_dir,
      treePL_out = file.path(tmp_dir, "out"),
      outgroup = "Portulaca_fulgens",
      num_bs = 0
    ),
    "Failed to extract optimal smoothing parameter"
  )
})

# ---------------------------------------------------------------------------------------------
# Added 2026-09-02. treePL's configuration keyword is `smooth`. A line reading `smoothing = X` is
# not recognised, is discarded without a message, and the run continues on the built-in default of
# 10. The wrapper wrote `smoothing`, so every chronogram this project produced was dated at 10 and
# the cross-validation that selects the value never reached the program. No output revealed it:
# the trees parsed, were ultrametric, and gave plausible ages. The only visible trace was in
# treePL's own log, which prints the smoothing it is using.
# ---------------------------------------------------------------------------------------------

test_that(".verify_treepl_smoothing accepts a run that used the smoothing it was given", {
  tmp <- withr::local_tempdir()
  cfg <- file.path(tmp, "cfg"); log <- file.path(tmp, "log")
  writeLines(c("treefile = t.tre", "numsites = 128", "smooth = 0.1"), cfg)
  writeLines(c("numparams:3042", "smoothing:0.1", "after opt calc: 100.0"), log)

  expect_true(PhyloCactus:::.verify_treepl_smoothing(log, cfg, "ok_run"))
})

test_that(".verify_treepl_smoothing stops when treePL silently used its default", {
  tmp <- withr::local_tempdir()
  cfg <- file.path(tmp, "cfg"); log <- file.path(tmp, "log")
  # Exactly the situation found on 2026-09-02: the configuration asks for 1e-04 and treePL
  # reports the default.
  writeLines(c("treefile = t.tre", "numsites = 128", "smooth = 0.0001"), cfg)
  writeLines(c("numparams:3042", "smoothing:10", "after opt calc: 100.0"), log)

  expect_error(PhyloCactus:::.verify_treepl_smoothing(log, cfg, "silent_default"),
               "asked for")
  # The message has to name the keyword, because that is the fix and it is not guessable.
  expect_error(PhyloCactus:::.verify_treepl_smoothing(log, cfg, "silent_default"),
               "`smooth`")
})

test_that(".verify_treepl_smoothing does not read the old keyword as a request", {
  tmp <- withr::local_tempdir()
  cfg <- file.path(tmp, "cfg"); log <- file.path(tmp, "log")
  # A pre-2026-09-02 configuration. `smoothing = 0.1` is not something treePL was ever asked to
  # honour, so comparing it against the default and failing would report a defect in the wrong
  # place. Such a run is unverifiable, not wrong, and the function has to say so by declining.
  writeLines(c("treefile = t.tre", "numsites = 128", "smoothing = 0.1"), cfg)
  writeLines(c("numparams:3042", "smoothing:10"), log)

  expect_true(is.na(PhyloCactus:::.verify_treepl_smoothing(log, cfg, "legacy")))
})

test_that(".verify_treepl_smoothing declines rather than fails when it cannot compare", {
  tmp <- withr::local_tempdir()
  cfg <- file.path(tmp, "cfg"); log <- file.path(tmp, "log")
  writeLines(c("treefile = t.tre", "smooth = 1"), cfg)
  # A log with no smoothing line: an older treePL build, or a run that died before reporting.
  writeLines(c("finished reading config file", "tiny branch length at X."), log)
  expect_true(is.na(PhyloCactus:::.verify_treepl_smoothing(log, cfg, "no_line")))

  expect_true(is.na(PhyloCactus:::.verify_treepl_smoothing(
    file.path(tmp, "absent.log"), cfg, "no_log")))
})

# ---------------------------------------------------------------------------------------------
# report_smoothing_sensitivity(). Cross-validation selected 1e-04 for this dataset, which is the
# lowest value on the tested grid and therefore the edge of the search rather than an optimum.
# The published justification rests on this table instead: the ages barely move between 1e-04 and
# 100, so how the value was chosen does not affect what is reported.
# ---------------------------------------------------------------------------------------------

test_that(".mrca_sets_from_cfg reads the calibrated nodes back out of a configuration", {
  tmp <- withr::local_tempdir()
  cfg <- file.path(tmp, "cfg")
  writeLines(c(
    "treefile = t.tre",
    "numsites = 128",
    "mrca = ACP_root A B C",
    "min = ACP_root 27.8",
    "max = ACP_root 53.4",
    "mrca = Cact_Anac A B",
    "smooth = 0.1",
    "outfile = out.tre"
  ), cfg)

  sets <- PhyloCactus:::.mrca_sets_from_cfg(cfg)
  expect_named(sets, c("ACP_root", "Cact_Anac"))
  expect_equal(sets$ACP_root, c("A", "B", "C"))
  expect_equal(sets$Cact_Anac, c("A", "B"))
  # min/max lines carry a node name too and must not be mistaken for node definitions.
  expect_length(sets, 2L)
})

test_that(".smoothing_cfg replaces both spellings and writes the value plainly", {
  # A configuration from before 2026-09-02 carries `smoothing`. Leaving it in place would put two
  # smoothing lines in the file, and treePL would read neither as intended.
  old_style <- c("treefile = t.tre", "numsites = 128", "smoothing = 0.0001",
                 "outfile = old.tre")
  out <- PhyloCactus:::.smoothing_cfg(old_style, 1e-4, "new.tre")

  expect_false(any(grepl("^smoothing *=", out)))
  expect_equal(sum(grepl("^smooth *=", out)), 1L)
  expect_true(any(out == "outfile = new.tre"))
  expect_false(any(out == "outfile = old.tre"))
  # Written as 0.0001, not 1e-04: treePL's parser does not take scientific notation, and a
  # smoothing quietly reset to the default is the failure this table exists to expose.
  expect_true(any(out == "smooth = 0.0001"))
  expect_false(any(grepl("e-0", out)))
  # Everything else survives untouched.
  expect_true(any(out == "numsites = 128"))
})

test_that("report_smoothing_sensitivity refuses inputs that cannot produce a sensitivity analysis", {
  tmp <- withr::local_tempdir()
  cfg <- file.path(tmp, "cfg")
  writeLines(c("treefile = t.tre", "mrca = root A B", "min = root 1"), cfg)

  expect_error(report_smoothing_sensitivity(file.path(tmp, "absent"), treepl_bin = "true"),
               "does not exist")
  # One value is not a sensitivity analysis.
  expect_error(report_smoothing_sensitivity(cfg, smoothing_values = 1, treepl_bin = "true"),
               "two or more")
  expect_error(report_smoothing_sensitivity(cfg, smoothing_values = c(0, 1), treepl_bin = "true"),
               "positive")

  # Without mrca lines there are no calibrated nodes and the table would have no rows.
  no_nodes <- file.path(tmp, "cfg_bare")
  writeLines(c("treefile = t.tre", "numsites = 128"), no_nodes)
  expect_error(report_smoothing_sensitivity(no_nodes, treepl_bin = "true"), "No `mrca` lines")
})

test_that("report_smoothing_sensitivity records a run that produced nothing instead of inventing ages", {
  skip_if(.Platform$OS.type != "unix", "needs a no-op executable on PATH")
  tmp <- withr::local_tempdir()
  cfg <- file.path(tmp, "cfg")
  writeLines(c("treefile = t.tre", "numsites = 128", "mrca = root A B C",
               "min = root 1", "max = root 5"), cfg)

  # `true` exits cleanly and writes no chronogram: the shape of a treePL run that failed without
  # saying so, which is what the converged column exists to record.
  out <- suppressWarnings(suppressMessages(
    report_smoothing_sensitivity(cfg, smoothing_values = c(0.1, 1),
                                 treepl_bin = "true", work_dir = file.path(tmp, "wd"),
                                 out_csv = file.path(tmp, "sens.csv"))
  ))

  expect_equal(nrow(out), 2L)
  expect_true(all(!out$converged))
  expect_true(all(is.na(out$age_ma)))
  expect_equal(unique(out$node), "root")
  expect_true(file.exists(file.path(tmp, "sens.csv")))
  # The requested values are recorded even when nothing came back, so the failed rows are
  # attributable.
  expect_equal(sort(unique(out$smoothing)), c(0.1, 1))
})
