# treePL exits 0 and writes an ultrametric tree whether or not its optimisation actually ran. The
# only place the difference shows is the log: `exit siman:` is the objective simulated annealing
# reached, `after opt calc2:` the one gradient descent reached from there. When the two are equal,
# the gradient phase moved nothing and the chronogram is the annealing endpoint.
#
# The numbers in these fixtures are the real ones from 2026-09-02, twenty-one runs of the same tree
# at seven smoothing values. The coherent runs improved by 0.90% and 0.71%; the ones that returned
# ACP_root at 188 Ma, or moved a node from its lower bound to its upper bound, improved by nothing.

log_with <- function(dir, label, siman, final, extra = character(0)) {
  f <- file.path(dir, paste0("treepl_run_", label, ".log"))
  writeLines(c("setting NLOPT: LD_LBFGS", "result: 3",
               paste0("exit siman: ", siman),
               paste0("after opt calc1: ", siman),
               paste0("after opt calc2: ", final),
               extra), f)
  f
}

test_that("a run whose optimisation improved on the annealing result passes quietly", {
  tmp <- withr::local_tempdir()
  # ambas, smoothing 1e-04: 0.90% improvement, and the ages were coherent.
  f <- log_with(tmp, "good", 5534.9608, 5485.1612)

  w <- character(0)
  imp <- withCallingHandlers(
    .check_treepl_convergence(f, "good"),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })

  expect_length(w, 0L)
  expect_equal(round(imp, 5), round((5534.9608 - 5485.1612) / 5534.9608, 5))
})

test_that("a run whose optimisation moved nothing is reported", {
  tmp <- withr::local_tempdir()
  # ambas, smoothing 1: identical to every printed digit. This is the run that put
  # Opuntioideae_mrca on its upper bound while every other smoothing left it on the lower one.
  f <- log_with(tmp, "stalled", 16430.919, 16430.919)

  expect_warning(.check_treepl_convergence(f, "stalled"), "nothing to do")
  expect_warning(.check_treepl_convergence(f, "stalled"), "report_smoothing_sensitivity")
  expect_equal(suppressWarnings(.check_treepl_convergence(f, "stalled")), 0)
})

test_that("an improvement far below the threshold is still reported", {
  tmp <- withr::local_tempdir()
  # Hernandez, smoothing 1: 0.0034%, the largest improvement among the incoherent runs. This is
  # the one that returned ACP_root at 188.13 Ma.
  f <- log_with(tmp, "barely", 16523.106, 16522.538)
  expect_warning(.check_treepl_convergence(f, "barely"), "gradient optimisation")
})

test_that("treePL giving up on a feasible start is an error, not a warning", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "treepl_run_infeasible.log")
  writeLines(c("attempting to get feasible start rates/dates.",
               "problem initializing. trying again.",
               "Failed setting feasible start rates/dates after 10 attempts. Aborting."), f)

  expect_error(.check_treepl_convergence(f, "infeasible"),
               "could not find feasible starting rates")
})

test_that("a log without both objectives returns NA rather than guessing", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "treepl_run_quiet.log")
  writeLines(c("smoothing:0.0001", "tiny branch length at Foo_bar. setting to 0.0078125"), f)

  expect_true(is.na(.check_treepl_convergence(f, "quiet")))
  expect_true(is.na(.check_treepl_convergence(file.path(tmp, "absent.log"), "absent")))
})

test_that("calc1 is used when the log carries no calc2", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "treepl_run_calc1.log")
  writeLines(c("exit siman: 1000", "after opt calc1: 900"), f)
  expect_equal(.check_treepl_convergence(f, "calc1"), 0.1)
})

# ---------------------------------------------------------------------------
# The seed treePL actually reads
# ---------------------------------------------------------------------------
# treePL has a `seed` keyword and seeds itself from the clock when it is absent. Both the
# cross-validation, which partitions sites, and the simulated annealing before every optimisation
# consume that randomness, so an unseeded run need not choose the same smoothing value or return
# the same ages twice. Until 2026-09-02 the `seed` argument of automate_treePL() seeded only R's
# choice of which replicates to date; nothing reached treePL.

test_that("the seed reaches the configuration treePL is handed", {
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()

  tr <- ape::read.tree(text = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);")
  ml <- file.path(tmp, "ml.tree"); ape::write.tree(tr, ml)
  bs <- file.path(tmp, "bs.tree"); ape::write.tree(c(tr, tr), bs)

  # Pre-seed the chronogram so Step 1 is skipped. Without it run_treePL_cv() primes n_prime times
  # against a real treePL installation, which is a hundred invocations per test for a stage none of
  # these tests is about: the configuration under inspection is written before Step 1 begins.
  ml_dir <- file.path(tmp, "res", "ML_tree"); dir.create(ml_dir, recursive = TRUE)
  ape::write.tree(tr, file.path(ml_dir, "treepl_ML_tree.tre"))
  writeLines(c("treefile = x", "numsites = 10", "smooth = 0.1"),
             file.path(ml_dir, "configure_smooth_ML_tree"))
  cfg <- file.path(tmp, "user.cfg")
  writeLines(c(paste0("treefile = ", ml), "numsites = 1000",
               "mrca = root A D", "min = root 1", "max = root 5"), cfg)

  # The run itself needs the treePL binary; the configuration is written before priming starts, so
  # a failed run still leaves the file to inspect.
  try(suppressWarnings(suppressMessages(automate_treePL(
    cfg_file = cfg, ml_tree_file = ml,
    bs_trees_file = bs, results_dir = file.path(tmp, "res"),
    treePL_out = file.path(tmp, "out"), num_bs = 1, numsites = 1000,
    outgroup = c("C", "D"), seed = 12345L))), silent = TRUE)

  scaled <- list.files(tmp, pattern = "calibrations_scaled\\.cfg$",
                       recursive = TRUE, full.names = TRUE)
  skip_if(length(scaled) == 0L, "automate_treePL() stopped before writing the configuration")
  lines <- readLines(scaled[[1]])
  expect_true(any(grepl("^\\s*seed\\s*=\\s*12345\\s*$", lines)))
  # Exactly one, or treePL would read whichever came last.
  expect_equal(sum(grepl("^\\s*seed\\s*=", lines)), 1L)
})

test_that("a seed already in the user configuration is replaced and not duplicated", {
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()

  tr <- ape::read.tree(text = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);")
  ml <- file.path(tmp, "ml.tree"); ape::write.tree(tr, ml)
  bs <- file.path(tmp, "bs.tree"); ape::write.tree(c(tr, tr), bs)

  # Pre-seed the chronogram so Step 1 is skipped. Without it run_treePL_cv() primes n_prime times
  # against a real treePL installation, which is a hundred invocations per test for a stage none of
  # these tests is about: the configuration under inspection is written before Step 1 begins.
  ml_dir <- file.path(tmp, "res", "ML_tree"); dir.create(ml_dir, recursive = TRUE)
  ape::write.tree(tr, file.path(ml_dir, "treepl_ML_tree.tre"))
  writeLines(c("treefile = x", "numsites = 10", "smooth = 0.1"),
             file.path(ml_dir, "configure_smooth_ML_tree"))
  cfg <- file.path(tmp, "user.cfg")
  writeLines(c(paste0("treefile = ", ml), "numsites = 1000", "seed = 999",
               "mrca = root A D", "min = root 1", "max = root 5"), cfg)

  try(suppressWarnings(suppressMessages(automate_treePL(
    cfg_file = cfg, ml_tree_file = ml,
    bs_trees_file = bs, results_dir = file.path(tmp, "res"),
    treePL_out = file.path(tmp, "out"), num_bs = 1, numsites = 1000,
    outgroup = c("C", "D"), seed = 777L))), silent = TRUE)

  scaled <- list.files(tmp, pattern = "calibrations_scaled\\.cfg$",
                       recursive = TRUE, full.names = TRUE)
  skip_if(length(scaled) == 0L, "automate_treePL() stopped before writing the configuration")
  lines <- readLines(scaled[[1]])
  expect_equal(sum(grepl("^\\s*seed\\s*=", lines)), 1L)
  expect_true(any(grepl("^\\s*seed\\s*=\\s*777\\s*$", lines)))
  expect_false(any(grepl("999", lines)))
})

test_that("a non-positive seed is refused", {
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  tr <- ape::read.tree(text = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);")
  ml <- file.path(tmp, "ml.tree"); ape::write.tree(tr, ml)
  bs <- file.path(tmp, "bs.tree"); ape::write.tree(c(tr, tr), bs)

  # Pre-seed the chronogram so Step 1 is skipped. Without it run_treePL_cv() primes n_prime times
  # against a real treePL installation, which is a hundred invocations per test for a stage none of
  # these tests is about: the configuration under inspection is written before Step 1 begins.
  ml_dir <- file.path(tmp, "res", "ML_tree"); dir.create(ml_dir, recursive = TRUE)
  ape::write.tree(tr, file.path(ml_dir, "treepl_ML_tree.tre"))
  writeLines(c("treefile = x", "numsites = 10", "smooth = 0.1"),
             file.path(ml_dir, "configure_smooth_ML_tree"))
  cfg <- file.path(tmp, "user.cfg")
  writeLines(c(paste0("treefile = ", ml), "numsites = 1000",
               "mrca = root A D", "min = root 1", "max = root 5"), cfg)

  expect_error(suppressMessages(automate_treePL(
    cfg_file = cfg, ml_tree_file = ml,
    bs_trees_file = bs, results_dir = file.path(tmp, "res"),
    treePL_out = file.path(tmp, "out"), num_bs = 1, numsites = 1000,
    outgroup = c("C", "D"), seed = 0L)), "positive integer")
})
