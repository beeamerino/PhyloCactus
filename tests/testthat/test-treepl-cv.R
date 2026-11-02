# run_treePL_cv() replaced a vendored shell script on 2026-09-02. These tests cover the three
# places where the script's behaviour was wrong or fragile, which are the reasons the replacement
# is not a straight port:
#
#   1. It read the priming block by taking the LAST CHARACTER of each line and reassembling them
#      positionally, so an absent `moredetail` flag shifted every field after it.
#   2. It took the MOST FREQUENT priming combination; Maurin (2020) instructs taking the lowest
#      `opt` and `optad`.
#   3. It returned the smoothing value with the lowest chi-square without noticing when that value
#      was the smallest one evaluated, which is the grid running out rather than a minimum. That is
#      what happened to this project in August 2026.

prime_output <- function(opt, optad, optcvad,
                         moredetail = TRUE, moredetailad = TRUE, moredetailcvad = TRUE) {
  c("Linear search failed", "Unable to progress",
    "PLACE THE LINES BELOW IN THE CONFIGURATION FILE",
    paste0("opt = ", opt),
    if (moredetail) "moredetail",
    paste0("optad = ", optad),
    if (moredetailad) "moredetailad",
    paste0("optcvad = ", optcvad),
    if (moredetailcvad) "moredetailcvad")
}

test_that("the priming block is read by keyword, with every flag present", {
  p <- .parse_treepl_prime(prime_output(2, 3, 5))
  expect_equal(p$opt, 2L)
  expect_equal(p$optad, 3L)
  expect_equal(p$optcvad, 5L)
  expect_true(p$moredetail && p$moredetailad && p$moredetailcvad)
})

test_that("an absent moredetail flag does not shift the values after it", {
  # The failure mode of the shell script: with `moredetail` gone, its positional reassembly read
  # optad's value as the moredetail flag and everything downstream moved by one.
  p <- .parse_treepl_prime(prime_output(2, 3, 5, moredetail = FALSE))
  expect_equal(p$opt, 2L)
  expect_equal(p$optad, 3L)
  expect_equal(p$optcvad, 5L)
  expect_false(p$moredetail)
  expect_true(p$moredetailad)

  p2 <- .parse_treepl_prime(prime_output(1, 4, 2, moredetailad = FALSE, moredetailcvad = FALSE))
  expect_equal(c(p2$opt, p2$optad, p2$optcvad), c(1L, 4L, 2L))
  expect_true(p2$moredetail)
  expect_false(p2$moredetailad)
  expect_false(p2$moredetailcvad)
})

test_that("`optad` is not read as `opt`, and `moredetail` is not matched by `moredetailad`", {
  # Both are prefixes of the other keyword, so an unanchored pattern silently swaps them.
  p <- .parse_treepl_prime(c("PLACE THE LINES BELOW IN THE CONFIGURATION FILE",
                             "optad = 7", "moredetailad", "opt = 1", "optcvad = 9"))
  expect_equal(p$opt, 1L)
  expect_equal(p$optad, 7L)
  expect_false(p$moredetail)
  expect_true(p$moredetailad)
})

test_that("output without the marker line yields NULL rather than a guess", {
  expect_null(.parse_treepl_prime(c("Linear search failed", "opt = 2")))
  expect_null(.parse_treepl_prime(character(0)))
  expect_null(.parse_treepl_prime(c("PLACE THE LINES BELOW IN THE CONFIGURATION FILE",
                                    "nothing usable here")))
})

test_that("the lowest rule takes the lowest opt and optad, and the modal rule the most frequent", {
  primes <- do.call(rbind, list(
    .parse_treepl_prime(prime_output(5, 5, 5)),
    .parse_treepl_prime(prime_output(5, 5, 5)),
    .parse_treepl_prime(prime_output(5, 5, 5)),
    .parse_treepl_prime(prime_output(1, 2, 3))   # rarest, and the lowest
  ))

  # Maurin (2020), Step 4: repeat priming and select the lines with the lowest opt and optad.
  lowest <- .prime_cfg_lines(primes, rule = "lowest")
  expect_true("opt = 1" %in% lowest)
  expect_true("optad = 2" %in% lowest)

  # The shell script's rule, kept so that a run made before 2026-09-02 can be reproduced.
  modal <- .prime_cfg_lines(primes, rule = "modal")
  expect_true("opt = 5" %in% modal)
  expect_true("optad = 5" %in% modal)
})

test_that("flags absent from the chosen run are absent from the configuration lines", {
  primes <- .parse_treepl_prime(prime_output(2, 3, 4, moredetailcvad = FALSE))
  lines <- .prime_cfg_lines(primes)
  expect_true("moredetail" %in% lines)
  expect_true("moredetailad" %in% lines)
  expect_false("moredetailcvad" %in% lines)
  # Order matters: treePL reads each flag as applying to the assignment above it.
  expect_equal(lines, c("opt = 2", "moredetail", "optad = 3", "moredetailad", "optcvad = 4"))
})

# ---------------------------------------------------------------------------
# Selecting the smoothing value
# ---------------------------------------------------------------------------

cv_file_with <- function(dir, smoothing, chisq) {
  f <- file.path(dir, "cv_out")
  writeLines(paste0("chisq: (", smoothing, ") ", chisq), f)
  f
}

test_that("the smoothing value with the lowest chi-square is selected", {
  tmp <- withr::local_tempdir()
  f <- cv_file_with(tmp, c(100, 10, 1, 0.1, 0.01), c(500, 400, 100, 300, 450))

  w <- character(0)
  sel <- withCallingHandlers(.select_cv_smoothing(f),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })

  expect_equal(sel$smoothing, 1)
  expect_false(sel$at_edge)
  expect_equal(sel$shape, "interior")
  expect_length(w, 0L)
  expect_equal(nrow(sel$table), 5L)
  # Returned in ascending smoothing order, whatever order treePL wrote them.
  expect_equal(sel$table$smoothing, sort(sel$table$smoothing))
})

test_that("a minimum on a still-falling floor of the grid is reported instead of returned silently", {
  tmp <- withr::local_tempdir()
  # The real August 2026 curve: monotone all the way down, minimum at the smallest value tried,
  # and still changing by 2.6% across the last three points.
  f <- cv_file_with(tmp, c(10000, 1000, 100, 10, 1, 0.1, 0.01, 0.001, 0.0001),
                    c(515633, 515568, 514954, 508216, 386526, 50695.5, 38602.4, 33893.8, 33027.9))

  expect_warning(.select_cv_smoothing(f), "edge of its own grid")
  expect_warning(.select_cv_smoothing(f), "Maurin")
  sel <- suppressWarnings(.select_cv_smoothing(f))
  expect_equal(sel$smoothing, 1e-4)
  expect_true(sel$at_edge)
  expect_equal(sel$shape, "edge")
})

test_that("a minimum on a flat floor is a plateau, and is not warned about", {
  tmp <- withr::local_tempdir()
  # The real 2026-09-02 curve, after moving to randomcv and extending the grid to 1e-08. The
  # minimum is again on the floor, but the last three values differ by 0.4%, and the node ages
  # moved by at most 1.25 Ma across the same four orders of magnitude. Extending the grid further
  # would change the number selected without changing any age.
  f <- cv_file_with(tmp, c(0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8),
                    c(539.316, 397.637, 337.907, 305.69, 301.773, 301.083, 301.056, 300.612))

  w <- character(0)
  sel <- withCallingHandlers(suppressMessages(.select_cv_smoothing(f)),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })

  expect_length(w, 0L)                     # a plateau is not a defect
  expect_equal(sel$smoothing, 1e-8)
  expect_true(sel$at_edge)
  expect_equal(sel$shape, "plateau")
  expect_message(.select_cv_smoothing(f), "reaches a plateau")
})

test_that("a flat low end with the minimum inside it is a plateau, not an interior optimum", {
  tmp <- withr::local_tempdir()
  # The reference Cactaceae curve (one thread, grid 1e+03 to 1e-14). The lowest chi-square is at
  # 1e-12, not at the floor, but the three lowest smoothing values differ by 0.06% and the minimum
  # lies 0.04 units below the value at the floor. The values from 10 upward are of the order of
  # 1e41, which is how treePL fails at high smoothing.
  f <- cv_file_with(tmp,
    c(1000, 100, 10, 1, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11,
      1e-12, 1e-13, 1e-14),
    c(9.27978e+40, 9.23935e+40, 1.12066e+41, 546.27, 722.437, 469.865, 288.245, 235.198, 227.567,
      224.911, 224.97, 223.22, 222.331, 222.663, 221.758, 221.122, 221.259, 221.164))

  w <- character(0)
  sel <- withCallingHandlers(suppressMessages(.select_cv_smoothing(f)),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })

  expect_length(w, 0L)
  expect_equal(sel$smoothing, 1e-12)
  expect_false(sel$at_edge)
  expect_equal(sel$shape, "plateau")
  expect_message(.select_cv_smoothing(f), "not an optimum")
})

test_that("a flat low end lying above an interior minimum is not a plateau", {
  tmp <- withr::local_tempdir()
  # The curve turns upward at 1e-3 and then levels off well above its minimum. The flat tail is not
  # where the data are best fitted, so the minimum stays an interior optimum.
  f <- cv_file_with(tmp, c(1, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6),
                    c(500, 300, 200, 150, 250, 251, 250.5))

  w <- character(0)
  sel <- withCallingHandlers(.select_cv_smoothing(f),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })

  expect_length(w, 0L)
  expect_equal(sel$smoothing, 0.001)
  expect_equal(sel$shape, "interior")
})

test_that("a minimum on the ceiling of the grid is reported too", {
  tmp <- withr::local_tempdir()
  f <- cv_file_with(tmp, c(0.01, 0.1, 1, 10), c(500, 400, 300, 100))
  expect_warning(.select_cv_smoothing(f), "edge of its own grid")
  expect_equal(suppressWarnings(.select_cv_smoothing(f))$smoothing, 10)
  expect_equal(suppressWarnings(.select_cv_smoothing(f))$shape, "edge")
})

test_that("an absent or unusable cross-validation file is an error", {
  tmp <- withr::local_tempdir()
  expect_error(.select_cv_smoothing(file.path(tmp, "absent")), "produced no output")

  empty <- file.path(tmp, "empty"); file.create(empty)
  expect_error(.select_cv_smoothing(empty), "produced no output")

  junk <- file.path(tmp, "junk")
  writeLines(c("some log line", "another one"), junk)
  expect_error(.select_cv_smoothing(junk), "No usable")
})

# ---------------------------------------------------------------------------
# The retired argument
# ---------------------------------------------------------------------------

test_that("automate_treePL() warns that wrapper_sh is retired instead of failing on it", {
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()

  tr <- ape::read.tree(text = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);")
  ml <- file.path(tmp, "ml.tree"); ape::write.tree(tr, ml)
  bs <- file.path(tmp, "bs.tree"); ape::write.tree(c(tr, tr), bs)
  cfg <- file.path(tmp, "user.cfg")
  writeLines(c("numsites = 1000", "mrca = root A D", "min = root 1", "max = root 5"), cfg)

  # Pre-seeded so Step 1 is skipped: the warning under test is raised before any stage runs, and
  # priming against a real treePL installation would cost n_prime invocations for nothing.
  ml_dir <- file.path(tmp, "res", "ML_tree"); dir.create(ml_dir, recursive = TRUE)
  ape::write.tree(tr, file.path(ml_dir, "treepl_ML_tree.tre"))
  writeLines(c("treefile = x", "numsites = 10", "smooth = 0.1"),
             file.path(ml_dir, "configure_smooth_ML_tree"))

  # A script written against the old signature must still run and say why the argument is gone,
  # rather than stopping on an absent file as it did between the two designs.
  expect_warning(
    try(suppressMessages(automate_treePL(
      cfg_file = cfg, ml_tree_file = ml, bs_trees_file = bs,
      results_dir = file.path(tmp, "res"), treePL_out = file.path(tmp, "out"),
      num_bs = 1, numsites = 1000, outgroup = c("C", "D"), seed = 1L, n_prime = 1L,
      wrapper_sh = "anything_at_all.sh")), silent = TRUE),
    "retired and ignored")
})

# ---------------------------------------------------------------------------
# Cross-validation threading and determinism (cv_nthreads)
# ---------------------------------------------------------------------------

test_that("run_treePL_cv() errors on invalid cv_nthreads", {
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  tr <- ape::read.tree(text = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);")
  tree_path <- file.path(tmp, "test.tree")
  ape::write.tree(tr, tree_path)
  cfg_path <- file.path(tmp, "test.cfg")
  writeLines(c("numsites = 1000", "min = root 1", "max = root 5"), cfg_path)

  expect_error(
    run_treePL_cv(cfg_path, tree_path, "invalid", cv_nthreads = 0L, work_dir = tmp),
    "positive integer"
  )
})


# ---------------------------------------------------------------------------
# Failed subsamples of the random cross-validation (decision D13 of the closure of 0.5.2)
# ---------------------------------------------------------------------------
# treePL writes, for each smoothing value, the running chi-square sum after each random subsample
# ("(chisq): v") and the total ("chisq: (s) v"); the last "smoothing:" block of the log is the fit
# at the selected value, with no subsamples. In the probes of the closure of 0.5.2 a failed
# subsample added 1e+36 to 1e+41. Seeds 20260919 and 8675309 failed at values other than the one
# they selected, which moved their selection without any message; the fixtures below reproduce
# that case, the case of a failure at every value, and the reference run (failures above 0.1 only).

# A cross-validation log and its cvoutfile from per-subsample chi-square, one row per smoothing
# value. The running sum is printed with eight significant digits, as treePL prints it.
cv_log_with <- function(dir, smoothing, per_subsample, selected = NULL) {
  log <- unlist(lapply(seq_along(smoothing), function(i) {
    run <- cumsum(per_subsample[i, ])
    c(paste0("smoothing:", smoothing[i]), "numparams:12",
      unlist(lapply(run, function(v) c("exit siman: 10.5", paste0("(chisq): ", signif(v, 8))))),
      paste0("chisq: (", smoothing[i], ") ", signif(run[length(run)], 8)), "--------------")
  }))
  if (!is.null(selected)) log <- c(log, paste0("smoothing:", selected), "after opt calc2: 5507.2184")
  writeLines(c("setting the random number seed to 8675309", "cvstart: 1", "cvstop: 1e-08", log),
             file.path(dir, "cv_run.log"))
  totals <- rowSums(per_subsample)
  writeLines(paste0("chisq: (", smoothing, ") ", signif(totals, 8)), file.path(dir, "cv_run"))
  c(log = file.path(dir, "cv_run.log"), cv = file.path(dir, "cv_run"))
}

grid_d13 <- c(1, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)

# Seed 8675309 of probe 3e: the seventh subsample fails at 1e-05 to 1e-03 and at 0.1 and 1; the
# values where it fails drop out, and the selection falls on 1e-08, an interior minimum.
fixture_away <- function() {
  m <- matrix(30 + rep(c(8, 6, 4, 3, 2, 1.5, 1, 0.5, 0, 0.4, 0.8), each = 10) + rep(0:9, 11) / 10,
              nrow = length(grid_d13), byrow = TRUE)
  m[c(1, 2, 4, 5, 6), 7] <- 3.2e39
  m
}

test_that("the failed subsamples are counted per smoothing value, and the final fit is not read as one", {
  tmp <- withr::local_tempdir()
  f <- cv_log_with(tmp, grid_d13, fixture_away(), selected = 1e-8)
  fl <- .read_cv_subsamples(f[["log"]])
  expect_equal(fl$smoothing, sort(grid_d13))
  expect_equal(unique(fl$subsamples), 10L)
  expect_equal(fl$failed[match(c(1, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6), fl$smoothing)],
               c(1L, 1L, 0L, 1L, 1L, 1L, 0L))
  # After a failure the running sum is printed at eight significant digits, so the subsamples that
  # follow add nothing visible; they are not counted as failures.
  expect_equal(sum(fl$failed), 5L)
})

test_that("a failed subsample away from the selected value is warned about, with each value named", {
  tmp <- withr::local_tempdir()
  f <- cv_log_with(tmp, grid_d13, fixture_away(), selected = 1e-8)
  sel <- .select_cv_smoothing(f[["cv"]])
  expect_equal(sel$smoothing, 1e-8)
  w <- character(0)
  chk <- withCallingHandlers(
    .check_cv_subsamples(.read_cv_subsamples(f[["log"]]), sel$table, sel$smoothing, "run"),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_length(w, 1L)
  # Values of 0.1 or below only: smoothing 1 failed too, and is left out.
  expect_match(w, "0.1 (1 of 10 subsamples)", fixed = TRUE)
  expect_match(w, "0.001 (1 of 10)", fixed = TRUE)
  expect_match(w, "1e-04 (1 of 10)", fixed = TRUE)
  expect_match(w, "1e-05 (1 of 10)", fixed = TRUE)
  expect_no_match(w, "(^|[^.0-9e-])1 \\(1 of 10")
  expect_match(w, "selected value, 1e-08, is not among them", fixed = TRUE)
  expect_match(w, "no subsample failed is 1e-08", fixed = TRUE)
  expect_equal(chk$smoothing_without_failures, 1e-8)
  expect_equal(chk$failures$failed[chk$failures$smoothing == 1e-5], 1L)
})

test_that("a failure at the selected value and at every other value is warned about too", {
  tmp <- withr::local_tempdir()
  # Scheme S4 of probe 3b: the same two subsamples fail at every value, so every total is of the
  # order of 1e+41 and the selected value is the lowest of failed values.
  m <- fixture_away()
  m[, 7] <- 30
  m[, c(2, 9)] <- 1.9e41 + abs(seq_len(nrow(m)) - 9L) * 1e39
  f <- cv_log_with(tmp, grid_d13, m, selected = 1e-8)
  sel <- .select_cv_smoothing(f[["cv"]])
  w <- character(0)
  chk <- withCallingHandlers(
    .check_cv_subsamples(.read_cv_subsamples(f[["log"]]), sel$table, sel$smoothing, "run"),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_length(w, 1L)
  expect_match(w, "1e-08 (2 of 10)", fixed = TRUE)
  expect_match(w, "is among them", fixed = TRUE)
  expect_match(w, "every value of 0.1 or below", fixed = TRUE)
  expect_true(is.na(chk$smoothing_without_failures))
})

test_that("failures above 0.1 only, as in the reference run, raise no warning", {
  tmp <- withr::local_tempdir()
  m <- fixture_away()
  m[, 7] <- 30
  m[1, c(3, 5)] <- 9.3e40
  f <- cv_log_with(tmp, grid_d13, m, selected = 1e-8)
  sel <- .select_cv_smoothing(f[["cv"]])
  w <- character(0)
  chk <- withCallingHandlers(
    .check_cv_subsamples(.read_cv_subsamples(f[["log"]]), sel$table, sel$smoothing, "run"),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_length(w, 0L)
  expect_equal(chk$failures$failed[chk$failures$smoothing == 1], 2L)
  expect_equal(chk$smoothing_without_failures, sel$smoothing)
})

test_that("a log without per-subsample lines is said so, not read as free of failures", {
  tmp <- withr::local_tempdir()
  f <- file.path(tmp, "cv.log")
  writeLines(c("smoothing:0.1", "chisq: (0.1) 300"), f)
  expect_null(.read_cv_subsamples(f))
  expect_null(.read_cv_subsamples(file.path(tmp, "absent.log")))
  expect_message(chk <- .check_cv_subsamples(NULL, data.frame(smoothing = 0.1, chisq = 300), 0.1, "run"),
                 "not checked")
  expect_true(is.na(chk$smoothing_without_failures))
})

test_that("run_treePL_cv() warns on the failed subsamples and returns them", {
  skip_if_not_installed("ape")
  skip_if(.Platform$OS.type != "unix", "requires a POSIX shell for the stand-in treePL")
  tmp <- withr::local_tempdir()
  fx <- cv_log_with(tmp, grid_d13, fixture_away(), selected = 1e-8)
  # Stand-in treePL: the priming block, then the fixture as the cross-validation, then an
  # ultrametric tree as the dating.
  bin <- file.path(tmp, "treePL_standin")
  writeLines(c("#!/bin/sh", "cfg=\"$1\"",
               "if grep -q '^prime$' \"$cfg\"; then",
               "  printf 'PLACE THE LINES BELOW IN THE CONFIGURATION FILE\\nopt = 1\\noptad = 1\\noptcvad = 1\\n'",
               "elif grep -q '^randomcv$' \"$cfg\"; then",
               paste0("  cp '", fx[["cv"]], "' \"$(sed -n 's/^cvoutfile = //p' \"$cfg\")\""),
               paste0("  cat '", fx[["log"]], "'"),
               "else",
               "  echo '((A:1,B:1):1,(C:1,D:1):1);' > \"$(sed -n 's/^outfile = //p' \"$cfg\")\"",
               "fi"), bin)
  Sys.chmod(bin, "0755")
  tr <- file.path(tmp, "t.tree")
  ape::write.tree(ape::read.tree(text = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);"), tr)
  cfg <- file.path(tmp, "t.cfg")
  writeLines(c("numsites = 1000", "mrca = root A D", "min = root 1", "max = root 5"), cfg)

  w <- character(0)
  res <- withCallingHandlers(
    run_treePL_cv(cfg, tr, "fx", n_prime = 1L, treepl_bin = bin,
                  work_dir = file.path(tmp, "wd"), quiet = TRUE),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("1e-05 (1 of 10)", w, fixed = TRUE)))
  expect_equal(res$smoothing, 1e-8)
  expect_equal(res$cv_smoothing_without_failures, 1e-8)
  expect_equal(sum(res$cv_failures$failed[res$cv_failures$smoothing <= 0.1]), 4L)
})

test_that("automate_treePL() repeats the check on a cross-validation read from disk", {
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  tr <- ape::read.tree(text = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);")
  ml <- file.path(tmp, "ml.tree"); ape::write.tree(tr, ml)
  bs <- file.path(tmp, "bs.tree"); ape::write.tree(c(tr, tr), bs)
  cfg <- file.path(tmp, "user.cfg")
  writeLines(c("numsites = 1000", "mrca = root A D", "min = root 1", "max = root 5"), cfg)
  ml_dir <- file.path(tmp, "res", "ML_tree"); dir.create(ml_dir, recursive = TRUE)
  ape::write.tree(ape::read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);"), file.path(ml_dir, "treepl_ML_tree.tre"))
  writeLines(c("treefile = x", "numsites = 10", "smooth = 1e-08"), file.path(ml_dir, "configure_smooth_ML_tree"))
  fx <- cv_log_with(tmp, grid_d13, fixture_away(), selected = 1e-8)
  file.copy(fx[["cv"]], file.path(ml_dir, "cv_ML_tree"))
  file.copy(fx[["log"]], file.path(ml_dir, "cv_ML_tree.log"))

  # The run stops later, at the dating of the replicates, without a treePL installation; the
  # warning under test is raised before that.
  w <- character(0)
  withCallingHandlers(
    try(suppressMessages(automate_treePL(
      cfg_file = cfg, ml_tree_file = ml, bs_trees_file = bs,
      results_dir = file.path(tmp, "res"), treePL_out = file.path(tmp, "out"),
      num_bs = 1, numsites = 1000, outgroup = c("C", "D"), seed = 1L, n_prime = 1L,
      treepl_bin = file.path(tmp, "absent_treePL"))), silent = TRUE),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("1e-04 (1 of 10)", w, fixed = TRUE)))
})

# Anomalous subsamples excluded from the selection (phase 3h of the closure of 0.5.2)
# ---------------------------------------------------------------------------
# In the re-dating of 30-10-2026, under four calibration schemes without internal bounds, one or two
# random subsamples of each cross-validation fitted badly from 0.1 down to 1e-07 (increments of 700
# to 2400, or a failure of 1e+38, against 20 to 80 for the other subsamples) and normally from 1e-08
# down. They moved every selection to 1e-09 to 1e-11; summed without them, each curve had its
# minimum at 1e-05. A subsample is anomalous when it fails (D13) or adds more than
# `cv_outlier_factor` times the median increment at any value of 0.1 or below. The values are
# compared on one set of subsamples: after a failure the running sum hides the subsamples that
# follow, so either the values with hidden subsamples or the hidden subsamples are left out,
# whichever keeps more subsample-by-value cells.

# Nine subsamples with a minimum at 1e-05; the fourth adds 950 from 0.1 to 1e-07 and 45 below, which
# moves the minimum of the plain sum to 1e-08 (the young tree of Zuntini et al. 2024 in 3h).
fixture_outlier <- function() {
  off <- c(8, 6, 4, 3, 1, 0, 0.5, 1, 1.5, 2, 2.5)
  m <- matrix(30 + rep(off, each = 10) + rep(0:9, 11) / 10, nrow = length(grid_d13), byrow = TRUE)
  m[1, 4] <- 80
  m[2:8, 4] <- 950
  m[9:11, 4] <- 45
  m
}

test_that("the increment of every subsample is read, and those hidden by a failure are marked lost", {
  tmp <- withr::local_tempdir()
  f <- cv_log_with(tmp, grid_d13, fixture_away(), selected = 1e-8)
  inc <- .read_cv_increments(f[["log"]])
  expect_setequal(names(inc), c("smoothing", "subsample", "increment", "failed", "lost"))
  expect_equal(nrow(inc), length(grid_d13) * 10L)
  at <- function(s, i) inc[abs(log10(inc$smoothing / s)) < 1e-9 & inc$subsample == i, ]
  expect_true(at(1e-5, 7)$failed)
  expect_true(all(vapply(8:10, function(i) at(1e-5, i)$lost, logical(1))))
  expect_true(all(is.na(vapply(8:10, function(i) at(1e-5, i)$increment, numeric(1)))))
  expect_false(at(1e-6, 8)$lost)
  expect_equal(at(1e-6, 8)$increment, 31.7, tolerance = 1e-6)
  expect_null(.read_cv_increments(file.path(tmp, "absent.log")))
})

test_that("an anomalous subsample is excluded and the selection falls on the minimum of the others", {
  tmp <- withr::local_tempdir()
  f <- cv_log_with(tmp, grid_d13, fixture_outlier(), selected = 1e-8)
  expect_equal(.select_cv_smoothing(f[["cv"]])$smoothing, 1e-8)
  r <- .robust_cv_curve(.read_cv_increments(f[["log"]]), factor = 5)
  expect_equal(r$excluded, 4L)
  expect_equal(r$smoothing, 1e-5)
  expect_equal(r$shape, "interior")
  expect_length(r$dropped_values, 0L)
  expect_length(r$dropped_subsamples, 0L)
  expect_true(all(r$table$smoothing <= 0.1))
  expect_null(.robust_cv_curve(.read_cv_increments(f[["log"]]), factor = NULL))
})

test_that("with a failure hiding later subsamples, the hidden subsamples are left out when that keeps more cells", {
  tmp <- withr::local_tempdir()
  # The seventh subsample fails from 0.1 to 1e-05 (Ramirez-Barahona et al. 2020 in 3h): leaving out
  # subsamples 8 to 10 keeps 6 x 10 cells; leaving out the five values keeps 9 x 5.
  m <- fixture_outlier(); m[, 4] <- 30 + c(8, 6, 4, 3, 1, 0, 0.5, 1, 1.5, 2, 2.5) + 0.3
  m[2:6, 7] <- 4.1e38
  f <- cv_log_with(tmp, grid_d13, m, selected = 1e-8)
  r <- .robust_cv_curve(.read_cv_increments(f[["log"]]), factor = 5)
  expect_equal(r$excluded, 7L)
  expect_equal(r$dropped_subsamples, 8:10)
  expect_length(r$dropped_values, 0L)
  expect_equal(r$smoothing, 1e-5)
})

test_that("with a failure at one value only, that value is left out when that keeps more cells", {
  tmp <- withr::local_tempdir()
  # The second subsample fails at 0.1 only (the old tree of Zuntini et al. 2024 in 3h): leaving out
  # the value keeps 9 x 9 cells; leaving out the hidden subsamples would keep 1 x 10.
  m <- fixture_outlier(); m[, 4] <- 30 + c(8, 6, 4, 3, 1, 0, 0.5, 1, 1.5, 2, 2.5) + 0.3
  m[2, 2] <- 3.7e40
  f <- cv_log_with(tmp, grid_d13, m, selected = 1e-5)
  r <- .robust_cv_curve(.read_cv_increments(f[["log"]]), factor = 5)
  expect_equal(r$excluded, 2L)
  expect_equal(r$dropped_values, 0.1)
  expect_length(r$dropped_subsamples, 0L)
  expect_equal(r$smoothing, 1e-5)
})

test_that("without an anomalous subsample nothing is excluded, so the selection is that of the plain sum", {
  tmp <- withr::local_tempdir()
  # As in the reference run: failures above 0.1 only.
  m <- fixture_away(); m[, 7] <- 30; m[1, c(3, 5)] <- 9.3e40
  f <- cv_log_with(tmp, grid_d13, m, selected = 1e-8)
  r <- .robust_cv_curve(.read_cv_increments(f[["log"]]), factor = 5)
  expect_length(r$excluded, 0L)
  expect_equal(r$smoothing, .select_cv_smoothing(f[["cv"]])$smoothing)
})

test_that("too few subsamples left is reported, and no curve is returned", {
  tmp <- withr::local_tempdir()
  # More than half of the subsamples cannot exceed five times the median at one value, so too few
  # remain only when subsamples fail: here subsamples 1 to 6 each fail at one value.
  m <- fixture_outlier()
  for (i in 1:6) m[i + 1L, i] <- 2e39
  f <- cv_log_with(tmp, grid_d13, m, selected = 1e-8)
  r <- .robust_cv_curve(.read_cv_increments(f[["log"]]), factor = 5)
  expect_null(r$table)
  expect_match(r$reason, "fewer than 5 subsamples", fixed = TRUE)
})

# A stand-in treePL that answers the priming, returns the fixture as the cross-validation and writes
# an ultrametric tree as the dating.
standin_treepl <- function(tmp, fx) {
  bin <- file.path(tmp, "treePL_standin")
  writeLines(c("#!/bin/sh", "cfg=\"$1\"",
               "if grep -q '^prime$' \"$cfg\"; then",
               "  printf 'PLACE THE LINES BELOW IN THE CONFIGURATION FILE\\nopt = 1\\noptad = 1\\noptcvad = 1\\n'",
               "elif grep -q '^randomcv$' \"$cfg\"; then",
               paste0("  cp '", fx[["cv"]], "' \"$(sed -n 's/^cvoutfile = //p' \"$cfg\")\""),
               paste0("  cat '", fx[["log"]], "'"),
               "else",
               "  echo '((A:1,B:1):1,(C:1,D:1):1);' > \"$(sed -n 's/^outfile = //p' \"$cfg\")\"",
               "fi"), bin)
  Sys.chmod(bin, "0755")
  bin
}

test_that("run_treePL_cv() selects without the anomalous subsample by default, says so and writes the curve", {
  skip_if_not_installed("ape")
  skip_if(.Platform$OS.type != "unix", "requires a POSIX shell for the stand-in treePL")
  tmp <- withr::local_tempdir()
  fx <- cv_log_with(tmp, grid_d13, fixture_outlier(), selected = 1e-5)
  bin <- standin_treepl(tmp, fx)
  tr <- file.path(tmp, "t.tree")
  ape::write.tree(ape::read.tree(text = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);"), tr)
  cfg <- file.path(tmp, "t.cfg")
  writeLines(c("numsites = 1000", "mrca = root A D", "min = root 1", "max = root 5"), cfg)
  wd <- file.path(tmp, "wd")
  w <- character(0)
  res <- withCallingHandlers(
    run_treePL_cv(cfg, tr, "fx", n_prime = 1L, treepl_bin = bin, work_dir = wd, quiet = TRUE),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_equal(res$smoothing, 1e-5)
  expect_equal(res$cv_smoothing_all_subsamples, 1e-8)
  expect_equal(res$cv_excluded_subsamples, 4L)
  expect_true(any(grepl("subsample 4", w, fixed = TRUE)))
  expect_true(file.exists(file.path(wd, "cv_fx_selection.csv")))
  expect_true(file.exists(file.path(wd, "cv_fx_selection.txt")))
  expect_true(any(grepl("^smooth = 0.00001$", readLines(file.path(wd, "configure_smooth_fx")))))

  res0 <- suppressWarnings(run_treePL_cv(cfg, tr, "fx", n_prime = 1L, treepl_bin = bin,
                                         work_dir = file.path(tmp, "wd0"), quiet = TRUE,
                                         cv_outlier_factor = NULL))
  expect_equal(res0$smoothing, 1e-8)
})

test_that("automate_treePL() passes cv_outlier_factor on, with 5 as its default", {
  expect_true("cv_outlier_factor" %in% names(formals(automate_treePL)))
  expect_equal(formals(automate_treePL)$cv_outlier_factor, 5)
  expect_equal(formals(run_treePL_cv)$cv_outlier_factor, 5)
})

# The threshold of the rule (review 4f of the closure of 0.5.2). Under one seed of the re-dating, two
# subsamples reached 6.1 and 6.4 times the median and a third 23.5 times, and the selection changed
# between a factor of 5 and 7. Here, from 0.1 to 1e-07, subsample 8 adds 25 times the median of the
# subsamples at that value, subsample 10 adds 5.5 times and subsample 4 adds 4 times; below 1e-07 all
# three are normal. With subsamples 4, 8 and 10 among them, the median is that of the other seven.
fixture_threshold <- function() {
  off <- c(8, 6, 4, 3, 1, 0, 0.5, 1, 1.5, 2, 2.5)
  m <- matrix(30 + rep(off, each = 10) + rep(0:9, 11) / 10, nrow = length(grid_d13), byrow = TRUE)
  med <- 30 + off + 0.55
  for (r in 2:8) {
    m[r, 8] <- 25 * med[r]
    m[r, 10] <- 5.5 * med[r]
    m[r, 4] <- 4 * med[r]
  }
  m
}

test_that("the factor decides which subsamples are left out, against the median at each value", {
  tmp <- withr::local_tempdir()
  f <- cv_log_with(tmp, grid_d13, fixture_threshold(), selected = 1e-8)
  inc <- .read_cv_increments(f[["log"]])
  expect_equal(.robust_cv_curve(inc, factor = 3)$excluded, c(4L, 8L, 10L))
  expect_equal(.robust_cv_curve(inc, factor = 5)$excluded, c(8L, 10L))
  expect_equal(.robust_cv_curve(inc, factor = 7)$excluded, 8L)
  expect_length(.robust_cv_curve(inc, factor = 30)$excluded, 0L)
  expect_equal(.robust_cv_curve(inc, factor = 3)$smoothing, 1e-5)
  expect_equal(.robust_cv_curve(inc, factor = 5)$smoothing, 1e-8)
  # The ratio that decides is written with one decimal, so a subsample left out for more than 5 times
  # the median does not read as 5 or 6.
  r5 <- .robust_cv_curve(inc, factor = 5)
  expect_true(any(grepl("subsample 10 added 5.5 to 5.5 times the median", r5$reasons, fixed = TRUE)))
})

test_that("when the subsample that fails at the selected value is left out, the D13 warning says so", {
  skip_if_not_installed("ape")
  skip_if(.Platform$OS.type != "unix", "requires a POSIX shell for the stand-in treePL")
  tmp <- withr::local_tempdir()
  # The seventh subsample fails from 0.1 to 1e-05 (the run of Ramirez-Barahona et al. 2020 in 3h);
  # without it the selection is 1e-05, a value where it failed.
  m <- fixture_outlier(); m[, 4] <- 30 + c(8, 6, 4, 3, 1, 0, 0.5, 1, 1.5, 2, 2.5) + 0.3
  m[2:6, 7] <- 4.1e38
  fx <- cv_log_with(tmp, grid_d13, m, selected = 1e-5)
  bin <- standin_treepl(tmp, fx)
  tr <- file.path(tmp, "t.tree")
  ape::write.tree(ape::read.tree(text = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);"), tr)
  cfg <- file.path(tmp, "t.cfg")
  writeLines(c("numsites = 1000", "mrca = root A D", "min = root 1", "max = root 5"), cfg)
  w <- character(0)
  res <- withCallingHandlers(
    run_treePL_cv(cfg, tr, "fx", n_prime = 1L, treepl_bin = bin, work_dir = file.path(tmp, "wd"), quiet = TRUE),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_equal(res$smoothing, 1e-5)
  d13 <- grep("failed in part", w, value = TRUE)
  expect_length(d13, 1L)
  expect_false(grepl("the selection falls on the values without failures", d13, fixed = TRUE))
  expect_match(d13, "left out of the selection", fixed = TRUE)
  expect_match(d13, "is among them", fixed = TRUE)
})

test_that("cv_outlier_factor must be NULL or one number above 1", {
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  tr <- file.path(tmp, "t.tree")
  ape::write.tree(ape::read.tree(text = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);"), tr)
  cfg <- file.path(tmp, "t.cfg")
  writeLines(c("numsites = 1000", "min = root 1", "max = root 5"), cfg)
  for (bad in list("5", 1, 0.5, c(5, 7), NA_real_, Inf)) {
    expect_error(run_treePL_cv(cfg, tr, "bad", cv_outlier_factor = bad, work_dir = tmp),
                 "cv_outlier_factor", fixed = TRUE)
  }
  expect_error(automate_treePL(cfg, tr, tr, results_dir = tmp, treePL_out = tmp, cv_outlier_factor = "5"),
               "cv_outlier_factor", fixed = TRUE)
})
