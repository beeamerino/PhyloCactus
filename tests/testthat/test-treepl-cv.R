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
  expect_false(sel$at_plateau)
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
  expect_true(sel$at_plateau)
  expect_message(.select_cv_smoothing(f), "plateau reached")
})

test_that("a minimum on the ceiling of the grid is reported too", {
  tmp <- withr::local_tempdir()
  f <- cv_file_with(tmp, c(0.01, 0.1, 1, 10), c(500, 400, 300, 100))
  expect_warning(.select_cv_smoothing(f), "edge of its own grid")
  expect_equal(suppressWarnings(.select_cv_smoothing(f))$smoothing, 10)
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
