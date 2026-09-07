test_that("run_local_bootstraps generates and reports a seed when none is given (real function call)", {
  tmp_dir <- withr::local_tempdir()
  aln <- file.path(tmp_dir, "aln.phy")
  part <- file.path(tmp_dir, "part.txt")
  cons <- file.path(tmp_dir, "constraint.nwk")
  file.create(aln, part, cons)

  # Seed generation happens before the binary check, so this exercises the real seed
  # logic without requiring RAxML-NG. A reproducible run depends on that seed being
  # emitted, because it is the only record of what an unseeded call actually used.
  msg <- NULL
  expect_error(
    withCallingHandlers(
      run_local_bootstraps(raxml_bin_path = "raxml-ng-binary-that-does-not-exist",
                           aln_file = aln, part_file = part, constraint_file = cons,
                           seed = NULL),
      message = function(m) {
        msg <<- c(msg, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    "not found in your system's PATH"
  )

  seed_line <- grep("Using generated seed:", msg, value = TRUE)
  expect_length(seed_line, 1)
  seed_val <- as.numeric(sub(".*Using generated seed: *", "", seed_line))
  expect_false(is.na(seed_val))
  expect_gt(seed_val, 0)
  expect_lte(seed_val, .Machine$integer.max)
})

test_that("run_local_bootstraps does not announce a seed when one is supplied (real function call)", {
  tmp_dir <- withr::local_tempdir()
  aln <- file.path(tmp_dir, "aln.phy")
  part <- file.path(tmp_dir, "part.txt")
  cons <- file.path(tmp_dir, "constraint.nwk")
  file.create(aln, part, cons)

  msg <- NULL
  expect_error(
    withCallingHandlers(
      run_local_bootstraps(raxml_bin_path = "raxml-ng-binary-that-does-not-exist",
                           aln_file = aln, part_file = part, constraint_file = cons,
                           seed = 4242L),
      message = function(m) {
        msg <<- c(msg, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    "not found in your system's PATH"
  )

  expect_length(grep("Using generated seed:", msg, value = TRUE), 0)
})

test_that("two unseeded calls do not reuse the same generated seed (real function call)", {
  tmp_dir <- withr::local_tempdir()
  aln <- file.path(tmp_dir, "aln.phy")
  part <- file.path(tmp_dir, "part.txt")
  cons <- file.path(tmp_dir, "constraint.nwk")
  file.create(aln, part, cons)

  grab_seed <- function() {
    msg <- NULL
    try(
      withCallingHandlers(
        run_local_bootstraps(raxml_bin_path = "raxml-ng-binary-that-does-not-exist",
                             aln_file = aln, part_file = part, constraint_file = cons,
                             seed = NULL),
        message = function(m) {
          msg <<- c(msg, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      ),
      silent = TRUE
    )
    as.numeric(sub(".*Using generated seed: *", "", grep("Using generated seed:", msg, value = TRUE)[[1]]))
  }

  seeds <- vapply(1:5, function(i) grab_seed(), numeric(1))
  expect_gt(length(unique(seeds)), 1)
})
