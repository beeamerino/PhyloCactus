test_that("run_treePL_direct restores the working directory after a failed run (real function call)", {
  skip_if(.Platform$OS.type != "unix", "requires a POSIX shell for system()")

  # run_treePL_direct() changes into cwd to run treePL, because treePL writes its outputs
  # relative to the working directory. The restoration is registered with on.exit(), so it
  # must survive the error path as well: a failed dating run that left the session in
  # another directory would silently misdirect every later write in the pipeline.
  tmp_dir <- withr::local_tempdir()
  cfg <- file.path(tmp_dir, "configure_broken")
  writeLines("this is not a valid treePL configuration", cfg)

  before <- normalizePath(getwd())

  expect_error(
    suppressWarnings(run_treePL_direct(cfg_file = basename(cfg), label = "broken", cwd = tmp_dir))
  )

  expect_equal(normalizePath(getwd()), before)
})

test_that("run_treePL_direct leaves the working directory untouched when cwd is NULL (real function call)", {
  skip_if(.Platform$OS.type != "unix", "requires a POSIX shell for system()")

  tmp_dir <- withr::local_tempdir()
  cfg <- file.path(tmp_dir, "configure_broken")
  writeLines("this is not a valid treePL configuration", cfg)

  before <- normalizePath(getwd())

  expect_error(
    suppressWarnings(run_treePL_direct(cfg_file = cfg, label = "broken_nocwd", cwd = NULL))
  )

  expect_equal(normalizePath(getwd()), before)
})

test_that("run_treePL_direct reports a non-zero exit status as an error (real function call)", {
  skip_if(.Platform$OS.type != "unix", "requires a POSIX shell for system()")

  # Whether treePL is installed or not, an invalid configuration must surface as an R
  # error naming the label, never as a silent success that yields no chronogram.
  tmp_dir <- withr::local_tempdir()
  cfg <- file.path(tmp_dir, "configure_broken")
  writeLines("this is not a valid treePL configuration", cfg)

  expect_error(
    suppressWarnings(run_treePL_direct(cfg_file = cfg, label = "labelled_run", cwd = NULL)),
    "labelled_run"
  )
})
