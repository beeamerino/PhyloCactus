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

  # With cwd = NULL the run writes treepl_run_<label>.log into the working directory, which for the
  # test suite is tests/testthat. Moving into the temporary directory first keeps the run's output
  # with the run instead of leaving it in the repository for someone to delete by hand.
  withr::local_dir(tmp_dir)

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

  withr::local_dir(tmp_dir)

  expect_error(
    suppressWarnings(run_treePL_direct(cfg_file = cfg, label = "labelled_run", cwd = NULL)),
    "labelled_run"
  )
})


# Added on 2026-09-25, closing the carry-over declared on 2026-09-19. run_treePL_direct() ran
# system2("treePL", ...) with the name written into the code, ignoring the treepl_bin that
# run_treePL_cv() accepts. On BMM's Mac it worked because treePL is installed; on GitHub Actions it
# failed with exit status 127 on macOS and cost two tests, deleted that day. The declared scope is
# these two functions; automate_treePL() still carries the fixed name inside and is not touched here.

test_that("run_treePL_direct keeps the fixed name as its default, so no existing caller changes", {
  expect_true("treepl_bin" %in% names(formals(run_treePL_direct)))
  expect_equal(eval(formals(run_treePL_direct)$treepl_bin), "treePL")
})

test_that("run_treePL_direct runs the binary it is handed, and not the one written in the code", {
  skip_if(.Platform$OS.type != "unix", "requires a POSIX shell for system()")
  tmp_dir <- withr::local_tempdir()
  marca <- file.path(tmp_dir, "lo_corrio.txt")
  falso <- file.path(tmp_dir, "mi_treepl")
  writeLines(c("#!/bin/sh", paste0("echo llamado > ", shQuote(marca)), "exit 0"), falso)
  Sys.chmod(falso, "0755")
  cfg <- file.path(tmp_dir, "configure_falso")
  writeLines("this is not a valid treePL configuration", cfg)

  # What happens after the call is not the point of this test: the log of a fake binary does not
  # satisfy the smoothing check. The point is that the binary handed in is the one executed
  try(suppressWarnings(run_treePL_direct(cfg_file = basename(cfg), label = "falso",
                                         cwd = tmp_dir, treepl_bin = falso)), silent = TRUE)

  expect_true(file.exists(marca))
})

test_that("run_treePL_cv hands its binary down to the smoothing stage", {
  # Structural check, and it is declared as such: exercising this path end to end would need a fake
  # binary that also satisfies the priming parser and the cross validation, which would test the
  # fake more than the code. What matters here is that the argument is forwarded at the call site.
  cuerpo <- paste(deparse(body(run_treePL_cv)), collapse = " ")
  expect_match(cuerpo, "run_treePL_direct\\(smooth_cfg, label, treepl_bin", fixed = FALSE)
})


# Added on 2026-09-25, at BMM's request, closing the other half of the same defect. The declared
# carry-over of 2026-09-19 named only run_treePL_direct() and run_treePL_cv(), but automate_treePL()
# reaches treePL through both of them and handed neither a binary, so the whole dating path still
# ran whatever was called treePL on the machine. On GitHub that is nothing.

test_that("automate_treePL and its worker keep the fixed name as their default", {
  expect_true("treepl_bin" %in% names(formals(automate_treePL)))
  expect_equal(eval(formals(automate_treePL)$treepl_bin), "treePL")
  expect_true("treepl_bin" %in% names(formals(PhyloCactus:::.automate_treePL_run)))
  expect_equal(eval(formals(PhyloCactus:::.automate_treePL_run)$treepl_bin), "treePL")
})

test_that("automate_treePL runs the binary it is handed, down in the bootstrap stage", {
  skip_if_not_installed("ape")
  skip_if(.Platform$OS.type != "unix", "requires a POSIX shell for system()")
  tmp <- withr::local_tempdir()

  tr <- ape::read.tree(text = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);")
  ml <- file.path(tmp, "ml.tree"); ape::write.tree(tr, ml)
  bs <- file.path(tmp, "bs.tree"); ape::write.tree(c(tr, tr), bs)

  # The chronogram is pre-seeded so the maximum-likelihood stage is skipped, exactly as the tests of
  # the seed do. What is left is the bootstrap stage, which is the one that reaches
  # run_treePL_direct(), and that is the path this test is about
  ml_dir <- file.path(tmp, "res", "ML_tree"); dir.create(ml_dir, recursive = TRUE)
  ape::write.tree(tr, file.path(ml_dir, "treepl_ML_tree.tre"))
  writeLines(c("treefile = x", "numsites = 10", "smooth = 0.1"),
             file.path(ml_dir, "configure_smooth_ML_tree"))
  cfg <- file.path(tmp, "user.cfg")
  writeLines(c(paste0("treefile = ", ml), "numsites = 1000",
               "mrca = root A D", "min = root 1", "max = root 5"), cfg)

  marca <- file.path(tmp, "lo_corrio.txt")
  falso <- file.path(tmp, "mi_treepl")
  writeLines(c("#!/bin/sh", paste0("echo llamado >> ", shQuote(marca)), "exit 0"), falso)
  Sys.chmod(falso, "0755")

  try(suppressWarnings(suppressMessages(automate_treePL(
    cfg_file = cfg, ml_tree_file = ml, bs_trees_file = bs,
    results_dir = file.path(tmp, "res"), treePL_out = file.path(tmp, "out"),
    num_bs = 1, numsites = 1000, outgroup = c("C", "D"), seed = 12345L,
    treepl_bin = falso))), silent = TRUE)

  expect_true(file.exists(marca))
})

test_that("the worker hands its binary to both stages", {
  # Structural, and declared as such: the cross validation stage needs a real treePL to get past
  # priming, so the test above exercises the bootstrap stage and this one pins the other call site.
  # deparse() breaks long calls across elements and they come back joined by several spaces, so the
  # patterns tolerate whitespace. Corrected on 2026-09-25 after the second of the two failed for
  # that reason while the behavioural test above, which is the decisive one, already passed
  cuerpo <- paste(deparse(body(PhyloCactus:::.automate_treePL_run)), collapse = " ")
  expect_match(cuerpo, "run_treePL_cv\\(.*treepl_bin\\s*=\\s*treepl_bin")
  expect_match(cuerpo, "run_treePL_direct\\(bs_cfg_abs, label, cwd = bs_folder,\\s*treepl_bin\\s*=\\s*treepl_bin")
})
