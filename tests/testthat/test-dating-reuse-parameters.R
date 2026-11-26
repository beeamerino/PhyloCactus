# AP-03 and AP-20 (closure review of 0.5.4). Until 0.5.4 a chronogram in results_dir was reused
# whenever the file existed, so a second call with other calibrations or another seed returned the
# old ages; leftover BS_* folders of a larger earlier call were collected; set.seed() was left in the
# session; and indented calibration lines were not copied into the replicate configurations.

make_inputs <- function(dir) {
  tr <- "(((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1):0.1,(Talinum_a:0.1,Talinum_b:0.1):0.1);"
  ml <- file.path(dir, "ml.tre"); writeLines(tr, ml)
  bs <- file.path(dir, "bs.tre"); writeLines(rep(tr, 3), bs)
  cfg <- file.path(dir, "cal.cfg")
  writeLines(c("numsites = 10000", "  mrca = AB A B", "  min = AB 1", "  max = AB 1"), cfg)
  list(ml = ml, bs = bs, cfg = cfg)
}

test_that("a chronogram without a record of its parameters is not reused", {
  dir <- withr::local_tempdir()
  x <- make_inputs(dir)
  res <- file.path(dir, "res")
  dir.create(file.path(res, "ML_tree"), recursive = TRUE)
  writeLines("((A:1,B:1):1,(C:1,D:1):1);", file.path(res, "ML_tree", "treepl_ML_tree.tre"))
  expect_error(
    suppressMessages(automate_treePL(x$cfg, x$ml, x$bs, res, file.path(dir, "out"), seed = 7,
                                     treepl_bin = "no-such-treepl")),
    "no record of the parameters"
  )
})

test_that("a chronogram produced with another seed or other calibrations is not reused", {
  dir <- withr::local_tempdir()
  x <- make_inputs(dir)
  path <- file.path(dir, "run_parameters_ML_tree.csv")
  sig <- list(seed = "7", cfg_file_md5 = "aaa", cvstart = 1000)
  PhyloCactus:::.write_treepl_signature(path, sig)
  expect_message(PhyloCactus:::.check_treepl_signature(path, sig, "the chronogram", dir), "Reusing")
  expect_error(PhyloCactus:::.check_treepl_signature(path, modifyList(sig, list(seed = "8")),
                                                      "the chronogram", dir),
               "seed: 7, now 8")
  expect_error(PhyloCactus:::.check_treepl_signature(path, modifyList(sig, list(cfg_file_md5 = "bbb")),
                                                      "the chronogram", dir),
               "cfg_file_md5")
})

test_that("the dating restores the session's random state and reads indented calibration lines", {
  cuerpo <- paste(deparse(body(PhyloCactus:::.automate_treePL_run)), collapse = " ")
  expect_match(cuerpo, ".bc_rng_state()", fixed = TRUE)
  expect_match(cuerpo, "(mrca|min|max)", fixed = TRUE)
  # Only the folders of the call are run and collected.
  expect_match(cuerpo, "bs_labels", fixed = TRUE)
})

test_that("reuse_unverified = TRUE reuses results without a record, with a warning, and never a mismatch", {
  dir <- withr::local_tempdir()
  path <- file.path(dir, "run_parameters_ML_tree.csv")
  expect_warning(PhyloCactus:::.check_treepl_signature(path, list(seed = "7"), "the chronogram", dir,
                                                        reuse_unverified = TRUE),
                 "reuse_unverified = TRUE")
  PhyloCactus:::.write_treepl_signature(path, list(seed = "7"))
  expect_error(PhyloCactus:::.check_treepl_signature(path, list(seed = "8"), "the chronogram", dir,
                                                      reuse_unverified = TRUE),
               "seed: 7, now 8")
  expect_false(eval(formals(automate_treePL)$reuse_unverified))
})
