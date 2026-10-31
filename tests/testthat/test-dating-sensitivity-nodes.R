# Closure of 0.5.2 into 0.5.3, phase 4e (BMM, 31-10): under one fixed-age calibration and no
# calibration inside Cactaceae, the only `mrca` line is the fixed node, so the sensitivity table
# held one constant age. report_smoothing_sensitivity() gains `nodes`, a named list of terminal
# sets whose ages are reported beside those of the calibrated nodes. Written before the change.

# A stand-in for treePL: it writes a fixed chronogram to the `outfile` of the configuration and
# logs the smoothing it was given, as treePL does.
.fake_treepl <- function(dir) {
  bin <- file.path(dir, "fake_treepl")
  writeLines(c(
    "#!/bin/sh",
    "out=$(grep '^outfile' \"$1\" | sed 's/.*= *//')",
    "lam=$(grep '^smooth' \"$1\" | sed 's/.*= *//')",
    "echo '(((A:1,B:1):2,C:3):4,(D:5,E:5):2);' > \"$out\"",
    "echo \"smoothing : $lam\""
  ), bin)
  Sys.chmod(bin, "0755")
  bin
}

.sens_cfg <- function(dir, mrca = TRUE) {
  cfg <- file.path(dir, "cfg")
  writeLines(c("treefile = t.tre", "numsites = 128",
               if (mrca) c("mrca = ABC A B C", "min = ABC 7", "max = ABC 7")), cfg)
  cfg
}

test_that("report_smoothing_sensitivity reports the ages of uncalibrated nodes given in `nodes`", {
  skip_if(.Platform$OS.type != "unix", "needs a shell script as treePL")
  tmp <- withr::local_tempdir()
  out <- suppressMessages(report_smoothing_sensitivity(
    .sens_cfg(tmp), smoothing_values = c(1e-5, 1e-2), treepl_bin = .fake_treepl(tmp),
    nodes = list(AB = c("A", "B"), DE = c("D", "E"), root = c("A", "B", "C", "D", "E")),
    work_dir = file.path(tmp, "wd"), out_csv = file.path(tmp, "sens.csv")))
  # The calibrated node first, then the given nodes in their order.
  expect_identical(unique(out$node), c("ABC", "AB", "DE", "root"))
  at <- function(n) unique(out$age_ma[out$node == n])
  expect_identical(at("ABC"), 3)
  expect_identical(at("AB"), 1)
  expect_identical(at("DE"), 5)
  expect_identical(at("root"), 7)
  expect_true(all(out$converged))
  expect_identical(nrow(utils::read.csv(file.path(tmp, "sens.csv"))), 8L)
})

test_that("with `nodes`, a configuration without `mrca` lines is accepted", {
  skip_if(.Platform$OS.type != "unix", "needs a shell script as treePL")
  tmp <- withr::local_tempdir()
  out <- suppressMessages(report_smoothing_sensitivity(
    .sens_cfg(tmp, mrca = FALSE), smoothing_values = c(1e-5, 1e-2), treepl_bin = .fake_treepl(tmp),
    nodes = list(AB = c("A", "B")), work_dir = file.path(tmp, "wd"), out_csv = file.path(tmp, "s.csv")))
  expect_identical(unique(out$node), "AB")
  expect_identical(unique(out$age_ma), 1)
})

test_that("`nodes` must be a named list of terminal sets, with names not taken by an mrca line", {
  tmp <- withr::local_tempdir()
  cfg <- .sens_cfg(tmp)
  expect_error(report_smoothing_sensitivity(cfg, treepl_bin = "true", nodes = c("A", "B")), "named list")
  expect_error(report_smoothing_sensitivity(cfg, treepl_bin = "true", nodes = list(c("A", "B"))), "named list")
  expect_error(report_smoothing_sensitivity(cfg, treepl_bin = "true", nodes = list(ABC = c("A", "B"))),
               "ABC")
  expect_error(report_smoothing_sensitivity(cfg, treepl_bin = "true", nodes = list(X = 1:2)), "character")
})

test_that("a node of `nodes` with fewer than two terminals in the chronogram gets no age and a warning", {
  skip_if(.Platform$OS.type != "unix", "needs a shell script as treePL")
  tmp <- withr::local_tempdir()
  expect_warning(
    out <- suppressMessages(report_smoothing_sensitivity(
      .sens_cfg(tmp), smoothing_values = c(1e-5, 1e-2), treepl_bin = .fake_treepl(tmp),
      nodes = list(AZ = c("A", "Z")), work_dir = file.path(tmp, "wd"), out_csv = file.path(tmp, "s.csv"))),
    "AZ")
  expect_true(all(is.na(out$age_ma[out$node == "AZ"])))
  expect_identical(unique(out$age_ma[out$node == "ABC"]), 3)
})
