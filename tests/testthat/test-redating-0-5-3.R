# Closure of 0.5.2 into 0.5.3, phase 4e (re-dating of 3h, BMM 30-10 and 31-10): Module 10 dates
# under three calibration schemes, each one fixed age on ACP_root and no calibration inside
# Cactaceae; scheme 1 (Zuntini et al. 2024, 53.22 Ma) is the distributed chronogram. PROVENANCE.csv
# records its scheme, smoothing and treePL seed. Written before the change.

.rd_extdata <- function(name) {
  f <- system.file("extdata", name, package = "PhyloCactus")
  skip_if(!nzchar(f), paste(name, "not installed"))
  f
}

test_that("PROVENANCE.csv records the scheme, smoothing and treePL seed of the chronogram, the active row of calibrations_bounds.csv", {
  pv <- utils::read.csv(.rd_extdata("PROVENANCE.csv"), stringsAsFactors = FALSE)
  row <- pv[pv$canonical == "phylocactus_chronogram_ml.tree", ]
  expect_equal(nrow(row), 1L)
  for (k in c("calibration_scheme", "smoothing_treepl", "seed_treepl")) expect_true(k %in% names(pv), info = k)
  cal <- utils::read.csv(.rd_extdata("calibrations_bounds.csv"), stringsAsFactors = FALSE)
  expect_identical(row$calibration_scheme, cal$reference[cal$used_in_analysis])
  expect_equal(as.numeric(row$smoothing_treepl), 1e-05)
})
