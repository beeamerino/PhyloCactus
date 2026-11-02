# Version 0.5.3 closes the audit of 26-10-2026 (closure plan, phase 5a). The package carries the
# version; PROVENANCE.csv records the seeds of the reference run read from its logs, the versions of
# the programs, the version of the package that published the run and the checklist it read. Written
# before the change.

.r53_extdata <- function(name) {
  f <- system.file("extdata", name, package = "PhyloCactus")
  skip_if(!nzchar(f), paste(name, "not installed"))
  f
}

test_that("the package is version 0.5.3", {
  expect_identical(as.character(utils::packageVersion("PhyloCactus")), "0.5.3")
})

test_that("PROVENANCE.csv records the seeds of the reference run read from its logs (M10)", {
  p <- utils::read.csv(.r53_extdata("PROVENANCE.csv"), stringsAsFactors = FALSE, colClasses = "character")
  seeds <- c(seed_ml_search = "1111", seed_bootstrap = "2150; 4300; 6450; 8600",
             seed_temporal_bootstrap = "87786", seed_treepl = "79992967")
  for (s in names(seeds)) {
    expect_true(s %in% names(p), info = s)
    if (s %in% names(p)) expect_identical(unique(p[[s]]), unname(seeds[s]), info = s)
  }
})

test_that("PROVENANCE.csv records the versions of the programs and of the package that published the run (M10)", {
  p <- utils::read.csv(.r53_extdata("PROVENANCE.csv"), stringsAsFactors = FALSE, colClasses = "character")
  versions <- c(mafft = "^7\\.526( \\(|$)", modeltest_ng = "^0\\.1\\.7( \\(|$)", phylotar = "^1\\.3\\.0( \\(|$)",
                blast = "^2\\.17\\.0( \\(|$)", decipher = "3\\.6\\.0.*3\\.8\\.1|3\\.8\\.1.*3\\.6\\.0",
                treeannotator = "^2\\.7\\.7$", run_package_version = "^0\\.5\\.0$")
  for (v in names(versions)) {
    expect_true(v %in% names(p), info = v)
    if (v %in% names(p)) {
      expect_length(unique(p[[v]]), 1L)
      expect_match(unique(p[[v]])[1], versions[[v]], info = v)
    }
  }
  # The versions already recorded before 0.5.3 stay.
  expect_true(all(c("raxml_ng", "treepl", "R_version", "package_version") %in% names(p)))
})

test_that("PROVENANCE.csv names the checklist the run read, with its md5 (Korotkova et al. 2021)", {
  p <- utils::read.csv(.r53_extdata("PROVENANCE.csv"), stringsAsFactors = FALSE, colClasses = "character")
  expect_true(all(c("checklist_file", "checklist_md5") %in% names(p)))
  skip_if(!all(c("checklist_file", "checklist_md5") %in% names(p)))
  expect_identical(unique(p$checklist_file), "CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx")
  expect_identical(unique(p$checklist_md5), "0ab7855cef6fe7db0e0105646d5fee08")
  expect_identical(unique(p$checklist_md5), PhyloCactus:::.phylocactus_checklist_md5)
})
