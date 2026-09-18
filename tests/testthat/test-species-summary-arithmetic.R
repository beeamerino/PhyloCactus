# Regression tests for defect M2 of the 2026-09-10 audit.
#
# TABLE_dataset_species_summary.csv is written by Stage 4 and revised by Stage 6, metric by metric.
# The revision reached "Accepted unique ingroup species retained (Cactaceae)" and not the recovery
# rate derived from it, so the published table asserted 1022 species retained and a rate of 53.59 %,
# which is 1053/1965, the Stage 4 value. The rate is now derived from the count it reports, and the
# three identities the table has to satisfy are checked before it is written.

summary_fixture <- function(...) {
  base <- c(
    "Total accepted species in Cactaceae checklist (Focal Ingroup)" = "1965",
    "Accepted unique ingroup species retained (Cactaceae)"          = "1022",
    "Cactaceae focal species recovery rate (%)"                     = "52.01%",
    "Unique Anacampserotaceae outgroup species retained"            = "15",
    "Unique Portulacaceae outgroup species retained"                = "17",
    "Unique Talinaceae outgroup species retained"                   = "6",
    "Total unique outgroup species retained"                        = "38",
    "Total unique species in final dataset (Joint)"                 = "1060"
  )
  overrides <- c(...)
  if (length(overrides) > 0L) base[names(overrides)] <- overrides
  data.frame(metric = names(base), value = unname(base),
             details = "", stringsAsFactors = FALSE)
}

test_that("a summary whose three identities hold passes without comment", {
  expect_silent(ok <- .check_species_summary_arithmetic(summary_fixture()))
  expect_true(ok)
})

test_that("the three outgroup families have to sum to the outgroup total", {
  df <- summary_fixture("Unique Talinaceae outgroup species retained" = "5")
  expect_warning(ok <- .check_species_summary_arithmetic(df), "outgroup families")
  expect_false(ok)
})

test_that("ingroup plus outgroup has to equal the joint total", {
  df <- summary_fixture("Total unique species in final dataset (Joint)" = "1061")
  expect_warning(ok <- .check_species_summary_arithmetic(df), "joint total")
  expect_false(ok)
})

test_that("the published recovery rate has to match the retained count", {
  # This is the published contradiction, reproduced exactly: 1022 retained of 1965, reported as the
  # Stage 4 figure of 53.59 %, per cent sign included, which is how Stage 4 writes it.
  df <- summary_fixture("Cactaceae focal species recovery rate (%)" = "53.59%")
  expect_warning(ok <- .check_species_summary_arithmetic(df), "recovery rate")
  expect_false(ok)

  # Rounding of the published figure is not a contradiction.
  expect_silent(.check_species_summary_arithmetic(
    summary_fixture("Cactaceae focal species recovery rate (%)" = "52.01%")
  ))
})

test_that("the per cent sign does not make the rate unreadable to the check", {
  # Stage 4 writes "53.59%" and a bare as.numeric() turns that into NA. Read that way, the identity
  # that this function exists to guard would be skipped without a word on exactly the table that
  # carries the defect.
  with_sign <- summary_fixture("Cactaceae focal species recovery rate (%)" = "53.59%")
  without_sign <- summary_fixture("Cactaceae focal species recovery rate (%)" = "53.59")

  expect_warning(.check_species_summary_arithmetic(with_sign), "recovery rate")
  expect_warning(.check_species_summary_arithmetic(without_sign), "recovery rate")
})

test_that("a missing row is not reported as a contradiction", {
  # Stage 4 writes the table before Stage 6 revises it, so a check that fired on absent rows would
  # warn on every run for reasons that are not defects.
  df <- summary_fixture()
  df <- df[df$metric != "Total unique species in final dataset (Joint)", , drop = FALSE]
  expect_silent(ok <- .check_species_summary_arithmetic(df))
  expect_true(ok)
})
