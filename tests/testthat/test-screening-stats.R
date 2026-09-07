test_that("get_alignment_stats decouples gaps from missing data (real function call)", {
  # Two sequences with known composition, 10 characters each:
  #   sp1: ACGT--NNNN  -> 4 bases, 2 gaps, 4 missing
  #   sp2: ACGTACGT--  -> 8 bases, 2 gaps, 0 missing
  # get_alignment_stats() reports the MEAN of the per-sequence fractions, not the
  # fraction over all matrix cells, so the expected values are:
  #   missing: mean(4/10, 0/10) = 0.20
  #   gaps   : mean(2/10, 2/10) = 0.20
  seqs <- c(sp1 = "ACGT--NNNN", sp2 = "ACGTACGT--")

  stats <- get_alignment_stats(seqs)

  expect_equal(stats$n_sequences, 2L)
  expect_equal(stats$mean_fraction_missing, 0.20)
  expect_equal(stats$mean_fraction_gaps, 0.20)
  # Gaps must not be counted inside the missing fraction, which was the double-count
  # this function was corrected for.
  expect_false(isTRUE(all.equal(stats$mean_fraction_missing, 0.40)))
})

test_that("get_alignment_stats counts IUPAC ambiguity separately from N and from gaps (real function call)", {
  # ACGT + the six two-fold ambiguity codes: 6 of 10 characters are ambiguous.
  seqs <- c(amb = "ACGTRYSWKM", clean = "ACGTACGTAC")

  stats <- get_alignment_stats(seqs)

  # mean(6/10, 0/10) = 0.30
  expect_equal(stats$mean_fraction_ambiguous, 0.30)
  expect_equal(stats$n_sites_ambiguous, 6L)
  # Ambiguity is neither missing nor gap.
  expect_equal(stats$mean_fraction_missing, 0)
  expect_equal(stats$mean_fraction_gaps, 0)
})

test_that("get_alignment_stats treats N and gaps as non-ambiguous (real function call)", {
  seqs <- c(only_n = "NNNNNNNNNN", only_gap = "----------")

  stats <- get_alignment_stats(seqs)

  expect_equal(stats$mean_fraction_ambiguous, 0)
  expect_equal(stats$n_sites_ambiguous, 0L)
  expect_equal(stats$mean_fraction_missing, 0.5)
  expect_equal(stats$mean_fraction_gaps, 0.5)
})

test_that("get_alignment_stats returns the full schema on empty input (real function call)", {
  # The empty branch must return the same columns as the populated branch, otherwise
  # the rbind that assembles the per-marker manifest fails on any empty marker.
  empty_stats <- get_alignment_stats(character(0))
  full_stats <- get_alignment_stats(c(a = "ACGT"))

  expect_equal(names(empty_stats), names(full_stats))
  expect_equal(empty_stats$n_sequences, 0L)
  expect_true(is.na(empty_stats$mean_fraction_ambiguous))
})
