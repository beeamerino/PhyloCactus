test_that("clean_ambiguous preserves IUPAC ambiguity codes by default (real function call)", {
  skip_if_not_installed("Biostrings")

  # The default must be preserve_iupac = TRUE. RAxML-NG and ModelTest-NG treat an
  # ambiguity code as a partial constraint on the state (R means A or G) whereas N
  # constrains nothing, so collapsing the codes discards information. This test is the
  # regression guard for that decision.
  x <- Biostrings::DNAStringSet(c(sp1 = "ACGTRYN-", sp2 = "ACGTSWKM"))

  out <- clean_ambiguous(x)
  chars <- unlist(strsplit(as.character(out), ""))

  expect_true(all(c("R", "Y", "S", "W", "K", "M") %in% chars))
  expect_true("-" %in% chars)
  expect_true("N" %in% chars)
  expect_equal(names(out), names(x))
  expect_equal(Biostrings::width(out), Biostrings::width(x))
})

test_that("clean_ambiguous collapses IUPAC codes to N when preserve_iupac = FALSE (real function call)", {
  skip_if_not_installed("Biostrings")

  x <- Biostrings::DNAStringSet(c(sp1 = "ACGTRYN-", sp2 = "ACGTSWKM"))

  out <- clean_ambiguous(x, preserve_iupac = FALSE)
  chars <- unlist(strsplit(as.character(out), ""))

  expect_false(any(c("R", "Y", "S", "W", "K", "M") %in% chars))
  # Gaps survive in both modes: only non-ACGT bases become N, never gaps.
  expect_true("-" %in% chars)
  expect_true(all(chars %in% c("A", "C", "G", "T", "N", "-")))
  expect_equal(Biostrings::width(out), Biostrings::width(x))
})

test_that("clean_ambiguous uppercases and never converts gaps to N (real function call)", {
  skip_if_not_installed("Biostrings")

  x <- Biostrings::DNAStringSet(c(lower = "acgtry-n"))

  preserved <- unlist(strsplit(as.character(clean_ambiguous(x)), ""))
  collapsed <- unlist(strsplit(as.character(clean_ambiguous(x, preserve_iupac = FALSE)), ""))

  expect_true(all(preserved %in% c("A", "C", "G", "T", "R", "Y", "N", "-")))
  expect_equal(sum(preserved == "-"), 1L)
  expect_equal(sum(collapsed == "-"), 1L)
  # In collapsed mode R and Y join the single N already present: three N in total.
  expect_equal(sum(collapsed == "N"), 3L)
})

test_that("get_alignment_stats and clean_ambiguous agree on what counts as ambiguous (real function calls)", {
  skip_if_not_installed("Biostrings")

  # The counter and the retainer share IUPAC_AMBIGUITY_CODES, so what is measured as
  # ambiguous must be exactly what preserve_iupac = TRUE keeps.
  x <- Biostrings::DNAStringSet(c(sp1 = "ACGTRYSWKMBDHV"))

  n_before <- get_alignment_stats(as.character(x))$n_sites_ambiguous
  n_after_preserved <- get_alignment_stats(as.character(clean_ambiguous(x)))$n_sites_ambiguous
  n_after_collapsed <- get_alignment_stats(as.character(clean_ambiguous(x, preserve_iupac = FALSE)))$n_sites_ambiguous

  expect_equal(n_before, 10L)
  expect_equal(n_after_preserved, 10L)
  expect_equal(n_after_collapsed, 0L)
})
