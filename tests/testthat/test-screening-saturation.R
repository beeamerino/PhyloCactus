test_that(".test_saturation_proxy() with the corrected method does not flag a saturation artifact for biologically realistic divergence (real function call)", {
  skip_if_not_installed("ape")

  set.seed(42)
  n_taxa <- 8
  n_sites <- 300
  ancestor <- sample(c("a", "c", "g", "t"), n_sites, replace = TRUE)
  seqs <- lapply(seq_len(n_taxa), function(i) {
    mut_pos <- sample(n_sites, size = round(0.12 * n_sites))
    child <- ancestor
    child[mut_pos] <- sample(c("a", "c", "g", "t"), length(mut_pos), replace = TRUE)
    child
  })
  names(seqs) <- paste0("sp_", seq_len(n_taxa))
  dna <- ape::as.DNAbin(seqs)

  res <- .test_saturation_proxy(dna, saturation_method = "corrected", saturation_flag_cutoff = 0.3)

  expect_equal(res$reason, "ok")
  expect_true(is.numeric(res$slope))
  # Under the corrected method (uncorrected p-distance regressed against a Gamma-corrected K80
  # distance), uncorrected distance grows more slowly than the model-corrected one, so the
  # regression slope must not exceed 1.0.
  expect_true(res$slope <= 1.0)
})

test_that(".test_saturation_proxy() with the legacy method reproduces the known slope >= 1.0 mathematical artifact (real function call)", {
  skip_if_not_installed("ape")

  set.seed(42)
  n_taxa <- 8
  n_sites <- 300
  ancestor <- sample(c("a", "c", "g", "t"), n_sites, replace = TRUE)
  seqs <- lapply(seq_len(n_taxa), function(i) {
    mut_pos <- sample(n_sites, size = round(0.12 * n_sites))
    child <- ancestor
    child[mut_pos] <- sample(c("a", "c", "g", "t"), length(mut_pos), replace = TRUE)
    child
  })
  names(seqs) <- paste0("sp_", seq_len(n_taxa))
  dna <- ape::as.DNAbin(seqs)

  res <- .test_saturation_proxy(dna, saturation_method = "legacy", saturation_flag_cutoff = 0.3)

  expect_equal(res$reason, "ok")
  # The legacy method regresses Gamma-corrected K80 against equal-rates K80; the Gamma-corrected
  # distance is always >= the equal-rates one, mechanically producing a slope >= 1.0 regardless
  # of true saturation. This is precisely why "corrected" is now the package default.
  expect_true(res$slope >= 1.0)
})

test_that(".test_saturation_proxy() reports insufficient_data for fewer than 5 pairwise comparisons (real function call)", {
  skip_if_not_installed("ape")

  # Two taxa yield exactly 1 pairwise comparison.
  dna <- ape::as.DNAbin(list(sp_1 = c("a", "c", "g", "t"), sp_2 = c("a", "c", "g", "t")))

  res <- .test_saturation_proxy(dna, saturation_method = "corrected")

  expect_equal(res$reason, "insufficient_data")
  expect_true(is.na(res$slope))
  expect_true(is.na(res$saturated))
})

test_that(".test_saturation_proxy() reports degenerate_fit instead of flagging a slope of machine zero", {
  skip_if_not_installed("ape")

  # The published screening table gave trnL-trnF as the only saturated locus of seventeen, on a
  # slope of -6.19e-17. That is not an eroded slope, it is a regression with no information in it,
  # read through isTRUE(slope < 0.3) as though it were one. Identical sequences reproduce the
  # condition: every pairwise distance is zero, so the fitted slope is zero to machine precision.
  n_sites <- 300
  one <- sample(c("a", "c", "g", "t"), n_sites, replace = TRUE)
  seqs <- rep(list(one), 8)
  names(seqs) <- paste0("sp_", seq_len(8))
  dna <- ape::as.DNAbin(seqs)

  res <- .test_saturation_proxy(dna, saturation_method = "corrected", saturation_flag_cutoff = 0.3)

  expect_equal(res$reason, "degenerate_fit")
  expect_true(is.na(res$saturated))
  expect_false(isTRUE(res$saturated))
})
