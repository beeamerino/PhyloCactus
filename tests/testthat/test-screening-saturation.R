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
