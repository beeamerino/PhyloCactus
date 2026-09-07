test_that("compute_pairwise_metric() rejects tree pairs sharing fewer than min_tips taxa (real function call)", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")

  # t1 vs t2 share exactly 3 tips (A, B, C); t1 vs t3 share all 4 tips (A, B, C, D).
  t1 <- ape::read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
  t2 <- ape::read.tree(text = "((A:1,B:1):1,C:2);")
  t3 <- ape::read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
  tree_list <- list(t1 = t1, t2 = t2, t3 = t3)

  res_strict <- compute_pairwise_metric(tree_list, phangorn::RF.dist, min_tips = 4L)
  expect_equal(res_strict$common_n["t1", "t2"], 3L)
  # Below min_tips = 4: distance must be reported as NA, not silently computed on 3 tips.
  expect_true(is.na(res_strict$distance["t1", "t2"]))
  expect_equal(res_strict$common_n["t1", "t3"], 4L)
  expect_false(is.na(res_strict$distance["t1", "t3"]))

  res_lenient <- compute_pairwise_metric(tree_list, phangorn::RF.dist, min_tips = 3L)
  # With min_tips = 3, the t1-t2 pair now clears the threshold and gets a real distance.
  expect_false(is.na(res_lenient$distance["t1", "t2"]))
})

test_that("safe_keep_tip() and compute_pairwise_metric() prune trees to an identical, correct common taxon set (real function calls)", {
  skip_if_not_installed("ape")
  skip_if_not_installed("phangorn")

  t1 <- ape::rtree(8, tip.label = paste0("sp", 1:8))
  t2 <- ape::rtree(6, tip.label = paste0("sp", c(1:4, 7, 9)))
  t3 <- ape::rtree(5, tip.label = paste0("sp", c(1:4, 10)))
  tree_set <- list(t1 = t1, t2 = t2, t3 = t3)

  common_tips <- Reduce(intersect, lapply(tree_set, function(tr) tr$tip.label))
  expect_equal(sort(common_tips), paste0("sp", 1:4))

  # safe_keep_tip() is the actual pruning helper used inside validate_phylogenies().
  pruned <- lapply(tree_set, function(tr) safe_keep_tip(tr, common_tips))
  for (p in pruned) {
    expect_equal(sort(p$tip.label), sort(common_tips))
  }

  # Every pairwise comparison among the pruned trees must now report all 4 tips as shared.
  res <- compute_pairwise_metric(pruned, phangorn::RF.dist, min_tips = 4L)
  expect_true(all(res$common_n[upper.tri(res$common_n)] == 4L))
})
