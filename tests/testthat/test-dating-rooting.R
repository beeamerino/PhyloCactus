# Rooting is imposed after the search, so these tests fix the two properties that decide where the
# root lands: the outgroup argument is a vector of terminals, and the maximum-likelihood tree and
# the bootstrap replicates are held to different standards of completeness.

test_that("automate_treePL() raises an error on an unrooted ML tree with a missing outgroup (real function call)", {
  skip_if_not_installed("ape")

  tmp_dir <- withr::local_tempdir()

  # cfg_file only needs to exist on disk (normalizePath(mustWork = TRUE) is evaluated first); its
  # content is irrelevant because automate_treePL() stops on the rooting-validation error before
  # ever reading it.
  cfg_file <- file.path(tmp_dir, "dummy.cfg")
  file.create(cfg_file)

  unrooted_tree_file <- file.path(tmp_dir, "unrooted.nwk")
  ape::write.tree(ape::read.tree(text = "(A:1,B:1,C:1,D:1);"), unrooted_tree_file)

  bs_trees_file <- file.path(tmp_dir, "bs.nwk") # never read; the function stops before this point
  file.create(bs_trees_file)

  expect_error(
    automate_treePL(
      cfg_file = cfg_file,
      ml_tree_file = unrooted_tree_file,
      bs_trees_file = bs_trees_file,
      results_dir = file.path(tmp_dir, "results"),
      treePL_out = file.path(tmp_dir, "out"),
      outgroup = "Nonexistent_outgroup"
    ),
    "strictly requires a rooted tree"
  )
})

test_that("automate_treePL() requires every rooting terminal to be present in the ML tree (real function call)", {
  skip_if_not_installed("ape")

  tmp_dir <- withr::local_tempdir()

  cfg_file <- file.path(tmp_dir, "dummy.cfg")
  file.create(cfg_file)

  # Only one of the two rooting terminals is in the tree. A partial rooting set must not be
  # silently accepted: rooting on the subset places the root on a different edge than the one
  # intended, and every calibration addressed by an MRCA shifts with it.
  unrooted_tree_file <- file.path(tmp_dir, "unrooted_partial.nwk")
  ape::write.tree(ape::read.tree(text = "(A:1,B:1,Portulaca_fulgens:1);"), unrooted_tree_file)

  bs_trees_file <- file.path(tmp_dir, "bs.nwk")
  file.create(bs_trees_file)

  expect_error(
    automate_treePL(
      cfg_file = cfg_file,
      ml_tree_file = unrooted_tree_file,
      bs_trees_file = bs_trees_file,
      results_dir = file.path(tmp_dir, "results"),
      treePL_out = file.path(tmp_dir, "out"),
      outgroup = c("Portulaca_fulgens", "Portulaca_grandiflora")
    ),
    "rooting terminal\\(s\\) are absent"
  )
})

test_that("automate_treePL() roots an unrooted ML tree on a clade of terminals, not just the first one (real function call)", {
  skip_if_not_installed("ape")
  skip_if(.Platform$OS.type != "unix", "requires a POSIX shell for system()")

  tmp_dir <- withr::local_tempdir()

  cfg_file <- file.path(tmp_dir, "dummy.cfg")
  writeLines("numsites = 100", cfg_file)

  # The configuration above carries no calibration, so priming produces nothing usable and
  # run_treePL_cv() stops there. That is what isolates this test to the rooting step: whatever
  # error surfaces must come from a stage after rooting, never from the "unrooted tree"
  # validation. Before 2026-09-02 the same isolation was obtained with a shell wrapper that
  # exited 1; the wrapper is gone and the priming stage now plays that role.

  # Fixed topology rather than a random coalescent tree, so the two Portulaca terminals are
  # guaranteed to form a clade and the root placement under test is deterministic.
  unrooted_with_outgroup <- file.path(tmp_dir, "unrooted_with_og.nwk")
  ape::write.tree(
    ape::read.tree(text = "((A:1,B:1):1,Portulaca_fulgens:1,Portulaca_grandiflora:1);"),
    unrooted_with_outgroup
  )

  bs_trees_file <- file.path(tmp_dir, "bs.nwk")
  file.create(bs_trees_file)

  rooting_set <- c("Portulaca_fulgens", "Portulaca_grandiflora")

  err_msg <- tryCatch({
    suppressMessages(automate_treePL(
      cfg_file = cfg_file,
      ml_tree_file = unrooted_with_outgroup,
      bs_trees_file = bs_trees_file,
      results_dir = file.path(tmp_dir, "results"),
      treePL_out = file.path(tmp_dir, "out"),
      outgroup = rooting_set,
      n_prime = 1L
    ))
    NULL
  }, error = function(e) conditionMessage(e))

  expect_false(is.null(err_msg))
  # Must fail downstream, never with the rooting-validation error, confirming root_on_clade()
  # succeeded with the whole vector rather than only its first element.
  expect_false(grepl("strictly requires a rooted tree", err_msg, fixed = TRUE))
  expect_true(grepl("PLACE THE LINES BELOW|treePL|error in running command", err_msg))

  # The rooted tree automate_treePL() wrote must carry the root on the Portulaca stem, which is
  # what keeps the clade monophyletic and keeps its MRCA off the root node.
  fixed_tree <- file.path(tmp_dir, "results", "ML_tree", "supportTree_raxml-ng_fixed.tree")
  expect_true(file.exists(fixed_tree))
  tr <- ape::read.tree(fixed_tree)
  expect_true(ape::is.rooted(tr))
  expect_true(ape::is.monophyletic(tr, rooting_set))
})

test_that("automate_treePL() fails loudly when the rooting terminals are not a clade (real function call)", {
  skip_if_not_installed("ape")

  tmp_dir <- withr::local_tempdir()

  cfg_file <- file.path(tmp_dir, "dummy.cfg")
  writeLines("numsites = 100", cfg_file)

  # The only internal edge separates {Portulaca_fulgens, A} from {Portulaca_grandiflora, B}, so no
  # edge roots the two Portulaca terminals as a clade. This cannot happen in the reference pipeline,
  # where the constraint scaffold forces Portulaca monophyletic, but a user supplying an arbitrary
  # rooting set must be stopped rather than handed an arbitrary root placement.
  interleaved <- file.path(tmp_dir, "interleaved.nwk")
  ape::write.tree(
    ape::read.tree(text = "((Portulaca_fulgens:1,A:1):1,Portulaca_grandiflora:1,B:1);"),
    interleaved
  )

  bs_trees_file <- file.path(tmp_dir, "bs.nwk")
  file.create(bs_trees_file)

  expect_error(
    suppressMessages(automate_treePL(
      cfg_file = cfg_file,
      ml_tree_file = interleaved,
      bs_trees_file = bs_trees_file,
      results_dir = file.path(tmp_dir, "results"),
      treePL_out = file.path(tmp_dir, "out"),
      outgroup = c("Portulaca_fulgens", "Portulaca_grandiflora")
    )),
    "not monophyletic"
  )
})

test_that("root_on_clade() roots on terminals split across a basal polytomy, where ape::root() alone does not", {
  skip_if_not_installed("ape")

  # RAxML-NG writes unrooted Newick with a basal trifurcation. Here two of the three rooting
  # terminals sit on separate basal branches, so ape::root() does not see them as a clade even
  # though {Portulaca | rest} is a valid bipartition of the unrooted topology.
  tr <- ape::read.tree(
    text = "(Portulaca_a:1,(Portulaca_b:1,Portulaca_c:1):1,((A:1,B:1):1,C:1):1);"
  )
  og <- c("Portulaca_a", "Portulaca_b", "Portulaca_c")

  rooted <- root_on_clade(tr, og)

  expect_true(ape::is.rooted(rooted))
  expect_true(ape::is.monophyletic(rooted, og))
  # The rooting clade must hang below the root, not span it.
  expect_false(ape::getMRCA(rooted, og) == ape::Ntip(rooted) + 1L)
})

test_that("root_on_clade() errors when the terminals are not a clade of the unrooted topology", {
  skip_if_not_installed("ape")

  # The only internal edge separates {Portulaca_a, A} from {Portulaca_b, B}, so no edge roots the
  # two Portulaca terminals together. Erroring is the right outcome: any root placed here would be
  # arbitrary, and every calibration addressed by an MRCA would inherit that arbitrariness.
  tr <- ape::read.tree(text = "((Portulaca_a:1,A:1):1,Portulaca_b:1,B:1);")

  expect_error(
    root_on_clade(tr, c("Portulaca_a", "Portulaca_b")),
    "not monophyletic"
  )
})

test_that("root_on_clade() uses only the terminals present in the tree", {
  skip_if_not_installed("ape")

  tr <- ape::read.tree(text = "((A:1,B:1):1,Portulaca_a:1,Portulaca_b:1);")

  rooted <- root_on_clade(tr, c("Portulaca_a", "Portulaca_b", "Portulaca_absent"))
  expect_true(ape::is.monophyletic(rooted, c("Portulaca_a", "Portulaca_b")))

  expect_error(root_on_clade(tr, "Nothing_here"), "None of the rooting terminals")
})
