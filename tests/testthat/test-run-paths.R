make_run <- function(dir, files = character()) {
  dir.create(file.path(dir, "7_Phylogenetics"), recursive = TRUE, showWarnings = FALSE)
  for (f in files) {
    p <- file.path(dir, f)
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    file.create(p)
  }
  invisible(dir)
}

test_that("every returned path is absolute even before the file exists", {
  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)
  make_run(tmp)

  paths <- resolve_run_paths("7_Phylogenetics", "cactus")

  # The whole point of the resolver is that resuming in a fresh session does not depend on the
  # working directory being what it was when the previous stage ran.
  path_fields <- names(paths$exists)
  for (k in path_fields) {
    expect_true(grepl("^(/|[A-Za-z]:[/\\\\])", paths[[k]]), info = k)
  }
})

test_that("analysed_phy follows the reduction that actually happened", {
  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)

  # No reduced matrix on disk: the analysed matrix is the supermatrix itself.
  make_run(tmp, "6_Concatenated/concatenated_alignments/ALIGNMENT_supermatrix.phy")
  p1 <- resolve_run_paths("7_Phylogenetics", "cactus")
  expect_false(p1$was_reduced)
  expect_match(p1$analysed_phy, "ALIGNMENT_supermatrix\\.phy$")

  # RAxML-NG collapsed identical terminals: the analysed matrix is the reduced one.
  file.create(file.path(tmp, "7_Phylogenetics", "cactus.raxml.reduced.phy"))
  p2 <- resolve_run_paths("7_Phylogenetics", "cactus")
  expect_true(p2$was_reduced)
  expect_match(p2$analysed_phy, "cactus\\.raxml\\.reduced\\.phy$")
})

test_that("best_tree follows the search to the cluster and back", {
  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)
  make_run(tmp)

  # Cluster run: generate_ml_search_script() writes into ml_search/.
  file.create2 <- function(p) {
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE); file.create(p)
  }
  file.create2(file.path(tmp, "7_Phylogenetics", "ml_search", "cactus_search.raxml.bestTree"))
  hpc <- resolve_run_paths("7_Phylogenetics", "cactus")
  expect_true(hpc$ml_on_cluster)
  expect_match(hpc$best_tree, "ml_search/cactus_search\\.raxml\\.bestTree$")

  # A local best tree in the run root takes precedence, because that is where a local run writes.
  file.create(file.path(tmp, "7_Phylogenetics", "cactus_search.raxml.bestTree"))
  local <- resolve_run_paths("7_Phylogenetics", "cactus")
  expect_false(local$ml_on_cluster)
  expect_false(grepl("ml_search", local$best_tree))
})

test_that("the temporal bootstrap path matches the prefix the function is given", {
  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)
  make_run(tmp)

  paths <- resolve_run_paths("7_Phylogenetics", "cactus")

  # calculate_temporal_bootstraps() names its output after the prefix, not after the directory.
  # Reading "<dir>/<dir>.raxml.bootstraps" points at a file that is never written.
  expect_match(paths$temporal_bs, "cactus_temporal_bs/cactus_temporal\\.raxml\\.bootstraps$")
  expect_equal(paths$temporal_bs_prefix, "cactus_temporal")
  expect_equal(basename(paths$temporal_bs),
               paste0(paths$temporal_bs_prefix, ".raxml.bootstraps"))
})

test_that("renaming the run renames every derived path", {
  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)
  make_run(tmp)

  paths <- resolve_run_paths("7_Phylogenetics", "opuntia")
  derived <- c("validated_part", "best_models", "constraint_tree", "best_tree",
               "ml_trees", "all_bootstraps", "support_tree", "temporal_bs")
  for (k in derived) {
    expect_match(basename(paths[[k]]), "^opuntia", info = k)
  }
})

test_that("require names the module that has not run yet", {
  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)
  make_run(tmp)

  # The failure this replaces is "object 'analysed_phy' not found", which says nothing about
  # which stage is missing or where its output was expected.
  expect_error(
    resolve_run_paths("7_Phylogenetics", "cactus", require = "best_tree"),
    "not on disk yet"
  )
  expect_error(
    resolve_run_paths("7_Phylogenetics", "cactus", require = "best_tree"),
    "Module 8"
  )
  expect_error(
    resolve_run_paths("7_Phylogenetics", "cactus", require = "best_models"),
    "Module 7"
  )
})

test_that("require rejects a name that is not part of the run", {
  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)
  make_run(tmp)

  expect_error(
    resolve_run_paths("7_Phylogenetics", "cactus", require = "besttree"),
    "Unknown path name"
  )
})

test_that("require passes once the files exist", {
  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)
  make_run(tmp, c("6_Concatenated/concatenated_alignments/ALIGNMENT_supermatrix.phy",
                  "7_Phylogenetics/cactus_partitions_validated.txt",
                  "7_Phylogenetics/cactus_modeltest.part.aicc"))

  paths <- resolve_run_paths("7_Phylogenetics", "cactus",
                             require = c("analysed_phy", "validated_part", "best_models"))
  expect_s3_class(paths, "cactus_run_paths")
  expect_true(paths$exists$best_models)
})

# ---------------------------------------------------------------------------------------------
# resolve_rooting_outgroup()
# ---------------------------------------------------------------------------------------------

test_that("resolve_rooting_outgroup() returns the whole matching clade, sorted and deduplicated", {
  tips <- c("Talinum_paniculatum", "Opuntia_ficus-indica", "Talinum_fruticosum",
            "Talinum_paniculatum", "Anacampseros_kurtzii")

  expect_equal(
    resolve_rooting_outgroup(tips),
    c("Talinum_fruticosum", "Talinum_paniculatum")
  )
})

test_that("the trailing underscore excludes Portulacaria from a Portulaca pattern", {
  # Portulacaria (Didiereaceae) is neither Portulacaceae nor part of the rooting sample. A bare
  # "^Portulaca" pattern would pull it in and place the root on an unrelated family. Portulacaceae
  # is not in the default pattern because it is not a clade in the distributed tree, so the guard
  # is exercised here with an explicit pattern.
  tips <- c("Portulaca_oleracea", "Portulacaria_afra", "Portulacaria_armiana")

  expect_equal(resolve_rooting_outgroup(tips, pattern = "^Portulaca_"), "Portulaca_oleracea")
})

test_that("resolve_rooting_outgroup() errors rather than returning an empty set", {
  # An empty rooting set would leave every downstream tree unrooted, and treePL would fail far
  # from the cause. The error has to surface here.
  expect_error(
    resolve_rooting_outgroup(c("Opuntia_ficus-indica", "Mammillaria_polyedra")),
    "No rooting terminal matched"
  )
  expect_error(resolve_rooting_outgroup(character(0)), "non-empty character vector")
})

test_that("resolve_rooting_outgroup() honours an alternative pattern for other datasets", {
  tips <- c("Portulaca_oleracea", "Talinum_paniculatum", "Talinum_fruticosum")

  expect_equal(
    resolve_rooting_outgroup(tips, pattern = "^Talinum_"),
    c("Talinum_fruticosum", "Talinum_paniculatum")
  )
})

test_that("resolve_rooting_outgroup() matches Talinaceae terminals by default", {
  tips <- c("Portulaca_oleracea", "Talinum_paniculatum", "Talinella_microphylla", "Opuntia_ficus-indica")

  expect_equal(
    resolve_rooting_outgroup(tips),
    c("Talinella_microphylla", "Talinum_paniculatum")
  )
})

test_that("the default rooting set is a clade of the distributed tree", {
  skip_if_not_installed("ape")

  # This is what the old default failed to be. It spanned Portulaca and Talinaceae, which Cactaceae
  # separates in this topology (Portulaca is sister to Anacampserotaceae), so root_on_clade() had no
  # monophyletic group to place the root on. Talinaceae does form one, and it is the family that
  # subtends the root.
  tree_path <- system.file("extdata", "phylocactus_ml_tree.tree", package = "PhyloCactus")
  skip_if(tree_path == "", "requires the distributed ML tree")

  tr <- ape::read.tree(tree_path)
  tips <- resolve_rooting_outgroup(tr$tip.label)

  expect_gt(length(tips), 1L)
  expect_true(all(tips %in% tr$tip.label))
  expect_true(ape::is.monophyletic(tr, tips))

  # And the whole family is taken, not a subset of it: a partial clade roots the tree in the
  # middle of Talinaceae.
  expect_setequal(tips, grep("^(Talinum|Talinella)_", tr$tip.label, value = TRUE))
})
