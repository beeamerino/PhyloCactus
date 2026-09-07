# Regression tests for the calibration defect audited on 2026-08-28.
#
# treePL addresses a calibrated node by the MRCA of the terminals declared for it, so a correct
# row label guarantees nothing about where the bound lands. Two failures followed from that:
# rooting on a single Portulaca terminal collapsed the Portulacaceae MRCA onto the tree root, and
# the three Ramirez-Barahona rows carried stem ages while being mapped to crown nodes, with two of
# them duplicating a single shared stem node across two different families.
#
# These tests hold the data contract of calibrations_bounds.csv and the node identities that the
# contract depends on, using the distributed maximum-likelihood topology rather than a toy tree.

read_calibrations <- function() {
  path <- system.file("extdata", "calibrations_bounds.csv", package = "PhyloCactus")
  skip_if(!nzchar(path), "calibrations_bounds.csv not installed")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

read_constraints <- function() {
  path <- system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus")
  skip_if(!nzchar(path), "cactus_constraints.csv not installed")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

# Mirrors the mapping used by the Module 10 calibration block: `value` is a semicolon-separated
# list of taxon-set names, required because a stem node is shared by the families it subtends and
# cannot be addressed by a single family name.
tips_for <- function(constraints, column, value, tip_labels) {
  vals <- trimws(strsplit(value, ";")[[1]])
  tips_all <- unique(constraints[constraints[[column]] %in% vals, "Specie_name"])
  tips_all[tips_all %in% tip_labels]
}

test_that("calibrations_bounds.csv declares valid structure and active calibrations", {
  calibs <- read_calibrations()

  expect_true("node_type" %in% names(calibs))
  expect_true(all(calibs$node_type %in% c("stem", "crown")))

  active <- calibs[calibs$used_in_analysis == TRUE, ]

  expect_true("ACP_root" %in% active$mrca)

  # The rows that applied a stem age to a crown node must not be active again.
  expect_false(any(c("Cactaceae_mrca", "Anacampserotaceae_mrca", "Portulacaceae_mrca") %in% active$mrca))

  # No active row may reuse a label, since treePL keys its mrca/min/max triplets by that name.
  expect_equal(anyDuplicated(active$mrca), 0L)
})

test_that("active ACP_root calibration addresses the multi-family clade", {
  calibs <- read_calibrations()
  active <- calibs[calibs$used_in_analysis == TRUE, ]

  acp <- active[active$mrca == "ACP_root", ]
  expect_equal(nrow(acp), 1L)

  set_of <- function(row) sort(trimws(strsplit(row$value, ";")[[1]]))
  expect_equal(set_of(acp), c("Anacampserotaceae", "Cactaceae", "Portulacaceae"))
})

test_that("every active calibration resolves to at least two terminals of the distributed topology", {
  skip_if_not_installed("ape")
  tree_path <- system.file("extdata", "phylocactus_ml_tree.tree", package = "PhyloCactus")
  skip_if(!nzchar(tree_path), "reference topology not installed")

  calibs <- read_calibrations()
  constraints <- read_constraints()
  tip_labels <- ape::read.tree(tree_path)$tip.label
  active <- calibs[calibs$used_in_analysis == TRUE, ]

  for (i in seq_len(nrow(active))) {
    row <- active[i, ]
    tips <- tips_for(constraints, row$column, row$value, tip_labels)
    # A calibration resolving to fewer than two terminals is skipped by the config builder, which
    # drops the bound without reporting it.
    expect_gte(length(tips), 2L)
  }
})

test_that("rooting on the outgroup clade keeps it monophyletic and places ACP_root correctly", {
  skip_if_not_installed("ape")
  tree_path <- system.file("extdata", "phylocactus_ml_tree.tree", package = "PhyloCactus")
  skip_if(!nzchar(tree_path), "reference topology not installed")

  constraints <- read_constraints()
  ml_tree <- ape::read.tree(tree_path)
  tip_labels <- ml_tree$tip.label

  # Which lineage roots the tree depends on what the installed topology samples, and the invariant
  # to check changes with it. Branching on the data rather than hard-coding one sampling keeps this
  # test valid across the republication of inst/extdata, instead of failing for the wrong reason.
  has_talinaceae <- any(grepl("^(Talinum|Talinella)_", tip_labels))
  rooting_pattern <- if (has_talinaceae) "^(Talinum|Talinella)_" else "^Portulaca_"

  rooting_set <- resolve_rooting_outgroup(tip_labels, pattern = rooting_pattern)
  # root_on_clade() rather than ape::root(): the distributed topology is written as an
  # unrooted trifurcation with outgroup terminals on more than one basal branch, which
  # ape::root() alone does not resolve into a clade.
  rooted <- root_on_clade(ml_tree, rooting_set)
  root_node <- ape::Ntip(rooted) + 1L

  # Rooting on the whole clade keeps the rooting lineage monophyletic, and its MRCA off the root.
  # Rooting on one of its terminals is what destroyed both properties in the original pipeline, and
  # this holds whichever lineage is used.
  expect_true(ape::is.monophyletic(rooted, rooting_set))
  expect_false(ape::getMRCA(rooted, rooting_set) == root_node)

  acp_tips <- tips_for(constraints, "Family", "Cactaceae;Anacampserotaceae;Portulacaceae", tip_labels)
  stem_tips <- tips_for(constraints, "Family", "Cactaceae;Anacampserotaceae", tip_labels)

  if (has_talinaceae) {
    # Talinaceae roots the tree, so the ACP crown is an internal node with a parent branch. That is
    # what makes it identifiable rather than a parameter parked at the deepest split, and it is the
    # reason the family was sampled.
    expect_false(ape::getMRCA(rooted, acp_tips) == root_node)
    expect_true(ape::is.monophyletic(rooted, acp_tips))
  } else {
    # Without Talinaceae the ACP crown is the deepest node of the tree, because the only outgroup
    # available sits inside the clade the calibration addresses.
    expect_equal(ape::getMRCA(rooted, acp_tips), root_node)
  }

  # The MRCA of Cactaceae and Anacampserotaceae is a node distinct from the root under either
  # sampling. Whether it is a clade is a result of the analysis, not a property asserted here:
  # Ramirez-Barahona et al. (2020) recover it, Zuntini et al. (2024) and Kew Tree of Life release
  # 4.0 do not, and the quartet support for the alternative is 0.47 with 0.33 on a third topology.
  expect_false(ape::getMRCA(rooted, stem_tips) == root_node)
})

test_that("rooting on a single Portulaca terminal reproduces the original defect", {
  skip_if_not_installed("ape")
  tree_path <- system.file("extdata", "phylocactus_ml_tree.tree", package = "PhyloCactus")
  skip_if(!nzchar(tree_path), "reference topology not installed")

  ml_tree <- ape::read.tree(tree_path)
  tip_labels <- ml_tree$tip.label
  portulaca <- resolve_rooting_outgroup(tip_labels)
  skip_if(length(portulaca) < 2L, "needs more than one Portulaca terminal")

  rooted_on_one <- ape::root(ml_tree, outgroup = portulaca[[1]], resolve.root = TRUE)
  root_node <- ape::Ntip(rooted_on_one) + 1L

  # Documented here so the contrast is explicit rather than folklore: a single-terminal root
  # leaves the genus paraphyletic and collapses its MRCA onto the root, which is exactly how the
  # Portulacaceae bound came to constrain the age of the whole tree.
  expect_false(ape::is.monophyletic(rooted_on_one, portulaca))
  expect_equal(ape::getMRCA(rooted_on_one, portulaca), root_node)
})

