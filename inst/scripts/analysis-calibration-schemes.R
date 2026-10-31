# -------------------------------------------------------------
# PhyloCactus: sensitivity of the chronogram to the choice of anchor on the ACP node
# [OPTIONAL / EXPLORATORY SENSITIVITY TEMPLATE]
# -------------------------------------------------------------
# NOTE: This script is an optional exploratory analysis for methodological research into
# secondary calibration anchor sensitivity. It is NOT executed during standard package
# builds or required for main manuscript production.
#
# Every published age for any node inside Cactaceae is a secondary calibration. This analysis makes
# that dependence explicit instead of hiding it behind a single number: it dates the same tree, with
# the same seed and the same replicates, under each defensible choice of anchor for the crown of the
# Cactaceae + Anacampserotaceae + Portulacaceae clade, and reports what changes.
#
# Why that node. Its identity does not depend on how the ACPT quartet resolves: the most recent
# common ancestor of the three families is the same node whether Anacampserotaceae groups with
# Cactaceae, with Portulacaceae, or with neither. Every other deep node in this tree is
# topology-dependent.
#
# Why a fixed value and not an interval. treePL takes `min` and `max`; a point age is expressed
# by setting them equal. The reason first given for the fixed value, that ACP_root returned the
# maximum of the 27.81-53.37 HPD interval of Ramirez-Barahona et al. (2020) in 84 of 100
# replicates, was measured on the 0.4.x tree, rooted on Portulaca without Talinaceae, and does not
# hold on the 0.5.x tree rooted on Talinaceae: there the interval places ACP_root on a bound in
# 11 % (first anchor probe of the closure of 0.5.2) and 15 % of the replicates (the probe at the
# common smoothing 1e-05). The reason that holds on this tree is another. 41.82 Ma is the estimate
# of the preferred analysis and method of the source (RC_complete, BEAST). Estimating the node
# inside that source's interval (S2) gives Cactaceae ages that the smoothing alone cannot tell apart
# from those under S1. The older anchors (67.01 and 79.97 Ma) push the crowns of Cactoideae,
# Cacteae and Opuntioideae onto the maxima of their internal bounds. S1 stays the main scheme and
# S2 is reported beside it as the scheme that carries the uncertainty of the source.
#
# The four anchors below all address the same biological node, verified in each source:
#   S1  41.82  Ramirez-Barahona (2020), RC_complete, Portulacaceae Stem_BEAST, the median height
#              of the node on the maximum clade credibility tree of that analysis. Their preferred
#              analysis and their preferred method. Read from Data4_Ages.xlsx.
#   S3  67.01  Ramirez-Barahona (2020), RC_complete, Portulacaceae Stem_treePL. Same data, same
#              node, method matched to this study. Their paper does not describe or validate these
#              columns, so this is a sensitivity value and not a recommendation.
#   S4  53.22  Zuntini et al. (2024), Suppl. Tab. 3, Cactaceae stem, young tree (maximum age of
#              154 Ma at the angiosperm crown). treePL, 200 primary fossil calibrations.
#   S5  79.97  Zuntini et al. (2024), same node, old tree (maximum age of 247 Ma).
#   S2  27.81-53.37  the HPD95 interval, retained so the fixed anchors can be compared against an
#              interval calibration.
#
# S1 and S3 come from the same six analyses of the same matrix and differ only in dating method, so
# the gap between them measures the penalized-likelihood against relaxed-clock offset directly, on
# this node. S4 and S5 come from a different matrix and a different fossil set but the same method,
# and the gap between them measures sensitivity to the maximum imposed at the angiosperm crown.
#
# Runtime: Module 10 of the reference run took about 37 min on the Mac for these steps under S1,
# and the cross-validation of S2 alone took 67 min in the anchor probe of the closure of 0.5.2, so
# budget an hour or more per scheme on a comparable machine. Cross-validation
# is repeated for every scheme because the optimal smoothing need not be the same under a different
# anchor.
# -------------------------------------------------------------
library(PhyloCactus)

tutorial_dir <- "~/Desktop/PhyloCactus_Tutorial"
setwd(tutorial_dir)

sweep_dir <- "8_Dating/anchor_sweep"
dir.create(sweep_dir, recursive = TRUE, showWarnings = FALSE)

# The seed is fixed across schemes so that the only difference between runs is the anchor. treePL
# seeds both the cross-validation site partitioning and the simulated annealing that precedes each
# optimization, so an unfixed seed would confound the comparison with run-to-run variation.
SEED <- 95054141

# Inputs are resolved as in Module 10 of Tutorial 2, from the run itself and not as literal paths or
# numbers. treePL reads a branch as edge.length * numsites expected substitutions, so a site count
# that differs from the analyzed matrix biases the rate smoothing without any warning.
output_dir       <- "7_Phylogenetics"
run_prefix       <- "cactus"
supermatrix_file <- "6_Concatenated/concatenated_alignments/ALIGNMENT_supermatrix.phy"

paths <- resolve_run_paths(output_dir, run_prefix, supermatrix_file,
                           require = c("analysed_phy", "best_tree", "temporal_bs"))

ml_tree_file  <- paths$best_tree
bs_trees_file <- paths$temporal_bs

phy_header <- strsplit(trimws(readLines(paths$analysed_phy, n = 1)), "\\s+")[[1]]
num_sites  <- as.integer(phy_header[2])
message("Matrix analysed: ", paths$analysed_phy, " - ", phy_header[1], " taxa x ", num_sites, " sites")

calibs <- read.csv(system.file("extdata", "calibrations_bounds.csv", package = "PhyloCactus"),
                   stringsAsFactors = FALSE)
constraints <- read.csv(system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus"),
                        stringsAsFactors = FALSE)

# The tree is rooted on Talinaceae. Resolved once here and reused by every
# scheme, so the rooting cannot drift between schemes.
rooting_outgroup <- resolve_rooting_outgroup(ape::read.tree(ml_tree_file)$tip.label,
                                             pattern = "^(Talinum|Talinella)_")

# Each scheme is identified by the substring that appears in the `reference` column of the anchor
# row it activates. Selecting by reference, not by row number, keeps the script valid if the
# file is reordered.
schemes <- list(
  S1 = "Portulacaceae_Stem_BEAST_point_41.82",
  S2 = "Portulacaceae_stem_BEAST_HPD95_interval",
  S3 = "Portulacaceae_Stem_treePL_point_67.01",
  S4 = "Zuntini_2024_SupplTab3_Cactaceae_stem_young_tree",
  S5 = "Zuntini_2024_SupplTab3_Cactaceae_stem_old_tree"
)

# The three internal bounds of Hernandez-Hernandez (2014) are carried unchanged through every
# scheme. Section 8y measured that removing them moves the ages by 0.12 Ma at the median and 0.51
# at most, so they are not what determines the result; they are retained as a declared consistency
# check against the only study that estimates Cactaceae ages from its own matrix.
internal_rows <- grepl("^Hernandez-Hernandez_2014_Table1", calibs$reference) &
  calibs$mrca %in% c("Opuntioideae_mrca", "Cactoideae_mrca", "Cacteae_mrca")

build_scheme <- function(anchor_pattern) {
  out <- calibs
  out$used_in_analysis <- FALSE
  out$used_in_analysis[internal_rows] <- TRUE
  hit <- which(out$mrca == "ACP_root" & grepl(anchor_pattern, out$reference, fixed = TRUE))
  if (length(hit) != 1L) {
    stop("Anchor pattern '", anchor_pattern, "' matched ", length(hit),
         " ACP_root rows; it must match exactly one.")
  }
  out$used_in_analysis[hit] <- TRUE
  out
}

results <- list()

for (tag in names(schemes)) {
  message("\n=== Scheme ", tag, ": ", schemes[[tag]], " ===")
  scheme_calibs <- build_scheme(schemes[[tag]])
  scheme_dir <- file.path(sweep_dir, tag)
  dir.create(scheme_dir, recursive = TRUE, showWarnings = FALSE)

  # Written out so that each scheme carries the exact calibration table it ran with, independent of the state of the packaged file at the time the sweep was executed.
  write.csv(scheme_calibs, file.path(scheme_dir, "calibrations_used.csv"), row.names = FALSE)

  # The calibration block is written the same way Tutorial 2 writes it: one mrca line naming the
  # terminals, then min and max. A point age is min and max set equal, which is how treePL expresses
  # a fixed node age; it has no separate keyword for one.
  tip_labels <- ape::read.tree(ml_tree_file)$tip.label
  active <- scheme_calibs[scheme_calibs$used_in_analysis == TRUE, , drop = FALSE]
  cfg_lines <- c()
  for (i in seq_len(nrow(active))) {
    row <- active[i, ]
    tips_in_tree <- calibration_tips(constraints, row$column, row$value, tip_labels)
    if (length(tips_in_tree) < 2L) {
      warning("Scheme ", tag, ": calibration '", row$mrca, "' resolves to ", length(tips_in_tree),
              " terminal(s) and is dropped.", call. = FALSE)
      next
    }
    cfg_lines <- c(cfg_lines,
                   paste("mrca =", row$mrca, paste(tips_in_tree, collapse = " ")),
                   sprintf("min = %s %f", row$mrca, row$min),
                   sprintf("max = %s %f", row$mrca, row$max))
  }
  cfg_file <- file.path(scheme_dir, paste0("configure_", tag))
  writeLines(c(paste0("numsites = ", num_sites), cfg_lines, "nthreads = 8", "thorough"), cfg_file)

  # `outgroup` is passed for the same reason Tutorial 2 passes it: the RAxML-NG tree is unrooted, and
  # treePL dates whatever rooting it is handed. Rooting on the Talinaceae clade, not on a
  # single terminal is what keeps ACP_root an internal node instead of collapsing the anchor onto
  # the deepest split. Omitting it here would have made every scheme in the sweep incomparable with
  # the main run for a reason unrelated to the anchor being tested.
  dating_out_dir <- file.path(scheme_dir, "dating_outputs")
  automate_treePL(
    cfg_file      = cfg_file,
    ml_tree_file  = ml_tree_file,
    bs_trees_file = bs_trees_file,
    results_dir   = file.path(scheme_dir, "auto_results"),
    treePL_out    = dating_out_dir,
    numsites      = num_sites,
    outgroup      = rooting_outgroup,
    seed          = SEED
  )

  chrono <- ape::read.tree(file.path(dating_out_dir, "BestTree_treePL.tree"))
  reps   <- ape::read.tree(file.path(dating_out_dir, "bsTree_treePL.tree"))

  adherence <- report_bound_adherence(chrono, scheme_calibs, constraints, bootstraps = reps)
  adherence$scheme <- tag
  results[[tag]] <- adherence
  write.csv(adherence, file.path(scheme_dir, "bound_adherence.csv"), row.names = FALSE)
}

sweep_table <- do.call(rbind, results)
write.csv(sweep_table, file.path(sweep_dir, "TABLE_anchor_sweep_bound_adherence.csv"),
          row.names = FALSE)

# Reported ages under each scheme, with the interval taken across the replicates. These are the
# figures that belong in the manuscript: one row per clade, one column per anchor, so that a reader
# can see how much of each age is data and how much is the choice of anchor.
clades <- list(
  ACP_crown          = list(column = "Family",
                            value = "Cactaceae;Anacampserotaceae;Portulacaceae"),
  Cactaceae_crown    = list(column = "Family", value = "Cactaceae"),
  Cactoideae_crown   = list(column = "Subfam", value = "Cactoideae"),
  Cacteae_crown      = list(column = "Clade",  value = "Cacteae"),
  Opuntioideae_crown = list(column = "Subfam", value = "Opuntioideae")
)

age_of <- function(tree, column, value) {
  tips <- calibration_tips(constraints, column, value, tree$tip.label)
  if (length(tips) < 2L) return(NA_real_)
  node <- ape::getMRCA(tree, tips)
  unname(ape::branching.times(tree)[as.character(node)])
}

age_rows <- list()
for (tag in names(schemes)) {
  dating_out_dir <- file.path(sweep_dir, tag, "dating_outputs")
  chrono <- ape::read.tree(file.path(dating_out_dir, "BestTree_treePL.tree"))
  reps   <- ape::read.tree(file.path(dating_out_dir, "bsTree_treePL.tree"))
  for (nm in names(clades)) {
    spec <- clades[[nm]]
    point <- age_of(chrono, spec$column, spec$value)
    rep_ages <- vapply(reps, function(t) age_of(t, spec$column, spec$value), numeric(1))
    rep_ages <- rep_ages[is.finite(rep_ages)]
    age_rows[[length(age_rows) + 1L]] <- data.frame(
      scheme = tag, anchor = schemes[[tag]], clade = nm,
      point = point,
      ci95_low  = if (length(rep_ages)) unname(stats::quantile(rep_ages, 0.025)) else NA_real_,
      ci95_high = if (length(rep_ages)) unname(stats::quantile(rep_ages, 0.975)) else NA_real_,
      n_replicates = length(rep_ages),
      stringsAsFactors = FALSE
    )
  }
}
age_table <- do.call(rbind, age_rows)
write.csv(age_table, file.path(sweep_dir, "TABLE_anchor_sweep_ages.csv"), row.names = FALSE)

print(age_table)

# What to read in the output.
#
# If the ages under the five schemes are close to proportional to their anchors, the chronogram is
# a rescaling of one external number and should be reported as such. Penalized likelihood is not
# exactly scale invariant, because the roughness penalty is not, so exact proportionality is not
# expected; the question is how far from it the result falls.
#
# The fossil pollen of Ramirez-Arriaga et al. (2026), from 15.6 Ma evaporites of the Tehuacan
# Formation, includes tricolpate grains that the authors compare with Cephalocereus and Mammillaria.
# Placed on the clade joining those genera, it implies a minimum for the common ancestor of
# Pachycereeae and Cacteae, which sits close to the crown of Cactoideae, and a scheme returning that
# node younger than 15.6 Ma would conflict with it. The placement is an assumption: the paper gives
# an affinity, not an assignment to a clade. It is not the only fossil of the family: the same paper
# cites Eopuntia douglassii, an Eocene macrofossil of the Green River Formation compared with Opuntia,
# which this script does not evaluate.
#   Ramirez-Arriaga et al. 2026. Cactaceae fossil pollen grains recovered from 15.6-Ma evaporites in
#   south-central Mexico. Review of Palaeobotany and Palynology 344: 105463.
#   doi:10.1016/j.revpalbo.2025.105463
