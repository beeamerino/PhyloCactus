# -------------------------------------------------------------
# PhyloCactus: the chronogram under the three calibration schemes of calibrations_bounds.csv
# [OPTIONAL ANALYSIS]
# -------------------------------------------------------------
# NOTE: This script is optional. It is not run during package builds, and Tutorial 2 already dates
# the scheme chosen in its Module 10 (scheme 1 by default, the distributed chronogram).
#
# Cactaceae lack a macrofossil record suitable for calibrating internal nodes, so the scale of the
# chronogram comes from one secondary calibration with a fixed age (min equal to max, how treePL
# expresses a fixed node age) on the crown of the clade of Anacampserotaceae, Cactaceae and
# Portulacaceae (ACP_root), and no node inside Cactaceae is calibrated. calibrations_bounds.csv
# holds three calibration schemes for that node, each identified by the suffix
# "_calibration_scheme_<n>" of its `reference`:
#   1  53.22 Ma  Zuntini et al. (2024), Supplementary Table 3, the stem of Cactaceae on the young tree
#                (maximum of 154 Ma at the angiosperm crown), dated with treePL. The main scheme.
#   2  41.82 Ma  Ramirez-Barahona et al. (2020), RC_complete, the median height of the stem node of
#                Portulacaceae on the maximum clade credibility tree of the BEAST analysis
#                (Data4_Ages.xlsx, column Stem_BEAST).
#   3  37.0 Ma   Arakaki et al. (2011), Table S4, applied as a fixed age and not as the interval of
#                31.90 to 42.10 Ma.
# Each age is that of the ACP crown in the source's own tree: on the topology of this tree, where
# Cactaceae is sister to Anacampserotaceae plus Portulacaceae, the stem of Cactaceae of Zuntini et al.
# (2024) is the ACP crown; in Ramirez-Barahona et al. (2020), where Cactaceae is sister to
# Anacampserotaceae, the stem of Portulacaceae is. Whatever the resolution of the three families,
# the most recent common ancestor of the three is the same node of this tree.
#
# Each scheme is dated as Module 10 of Tutorial 2 dates one: the same configuration with only the
# calibration line changed, automate_treePL() with priming, cross-validation over 1e+03 to 1e-14
# (anomalous random subsamples left out of the selection, cv_outlier_factor = 5), the
# maximum-likelihood tree and 100 temporal replicates at the selected smoothing, and
# report_smoothing_sensitivity() over eleven values, reading the root and the crowns of Cactaceae,
# Cactoideae, Cacteae and Opuntioideae. Each scheme is written to the folder Tutorial 2 uses for it,
# 8_Dating/<scheme folder>; a scheme whose folder already holds BestTree_treePL.tree is not dated
# again. The seed is that of the reference run, so scheme 1 reproduces the distributed chronogram.
# The re-dating of the reference run took 3.1 h on the Mac for four schemes, three at a time.
# -------------------------------------------------------------
library(PhyloCactus)

tutorial_dir <- "~/Desktop/PhyloCactus_Tutorial"
setwd(tutorial_dir)

SEED <- 79992967
scheme_dirs <- c("Zuntini2024_young", "RamirezBarahona2020", "Arakaki2011")
treepl_path <- Sys.getenv("PATH_TREEPL", "treePL")

# Inputs are resolved as in Module 10 of Tutorial 2, from the run itself and not as literal paths or
# numbers. treePL reads a branch as edge.length * numsites expected substitutions, so a site count
# that differs from the analyzed matrix biases the rate smoothing without any warning.
output_dir       <- "7_Phylogenetics"
run_prefix       <- "cactus"
supermatrix_file <- "6_Concatenated/concatenated_alignments/ALIGNMENT_supermatrix.phy"

paths <- resolve_run_paths(output_dir, run_prefix, supermatrix_file,
                           require = c("analysed_phy", "best_tree", "temporal_bs"))

phy_header <- strsplit(trimws(readLines(paths$analysed_phy, n = 1)), "\\s+")[[1]]
num_sites  <- as.integer(phy_header[2])
message("Matrix analysed: ", paths$analysed_phy, " - ", phy_header[1], " taxa x ", num_sites, " sites")

calibs_all <- read.csv(system.file("extdata", "calibrations_bounds.csv", package = "PhyloCactus"),
                       stringsAsFactors = FALSE)
constraints <- read.csv(system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus"),
                        stringsAsFactors = FALSE)
tip_labels <- ape::read.tree(paths$best_tree)$tip.label
tips_for <- function(column, value) calibration_tips(constraints, column, value, tip_labels)

# The tree is rooted on Talinaceae, resolved once and reused by every scheme.
rooting_outgroup <- resolve_rooting_outgroup(tip_labels, pattern = "^(Talinum|Talinella)_")

# The nodes read in every scheme: the root (the crown of the ACPT clade) and four crowns.
nodes <- list(
  root = tip_labels,
  Cactaceae_crown = tips_for("Family", "Cactaceae"),
  Cactoideae_crown = tips_for("Subfam", "Cactoideae"),
  Cacteae_crown = tips_for("Clade", "Cacteae"),
  Opuntioideae_crown = tips_for("Subfam", "Opuntioideae")
)
age_of <- function(tree, tips) {
  d <- ape::node.depth.edgelength(tree)
  max(d[seq_len(ape::Ntip(tree))]) - d[ape::getMRCA(tree, intersect(tips, tree$tip.label))]
}

age_rows <- list()
for (k in seq_along(scheme_dirs)) {
  calibs <- calibs_all[grepl(paste0("_calibration_scheme_", k, "$"), calibs_all$reference), ]
  if (nrow(calibs) != 1L) stop("Scheme ", k, " matches ", nrow(calibs), " rows of calibrations_bounds.csv, not 1")
  dating_dir <- file.path("8_Dating", scheme_dirs[k])
  dir.create(dating_dir, recursive = TRUE, showWarnings = FALSE)
  message("\n=== Scheme ", k, ": ", calibs$reference, " ===")
  write.csv(calibs, file.path(dating_dir, "calibrations_used.csv"), row.names = FALSE)

  if (!file.exists(file.path(dating_dir, "BestTree_treePL.tree"))) {
    cfg_file <- file.path(dating_dir, "calibrations_treePL_fulltips.cfg")
    writeLines(c(paste0("numsites = ", num_sites),
                 paste("mrca =", calibs$mrca, paste(tips_for(calibs$column, calibs$value), collapse = " ")),
                 sprintf("min = %s %f", calibs$mrca, calibs$min),
                 sprintf("max = %s %f", calibs$mrca, calibs$max),
                 "nthreads = 8", "thorough"), cfg_file)
    automate_treePL(
      cfg_file      = cfg_file,
      ml_tree_file  = paths$best_tree,
      bs_trees_file = paths$temporal_bs,
      numsites      = num_sites,
      outgroup      = rooting_outgroup,
      results_dir   = file.path(dating_dir, "auto_results"),
      treePL_out    = dating_dir,
      num_bs        = 100,
      seed          = SEED
    )
    report_smoothing_sensitivity(
      cfg_file = file.path(dating_dir, "auto_results", "ML_tree", "configure_smooth_ML_tree"),
      smoothing_values = c(1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-5, 1e-4, 1e-2, 1, 10, 100),
      treepl_bin = treepl_path,
      nodes = nodes
    )
  } else {
    message("8_Dating/", scheme_dirs[k], " already holds BestTree_treePL.tree: not dated again")
  }

  # Ages on the maximum-likelihood chronogram, with the 2.5 and 97.5 percentiles over the replicates.
  chrono <- ape::read.tree(file.path(dating_dir, "BestTree_treePL.tree"))
  reps   <- ape::read.tree(file.path(dating_dir, "bsTree_treePL.tree"))
  for (nm in names(nodes)) {
    rep_ages <- vapply(reps, age_of, numeric(1), tips = nodes[[nm]])
    age_rows[[length(age_rows) + 1L]] <- data.frame(
      scheme = k, fixed_age = calibs$min, node = nm, age_ml = age_of(chrono, nodes[[nm]]),
      lo95 = unname(stats::quantile(rep_ages, 0.025)), hi95 = unname(stats::quantile(rep_ages, 0.975)),
      n_replicates = length(rep_ages), stringsAsFactors = FALSE)
  }
}
age_table <- do.call(rbind, age_rows)
write.csv(age_table, file.path("8_Dating", "TABLE_ages_calibration_schemes.csv"), row.names = FALSE)
print(age_table)

# What to read in the output.
#
# Each scheme fixes one age, so the three chronograms are close to rescalings of one another, and
# every age is conditional on the fixed age of its scheme. Penalized likelihood is not exactly scale
# invariant, because the roughness penalty is not, so exact proportionality is not expected.
#
# The fossil pollen of Ramirez-Arriaga et al. (2026), from 15.6 Ma evaporites of the Tehuacan
# Formation, includes tricolpate grains that the authors compare with Cephalocereus and Mammillaria.
# Placed on the clade joining those genera, it implies a minimum for the common ancestor of
# Pachycereeae and Cacteae, which on this tree is the crown of Cactoideae, and a scheme returning that
# node younger than 15.6 Ma would conflict with it. The placement is an assumption: the paper gives
# an affinity, not an assignment to a clade. It is not the only fossil of the family: the same paper
# cites Eopuntia douglassii, an Eocene macrofossil of the Green River Formation compared with Opuntia,
# which this script does not evaluate.
#   Ramirez-Arriaga et al. 2026. Cactaceae fossil pollen grains recovered from 15.6-Ma evaporites in
#   south-central Mexico. Review of Palaeobotany and Palynology 344: 105463.
#   doi:10.1016/j.revpalbo.2025.105463
