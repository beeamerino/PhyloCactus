# -------------------------------------------------------------
# PhyloCactus: Tutorial 2 - Phylogenetics Pipeline: Inference & Dating
# -------------------------------------------------------------
# This script covers Stages 7 to 10 of the phylogenetic pipeline:
# Substitution Models, ML search, Bootstraps, and treePL dating.
# -------------------------------------------------------------
library(PhyloCactus)

# -------------------------------------------------------------
# Setup
# -------------------------------------------------------------
tutorial_dir <- "~/Desktop/PhyloCactus_Tutorial"
setwd(tutorial_dir)

output_dir <- "7_Phylogenetics"
dir.create(output_dir, showWarnings = FALSE)

# Get model test path. Note: Configure these in your .Renviron file
modeltest_path <- Sys.getenv("PATH_MODELTEST_NG", "modeltest-ng")
raxml_path     <- Sys.getenv("PATH_RAXML_NG", "raxml-ng")
treepl_path    <- Sys.getenv("PATH_TREEPL", "treePL")

# -------------------------------------------------------------
# Stage 7: Statistical Control of Mutational Heterogeneity
# -------------------------------------------------------------
# Preprocess partitions before substitution checks
preprocess_partitions(
  phy_matrix = "6_Concatenated/concatenated_alignments/ALIGNMENT_supermatrix.phy",
  part_file = "6_Concatenated/concatenated_alignments/PARTITION_iqtree.part",
  raxml_path = raxml_path,
  output_dir = output_dir
)

# Run ModelTest-NG to select substitution models for each partition
# Use alternative working directory to avoid modeltest error.
setwd(output_dir)
best_models <- run_modeltest_ng(
  modeltest_exec_path = modeltest_path,
  aln_file = "cactus_check.raxml.reduced.phy",
  part_file = "cactus_check.raxml.reduced.clean.partition",
  prefix = "cactus_phylo_modeltest",
  threads = 8
)
setwd(tutorial_dir) # return to original working directory

# -------------------------------------------------------------
# Stage 8: Enforcing Topological Constraints and Inferring Maximum-Likelihood Hypotheses
# -------------------------------------------------------------
# Assemble taxonomic classifications into a constraint scaffold
constraint_tree <- build_constraint_scaffold(
  alignment_path = file.path(output_dir, "cactus_check.raxml.reduced.phy"),
  constraints_csv_path = system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus"),
  output_dir = output_dir
)

# Execute maximum-likelihood search with RAxML-NG
ml_results <- calculate_ml_tree(
  raxml_bin_path = raxml_path,
  aln_file = file.path(output_dir, "cactus_check.raxml.reduced.phy"),
  part_file = file.path(output_dir, "cactus_phylo_modeltest.part.aicc"),
  constraint_file = file.path(output_dir, "cactus_constraints.tree"),
  outgroup = "Portulaca_fulgens",
  n_init_trees = "rand{25},pars{25}",
  seed = 1111,
  n_workers = 1,
  threads = 8, 
  output_dir = output_dir,
  prefix = "cactus_search"
)

# Calculate Robinson-Foulds distance among maximum-likelihood trees
cat("Calculating distances among maximum-likelihood trees...\n")
rf_dist <- calculate_rf_distances(
  raxml_bin_path = raxml_path,
  ml_trees_file = file.path(output_dir, "cactus_search.raxml.mlTrees"),
  output_dir = output_dir
)

# -------------------------------------------------------------
# Stage 9: Estimating Statistical Robustness via Bootstrap Resampling
# -------------------------------------------------------------
# RAxML-NG supports running Bootstraps either locally or via HPC array chunks.
# You can set run_local_bs to TRUE to run them now on your machine, 
# or FALSE to generate an HPC bash script to run them remotely.
run_local_bs <- FALSE

if (run_local_bs) {
  cat("\n1. Estimating bootstraps locally...\n")
  bs_dir <- file.path(output_dir, "local_bs")
  local_bs <- run_local_bootstraps(
    raxml_bin_path = raxml_path,
    aln_file = file.path(output_dir, "cactus_check.raxml.reduced.phy"),
    part_file = file.path(output_dir, "cactus_phylo_modeltest.part.aicc"),
    constraint_file = file.path(output_dir, "cactus_constraints.tree"),
    bs_trees = 700,
    outgroup = "Portulaca_fulgens",
    threads = 8,
    workers = 1,
    output_dir = bs_dir
  )
} else {
  cat("\n1. Generating BS script for HPC chunks...\n")
  bs_dir <- file.path(output_dir, "bs_chunks")
  bs_script <- generate_bootstrap_script(
    alignment_file = file.path(output_dir, "cactus_check.raxml.reduced.phy"),
    partition_file = file.path(output_dir, "cactus_phylo_modeltest.part.aicc"),
    constraint_file = file.path(output_dir, "cactus_constraints.tree"),
    outgroup = "Portulaca_fulgens",
    bs_per_rep = 500,
    max_reps = 2,   # E.g. to reach 1000 replicates using 2 scripts/jobs 
    threads = 120,
    workers = 20,
    output_dir = file.path(output_dir)
  )
}

# -------------------------------------------------------------------------
# Collect Bootstrap Replicates
# -------------------------------------------------------------------------
# Before running this step, ensure that your bootstrap replicates have finished.
# If you ran them on HPC, make sure all chunks completed successfully and
# optionally sync the files back to your local repository.

cat("\n2. Collecting bootstrap replicates...\n")
# bs_dir was defined in Step 1 depending on whether you ran locally or via HPC chunks
bs_concat_file <- collect_bootstraps(
  bs_dir = bs_dir,
  output_dir = output_dir
)
cat("   Bootstraps concatenated to:", bs_concat_file, "\n")

# -------------------------------------------------------------------------
# Check Bootstrap Convergence
# -------------------------------------------------------------------------
cat("\n3. Checking bootstrap convergence...\n")
converge_log <- check_bs_convergence(
  raxml_bin_path = raxml_path,
  bs_trees_file = file.path(output_dir, "cactus_ALL_bootstraps.tree"),
  bs_cutoff = 0.03,
  seed = 111,
  threads = 8,
  output_dir = output_dir,
  prefix = "cactus_bs_convergence"
)
cat("   Convergence log saved to:", converge_log, "\n")

# -------------------------------------------------------------------------
# Map TBE Supports onto Best Tree
# -------------------------------------------------------------------------
cat("\n4. Mapping TBE supports onto the best tree...\n")
support_tree <- map_branch_supports(
  raxml_bin = raxml_path,
  best_tree = file.path(output_dir, "cactus_search.raxml.bestTree"),
  bootstraps_file = file.path(output_dir, "cactus_ALL_bootstraps.tree"),
  metric = "tbe",
  output_dir = output_dir,
  prefix = "cactus_support"
)
cat("   Support tree generated at:", support_tree, "\n")

# Note: for Dating, run sequential temporal bootstraps constrained over the best tree.
cat("Temporal Bootstraps for date intervals...\n")
temporal_bs <- calculate_temporal_bootstraps(
  raxml_bin_path = raxml_path,
  aln_file = file.path(output_dir, "cactus_check.raxml.reduced.phy"),
  part_file = file.path(output_dir, "cactus_phylo_modeltest.part.aicc"),
  best_tree_file = file.path(output_dir, "cactus_search.raxml.bestTree"),
  outgroup = "Portulaca_fulgens",
  bs_trees = 100, # For demo speed. Adjust to 500 or 1000.
  threads = 8,
  output_dir = file.path(output_dir, "cactus_temporal_bs")
)

# -------------------------------------------------------------
# Stage 10: Accommodating Evolutionary Rate Heterogeneity via Penalized Likelihood (treePL)
# -------------------------------------------------------------
cat("\n=======================================================\n")
cat("Stage 10: Accommodating Evolutionary Rate Heterogeneity via treePL\n")
cat("=======================================================\n")

dating_dir <- "8_Dating"
dir.create(dating_dir, showWarnings = FALSE)

# 10.1 Prepare Calibrations
cat("Creating treePL calibrations configurations...\n")
calibs_all <- read.csv(system.file("extdata", "calibrations_bounds.csv", package = "PhyloCactus"))
calibs <- calibs_all[calibs_all$used_in_analysis == TRUE, ]

# Load constraints and phylogeny tips to map calibrations
constraints <- read.csv(system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus"))
ml_tree <- ape::read.tree(file.path(output_dir, "cactus_search.raxml.bestTree"))
tip_labels <- ml_tree$tip.label
num_sites <- 12700

cfg_lines <- c()
for (i in seq_len(nrow(calibs))) {
  row <- calibs[i, ]
  tips_all <- unique(constraints[constraints[[row$column]] == row$value, "Specie_name"])
  tips_in_tree <- tips_all[tips_all %in% tip_labels]
  
  if (length(tips_in_tree) < 2) next
  
  mrca_line <- paste("mrca =", row$mrca, paste(tips_in_tree, collapse = " "))
  min_line  <- sprintf("min = %s %f", row$mrca, row$min)
  max_line  <- sprintf("max = %s %f", row$mrca, row$max)
  
  cfg_lines <- c(cfg_lines, mrca_line, min_line, max_line)
}

treepl_cfg <- c(
  paste0("numsites = ", num_sites),
  cfg_lines,
  "nthreads = 8",
  "thorough"
)

calibrations_cfg_path <- file.path(dating_dir, "calibrations_treePL_fulltips.cfg")
writeLines(treepl_cfg, calibrations_cfg_path)
cat("   Calibrations compiled successfully to:", calibrations_cfg_path, "\n")

# 10.2 Run Fast Automated treePL Dating Pipeline
cat("\nRunning automated treePL wrapper script over maximum-likelihood tree and bootstrap replicates...\n")

# Provide the wrapper shell script provided in your extdata
# See reference: https://github.com/tongjial/treepl_wrapper
wrapper_sh_path <- system.file("extdata", "treepl_wrapper_v1.sh", package = "PhyloCactus")

# Note: For tutorial purposes, you can limit the number of bootstrap trees to process
# by setting `num_bs = 100` (or any other number). If not provided, it will process all
# available bootstrap trees. Here we set it to 100 for faster tutorial execution.
automate_treePL(
  cfg_file = calibrations_cfg_path,
  wrapper_sh = wrapper_sh_path,
  ml_tree_file = file.path(output_dir, "cactus_search.raxml.bestTree"),
  bs_trees_file = file.path(output_dir, "cactus_temporal_bs", "cactus_temporal_bs.raxml.bootstraps"),
  results_dir = file.path(dating_dir, "auto_results"),
  treePL_out = dating_dir,
  num_bs = 100
)

cat("\n--- Chronological Dating Results Summary ---\n")
if (file.exists(file.path(dating_dir, "BestTree_treePL.tree"))) {
  ml_chronogram <- ape::read.tree(file.path(dating_dir, "BestTree_treePL.tree"))
  cat("   Best ML Chronogram:\n")
  cat("     - File site:", file.path(dating_dir, "BestTree_treePL.tree"), "\n")
  cat("     - Root age:", max(ape::node.depth.edgelength(ml_chronogram)), "Mya\n")
  cat("     - Number of tips:", length(ml_chronogram$tip.label), "\n")
}
cat("--------------------------------------------\n")
