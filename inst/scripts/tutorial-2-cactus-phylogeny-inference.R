# -------------------------------------------------------------
# PhyloCactus: Tutorial 2 - Phylogenetics Pipeline: Inference & Dating
# -------------------------------------------------------------
# This script covers Stages 7 to 10 of the phylogenetic pipeline:
# Substitution Models, ML search, Bootstraps, and treePL dating.
#
# RESUMING AT A LATER STAGE. The stages are long, and a full run spans several sessions. Every
# stage below therefore begins by calling resolve_run_paths(), which reconstructs the absolute
# path of every file of the run from output_dir and run_prefix instead of relying on variables
# left in the session. Run this SETUP block, then jump to whichever stage you need: nothing
# carries over in memory, so nothing breaks with "object not found" after a restart.
# -------------------------------------------------------------
library(PhyloCactus)

# -------------------------------------------------------------
# SETUP (run this block first, in every session)
# -------------------------------------------------------------
tutorial_dir <- "~/Desktop/PhyloCactus_Tutorial"
setwd(tutorial_dir)

output_dir <- "7_Phylogenetics"
dir.create(output_dir, showWarnings = FALSE)

# Every file of the run carries this prefix, so the run is renamed by changing one value.
run_prefix <- "cactus"

# Where Module 6 left the concatenated supermatrix.
supermatrix_file <- "6_Concatenated/concatenated_alignments/ALIGNMENT_supermatrix.phy"
partition_file   <- "6_Concatenated/concatenated_alignments/PARTITION_raxml_style.txt"

# Get model test path. Note: Configure these in your .Renviron file
modeltest_path <- Sys.getenv("PATH_MODELTEST_NG", "modeltest-ng")
raxml_path     <- Sys.getenv("PATH_RAXML_NG", "raxml-ng")
treepl_path    <- Sys.getenv("PATH_TREEPL", "treePL")

# -------------------------------------------------------------
# Rooting terminals, declared once and reused by every step that needs them.
# -------------------------------------------------------------
# Rooting is imposed after the search, not during it. RAxML-NG returns an unrooted topology and
# --outgroup only places the named terminals first in the output; the root is set later by
# automate_treePL() with root_on_clade(). What the choice controls is which edge the root sits on,
# and that governs which groups are monophyletic, which nodes exist, and where each calibration
# lands, since treePL addresses nodes by the MRCA of the terminals declared for them.
#
# Root on the whole outgroup clade, not on one of its terminals. Naming a single terminal places
# the root INSIDE the clade: the remaining terminals of that lineage fall on the ingroup side, the
# lineage is left paraphyletic, and its crown node collapses onto the root. Any calibration
# addressed by the MRCA of that lineage then lands on the root instead of its intended node.
#
# Occupancy no longer selects a terminal, since none is selected, but it is still worth checking:
# a rooting clade whose terminals are all sparsely sampled yields a poorly supported root.
#
#   s <- read.csv("6_Concatenated/final_tables/TABLE_final_species_alignment_summary.csv")
#   s <- s[s$species_class == "outgroup", c("species", "retained_markers", "pct_markers")]
#   head(s[order(-s$retained_markers), ])
#
# WHY THE ROOT MOVED TO TALINACEAE (2026-09-04).
#
# Rooting on Portulacaceae, as earlier versions of this script did, asserted the topology
# (Portulacaceae,(Cactaceae,Anacampserotaceae)). With only those three families sampled that was
# not a hypothesis the data could test: once Cactaceae is constrained to be monophyletic and the
# root is placed in Portulacaceae, the grouping of Cactaceae with Anacampserotaceae follows
# arithmetically, because no fourth lineage is present that would allow any other resolution. The
# node so produced is the one the Cactaceae_Anacampserotaceae_stem calibration addressed, so that
# calibration rested on a relationship the matrix could not evaluate.
#
# The relationship is genuinely open. Ramirez-Barahona et al. (2020) recover Cactaceae with
# Anacampserotaceae, on a tree whose interfamilial relationships were left unconstrained.
# Zuntini et al. (2024) and Kew Tree of Life release 4.0 recover Anacampserotaceae with
# Portulacaceae, but the quartet support for that node is 0.47 with 0.33 on an alternative, so
# neither resolution is established. de Vos et al. (2025) describe the region as a zone of gene
# tree conflict.
#
# Sampling Talinaceae supplies the fourth lineage, moves the root outside the ACP clade, and leaves
# all three resolutions of the quartet available to the likelihood. It also gives ACP_root a parent
# branch: the node stops being the root of the tree and becomes identifiable, which is the point
# raised in section 8o of the decision record.
#
# The scaffold written by build_constraint_scaffold() places the four families in a basal polytomy
# for the same reason: a nested outgroup would impose one of the three resolutions.
#
# Sensitivity check (optional): re-root on Portulacaceae, repeat the dating, and compare the deep
# topology and the root age against the Talinaceae-rooted result.
#
#   por_tips <- grep("^Portulaca_", tr$tip.label, value = TRUE)
#   tr_alt   <- root_on_clade(tr, por_tips)
#
# Taxon names are read from the first field of each PHYLIP record. Adapt `pattern` to the outgroup
# lineage sampled in your own dataset.
supermatrix_taxa <- sub("\\s.*$", "", readLines(supermatrix_file)[-1])
rooting_outgroup <- resolve_rooting_outgroup(supermatrix_taxa[nzchar(supermatrix_taxa)],
                                             pattern = "^(Talinum|Talinella)_")
cat("Rooting on", length(rooting_outgroup), "terminals:",
    paste(rooting_outgroup, collapse = ", "), "\n")

# -------------------------------------------------------------
# Stage 7: Statistical Control of Mutational Heterogeneity
# -------------------------------------------------------------
# Preprocess partitions before substitution checks.
# model_handling = "force_dna" writes the datatype token DNA into the model field. RAxML-NG expands
# a bare DNA token into its own default (GTR+FC+G4m+B) and writes that expansion into the reduced
# partition file; handing that to ModelTest-NG would both pre-empt the model selection ModelTest-NG
# exists to perform and crash its partition parser.
preprocess_partitions(
  phy_matrix = supermatrix_file,
  part_file = partition_file,
  raxml_path = raxml_path,
  output_dir = output_dir,
  prefix = run_prefix,
  force_check = FALSE,
  model_handling = "force_dna"
)

# Resolve every path of the run from the naming convention. RAxML-NG --check collapses terminals
# whose concatenated sequences are identical and writes a reduced matrix; those terminals are not
# analysed and will be absent from the tree, so paths$analysed_phy points at the reduced matrix when
# one was written and at the supermatrix when it was not.
paths <- resolve_run_paths(output_dir, run_prefix, supermatrix_file)
print(paths)

dup_groups_file <- "6_Concatenated/logs_and_qc/SUPP_TABLE_identical_sequence_groups.csv"
if (file.exists(dup_groups_file)) {
  dup_groups <- utils::read.csv(dup_groups_file)
  cat("Identical-sequence groups exported by Module 6:", length(unique(dup_groups$group)),
      "covering", nrow(dup_groups), "terminals;",
      nrow(dup_groups) - length(unique(dup_groups$group)), "will be collapsed.\n")
}
if (isTRUE(paths$was_reduced)) {
  cat("Reduced matrix header (taxa sites):", readLines(paths$analysed_phy, n = 1), "\n")
} else {
  cat("No reduced matrix written: RAxML-NG kept every terminal of the supermatrix.\n")
}

# Run ModelTest-NG to select substitution models for each partition.
# ModelTest-NG resolves its output prefix relative to the working directory, so the call is made
# from inside output_dir. resolve_run_paths() already returns absolute paths, which matters because
# when RAxML-NG applied no reduction the analysed matrix still lives in 6_Concatenated.
setwd(output_dir)
run_modeltest_ng(
  modeltest_exec_path = modeltest_path,
  aln_file = paths$analysed_phy,
  part_file = paths$validated_part,
  prefix = paste0(run_prefix, "_modeltest"),
  threads = 8
)
setwd(tutorial_dir) # return to original working directory

# -------------------------------------------------------------
# Stage 8: Enforcing Topological Constraints and Inferring Maximum-Likelihood Hypotheses
# -------------------------------------------------------------
# Stage 8 needs the outputs of Stage 7. Asking for them by name turns a missing input into an
# error that says which module has not run yet.
paths <- resolve_run_paths(output_dir, run_prefix, supermatrix_file,
                           require = c("analysed_phy", "best_models"))

# Assemble taxonomic classifications into a constraint scaffold
build_constraint_scaffold(
  alignment_path = paths$analysed_phy,
  constraints_csv_path = system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus"),
  output_dir = output_dir,
  prefix = run_prefix
)

# -------------------------------------------------------------
# Execute the constrained maximum-likelihood search
# -------------------------------------------------------------
# RAxML-NG can run the search either on this machine or as a SLURM job on an HPC allocation.
# Set run_local_ml to TRUE to search now, or FALSE to write a batch script to submit on the cluster.
run_local_ml <- FALSE

# Starting trees, declared once so the local run and the cluster script search the same space.
ml_start_trees <- "rand{25},pars{25}"

if (run_local_ml) {
  # n_workers = NULL derives the worker count from threads and the starting tree count. Wall-clock
  # is governed by how many sequential rounds each worker runs, ceiling(n_trees / workers), not by
  # the thread count alone: leaving n_workers at 1 puts each of the 50 starting trees in its own
  # round, which is the single largest avoidable cost in this stage.
  calculate_ml_tree(
    raxml_bin_path = raxml_path,
    aln_file = paths$analysed_phy,
    part_file = paths$best_models,
    constraint_file = paths$constraint_tree,
    outgroup = rooting_outgroup,
    n_init_trees = ml_start_trees,
    seed = 1111,
    n_workers = NULL,
    min_threads_per_worker = 4L,
    threads = 8,
    output_dir = output_dir,
    prefix = paste0(run_prefix, "_search")
  )
} else {
  ml_script <- generate_ml_search_script(
    alignment_file = paths$analysed_phy,
    partition_file = paths$best_models,
    constraint_file = paths$constraint_tree,
    outgroup = rooting_outgroup,
    n_init_trees = ml_start_trees,
    seed = 1111,
    # --- SLURM & HPC Resource Configuration ---
    cluster_job_name = "cactus_ml",
    cluster_partition = "main",                        # Partition name (e.g., "main", "standard", "general")
    cluster_nodes = 1L,                               # Number of compute nodes requested
    threads = 75,                                     # CPU cores requested (--cpus-per-task)
    cluster_mem = "16G",                              # RAM memory allocation (--mem, ~11 GB max used)
    cluster_time = "02:00:00",                        # Walltime limit (--time HH:MM:SS, ~41 min real runtime)
    cluster_mail_user = Sys.getenv("MY_EMAIL", ""),   # Notification email address (reads MY_EMAIL from .Renviron)
    load_module = c("gcc/14.2.0-nlhpc", "openmpi/5.0.3-o", "raxml-ng/1.1.0-mpi-zen4-n"), # Cluster environment modules
    raxml_exec = "raxml-ng-mpi",                      # Executable binary in $PATH
    # --- Search Tuning ---
    workers = NULL,
    min_threads_per_worker = 3L,
    preparse = TRUE,
    output_dir = output_dir,
    prefix = paste0(run_prefix, "_search")
  )
  cat("Submit on cluster with: cd ml_search && sbatch run_ml_search.sh\n")
  cat("When the job finishes, sync results back and re-run from the Stage 9 block below.\n")
}

# Re-resolve now that Stage 8 has written its outputs. resolve_run_paths() finds the best tree in
# ml_search/ when the search ran on the cluster and in the run root when it ran locally, so the
# steps below are identical either way.
paths <- resolve_run_paths(output_dir, run_prefix, supermatrix_file)

# Calculate Robinson-Foulds distance among maximum-likelihood trees
cat("Calculating distances among maximum-likelihood trees...\n")
rf_dist <- calculate_rf_distances(
  raxml_bin_path = raxml_path,
  ml_trees_file = paths$ml_trees,
  output_dir = output_dir
)

# -------------------------------------------------------------
# Stage 9: Estimating Statistical Robustness via Bootstrap Resampling
# -------------------------------------------------------------
# Resuming here after a restart: run the SETUP block above, then continue from this line.
paths <- resolve_run_paths(output_dir, run_prefix, supermatrix_file,
                           require = c("analysed_phy", "best_models",
                                       "constraint_tree"))

# RAxML-NG supports running Bootstraps either locally or via HPC array chunks.
# You can set run_local_bs to TRUE to run them now on your machine,
# or FALSE to generate an HPC bash script to run them remotely.
run_local_bs <- FALSE

if (run_local_bs) {
  cat("\n1. Estimating bootstraps locally...\n")
  bs_dir <- file.path(output_dir, "local_bs")
  local_bs <- run_local_bootstraps(
    raxml_bin_path = raxml_path,
    aln_file = paths$analysed_phy,
    part_file = paths$best_models,
    constraint_file = paths$constraint_tree,
    bs_trees = 700,
    outgroup = rooting_outgroup,
    threads = 8,
    workers = 1,
    output_dir = bs_dir
  )
} else {
  cat("\n1. Generating BS script for HPC chunks...\n")
  bs_dir <- file.path(output_dir, "bs_chunks")
  bs_script <- generate_bootstrap_script(
    alignment_file = paths$analysed_phy,
    partition_file = paths$best_models,
    constraint_file = paths$constraint_tree,
    outgroup = rooting_outgroup,
    # --- SLURM & HPC Resource Configuration ---
    cluster_job_name = "cactus_bs",
    cluster_partition = "main",                        # Partition name (e.g., "main", "standard", "general")
    cluster_nodes = 1L,                               # Number of compute nodes requested
    threads = 40,                                     # CPU cores per array task (--cpus-per-task)
    cluster_mem = "10G",                              # RAM memory allocation per task (--mem, ~5.3 GB max used)
    cluster_time = "24:00:00",                        # Walltime limit per task (--time HH:MM:SS, ~17 h real runtime)
    cluster_mail_user = Sys.getenv("MY_EMAIL", ""),   # Notification email address (reads MY_EMAIL from .Renviron)
    load_module = c("gcc/14.2.0-nlhpc", "openmpi/5.0.3-o", "raxml-ng/1.1.0-mpi-zen4-n"), # Cluster environment modules
    raxml_exec = "raxml-ng-mpi",                      # Executable binary in $PATH
    # --- Bootstrap Array Settings ---
    bs_per_rep = 500,                                 # Bootstrap replicates per array task
    max_reps = 2,                                     # Number of array tasks (2 x 500 = 1000 total replicates)
    workers = 10,
    output_dir = file.path(output_dir)
  )
  cat("Submit on cluster with: cd bs_chunks && sbatch run_bs_chunks.sh\n")
  cat("When all array tasks finish, sync results back and continue with collect_bootstraps().\n")
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
  output_dir = output_dir,
  prefix = paste0(run_prefix, "_ALL_bootstraps")
)
cat("   Bootstraps concatenated to:", bs_concat_file, "\n")

# Pick up the file collect_bootstraps() just wrote.
paths <- resolve_run_paths(output_dir, run_prefix, supermatrix_file,
                           require = c("all_bootstraps", "best_tree"))

# -------------------------------------------------------------------------
# Check Bootstrap Convergence
# -------------------------------------------------------------------------
# autoMRE bootstopping criterion evaluated under FBP (Felsenstein's Bootstrap Proportions),
# the standard metric comparable across the systematic literature of Cactaceae.
cat("\n3. Checking bootstrap convergence (FBP autoMRE)...\n")
converge_log <- check_bs_convergence(
  raxml_bin_path = raxml_path,
  bs_trees_file = paths$all_bootstraps,
  bs_metric = "fbp",
  bs_cutoff = 0.03,
  seed = 111,
  threads = 8,
  output_dir = output_dir,
  prefix = paste0(run_prefix, "_bs_convergence")
)
print(converge_log)

# -------------------------------------------------------------------------
# Map Bootstrap Supports onto Best Tree (FBP Primary, TBE Secondary)
# -------------------------------------------------------------------------
# EDITORIAL & METHODOLOGICAL NOTE:
# 1. Primary metric: Felsenstein's Bootstrap Proportions (FBP). This is the canonical metric
#    reported across published Cactaceae phylogenies (e.g. Hernandez-Hernandez et al. 2014,
#    Arakaki et al. 2011). TBE systematically inflates support on large clades (median FBP 0.468
#    vs median TBE 0.783 in unconstrained nodes).
# 2. Secondary metric: Transfer Bootstrap Expectation (TBE; Lemoine et al. 2018), retained as a
#    complementary assessment of stability in supermatrices with missing data.
# 3. Constrained nodes warning: The 34 nodes enforced by cactus_constraints.tree exhibit 1.000
#    support by construction. Do not report empirical bootstrap values for constrained nodes.
cat("\n4. Mapping bootstrap supports onto the best tree...\n")

# Primary: FBP support (written to paths$support_tree: cactus_support.raxml.support)
support_tree <- map_branch_supports(
  raxml_bin = raxml_path,
  best_tree = paths$best_tree,
  bootstraps_file = paths$all_bootstraps,
  metric = "fbp",
  output_dir = output_dir,
  prefix = paste0(run_prefix, "_support")
)
cat("   Primary FBP support tree generated at:", support_tree, "\n")

# Secondary: TBE support (written to cactus_support_tbe.raxml.support)
support_tree_tbe <- map_branch_supports(
  raxml_bin = raxml_path,
  best_tree = paths$best_tree,
  bootstraps_file = paths$all_bootstraps,
  metric = "tbe",
  output_dir = output_dir,
  prefix = paste0(run_prefix, "_support_tbe")
)
cat("   Secondary TBE support tree generated at:", support_tree_tbe, "\n")

# Note: for Dating, run sequential temporal bootstraps constrained over the best tree.
run_local_temporal_bs <- FALSE  # Set to TRUE to execute locally on your machine

if (run_local_temporal_bs) {
  cat("Running temporal bootstraps locally on machine...\n")
  temporal_bs <- calculate_temporal_bootstraps(
    raxml_bin_path = raxml_path,
    aln_file = paths$analysed_phy,
    part_file = paths$best_models,
    best_tree_file = paths$best_tree,
    outgroup = rooting_outgroup,
    bs_trees = 100, # Adjust to 500 or 1000 for production.
    threads = 8,
    workers = 2,
    output_dir = paths$temporal_bs_dir,
    prefix = paths$temporal_bs_prefix
  )
} else {
  cat("Generating HPC SLURM submission script for temporal bootstraps...\n")
  hpc_temp_bs_script <- generate_temporal_bootstrap_script(
    alignment_file = paths$analysed_phy,
    partition_file = paths$best_models,
    best_tree_file = paths$best_tree,
    outgroup = rooting_outgroup,
    bs_trees = 500,
    threads = 40,
    workers = 10,
    threads_per_worker = 4L,
    preparse = TRUE,
    output_dir = output_dir,
    script_name = "run_temporal_bs.sh",
    cluster_job_name = "cactus_temp_bs",
    cluster_partition = "main",
    cluster_nodes = 1L,
    cluster_mem = "16G",
    cluster_time = "04:00:00",
    cluster_mail_user = Sys.getenv("MY_EMAIL", "")
  )
  cat("   Temporal bootstrap HPC script generated at:", hpc_temp_bs_script, "\n")
  cat("   Submit on cluster: cd cactus_temporal_bs && sbatch run_temporal_bs.sh\n")
}

# -------------------------------------------------------------
# Stage 10: Accommodating Evolutionary Rate Heterogeneity via Penalized Likelihood (treePL)
# -------------------------------------------------------------
cat("\n=======================================================\n")
cat("Stage 10: Accommodating Evolutionary Rate Heterogeneity via treePL\n")
cat("=======================================================\n")

# Resuming here after a restart: run the SETUP block above, then continue from this line.
paths <- resolve_run_paths(output_dir, run_prefix, supermatrix_file,
                           require = c("analysed_phy", "best_tree", "temporal_bs"))

dating_dir <- "8_Dating"
dir.create(dating_dir, showWarnings = FALSE)

# 10.1 Prepare Calibrations
cat("Creating treePL calibrations configurations...\n")
calibs_all <- read.csv(system.file("extdata", "calibrations_bounds.csv", package = "PhyloCactus"))
calibs <- calibs_all[calibs_all$used_in_analysis == TRUE, ]

# Load constraints and phylogeny tips to map calibrations
constraints <- read.csv(system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus"))
ml_tree <- ape::read.tree(paths$best_tree)
tip_labels <- ml_tree$tip.label

# numsites must be the length of the matrix RAxML-NG actually analysed, not a rounded figure.
# treePL uses it to convert branch lengths into expected substitution counts, so it is read from
# the header of the matrix Stage 7 reported as analysed. paths$analysed_phy already resolves to the
# reduced PHYLIP when terminals were collapsed and to the supermatrix when they were not.
phy_header <- strsplit(trimws(readLines(paths$analysed_phy, n = 1)), "\\s+")[[1]]
num_sites <- as.integer(phy_header[2])
cat("   Matrix analysed:", paths$analysed_phy, "-", phy_header[1], "taxa x", num_sites, "sites\n")

# Rooting terminals, declared once in Module 8 and re-derived here so this block also runs
# standalone after a restart.
if (!exists("rooting_outgroup")) rooting_outgroup <- resolve_rooting_outgroup(tip_labels, pattern = "^(Talinum|Talinella)_")

# treePL resolves every mrca entry by the MRCA of the listed terminals, so a correct label is
# no guarantee of a correct node. The rooted topology is reconstructed here to verify node
# identity before any bound is written.
rooted_ml <- root_on_clade(ml_tree, rooting_outgroup)
root_node <- ape::Ntip(rooted_ml) + 1L

# value accepts a semicolon-separated list, required by the two stem calibrations: a stem node
# is shared by the families it subtends and cannot be addressed by a single family name.
tips_for <- function(column, value) calibration_tips(constraints, column, value, tip_labels)

# Bounds are declared per node, but the nodes are nested and an ultrametric tree forces every
# ancestor to be older than each of its descendants. A set that is defensible row by row can still
# describe a tree that cannot exist, and treePL does not report the violation: it pins the
# offending node to a bound and returns a chronogram that reads like an estimate.
check_calibration_consistency(calibs, rooted_ml, constraints)

cfg_lines <- c()
for (i in seq_len(nrow(calibs))) {
  row <- calibs[i, ]
  tips_in_tree <- tips_for(row$column, row$value)

  if (length(tips_in_tree) < 2) next

  mrca_line <- paste("mrca =", row$mrca, paste(tips_in_tree, collapse = " "))
  min_line  <- sprintf("min = %s %f", row$mrca, row$min)
  max_line  <- sprintf("max = %s %f", row$mrca, row$max)

  cfg_lines <- c(cfg_lines, mrca_line, min_line, max_line)
}

# Node-identity assertions, rewritten on 2026-09-04 when Talinaceae entered the sampling.
#
# With Talinaceae rooting the tree, the root is the crown of the ACPT clade and ACP_root is an
# internal node with a parent branch, which is what makes it identifiable rather than a parameter
# parked at the end of the tree. The earlier assertion required ACP_root to BE the root and would
# now fail, correctly: it encoded the previous sampling.
#
# What still has to hold, and what these lines check:
#   1. The rooting clade is monophyletic. Rooting on a single terminal instead of the clade leaves
#      the lineage paraphyletic and collapses its crown onto the root, which silently reassigns
#      every bound addressed by that lineage.
#   2. ACP_root resolves to a node that is NOT the root, so the anchor is placed on an internal
#      node rather than on the deepest split of the tree.
#   3. ACP_root is monophyletic. If the likelihood placed Talinaceae inside the ACP clade, the
#      anchor would be addressing something other than the node it was written for.
stopifnot(ape::is.monophyletic(rooted_ml, grep("^(Talinum|Talinella)_", tip_labels, value = TRUE)))
acp_tips <- tips_for("Family", "Cactaceae;Anacampserotaceae;Portulacaceae")
stopifnot(ape::getMRCA(rooted_ml, acp_tips) != root_node)
stopifnot(ape::is.monophyletic(rooted_ml, acp_tips))

# The resolution of the ACPT quartet is a result of this run, not an assumption, so it is reported
# rather than asserted. Whether Cactaceae_Anacampserotaceae_stem can be reactivated in
# calibrations_bounds.csv depends on what this prints.
cact_anac <- tips_for("Family", "Cactaceae;Anacampserotaceae")
anac_por  <- tips_for("Family", "Anacampserotaceae;Portulacaceae")
cact_por  <- tips_for("Family", "Cactaceae;Portulacaceae")
cat("ACPT quartet resolution recovered by this run:\n")
cat("   Cactaceae + Anacampserotaceae monophyletic: ",
    ape::is.monophyletic(rooted_ml, cact_anac), "\n")
cat("   Anacampserotaceae + Portulacaceae monophyletic: ",
    ape::is.monophyletic(rooted_ml, anac_por), "\n")
cat("   Cactaceae + Portulacaceae monophyletic: ",
    ape::is.monophyletic(rooted_ml, cact_por), "\n")

# num_sites here is the true length of the analysed matrix, which is what belongs in a file
# meant to document the analysis. automate_treePL() rescales branch lengths and divides this
# value by the same factor before handing a copy to treePL, so the substitution count per
# branch stays correct. The copy it actually runs is written to auto_results/ML_tree/.
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
cat("\nRunning treePL priming, cross-validation and dating over maximum-likelihood tree and bootstrap replicates...\n")

# Priming, cross-validation and dating follow the empirical protocol of Maurin (2020,
# arXiv:2008.07054, CC BY 4.0) and run through run_treePL_cv(), in R. Until 2026-09-02 they were
# driven by a copy of the shell script at https://github.com/tongjial/treepl_wrapper, which was an
# earlier implementation of the same protocol; that script carries no licence, so it could not be
# redistributed inside this GPL-3 package, and it wrote `smoothing = ` where treePL reads
# `smooth = `, silently dating every chronogram at the built-in default of 10.
#
# Three defaults follow Maurin where the shell script did not: the priming parameters are the
# lowest rather than the most frequent, cross-validation uses `randomcv` rather than leave-one-out
# `cv`, and the smoothing grid reaches 1e-08 rather than stopping at 1e-04. Maurin reports optimal
# smoothing between 1e-06 and 1e-08 for a tree whose branch lengths were rescaled as they are here.

# Note: For tutorial purposes, you can limit the number of bootstrap trees to process
# by setting `num_bs = 100` (or any other number). If not provided, it will process all
# available bootstrap trees. Here we set it to 100 for faster tutorial execution.
automate_treePL(
  cfg_file = calibrations_cfg_path,
  ml_tree_file = paths$best_tree,
  bs_trees_file = paths$temporal_bs,
  numsites = num_sites,
  outgroup = rooting_outgroup,
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

  # A node whose age equals one of its own bounds was not estimated: penalized likelihood returned
  # the constraint. In the August 2026 run all five calibrated nodes came back on a bound and the
  # chronogram gave no sign of it.
  # The bootstrap chronograms are passed deliberately. A single tree cannot tell a node that was
  # estimated from one whose age the data cannot identify: on 2026-09-02 the maximum-likelihood
  # tree put ACP_root at 52.96 Ma, 0.41 Ma inside its upper bound, and this check called it
  # interior, while 96 of the 100 replicates returned the bound exactly.
  bs_chronograms <- file.path(dating_dir, "bsTree_treePL.tree")
  adherence <- report_bound_adherence(
    ml_chronogram, calibs, constraints,
    bootstraps = if (file.exists(bs_chronograms)) bs_chronograms else NULL
  )
  cat("\n   Calibrated node ages against their bounds:\n")
  print(adherence)
}
cat("--------------------------------------------\n")

# -------------------------------------------------------------
# Smoothing sensitivity
# -------------------------------------------------------------
# Cross-validation selects the rate-smoothing parameter, and for this dataset it selects 1e-04,
# the lowest value on the tested grid. A minimum on the edge of a grid is the boundary of the
# search rather than an optimum, so the analysis cannot rest on the selection alone and a reader
# is entitled to ask whether the ages are an artefact of it.
#
# This answers the question rather than arguing about it: the same tree, under the same
# calibrations, dated at five smoothing values spanning six orders of magnitude. Measured on
# 2026-09-02, the Cactaceae/Anacampserotaceae divergence moved 1.71 Ma across that range
# (40.51 to 42.22), and ACP_root returned its upper bound of 53.37 at every value, which is a
# property of the sampling rather than of the smoothing: nothing outside the three sampled
# families constrains the root.
#
# The methods statement is therefore short. Cross-validation selected 1e-04; a sensitivity
# analysis over 1e-04 to 100 moved the reported ages by less than 2 Ma; the table is supplementary.
#
# Roughly a minute per value on 1000 terminals. Each run is verified against treePL's own log:
# the keyword treePL reads is `smooth`, and a configuration line it does not recognise is
# discarded in silence, which is how every chronogram this project produced before 2026-09-02 came
# to be dated at the built-in default of 10.
sensitivity <- report_smoothing_sensitivity(
  cfg_file = file.path(dating_dir, "auto_results", "ML_tree", "configure_smooth_ML_tree"),
  smoothing_values = c(1e-4, 1e-2, 1, 10, 100)
)
print(sensitivity)
cat("--------------------------------------------\n")
