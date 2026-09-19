# Tutorial 2: Phylogenetic Inference and Divergence Time Estimation

## Abstract

This tutorial covers the second stage of the `PhyloCactus` workflow
(Modules 7 to 10): substitution model selection per partition,
constrained maximum-likelihood inference, bootstrap support, and
divergence time estimation by penalized likelihood, starting from the
supermatrix assembled in [Tutorial
1](https://beeamerino.github.io/PhyloCactus/articles/tutorial-1-cactus-phylogeny-prep.html).

Taxon sampling and sequence availability are uneven across a matrix of
this size, while relationships among the major lineages of **Cactaceae**
are established by previous phylogenetic studies. The maximum-likelihood
search is therefore run under a topological constraint derived from
`cactus_constraints.csv`, which assigns the sampled species to higher
clades and subfamilies. The constraint fixes the relationships among
those clades and leaves the relationships within them to be estimated
from the data.

**Cactaceae** lack fossils suitable for calibrating the family directly.
Divergence times are estimated with secondary calibrations
(`calibrations_bounds.csv`): crown ages of major clades from
Hernández-Hernández *et al*. (2014), and an age for the ACP clade from
the angiosperm-wide chronogram of Ramírez-Barahona *et al*. (2020). The
resulting chronogram is conditional on those studies and on the models
selected here.

## Complete Pipeline Execution Workflow

### Setup: Creating a Clean Workspace

The external programs `ModelTest-NG`, `RAxML-NG` and `treePL` must be
installed and available on the system `PATH` before Module 7.

*On macOS, several of them can be installed with Homebrew
(`brew install raxml-ng`). If R reports “command not found”, pass the
absolute path of the executable in the corresponding argument (for
example, `raxml_path = "/usr/local/bin/raxml-ng"`).*

The complete script of this tutorial can be copied to the working
directory:

``` r

tutorial_dir <- "~/Desktop/PhyloCactus_Tutorial"
setwd(tutorial_dir)

# Copy the entire tutorial script for easy execution
file.copy(
  system.file("scripts", "tutorial-2-cactus-phylogeny-inference.R", package = "PhyloCactus"),
  file.path(tutorial_dir, "tutorial-2-cactus-phylogeny-inference.R")
)
```

### Module 7: Preprocess and Substitution Models

Loci differ in substitution pattern, base composition and rate, and a
single model applied to the whole matrix can describe them poorly and
bias likelihood-based inference. Each partition of the supermatrix is
therefore assigned its own substitution model.

The
**[`preprocess_partitions()`](https://beeamerino.github.io/PhyloCactus/reference/preprocess_partitions.md)**
function reads the supermatrix and partition file produced in Module 6,
checks their compatibility with the `RAxML-NG` parser, and writes the
input files for model selection and maximum-likelihood inference.

The modules of this tutorial run for hours, and a complete analysis
spans several R sessions, so each module has to run on its own. Passing
results between modules in variables fails after a restart, and typing
the file names at each call site lets them drift out of step with the
run prefix.

The
**[`resolve_run_paths()`](https://beeamerino.github.io/PhyloCactus/reference/resolve_run_paths.md)**
function reconstructs the path of every file of the run from the output
directory, the run prefix and the supermatrix location, and each module
below starts by calling it. Two entries depend on what happened in the
run: `analysed_phy` points to the reduced matrix when `RAxML-NG`
collapsed identical terminals and to the supermatrix otherwise, and
`best_tree` points to `ml_search/` when the maximum-likelihood search
ran on a cluster and to the run root when it ran locally. The `require`
argument names the entries a module needs, so that a missing input
raises an error naming the module that produces it.

The
**[`run_modeltest_ng()`](https://beeamerino.github.io/PhyloCactus/reference/run_modeltest_ng.md)**
function selects a substitution model for each partition with
`ModelTest-NG` (Darriba *et al*., 2020; Flouri *et al*., 2015) under the
corrected Akaike information criterion (AICc; Akaike, 1974; Hurvich &
Tsai, 1989). The selected scheme is the model used in the
maximum-likelihood search.

``` r

library(PhyloCactus)

# --- SETUP: run this block first, in every session -----------------------------------------
output_dir <- "7_Phylogenetics"
dir.create(output_dir, showWarnings = FALSE)

# Get model test path. Note: Configure these in your .Renviron file
modeltest_path <- Sys.getenv("PATH_MODELTEST_NG", "modeltest-ng")
raxml_path     <- Sys.getenv("PATH_RAXML_NG", "raxml-ng")
treepl_path    <- Sys.getenv("PATH_TREEPL", "treePL")

# Email notification for the long steps that run on this machine. Module 10 dates the tree and its
# bootstrap replicates locally, for hours, with no scheduler to mail on its behalf. Declared here,
# once, so the whole script is switched from one line.
#
# FALSE sends nothing and needs no configuration. TRUE needs two things, both outside this
# repository: MY_EMAIL in your .Renviron, and a blastula credentials file created once with
# blastula::create_smtp_creds_file() (use an application password if the account is Gmail). The
# package reads neither: it hands the path to blastula. A notification that cannot be sent is
# reported and ignored, so it never fails a run. See ?send_run_notification.
notify_email <- FALSE

# Every file of the run carries this prefix, so the run is renamed by changing one value.
run_prefix <- "cactus"

# Where Module 6 left the concatenated supermatrix.
supermatrix_file <- "6_Concatenated/concatenated_alignments/ALIGNMENT_supermatrix.phy"
partition_file   <- "6_Concatenated/concatenated_alignments/PARTITION_raxml_style.txt"

# Rooting terminals, declared once and reused by every step that needs them. The whole Talinaceae
# sample is used, not one terminal of it; see the note below for why. Taxon names are read from the
# first field of each PHYLIP record. Adapt `pattern` to the outgroup lineage sampled in your dataset.
supermatrix_taxa <- sub("\\s.*$", "", readLines(supermatrix_file)[-1])
rooting_outgroup <- resolve_rooting_outgroup(supermatrix_taxa[nzchar(supermatrix_taxa)],
                                             pattern = "^(Talinum|Talinella)_")
cat("Rooting on", length(rooting_outgroup), "terminals:",
    paste(rooting_outgroup, collapse = ", "), "\n")
# --- end of SETUP --------------------------------------------------------------------------

# model_handling = "force_dna" writes the datatype token DNA into the model field. RAxML-NG
# expands a bare DNA token into its own default (GTR+FC+G4m+B) and writes that expansion into the
# reduced partition file; handing that to ModelTest-NG would both pre-empt the model selection and
# crash its partition parser.
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
# whose concatenated sequences are identical and writes a reduced matrix; those terminals are
# absent from the tree, so `paths$analysed_phy` points at the reduced matrix when one was written
# and at the supermatrix when it was not.
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

# Run `ModelTest-NG` to select substitution models for each partition.
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
setwd("..")
```

### Module 8: Constraint Trees and Maximum Likelihood Search

The maximum-likelihood tree is the topology that maximizes the
likelihood of the alignment under the selected models and the
topological constraint. Taxon sampling is incomplete and uneven, and
previously established relationships among major lineages are used to
restrict the search to the part of tree space that remains unresolved.

The
**[`build_constraint_scaffold()`](https://beeamerino.github.io/PhyloCactus/reference/build_constraint_scaffold.md)**
function builds the constraint tree from `cactus_constraints.csv`, which
assigns the sampled species to subfamilies and major clades of
**Cactaceae** following published phylogenies, principally Guerrero *et
al*. (2019).

The constraint does not impose a complete topology. In the reference
dataset it fixes 34 of the 1021 internal bipartitions, the relationships
among 21 named clades, and leaves the other 987 to be estimated from the
data, including the placement of every species within its clade.

The
**[`calculate_ml_tree()`](https://beeamerino.github.io/PhyloCactus/reference/calculate_ml_tree.md)**
function runs the constrained search in `RAxML-NG` (Kozlov *et al*.,
2019) with the partition models selected in Module 7.

The maximum-likelihood tree is the reference topology for bootstrap
support and divergence time estimation.

The search runs on a local workstation with
[`calculate_ml_tree()`](https://beeamerino.github.io/PhyloCactus/reference/calculate_ml_tree.md),
or on a cluster with
**[`generate_ml_search_script()`](https://beeamerino.github.io/PhyloCactus/reference/generate_ml_search_script.md)**,
which writes a SLURM batch script for the same constrained search. The
`run_local_ml` flag below selects between them, as `run_local_bs` does
for the bootstrap replicates of Module 9. The search is a single job:
`RAxML-NG` distributes the starting trees among workers within the job
and writes one best tree.

The wall-clock time of this stage depends mainly on `n_workers`.
`RAxML-NG` gives each worker `threads / workers` threads and one
starting tree at a time, so the search takes
`ceiling(n_init_trees / workers)` rounds: 50 starting trees take 50
rounds on one worker and two rounds on 25 workers. On an Apple M2 Pro
with one worker and 8 threads, 50 starting trees on an earlier version
of the supermatrix took 47115 s (about 938 s per tree). On an AMD EPYC
9754 node with 25 workers and 75 threads (3 per worker), 50 starting
trees on the reference supermatrix (1024 taxa, 12809 sites, 5959
patterns, 11 partitions) took 2739 s (about 46 minutes).

With `n_workers = NULL` the worker count is chosen as a divisor of the
number of starting trees, so that no worker is idle in the last round;
16 workers on 50 trees run four rounds with two workers busy in the last
one, and finish no sooner than 25 workers on two full rounds. Each
worker should keep several hundred site patterns per thread, below which
parallel efficiency within a worker saturates. `min_threads_per_worker`
sets that floor; on unfamiliar hardware it can be checked with a short
trial search (`n_init_trees = "rand{2}"`).

The
**[`calculate_rf_distances()`](https://beeamerino.github.io/PhyloCactus/reference/calculate_rf_distances.md)**
function computes Robinson-Foulds distances (Robinson & Foulds, 1981)
among the trees of independent searches, as a measure of their
convergence.

The root is placed after the search. Under time-reversible substitution
models the likelihood of a topology is the same for all its rootings
(Felsenstein, 1981), so the data carry no information on the position of
the root, and constraint trees, like the search, are unrooted. The
`--outgroup` option of `RAxML-NG` only orders the output. The root is
placed by
[`automate_treePL()`](https://beeamerino.github.io/PhyloCactus/reference/automate_treePL.md)
through
[`root_on_clade()`](https://beeamerino.github.io/PhyloCactus/reference/root_on_clade.md),
and its position determines which groups are monophyletic, which nodes
exist and where each calibration falls, because `treePL` addresses nodes
by the MRCA of the terminals declared for them.

The root is placed on the whole outgroup clade. Rooting on a single
terminal places the root inside the clade: the other terminals of that
lineage fall on the ingroup side, the lineage becomes paraphyletic, and
its crown node collapses onto the root, so a calibration addressed by
the MRCA of that lineage falls on the root.
[`resolve_rooting_outgroup()`](https://beeamerino.github.io/PhyloCactus/reference/resolve_rooting_outgroup.md)
returns the full set, and
[`root_on_clade()`](https://beeamerino.github.io/PhyloCactus/reference/root_on_clade.md)
places the root on its stem edge and stops with an error when the
terminals are not a clade of the unrooted topology.
[`ape::root()`](https://rdrr.io/pkg/ape/man/root.html) is not used
directly because `RAxML-NG` writes the tree as an unrooted trifurcation,
and when the rooting terminals fall on more than one of its basal
branches, as in the reference dataset,
[`ape::root()`](https://rdrr.io/pkg/ape/man/root.html) does not
recognize them as a clade.

In the reference dataset the rooting set comprises the 8 sampled
terminals of Talinaceae (*Talinum* and *Talinella*); their occupancy is
reported in
`6_Concatenated/final_tables/TABLE_final_species_alignment_summary.csv`.
Rooting on Talinaceae places the root on the stem of the ACP clade
(Anacampserotaceae, Cactaceae, Portulacaceae). Rooting on Portulacaceae,
without a fourth lineage, would impose a sister relationship between
Cactaceae and Anacampserotaceae.

With this rooting the crown of the ACP clade is an internal node, and
the likelihood can choose among the three resolutions of the three
families (Ramírez-Barahona *et al.*, 2020; Zuntini *et al.*, 2024; de
Vos *et al.*, 2025). Relationships among these families are in conflict
across published analyses, and the resolution in the tree is an estimate
from this supermatrix.

``` r

# Module 8 needs the outputs of Module 7. Asking for them by name turns a missing input into an
# error that says which module has not run yet, instead of a missing object after a restart.
paths <- resolve_run_paths(output_dir, run_prefix, supermatrix_file,
                           require = c("analysed_phy", "best_models"))

# Assembles taxonomic classifications into a constraint scaffold
build_constraint_scaffold(
  alignment_path = paths$analysed_phy,
  constraints_csv_path = system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus"),
  output_dir = output_dir,
  prefix = run_prefix
)

# Execute the constrained maximum-likelihood search with `RAxML-NG`, either on this machine or as a
# SLURM job on an HPC allocation. Set run_local_ml to TRUE to search now, or FALSE to write a batch
# script to submit on the cluster.
run_local_ml <- FALSE

# Starting trees, declared once so that the local run and the cluster script search the same space.
ml_start_trees <- "rand{25},pars{25}"

if (run_local_ml) {
  # n_workers = NULL derives the worker count from threads and the starting tree count. Leaving it
  # at 1 puts every starting tree in its own sequential round, which is the single largest avoidable
  # cost in this stage.
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
  cat("When the job finishes, sync results back and resume from the Module 9 block below.\n")
}

# Re-resolve now that Module 8 has written its outputs. `resolve_run_paths()` finds the best tree in
# `ml_search/` when the search ran on the cluster and in the run root when it ran locally, so every
# step below is identical either way.
paths <- resolve_run_paths(output_dir, run_prefix, supermatrix_file)

# Calculate Robinson-Foulds distance among maximum-likelihood trees
cat("Calculating distances among maximum-likelihood trees...\n")
rf_dist <- calculate_rf_distances(
  raxml_bin_path = raxml_path,
  ml_trees_file = paths$ml_trees,
  output_dir = output_dir
)
```

### Module 9: Bootstrap Support

Bootstrap analysis measures how consistently each clade of the
maximum-likelihood tree is recovered when the alignment is resampled.

The
**[`generate_bootstrap_script()`](https://beeamerino.github.io/PhyloCactus/reference/generate_bootstrap_script.md)**
and
**[`run_local_bootstraps()`](https://beeamerino.github.io/PhyloCactus/reference/run_local_bootstraps.md)**
functions run the bootstrap searches in `RAxML-NG` with the same models,
constraint and partition scheme as the maximum-likelihood search.

The
**[`collect_bootstraps()`](https://beeamerino.github.io/PhyloCactus/reference/collect_bootstraps.md)**
function gathers the replicates, and
**[`check_bs_convergence()`](https://beeamerino.github.io/PhyloCactus/reference/check_bs_convergence.md)**
tests whether their number is sufficient with the autoMRE bootstopping
criterion (Pattengale *et al*., 2010) on Felsenstein bootstrap
proportions (FBP), with a cutoff of 0.03. In the reference run the
criterion was met after 650 of 1000 replicates.

The
**[`map_branch_supports()`](https://beeamerino.github.io/PhyloCactus/reference/map_branch_supports.md)**
function maps the bootstrap support onto the maximum-likelihood tree.
FBP is the primary measure, as in published phylogenies of Cactaceae,
which report bipartition frequencies or posterior probabilities (Arakaki
*et al*., 2011; Hernández-Hernández *et al*., 2014). The transfer
bootstrap expectation (TBE; Lemoine *et al*., 2018) is computed as a
secondary measure. TBE is less sensitive to unstable terminals in large
matrices with missing data, but it is higher than FBP on large clades:
on the 987 unconstrained nodes of the reference tree, clades of 51 to
200 terminals have a median FBP of 0.220 and a median TBE of 0.952. FBP
is reported for comparability with published studies, and TBE as a
measure of the stability of taxon placement.

Nodes imposed by the constraint (`cactus_constraints.tree`) are
recovered in every bootstrap replicate by construction and receive a
support of 1.000. That value reproduces the constraint and does not
measure support from the data; these nodes are reported as constrained
in publication tables.

#### Temporal Bootstrap Replicates for Divergence Time Estimation

Following Maurin (2020), a second set of bootstrap replicates is
generated for dating, to measure the uncertainty in branch lengths.

The
**[`calculate_temporal_bootstraps()`](https://beeamerino.github.io/PhyloCactus/reference/calculate_temporal_bootstraps.md)**
function resamples the alignment and estimates branch lengths for each
replicate with the topology fixed to the maximum-likelihood tree.

Because the topology is fixed, the dated replicates measure uncertainty
in branch lengths only. The confidence intervals of node ages do not
include topological uncertainty, which is measured separately by the
bootstrap support of Module 9.

Following Maurin (2020), the penalized likelihood parameters and the
smoothing value are selected on the maximum-likelihood tree with
`treePL` (Smith & O’Meara, 2012) and then applied to the temporal
replicates, whose dated trees give the confidence intervals of node
ages.

The two bootstrap procedures therefore measure different things: the
standard bootstrap measures topological support, and the temporal
bootstrap measures uncertainty in branch lengths for dating.

``` r

# -------------------------------------------------------------
# Estimate Bootstrap Replicates
# -------------------------------------------------------------
# Resuming here after a restart: run the SETUP block of Module 7, then continue from this line.
paths <- resolve_run_paths(output_dir, run_prefix, supermatrix_file,
                           require = c("analysed_phy", "best_models",
                                       "constraint_tree"))

# `RAxML-NG` supports running Bootstraps either locally or via HPC array chunks.
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
    bs_trees = 500,
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

# -------------------------------------------------------------
# Collect Bootstrap Replicates
# -------------------------------------------------------------
# Run this step after all bootstrap replicates have finished.
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

# -------------------------------------------------------------
# Check Bootstrap Convergence
# -------------------------------------------------------------
# autoMRE bootstopping criterion evaluated under FBP (Felsenstein's Bootstrap Proportions)
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

# -------------------------------------------------------------
# Map Bootstrap Supports onto Best Tree (FBP Primary, TBE Secondary)
# -------------------------------------------------------------
cat("\nMapping bootstrap supports onto the best tree...\n")

# Primary metric: FBP support (comparable with published literature)
support_tree <- map_branch_supports(
  raxml_bin = raxml_path,
  best_tree = paths$best_tree,
  bootstraps_file = paths$all_bootstraps,
  metric = "fbp",
  output_dir = output_dir,
  prefix = paste0(run_prefix, "_support")
)
cat("   Primary FBP support tree generated at:", support_tree, "\n")

# Secondary metric: TBE support (complementary missing-data assessment)
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
```

### Module 10: Estimating Divergence Times via Penalized Likelihood

Branch lengths of the maximum-likelihood tree are in substitutions per
site. Converting them to time requires calibrations, which here come
from previous dating studies.

Divergence times are estimated by penalized likelihood in `treePL`
(Sanderson, 2002; Smith & O’Meara, 2012), which allows rates to vary
among lineages and is computationally feasible for trees of thousands of
terminals.

Temporal constraints are provided through the `calibrations_bounds.csv`
reference file distributed with the package. Each row declares the node
it constrains as a set of taxa (`column` and `value`), the kind of node
it represents (`node_type`, either `stem` or `crown`), the minimum and
maximum age bounds, and the study the bounds were taken from. `treePL`
addresses nodes by the MRCA of the terminals listed for them, so the
taxon set, not the row label, is what determines where a bound lands.

The distinction between stem and crown nodes governs which bound may be
assigned to which node. A crown node is the first divergence within a
clade; a stem node is the divergence separating that clade from its
sister, and is therefore always older. Assigning a stem age to a crown
node inflates the node by construction, and propagates that inflation to
every node below it.

Internal calibrations are crown ages of Opuntioideae, Cactoideae and
Cacteae from Hernández-Hernández *et al.* (2014), applied to the
corresponding crown nodes.

Cactaceae lack a macrofossil record suitable for calibrating internal
nodes. The deepest calibration is therefore a secondary age from the
angiosperm-wide chronogram of Ramírez-Barahona *et al.* (2020),
`RC_complete` scheme. With Talinaceae sampled to root the tree, the
calibrated node is:

- `ACP_root`: Constrains the crown divergence of the Anacampserotaceae +
  Portulacaceae + **Cactaceae** clade (corresponding to the stem node of
  Portulacaceae in Ramírez-Barahona *et al.*, 2020). The primary
  calibration scheme (Scheme S1) applies a fixed point of 41.82 Ma (the
  median stem age under BEAST relaxed clock), with alternative
  sensitivity schemes available in `calibrations_bounds.csv` (Scheme S2:
  95% HPD interval \[27.81, 53.37\] Ma; Scheme S4: 53.22 Ma from Zuntini
  *et al.*, 2024 young tree).
- `Cactaceae_Anacampserotaceae_stem` (\[21.95, 48.51\] Ma) is retained
  in `calibrations_bounds.csv` with `used_in_analysis = FALSE`, because
  it assumes that Cactaceae and Anacampserotaceae are sisters. The
  maximum-likelihood tree recovers Anacampserotaceae and Portulacaceae
  as sisters (FBP 45; the node is not constrained), as do Zuntini *et
  al.* (2024).
- `ACPT_root` (\[33.28, 59.86\] Ma, Ramírez-Barahona *et al.*, 2020) is
  retained with `used_in_analysis = FALSE`. The root of the chronogram
  is the crown of the ACPT clade, so the deepest node of the tree is not
  calibrated: its age of 45.61 Ma is extrapolated above the fixed
  calibration and depends on the smoothing value. It falls inside the
  published interval, which is a consistency check and not a constraint.
  Activating the bound would add a second calibration on the deepest
  node and would be a different analysis. Methods should report the root
  age as uncalibrated.

The crown intervals of Ramírez-Barahona *et al.* (2020) for these
families are also listed with `used_in_analysis = FALSE`, for reference.
That study samples two terminals per family, so its crown ages are the
divergence between the two sampled species and underestimate the crown
of each family; its stem ages do not depend on within-family sampling
and are the ones used.

Secondary calibrations do not carry the uncertainty of the studies they
come from, so the node ages are conditional on those studies. The bounds
of Hernández-Hernández *et al.* (2014) and Ramírez-Barahona *et al.*
(2020) also come from analyses with different calibration schemes. Both
limitations should be stated in Methods.

The
**[`automate_treePL()`](https://beeamerino.github.io/PhyloCactus/reference/automate_treePL.md)**
function runs the dating: it prepares the maximum-likelihood tree and
the calibrations for `treePL`, selects the optimization parameters and
the smoothing value, and dates the maximum-likelihood tree and the
temporal replicates.

Following Maurin (2020), the optimization parameters and the smoothing
value are selected on the maximum-likelihood tree and applied unchanged
to every temporal replicate. Repeating the cross-validation for each
replicate is computationally prohibitive; as a consequence, the
confidence intervals do not include uncertainty in the smoothing value.

The temporal replicates of Module 9 are then dated with those
parameters.

The dated replicates give the distribution of each node age, from which
the confidence intervals are obtained.

This is the longest step that runs on the local machine, and it has no
scheduler to send a notification. With `notify_email <- TRUE` in the
SETUP block,
[`automate_treePL()`](https://beeamerino.github.io/PhyloCactus/reference/automate_treePL.md)
sends an email when the run ends, with its status, start and end times,
elapsed time, host, output paths and, if the run failed, the R error.
This requires `MY_EMAIL` in `.Renviron` and a `blastula` credentials
file created with
[`blastula::create_smtp_creds_file()`](https://rstudio.github.io/blastula/reference/create_smtp_creds_file.html);
`PhyloCactus` passes the path of that file to `blastula` and does not
read it. A notification that cannot be sent is reported and does not
fail the run. See
[`?send_run_notification`](https://beeamerino.github.io/PhyloCactus/reference/send_run_notification.md).

#### Summarizing Dated Temporal Bootstrap Replicates

The dated replicates are summarized with `TreeAnnotator` (Helfrich *et
al.*, 2018), following Maurin (2020).

All replicates share the maximum-likelihood topology, so the summary
estimates the distribution of node ages and does not compare topologies.

`TreeAnnotator` is run with a maximum sum of clade credibilities target
tree, mean node heights, and 0% burn-in. It computes the mean age of
each node and its interval across the dated replicates.

`TreeAnnotator` labels these intervals as highest posterior density
(HPD) intervals because it is written for Bayesian analyses. Here they
are confidence intervals from independently dated bootstrap replicates,
not posterior distributions.

#### Selecting the Rate-Smoothing Parameter

In penalized likelihood dating, the rate-smoothing parameter
($`\lambda`$) governs the roughness penalty that balances adherence to
observed branch lengths against the enforcement of rate autocorrelation
across ancestral-descendant lineages (Sanderson, 2002). When $`\lambda`$
approaches zero, each branch is permitted its own independent
substitution rate (approaching a clockless model); as $`\lambda`$
increases toward infinity, rate variation across lineages is penalized
until rates become identical across the entire tree (enforcing a strict
molecular clock).

Cross-validation (`randomcv`) removes subsets of the data, re-estimates
rates and times without them, and scores the prediction error of each
smoothing value as a chi-square ($`\chi^2`$); the value with the lowest
$`\chi^2`$ is selected (Maurin, 2020).
[`automate_treePL()`](https://beeamerino.github.io/PhyloCactus/reference/automate_treePL.md)
evaluates one value per order of magnitude from $`10^3`$ to
$`10^{-14}`$. The meaning of the selected value depends on the shape of
the curve, which
[`run_treePL_cv()`](https://beeamerino.github.io/PhyloCactus/reference/run_treePL_cv.md)
classifies and reports:

1.  **Interior minimum.** The curve decreases and then increases within
    the grid. The lowest value is an optimum, and the chronogram is
    dated at it.
2.  **Edge.** The lowest value lies at the end of the grid, and the
    curve is still decreasing there. The optimum lies beyond the grid;
    `cvstop` has to be extended and the analysis repeated before any age
    is interpreted (Maurin, 2020).
3.  **Plateau.** The curve decreases and then levels off, so that
    smoothing values below a threshold fit the data about equally well.
    The lowest $`\chi^2`$ is one of several equivalent values, and node
    ages have to be reported across the plateau.

The reference dataset produces a plateau. Between $`10^{-4}`$ and
$`10^{-6}`$ the $`\chi^2`$ decreases by 10.3 units; across the eight
orders of magnitude from $`10^{-6}`$ to $`10^{-14}`$ it decreases by 3.7
units, with three reversals. The lowest value, 221.12 at $`10^{-12}`$,
is 0.04 units below the value at $`10^{-14}`$. Cross-validation
therefore supports smoothing values below approximately $`10^{-6}`$
without distinguishing among them. The chronogram is dated at
$`10^{-12}`$, and
[`report_smoothing_sensitivity()`](https://beeamerino.github.io/PhyloCactus/reference/report_smoothing_sensitivity.md)
dates the same tree at values across and beyond the plateau. Across the
plateau ($`10^{-14}`$ to $`10^{-6}`$), the age of Opuntioideae varies by
0.15 Ma, that of Cacteae by 0.25 Ma, that of Cactoideae by 1.09 Ma, and
that of the uncalibrated root by 0.89 Ma. Methods should report the
plateau and these ranges.

The code for Module 10:

``` r

# Resuming here after a restart: run the SETUP block of Module 7, then continue from this line.
paths <- resolve_run_paths(output_dir, run_prefix, supermatrix_file,
                           require = c("analysed_phy", "best_tree", "temporal_bs"))

dating_dir <- "8_Dating"
dir.create(dating_dir, showWarnings = FALSE)

# 1. Prepare Calibrations
cat("Creating `treePL` calibrations configurations...\n")
calibs_all <- read.csv(system.file("extdata", "calibrations_bounds.csv", package = "PhyloCactus"))
calibs <- calibs_all[calibs_all$used_in_analysis == TRUE, ]

constraints <- read.csv(system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus"))
ml_tree <- ape::read.tree(paths$best_tree)
tip_labels <- ml_tree$tip.label

# `numsites` must be the length of the matrix RAxML-NG actually analysed, not a rounded figure.
# treePL uses it to convert branch lengths into expected substitution counts, so it is read from
# the header of the matrix Module 7 reported as analysed. `paths$analysed_phy` already resolves to the
# reduced PHYLIP when terminals were collapsed and to the supermatrix when they were not, so no
# second guess at the filename is needed here.
phy_header <- strsplit(trimws(readLines(paths$analysed_phy, n = 1)), "\\s+")[[1]]
num_sites <- as.integer(phy_header[2])
cat("   Matrix analysed:", paths$analysed_phy, "-", phy_header[1], "taxa x", num_sites, "sites\n")

# Rooting terminals, declared once in Module 8 and re-derived here so this block also runs
# standalone after a restart.
if (!exists("rooting_outgroup")) rooting_outgroup <- resolve_rooting_outgroup(tip_labels, pattern = "^(Talinum|Talinella)_")

# `treePL` resolves every `mrca` entry by the MRCA of the listed terminals, so a correct
# label is no guarantee of a correct node. The rooted topology is reconstructed here to
# verify node identity before any bound is written.
rooted_ml <- root_on_clade(ml_tree, rooting_outgroup)
root_node <- ape::Ntip(rooted_ml) + 1L

# `value` accepts a semicolon-separated list, which is required by the two stem
# calibrations: a stem node is shared by the families it subtends and cannot be addressed
# by a single family name.
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

# Node-identity assertions. With Talinaceae rooting the tree, the root is the crown of the ACPT
# clade and ACP_root is an internal node with a parent branch, which makes its age identifiable.
#
# What still has to hold, and what these lines check:
#   1. The rooting clade is monophyletic. Rooting on a single terminal instead of the clade leaves
#      the lineage paraphyletic and collapses its crown onto the root, which silently reassigns
#      every bound addressed by that lineage.
#   2. ACP_root resolves to a node that is NOT the root, so the anchor is placed on an internal
#      node and not on the deepest split of the tree.
#   3. ACP_root is monophyletic. If the likelihood placed Talinaceae inside the ACP clade, the
#      anchor would be addressing something other than the node it was written for.
stopifnot(ape::is.monophyletic(rooted_ml, grep("^(Talinum|Talinella)_", tip_labels, value = TRUE)))
acp_tips <- tips_for("Family", "Cactaceae;Anacampserotaceae;Portulacaceae")
stopifnot(ape::getMRCA(rooted_ml, acp_tips) != root_node)
stopifnot(ape::is.monophyletic(rooted_ml, acp_tips))

# The resolution of the ACPT quartet is estimated by this run, so it is reported and not asserted. Whether Cactaceae_Anacampserotaceae_stem can be reactivated in
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

# `num_sites` here is the true length of the analysed matrix, which is what belongs in a file
# meant to document the analysis. `automate_treePL()` rescales branch lengths and divides this
# value by the same factor before handing a copy to `treePL`, so the substitution count per
# branch stays correct. The copy it actually runs is written to `auto_results/ML_tree/`.
treepl_cfg <- c(
  paste0("numsites = ", num_sites),
  cfg_lines,
  "nthreads = 8",
  "thorough"
)

calibrations_cfg_path <- file.path(dating_dir, "calibrations_treePL_fulltips.cfg")
writeLines(treepl_cfg, calibrations_cfg_path)
cat("   Calibrations compiled successfully to:", calibrations_cfg_path, "\n")

# 2. Run Fast Automated `treePL` Dating Pipeline
cat("\nRunning automated `treePL` wrapper script over maximum-likelihood tree and bootstrap replicates...\n")

# Priming, cross-validation and dating follow the protocol of Maurin (2020, arXiv:2008.07054) and
# run through run_treePL_cv(). The priming parameters are the lowest of the repeats,
# cross-validation uses `randomcv`, and the smoothing grid runs from 1e+03 to 1e-14, one value per
# order of magnitude. The grid extends well below the 1e-06 to 1e-08 range reported by Maurin, so
# that a plateau can be distinguished from a search that stopped at the edge of the grid.
#
# The seed fixes every stochastic step of the dating. 79992967 is the seed of the reference run:
# keep it to reproduce the published chronogram, or change it for an independent run.
#
# num_bs = 100 dates 100 temporal bootstrap replicates; num_bs = NULL dates all of them.
automate_treePL(
  cfg_file = calibrations_cfg_path,
  ml_tree_file = paths$best_tree,
  bs_trees_file = paths$temporal_bs,
  numsites = num_sites,
  outgroup = rooting_outgroup,
  results_dir = file.path(dating_dir, "auto_results"),
  treePL_out = dating_dir,
  num_bs = 100,
  seed = 79992967,
  # Declared in the SETUP block. Reports the outcome, the timings and the output paths by email.
  notify = notify_email
)

# 3. View and evaluate dating metrics cleanly in R
cat("\n--- Chronological Dating Results Summary ---\n")
if (file.exists(file.path(dating_dir, "BestTree_treePL.tree"))) {
  ml_chronogram <- ape::read.tree(file.path(dating_dir, "BestTree_treePL.tree"))
  cat("   Best ML Chronogram:\n")
  cat("     - File path:", file.path(dating_dir, "BestTree_treePL.tree"), "\n")
  cat("     - Root age:", max(ape::node.depth.edgelength(ml_chronogram)), "Mya\n")
  cat("     - Number of tips:", length(ml_chronogram$tip.label), "\n")

  # A node whose age equals one of its own bounds was not estimated: penalized likelihood returned
  # the constraint, and the chronogram gives no sign of it. ACP_root is fixed (min = max), so it
  # sits on its bound by construction; the other three calibrations are expected to be estimated.
  # The bootstrap chronograms are passed because a single point estimate can land just inside
  # a bound while the underlying age is unidentifiable, and only the replicates separate the two.
  bs_chronograms <- file.path(dating_dir, "bsTree_treePL.tree")
  adherence <- report_bound_adherence(
    ml_chronogram, calibs, constraints,
    bootstraps = if (file.exists(bs_chronograms)) bs_chronograms else NULL
  )
  cat("\n   Calibrated node ages against their bounds:\n")
  print(adherence)
}
cat("--------------------------------------------\n")

# 4. The cross-validation curve
# run_treePL_cv() reported the shape of this curve when it ran (interior minimum, edge or plateau);
# the table shows the curve itself. Values from 10 upward are of the order of 1e41: treePL fails
# numerically at high smoothing, and those rows carry no information.
cv_file <- file.path(dating_dir, "auto_results", "ML_tree", "cv_ML_tree")
if (file.exists(cv_file)) {
  lines_cv <- grep("chisq", readLines(cv_file, warn = FALSE), value = TRUE)
  cv_tab <- data.frame(
    smoothing = as.numeric(gsub(".*\\(([^)]*)\\).*", "\\1", lines_cv)),
    chisq = as.numeric(sub(".*\\)\\s*", "", lines_cv))
  )
  print(cv_tab[order(cv_tab$smoothing), ], row.names = FALSE)
}

# 5. Smoothing sensitivity
# The same tree dated across the plateau and beyond it. Read the range of each node across 1e-14
# to 1e-06, the plateau of the reference run. At 1 and 10 the calibrated nodes sit on their upper
# bounds; that end of the table lies outside the range supported by cross-validation.

sensitivity_cfg <- file.path(dating_dir, "auto_results", "ML_tree", "configure_smooth_ML_tree")
if (file.exists(sensitivity_cfg)) {
  cat("\nRunning smoothing sensitivity analysis across 1e-14 to 100...\n")
  sensitivity <- report_smoothing_sensitivity(
    cfg_file = sensitivity_cfg,
    smoothing_values = c(1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-5, 1e-4, 1e-2, 1, 10, 100),
    treepl_bin = treepl_path
  )
  cat("\n--- Calibrated Node Age Sensitivity Table ---\n")
  print(sensitivity)
  cat("--------------------------------------------\n")
}
```

## Conclusion

Stage 2 ends with the maximum-likelihood tree, its bootstrap support,
and the dated chronogram with confidence intervals for node ages.

Tutorial 3 adds IUCN Red List assessments and produces the tree,
chronogram and conservation figures.

[Continue to Tutorial 3: Visualization and Metadata
Integration](https://beeamerino.github.io/PhyloCactus/articles/tutorial-3-cactus-phylogeny-visualization.html)

## References

- Akaike, H. 1974. A new look at the statistical model identification.
  *IEEE Transactions on Automatic Control*, 19(6), 716-723.
  <https://doi.org/10.1109/TAC.1974.1100705>
- Arakaki *et al*. 2011. Contemporaneous and recent radiations of the
  world’s major succulent plant lineages. *Proceedings of the National
  Academy of Sciences of the United States of America*, 108(20),
  8379–8384. <https://doi.org/10.1073/pnas.1100628108>
- Darriba *et al*. 2020. ModelTest-NG: a new and scalable tool for the
  selection of DNA and protein evolutionary models. *Molecular Biology
  and Evolution*, 37(1), 291-294.
  <https://doi.org/10.1093/molbev/msz189>
- Felsenstein, J. 1981. Evolutionary trees from DNA sequences: A maximum
  likelihood approach. *Journal of Molecular Evolution*, 17(6), 368–376.
  <https://doi.org/10.1007/BF01734359>
- Flouri *et al*. 2014. The Phylogenetic Likelihood Library. *Systematic
  Biology*, 64(2), 356-362. <https://doi.org/10.1093/sysbio/syu084>
- Helfrich *et al*. 2018. TreeAnnotator: versatile visual annotation of
  hierarchical text relations. *Proceedings of the Eleventh
  International Conference on Language Resources and Evaluation*.
  <https://lrec.elra.info/lrec2018-main-308>
- Hernández-Hernández *et al*. 2014. Beyond aridification: Multiple
  explanations for the elevated diversification of cacti in the New
  World Succulent Biome. *New Phytologist*, 202(4), 1382–1397.
  <https://doi.org/10.1111/nph.12752>
- Kozlov *et al*. 2019. RAxML-NG: A fast, scalable and user-friendly
  tool for maximum likelihood phylogenetic inference. *Bioinformatics*,
  35(21), 4453–4455. <https://doi.org/10.1093/bioinformatics/btz305>
- Lemoine *et al*. 2018. Renewing Felsenstein’s phylogenetic bootstrap
  in the era of big data. *Nature*, 556(7702), 452-456.
  <https://doi.org/10.1038/s41586-018-0043-0>
- Maurin, K. J. 2020. An empirical guide for producing a dated phylogeny
  with treePL in a maximum likelihood framework. *arXiv preprint
  arXiv:2008.07054*. <https://doi.org/10.48550/arXiv.2008.07054>
- Pattengale *et al*. 2010. How many bootstrap replicates are
  necessary?. *Journal of Computational Biology*, 17(3), 337–354.
  <https://doi.org/10.1089/cmb.2009.0179>
- Ramírez-Barahona *et al*. 2020. The delayed and geographically
  heterogeneous diversification of flowering plant families. *Nature
  Ecology and Evolution*, 4(9), 1232–1238.
  <https://doi.org/10.1038/s41559-020-1241-3>
- Robinson, D. F., & Foulds, L. R. 1981. Comparison of phylogenetic
  trees. *Mathematical Biosciences*, 53(1-2), 131-147.
  <https://doi.org/10.1016/0025-5564(81)90043-2>
- Sanderson, M. J. 2002. Estimating Absolute Rates of Molecular
  Evolution and Divergence Times: A Penalized Likelihood Approach.
  *Molecular Biology and Evolution*, 19(1), 101–109.
  <https://doi.org/10.1093/oxfordjournals.molbev.a003974>
- Smith, S. A., & O’Meara, B. C. 2012. TreePL: Divergence time
  estimation using penalized likelihood for large phylogenies.
  *Bioinformatics*, 28(20), 2689–2690.
  <https://doi.org/10.1093/bioinformatics/bts492>
