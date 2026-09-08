# Generate a SLURM Batch Script for the Constrained Maximum-Likelihood Search

Writes a single SLURM batch script running the same constrained
`RAxML-NG` search as
[`calculate_ml_tree()`](https://beeamerino.github.io/PhyloCactus/reference/calculate_ml_tree.md),
sized for a compute node instead of a workstation. The search is one job
rather than a job array: independent starting-tree searches are
distributed across `--workers` inside the job, and RAxML-NG writes one
`.raxml.bestTree` directly, so no collection step is needed to compare
log-likelihoods across tasks.

## Usage

``` r
generate_ml_search_script(
  alignment_file,
  partition_file,
  constraint_file,
  outgroup = NULL,
  n_init_trees = "rand{25},pars{25}",
  seed = NULL,
  threads = 75,
  workers = NULL,
  min_threads_per_worker = 3L,
  preparse = TRUE,
  output_dir = getwd(),
  script_name = "run_ml_search.sh",
  prefix = "cactus_search",
  cluster_job_name = "cactus_ml",
  cluster_partition = "main",
  cluster_nodes = 1L,
  cluster_mem = "16G",
  cluster_time = "02:00:00",
  cluster_queue = NULL,
  cluster_mail_user = Sys.getenv("MY_EMAIL", ""),
  load_module = c("gcc/14.2.0-nlhpc", "openmpi/5.0.3-o", "raxml-ng/1.1.0-mpi-zen4-n"),
  raxml_exec = "raxml-ng-mpi"
)
```

## Arguments

- alignment_file:

  Character. Path to input PHYLIP supermatrix alignment file.

- partition_file:

  Character. Path to partition file specifying substitution models per
  partition.

- constraint_file:

  Character. Path to Newick topological constraint scaffold file.

- outgroup:

  Character vector of terminals passed to `RAxML-NG --outgroup`, or
  `NULL` (default); multiple terminals are joined with commas.
  `RAxML-NG` writes an unrooted topology with these terminals placed
  first, so this argument orders the output rather than rooting the
  tree: the root is imposed downstream by
  [`automate_treePL()`](https://beeamerino.github.io/PhyloCactus/reference/automate_treePL.md)
  via `ape::root(..., resolve.root = TRUE)`. Declaring the same set at
  every stage keeps the output ordering consistent across the
  maximum-likelihood search, the bootstrap replicates and the temporal
  bootstraps. Derive it with
  [`resolve_rooting_outgroup()`](https://beeamerino.github.io/PhyloCactus/reference/resolve_rooting_outgroup.md)
  rather than naming a terminal by hand.

- n_init_trees:

  Character. Starting tree specification passed to `--tree`. Defaults to
  `"rand{25},pars{25}"`.

- seed:

  Integer. Random seed for reproducible tree search initialization.
  Defaults to `NULL` (random).

- threads:

  Integer. CPU cores requested per SLURM task (`--cpus-per-task`).
  Defaults to `75`.

- workers:

  Integer. Worker process count. If `NULL` (default), derived from
  `threads`, `n_init_trees` and `min_threads_per_worker` so that workers
  divide the starting tree count evenly.

- min_threads_per_worker:

  Integer. Lower bound on threads per worker when `workers` is derived.
  Defaults to `3L`.

- preparse:

  Logical. Pre-parse the alignment into compressed binary `.rba` format
  before the search? Defaults to `TRUE`.

- output_dir:

  Character. Destination directory for the script and its outputs.
  Defaults to [`getwd()`](https://rdrr.io/r/base/getwd.html).

- script_name:

  Character. Name of the generated Bash script. Defaults to
  `"run_ml_search.sh"`.

- prefix:

  Character. RAxML-NG output prefix inside the job. Defaults to
  `"cactus_search"`.

- cluster_job_name:

  Character. SLURM job name. Defaults to `"cactus_ml"`.

- cluster_partition:

  Character. SLURM partition. Defaults to `"main"`.

- cluster_nodes:

  Integer. Number of compute nodes requested (`--nodes`). Defaults to
  `1L`.

- cluster_mem:

  Character. Memory allocation (`--mem`). Defaults to `"16G"`.

- cluster_time:

  Character. Time limit (`--time`). Defaults to `"02:00:00"`.

- cluster_queue:

  Character. Optional SLURM queue / QoS. Defaults to `NULL`.

- cluster_mail_user:

  Character. Notification recipient email address for SLURM
  (`--mail-user`). Defaults to `Sys.getenv("MY_EMAIL", "")`.

- load_module:

  Character vector. Environment modules to load. Defaults to the NLHPC
  Leftraru `RAxML-NG` toolchain.

- raxml_exec:

  Character. `RAxML-NG` executable invoked inside the job. Defaults to
  `"raxml-ng-mpi"`.

## Value

Character path to the generated SLURM batch script, carrying the
resolved parallel plan as the `ml_plan` attribute.

## Details

Wall-clock is governed by the number of sequential rounds each worker
must run, `ceiling(n_trees / workers)`, not by the thread count alone.
Total core-work is fixed, so adding threads to a single worker shortens
each tree while adding workers shortens the number of rounds, and the
second lever is the one that matters once a search is already using an
efficient thread count per worker. Running 50 starting trees on one
worker executes 50 rounds; the same 50 trees across 25 workers execute
2.

Reference timing illustrating worker parallelization: on an Apple M2 Pro
(8 threads, `--workers 1`, `RAxML-NG` 1.2.2, SSE3 kernels), searching 50
starting trees sequentially required 47115 s (~938 s per tree). In
contrast, the production run on an HPC cluster node (AMD EPYC 9754, 75
threads, `--workers 25`) over the full supermatrix (1023 terminals,
12806 sites, 5954 patterns, 11 partitions) completed 50 starting trees
in 2393 s (~40 minutes), executing two parallel rounds.

`workers` is derived automatically and constrained to a divisor of the
starting tree count, so that no worker sits idle in the final round.

## References

Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A.
(2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum
likelihood phylogenetic inference. *Bioinformatics*, 35(21), 4453-4455.
[doi:10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

## See also

[`calculate_ml_tree()`](https://beeamerino.github.io/PhyloCactus/reference/calculate_ml_tree.md)
for the equivalent local run.

## Examples

``` r
if (FALSE) { # \dontrun{
generate_ml_search_script(
  alignment_file = "cactus.raxml.reduced.phy",
  partition_file = "cactus_modeltest.part.aicc",
  constraint_file = "cactus_constraints.tree",
  outgroup = resolve_rooting_outgroup(ape::read.tree("bestTree.tree")$tip.label),
  output_dir = "7_Phylogenetics"
)
} # }
```
