# Generate High-Performance Computing (HPC/SLURM) Batch Script for Temporal Bootstrapping

Writes a self-contained SLURM batch script for running high-throughput
temporal bootstrap replicate inference using MPI-enabled `RAxML-NG`
(`raxml-ng-mpi`) on a compute cluster. Temporal bootstraps constrain the
maximum-likelihood tree topology and re-estimate branch lengths across
bootstrap matrices, generating empirical distributions of branch lengths
for downstream divergence time estimation (`treePL`) without introducing
topological discordance.

## Usage

``` r
generate_temporal_bootstrap_script(
  alignment_file,
  partition_file,
  best_tree_file,
  outgroup = NULL,
  bs_trees = 500,
  base_seed = NULL,
  threads = 40,
  workers = NULL,
  threads_per_worker = 4L,
  preparse = TRUE,
  output_dir = getwd(),
  script_name = "run_temporal_bs.sh",
  cluster_job_name = "cactus_temp_bs",
  cluster_partition = "main",
  cluster_nodes = 1L,
  cluster_mem = "16G",
  cluster_time = "04:00:00",
  cluster_queue = NULL,
  cluster_mail_user = Sys.getenv("MY_EMAIL", ""),
  load_module = c("gcc/14.2.0-nlhpc", "openmpi/5.0.3-o", "raxml-ng/1.1.0-mpi-zen4-n"),
  raxml_exec = "raxml-ng-mpi"
)
```

## Arguments

- alignment_file:

  Character. Path to input PHYLIP alignment file.

- partition_file:

  Character. Path to partition file specifying substitution models.

- best_tree_file:

  Character. Path to best scoring maximum-likelihood tree topology file
  (used as constraint).

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

- bs_trees:

  Integer. Total number of temporal bootstrap trees to generate.
  Defaults to `500`.

- base_seed:

  Integer. Random seed for reproducible initialization. Defaults to
  `NULL` (random).

- threads:

  Integer. Number of CPU cores requested per SLURM task
  (`--cpus-per-task`). Defaults to `40`.

- workers:

  Integer or NULL. Number of RAxML-NG worker processes. If `NULL`
  (default), calculated dynamically as
  `max(1L, as.integer(threads / threads_per_worker))`.

- threads_per_worker:

  Integer. Thread-to-worker ratio used to derive `workers` when
  `workers = NULL`. Defaults to `4L`.

- preparse:

  Logical. Generate compressed binary `.rba` format with
  `raxml-ng --parse` prior to execution to minimize disk I/O contention?
  Defaults to `TRUE`.

- output_dir:

  Character. Destination directory for script and outputs. Defaults to
  [`getwd()`](https://rdrr.io/r/base/getwd.html).

- script_name:

  Character. Name of output Bash script file. Defaults to
  `"run_temporal_bs.sh"`.

- cluster_job_name:

  Character. SLURM job name identifier. Defaults to `"cactus_temp_bs"`.

- cluster_partition:

  Character. SLURM partition name. Defaults to `"main"`.

- cluster_nodes:

  Integer. Number of compute nodes requested (`--nodes`). Defaults to
  `1L`.

- cluster_mem:

  Character. Memory allocation string for SLURM (`--mem`). Defaults to
  `"16G"`.

- cluster_time:

  Character. Time limit allocation string for SLURM (`--time`). Defaults
  to `"04:00:00"`.

- cluster_queue:

  Character. Optional SLURM queue / QoS name. Defaults to `NULL`.

- cluster_mail_user:

  Character. Notification recipient email address for SLURM
  (`--mail-user`). Defaults to `Sys.getenv("MY_EMAIL", "")`.

- load_module:

  Character vector or string. Environment module(s) to load prior to
  execution. Defaults to
  `c("gcc/14.2.0-nlhpc", "openmpi/5.0.3-o", "raxml-ng/1.1.0-mpi-zen4-n")`.

- raxml_exec:

  Character. Executable `RAxML-NG` binary command. Defaults to
  `"raxml-ng-mpi"`.

## Value

Character path to the generated SLURM batch script file.
