# Generate HPC SLURM Batch Script for Parallel Bootstrapping

Generates an executable Bash script with SLURM scheduler directives to
parallelize non-parametric bootstrapping across HPC compute nodes
(optimized for NLHPC Leftraru Epu and general Slurm clusters). Uses
coarse-grained parallelization via Slurm Job Arrays
(`#SBATCH --array=1-N`), pre-parsing to compressed binary `.rba` format
to optimize disk I/O, dynamic per-task seed multiplication for
statistical independence, and a configurable thread-to-worker ratio.

## Usage

``` r
generate_bootstrap_script(
  alignment_file,
  partition_file,
  constraint_file,
  outgroup = NULL,
  bs_per_rep = 500,
  max_reps = 2,
  base_seed = NULL,
  threads = 40,
  workers = NULL,
  threads_per_worker = 4L,
  preparse = TRUE,
  output_dir = getwd(),
  script_name = "run_bs_chunks.sh",
  cluster_job_name = "cactus_bs",
  cluster_partition = "main",
  cluster_nodes = 1L,
  cluster_mem = "10G",
  cluster_time = "24:00:00",
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

  Character. Path to partition file.

- constraint_file:

  Character. Path to constraint scaffold tree file.

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

- bs_per_rep:

  Integer. Number of bootstrap trees generated per chunk replicate.
  Defaults to `500`.

- max_reps:

  Integer. Total number of parallel chunk replicates (tasks) to spawn in
  the SLURM array (`1-max_reps`). Defaults to `2`, giving a total array
  target of `bs_per_rep * max_reps = 1000` bootstrap trees. This
  provides a safety margin over the autoMRE convergence point observed
  in a production run of this exact dataset (1023 taxa, 11 partitions),
  which converged (`bs-cutoff = 0.03`, FBP) after 600 of 1000 collected
  trees, i.e. convergence is not guaranteed at a fixed replicate count
  for every dataset or taxon sampling scheme. Always confirm convergence
  with
  [`check_bs_convergence()`](https://beeamerino.github.io/PhyloCactus/reference/check_bs_convergence.md)
  on the collected trees (via
  [`collect_bootstraps()`](https://beeamerino.github.io/PhyloCactus/reference/collect_bootstraps.md))
  rather than assuming `bs_per_rep * max_reps` is sufficient; increase
  `max_reps` and re-run
  [`collect_bootstraps()`](https://beeamerino.github.io/PhyloCactus/reference/collect_bootstraps.md)
  if it is not.

- base_seed:

  Integer. Base random seed for dynamic seed calculation
  (`SEED=$(( SLURM_ARRAY_TASK_ID * base_seed ))`). Defaults to `NULL`
  (random).

- threads:

  Integer. Number of CPU cores requested per SLURM task
  (`--cpus-per-task`). Defaults to `40`.

- workers:

  Integer. Number of RAxML-NG worker processes. If `NULL` (default),
  calculated dynamically as
  `max(1L, as.integer(threads / threads_per_worker))`.

- threads_per_worker:

  Integer. Thread-to-worker ratio used to derive `workers` when
  `workers = NULL`. Defaults to `4L` (empirically validated for this
  workload on an AMD EPYC node; see Details). Ignored if `workers` is
  set explicitly.

- preparse:

  Logical. Generate compressed binary `.rba` format with
  `raxml-ng --parse` prior to array execution to minimize disk I/O
  contention? Defaults to `TRUE`.

- output_dir:

  Character. Destination directory for script and chunk logs. Defaults
  to [`getwd()`](https://rdrr.io/r/base/getwd.html).

- script_name:

  Character. Name of output Bash script file. Defaults to
  `"run_bs_chunks.sh"`.

- cluster_job_name:

  Character. SLURM job name identifier. Defaults to `"cactus_bs"`.

- cluster_partition:

  Character. SLURM partition name. Defaults to `"main"` (Leftraru Epu
  AMD EPYC 9754 partition).

- cluster_nodes:

  Integer. Number of compute nodes requested (`--nodes`). Defaults to
  `1L`.

- cluster_mem:

  Character. Memory allocation string for SLURM (`--mem`). Defaults to
  `"10G"`.

- cluster_time:

  Character. Time limit allocation string for SLURM (`--time`). Defaults
  to `"24:00:00"`.

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

## Details

The default `threads_per_worker = 4L` is configured to optimize
throughput on multi-core compute nodes. In a production run of this
workload (1023 taxa, 11 partitions), using `--threads 40 --workers 10`
(4 threads per worker) on an AMD EPYC node provided efficient per-worker
memory and CPU allocation. This ratio is dataset- and hardware-dependent
(it trades per-worker single-tree search speed against the number of
trees searched in parallel); if you migrate to different node hardware
or a markedly different supermatrix size, re-validate it empirically
with a short trial run before committing a full job array to it.
