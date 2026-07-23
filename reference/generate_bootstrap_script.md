# Generate HPC SLURM Batch Script for Parallel Bootstrapping

Generates an executable Bash script with SLURM scheduler directives to
parallelize non-parametric bootstrapping across HPC compute nodes.
Chunking bootstrap searches into parallel sub-jobs accelerates support
estimation for large supermatrices.

## Usage

``` r
generate_bootstrap_script(
  alignment_file,
  partition_file,
  constraint_file,
  outgroup = NULL,
  bs_per_rep = 500,
  max_reps = 4,
  base_seed = 111,
  seed_step = 1000,
  threads = 120,
  workers = 20,
  output_dir = getwd(),
  script_name = "run_bs_chunks.sh",
  cluster_job_name = "C-raxml.bs",
  cluster_mem = "32G",
  cluster_time = "7-00:00",
  cluster_partition = "general",
  cluster_queue = "public",
  load_module = "raxml-ng-1.1.0-gcc-11.2.0",
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

  Character. Optional outgroup taxon binomial. Defaults to \`NULL\`.

- bs_per_rep:

  Integer. Number of bootstrap trees generated per chunk replicate.
  Defaults to \`500\`.

- max_reps:

  Integer. Total number of parallel chunk replicates to spawn. Defaults
  to \`4\`.

- base_seed:

  Integer. Base random seed. Defaults to \`111\`.

- seed_step:

  Integer. Seed increment value between consecutive chunks. Defaults to
  \`1000\`.

- threads:

  Integer. Number of CPU cores requested per SLURM task. Defaults to
  \`120\`.

- workers:

  Integer. Number of RAxML-NG worker processes. Defaults to \`20\`.

- output_dir:

  Character. Destination directory for script and chunk logs. Defaults
  to \`getwd()\`.

- script_name:

  Character. Name of output Bash script file. Defaults to
  \`"run_bs_chunks.sh"\`.

- cluster_job_name:

  Character. SLURM job name identifier. Defaults to \`"C-raxml.bs"\`.

- cluster_mem:

  Character. Memory allocation string for SLURM. Defaults to \`"32G"\`.

- cluster_time:

  Character. Time limit allocation string for SLURM. Defaults to
  \`"7-00:00"\`.

- cluster_partition:

  Character. SLURM partition name. Defaults to \`"general"\`.

- cluster_queue:

  Character. SLURM queue name. Defaults to \`"public"\`.

- load_module:

  Character. Environment module command to load prior to execution.
  Defaults to \`"raxml-ng-1.1.0-gcc-11.2.0"\`.

- raxml_exec:

  Character. Executable \`RAxML-NG\` binary name. Defaults to
  \`"raxml-ng-mpi"\`.

## Value

Character path to the generated SLURM batch script file.
