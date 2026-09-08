# Generate Non-Parametric Bootstrap Trees Locally

Performs non-parametric bootstrap resampling over supermatrix site
columns to infer a distribution of bootstrap tree topologies
(`RAxML-NG`). Evaluates topological variation under non-parametric
resampling to quantify node support via Felsenstein Bootstrap
Proportions (FBP) or Transfer Bootstrap Expectation (TBE). Optimizes
multi-threading for local workstations (e.g., Apple Silicon) by
allocating threads to Performance cores (P-cores) and configuring
parallel workers according to a configurable thread-to-worker ratio,
avoiding thread contention with Efficiency cores (E-cores).

## Usage

``` r
run_local_bootstraps(
  raxml_bin_path,
  aln_file,
  part_file = NULL,
  constraint_file,
  bs_trees = 500,
  outgroup = NULL,
  seed = NULL,
  threads = 8,
  workers = NULL,
  threads_per_worker = 4L,
  output_dir = dirname(aln_file),
  prefix = "cactus_bs"
)
```

## Arguments

- raxml_bin_path:

  Character. System command or full path to executable `RAxML-NG`
  binary.

- aln_file:

  Character. Path to input PHYLIP supermatrix alignment file or
  pre-parsed `.rba` binary file.

- part_file:

  Character. Path to partition file specifying substitution models.
  Optional if `aln_file` is an `.rba` binary file. Defaults to `NULL`.

- constraint_file:

  Character. Path to Newick topological constraint scaffold file.

- bs_trees:

  Integer. Total number of non-parametric bootstrap trees to generate.
  Defaults to `500`.

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

- seed:

  Integer. Random seed for reproducible bootstrap initialization.
  Defaults to `NULL` (random).

- threads:

  Integer. Number of CPU threads. Defaults to `8` (or detected P-cores).

- workers:

  Integer. Parallel worker process count. If `NULL` (default),
  automatically calculated as
  `max(1L, as.integer(threads / threads_per_worker))`.

- threads_per_worker:

  Integer. Thread-to-worker ratio used to derive `workers` when
  `workers = NULL`. Defaults to `4L` (unvalidated heuristic for local
  bootstrap workloads; see Details). Ignored if `workers` is set
  explicitly.

- output_dir:

  Character. Output directory for generated bootstrap trees. Defaults to
  `dirname(aln_file)`.

- prefix:

  Character. Output file prefix. Defaults to `"cactus_bs"`.

## Value

Character path to the output bootstrap trees file (`.raxml.bootstraps`).

## Details

Unlike
[`generate_bootstrap_script()`](https://beeamerino.github.io/PhyloCactus/reference/generate_bootstrap_script.md)'s
HPC-oriented default (empirically calibrated on a 128-core AMD EPYC
node), the `threads_per_worker = 4L` default here has not been
empirically validated for bootstrap workloads on Apple Silicon
specifically (only single-threaded, non-coarse-grained ML tree search
timings are currently available for that hardware). Treat this default
as a reasonable starting heuristic and re-validate empirically (e.g.,
time a short run with `bs_trees` set low, at a couple of
`threads_per_worker` values) before committing to a long local run.

## References

Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A.
(2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum
likelihood phylogenetic inference. *Bioinformatics*, 35(21), 4453-4455.
[doi:10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)
