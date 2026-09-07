# Estimate Temporal Bootstrap Replicates Constrained to Best ML Topology

Generates non-parametric bootstrap tree replicates where branch lengths
are re-estimated while holding the focal maximum-likelihood topology
constrained. Temporal bootstraps propagate branch length uncertainty
into downstream penalized likelihood divergence time estimation
(`treePL`), providing empirical confidence intervals for node ages
without introducing topological variance.

## Usage

``` r
calculate_temporal_bootstraps(
  raxml_bin_path,
  aln_file,
  part_file,
  best_tree_file,
  bs_trees = 500,
  outgroup = NULL,
  seed = NULL,
  threads = 4,
  workers = NULL,
  blopt = "nr_safe",
  output_dir = dirname(aln_file),
  prefix = "cactus_temporal"
)
```

## Arguments

- raxml_bin_path:

  Character. System command or full path to executable `RAxML-NG`
  binary.

- aln_file:

  Character. Path to input PHYLIP supermatrix alignment file.

- part_file:

  Character. Path to partition file specifying substitution models.

- best_tree_file:

  Character. Path to reference maximum-likelihood tree topology file
  (used as constraint).

- bs_trees:

  Integer. Total number of temporal bootstrap trees to generate.
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

  Integer. Random seed for reproducible temporal bootstrap
  initialization. Defaults to `NULL` (random).

- threads:

  Integer. Number of CPU threads. Defaults to `4`.

- workers:

  Integer or NULL. Number of parallel tree search workers. If `NULL`,
  auto-configured. Defaults to `NULL`.

- blopt:

  Character. Branch length optimization method (`"nr_safe"` or
  `"nr_fast"`). Defaults to `"nr_safe"` for robust numerical convergence
  across partitioned alignments.

- output_dir:

  Character. Directory path to save output temporal bootstrap trees.
  Defaults to `dirname(aln_file)`.

- prefix:

  Character. Output file prefix. Defaults to `"cactus_temporal"`.

## Value

Character path to the resulting temporal bootstrap trees file
(`.raxml.bootstraps`).
