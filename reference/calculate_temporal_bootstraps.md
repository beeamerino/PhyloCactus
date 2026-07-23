# Estimate Temporal Bootstrap Replicates Constrained to Best ML Topology

Generates non-parametric bootstrap tree replicates where branch lengths
are re-estimated while holding the focal maximum-likelihood topology
constrained. Temporal bootstraps propagate branch length uncertainty
into downstream penalized likelihood divergence time estimation
(\`treePL\`), providing empirical confidence intervals for node ages
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
  seed = 111,
  threads = 4,
  output_dir = dirname(aln_file),
  prefix = "cactus_temporal"
)
```

## Arguments

- raxml_bin_path:

  Character. System command or full path to executable \`RAxML-NG\`
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
  Defaults to \`500\`.

- outgroup:

  Character. Optional outgroup taxon binomial. Defaults to \`NULL\`.

- seed:

  Integer. Random seed for reproducible temporal bootstrap
  initialization. Defaults to \`111\`.

- threads:

  Integer. Number of CPU threads. Defaults to \`4\`.

- output_dir:

  Character. Directory path to save output temporal bootstrap trees.
  Defaults to \`dirname(aln_file)\`.

- prefix:

  Character. Output file prefix. Defaults to \`"cactus_temporal"\`.

## Value

Character path to the resulting temporal bootstrap trees file
(\`.raxml.bootstraps\`).
