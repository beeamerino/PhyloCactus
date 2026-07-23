# Compute Robinson-Foulds Distances Across Maximum-Likelihood Trees

Computes pairwise Robinson-Foulds (RF) topological distances across tree
topologies generated during maximum-likelihood search in \`RAxML-NG\`.
Quantifying topological variance evaluates whether independent search
runs converged on identical tree topologies.

## Usage

``` r
calculate_rf_distances(
  raxml_bin_path,
  ml_trees_file,
  output_dir = dirname(ml_trees_file),
  prefix = "cactus_RF"
)
```

## Arguments

- raxml_bin_path:

  Character. System command or full path to executable \`RAxML-NG\`
  binary.

- ml_trees_file:

  Character. Path to input \`.raxml.mlTrees\` file containing multiple
  ML tree search replicates.

- output_dir:

  Character. Output directory for RF distance calculations. Defaults to
  \`dirname(ml_trees_file)\`.

- prefix:

  Character. Output file prefix. Defaults to \`"cactus_RF"\`.

## Value

Character path to the resulting RF distance output file
(\`.raxml.rfdist\`).
