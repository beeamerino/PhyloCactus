# Compute Most Recent Common Ancestor (MRCA) Nodes for Taxon Groups

Identifies internal MRCA node numbers and evaluates monophyly for
specified taxonomic groups within a phylogenetic tree.

## Usage

``` r
compute_group_nodes(
  tree,
  tip_tbl,
  group_col,
  group_values = NULL,
  min_tips = 2L
)
```

## Arguments

- tree:

  Object of class \`phylo\`.

- tip_tbl:

  Data frame containing tip annotations.

- group_col:

  Character. Column name in \`tip_tbl\` specifying group assignments.

- group_values:

  Character vector. Optional subset of group values to evaluate.
  Defaults to \`NULL\`.

- min_tips:

  Integer. Minimum required tips per group. Defaults to \`2L\`.

## Value

A \`tibble\` containing group names, MRCA node numbers, tip counts, and
monophyly logical indicators.
