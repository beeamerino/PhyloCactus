# Prepare Tip Annotations Matched to Tree Leaves

Matches and standardizes taxonomic metadata from a classification table
against tree tip labels.

## Usage

``` r
prepare_tip_annotation(tree, constraints_tbl)
```

## Arguments

- tree:

  Object of class \`phylo\` representing a phylogenetic tree.

- constraints_tbl:

  Data frame containing taxonomic constraint annotations.

## Value

Data frame of tip annotations filtered to matching tree tips.
