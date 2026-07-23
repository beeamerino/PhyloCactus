# Extract Monophyletic Internal Annotation Nodes

Retrieves internal node numbers for monophyletic taxonomic groups across
a phylogenetic tree.

## Usage

``` r
get_annotation_nodes(
  tree,
  constraints_tbl,
  group_col,
  group_values = NULL,
  min_tips = 2L,
  require_monophyly = TRUE
)
```

## Arguments

- tree:

  Object of class \`phylo\`.

- constraints_tbl:

  Data frame containing taxonomic classification.

- group_col:

  Character. Column name defining group annotations.

- group_values:

  Character vector. Optional subset of group values. Defaults to
  \`NULL\`.

- min_tips:

  Integer. Minimum tip count. Defaults to \`2L\`.

- require_monophyly:

  Logical. Exclude non-monophyletic groups? Defaults to \`TRUE\`.

## Value

A \`tibble\` of node annotations.
