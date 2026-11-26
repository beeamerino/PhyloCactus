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

  Object of class `phylo`.

- constraints_tbl:

  Data frame containing taxonomic classification.

- group_col:

  Character. Column name defining group annotations.

- group_values:

  Character vector. Optional subset of group values. Defaults to `NULL`.

- min_tips:

  Integer. Minimum tip count. Defaults to `2L`.

- require_monophyly:

  Logical. Exclude non-monophyletic groups? Defaults to `TRUE`.

## Value

A `tibble` of node annotations.

## Examples

``` r
# A pruned copy of the distributed chronogram and the constraint table
ext <- function(f) system.file("extdata", f, package = "PhyloCactus")
tree <- ape::read.tree(ext("phylocactus_chronogram_ml.tree"))
genera <- c("Talinum", "Portulaca", "Leuenbergeria", "Pereskia", "Opuntia",
            "Cylindropuntia", "Tephrocactus", "Maihuenia")
tree <- ape::keep.tip(tree, tree$tip.label[sub("_.*", "", tree$tip.label) %in% genera])
constraints_tbl <- utils::read.csv(ext("cactus_constraints.csv"))
get_annotation_nodes(tree, constraints_tbl, "Annot_level_2")
#> # A tibble: 1 × 5
#>   label         node n_tips_target n_tips_sub monophyletic
#>   <chr>        <int>         <int>      <int> <lgl>       
#> 1 Opuntioideae   198           155        155 TRUE        
```
