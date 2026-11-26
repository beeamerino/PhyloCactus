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

  Object of class `phylo`.

- tip_tbl:

  Data frame containing tip annotations.

- group_col:

  Character. Column name in `tip_tbl` specifying group assignments.

- group_values:

  Character vector. Optional subset of group values to evaluate.
  Defaults to `NULL`.

- min_tips:

  Integer. Minimum required tips per group. Defaults to `2L`.

## Value

A `tibble` containing group names, MRCA node numbers, tip counts, and
monophyly logical indicators.

## Examples

``` r
# A pruned copy of the distributed chronogram and the constraint table
ext <- function(f) system.file("extdata", f, package = "PhyloCactus")
tree <- ape::read.tree(ext("phylocactus_chronogram_ml.tree"))
genera <- c("Talinum", "Portulaca", "Leuenbergeria", "Pereskia", "Opuntia",
            "Cylindropuntia", "Tephrocactus", "Maihuenia")
tree <- ape::keep.tip(tree, tree$tip.label[sub("_.*", "", tree$tip.label) %in% genera])
constraints_tbl <- utils::read.csv(ext("cactus_constraints.csv"))
tips <- prepare_tip_annotation(tree, constraints_tbl)
compute_group_nodes(tree, tips, "Annot_level_4")
#> # A tibble: 8 × 5
#>   label             node n_tips_target n_tips_sub monophyletic
#>   <chr>            <int>         <int>      <int> <lgl>       
#> 1 Cylindropuntieae   199            34         34 TRUE        
#> 2 Leuenbergeria      359             8          8 TRUE        
#> 3 Maihuenia          352             2          2 TRUE        
#> 4 Opuntieae          233           109        109 TRUE        
#> 5 Pereskia           353             7          7 TRUE        
#> 6 Portulaca          366             9          9 TRUE        
#> 7 Talinum            189             6          6 TRUE        
#> 8 Tephrocacteae      341            12         12 TRUE        
```
