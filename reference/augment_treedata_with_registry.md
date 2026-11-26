# Augment treedata Object with Registry Metadata

Joins multi-level taxonomic clade annotations from a registry onto a
`treedata` object.

## Usage

``` r
augment_treedata_with_registry(treedata_obj, registry)
```

## Arguments

- treedata_obj:

  An object of class `treedata` or `phylo`.

- registry:

  A named list of annotation registry data frames.

## Value

An augmented `treedata` object.

## Examples

``` r
# A pruned copy of the distributed chronogram and the constraint table
ext <- function(f) system.file("extdata", f, package = "PhyloCactus")
tree <- ape::read.tree(ext("phylocactus_chronogram_ml.tree"))
genera <- c("Talinum", "Portulaca", "Leuenbergeria", "Pereskia", "Opuntia",
            "Cylindropuntia", "Tephrocactus", "Maihuenia")
tree <- ape::keep.tip(tree, tree$tip.label[sub("_.*", "", tree$tip.label) %in% genera])
constraints_tbl <- utils::read.csv(ext("cactus_constraints.csv"))
registry <- build_annotation_registry(tree, constraints_tbl, supp_level4_min_tips = 3L)
annotated <- augment_treedata_with_registry(tree, registry)
head(tidytree::as_tibble(annotated))
#> ! # Invaild edge matrix for <phylo>. A <tbl_df> is returned.
#> # A tibble: 6 × 8
#>   parent  node branch.length label     level_1 level_2 level_4_main level_4_supp
#>    <int> <int>         <dbl> <chr>     <chr>   <chr>   <chr>        <chr>       
#> 1    190     1         20.6  Talinum_… NA      NA      NA           NA          
#> 2    191     2         18.7  Talinum_… NA      NA      NA           NA          
#> 3    192     3          5.63 Talinum_… NA      NA      NA           NA          
#> 4    192     4          5.63 Talinum_… NA      NA      NA           NA          
#> 5    193     5         15.8  Talinum_… NA      NA      NA           NA          
#> 6    193     6         15.8  Talinum_… NA      NA      NA           NA          
```
