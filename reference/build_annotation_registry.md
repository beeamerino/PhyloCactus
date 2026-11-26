# Build Taxonomic Annotation Registry across Tree Levels

Builds a registry mapping internal nodes to multi-level taxonomic
annotations (subfamilies, tribes, subtribes, genera).

## Usage

``` r
build_annotation_registry(
  tree,
  constraints_tbl,
  main_level4_values = NULL,
  supp_level4_min_tips = 10L,
  require_monophyly = TRUE
)
```

## Arguments

- tree:

  Object of class `phylo`.

- constraints_tbl:

  Data frame containing taxonomic classifications.

- main_level4_values:

  Character vector of main level 4 (genus) annotations. Defaults to
  `NULL`.

- supp_level4_min_tips:

  Integer. Minimum tip threshold for secondary level 4 annotations.
  Defaults to `10L`.

- require_monophyly:

  Logical. Enforce monophyly constraint? Defaults to `TRUE`.

## Value

A named list of annotation node data frames organized by taxonomic
level.

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
registry$level_4_supp
#> # A tibble: 7 × 5
#>   label             node n_tips_target n_tips_sub monophyletic
#>   <chr>            <int>         <int>      <int> <lgl>       
#> 1 Cylindropuntieae   199            34         34 TRUE        
#> 2 Leuenbergeria      359             8          8 TRUE        
#> 3 Opuntieae          233           109        109 TRUE        
#> 4 Pereskia           353             7          7 TRUE        
#> 5 Portulaca          366             9          9 TRUE        
#> 6 Talinum            189             6          6 TRUE        
#> 7 Tephrocacteae      341            12         12 TRUE        
```
