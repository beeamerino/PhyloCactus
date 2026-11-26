# Compute Maximum Depth of a Dated Chronogram Tree

Calculates the maximum root-to-tip temporal distance of an ultrametric
tree in Millions of Years ago (Ma).

## Usage

``` r
get_tree_max_depth(tree)
```

## Arguments

- tree:

  Object of class `phylo`.

## Value

Numeric scalar representing maximum tree depth.

## Examples

``` r
ext <- function(f) system.file("extdata", f, package = "PhyloCactus")
tree <- ape::read.tree(ext("phylocactus_chronogram_ml.tree"))
get_tree_max_depth(tree)   # the root age in Ma
#> [1] 58.10159
```
