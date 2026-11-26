# Prepare Tip Annotations Matched to Tree Leaves

Matches and standardizes taxonomic metadata from a classification table
against tree tip labels.

## Usage

``` r
prepare_tip_annotation(tree, constraints_tbl)
```

## Arguments

- tree:

  Object of class `phylo` representing a phylogenetic tree.

- constraints_tbl:

  Data frame containing taxonomic constraint annotations.

## Value

Data frame of tip annotations filtered to matching tree tips.

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
table(tips$Annot_level_2)
#> 
#>              Opuntioideae 
#>           32          155 
```
