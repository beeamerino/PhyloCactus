# Apply Multiple Taxonomic Clade Label Layers onto ggtree Plot

Sequentially overlays hierarchical taxonomic clade annotations (e.g.,
subfamilies, tribes, genera) onto a `ggtree` plot.

## Usage

``` r
apply_clade_label_layers(p, registry, layer_specs)
```

## Arguments

- p:

  A `ggtree` plot object.

- registry:

  Named list of node annotation data frames.

- layer_specs:

  Data frame defining visual style specifications per annotation layer.
  An optional `hjust` column sets the justification of each layer;
  without it, rotated labels are centred on their bars (see
  [`add_clade_labels_by_level()`](https://beeamerino.github.io/PhyloCactus/reference/add_clade_labels_by_level.md)).

## Value

Updated `ggtree` plot object.

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
specs <- data.frame(reg_name = c("level_4_supp", "level_2"),
                    fontsize = c(2.5, 3), barsize = c(0.3, 0.5),
                    offset = c(0.01, 0.15), offset_text = c(0.004, 0.006),
                    fontface = c("plain", "bold"), sort_desc = TRUE,
                    angle = c(0, 270), align = TRUE)
p <- apply_clade_label_layers(ggtree::ggtree(tree), registry, specs)
#> ℹ invalid tbl_tree object. Missing column: parent,label.
#> ℹ invalid tbl_tree object. Missing column: parent,label.
#> ℹ invalid tbl_tree object. Missing column: parent,label.
#> ℹ invalid tbl_tree object. Missing column: parent,label.
```
