# Add Geological Time Scale Axis to a Chronogram ggtree Plot

Draws a standardized geological timeline axis (Ma before present)
beneath an ultrametric chronogram `ggtree` plot.

## Usage

``` r
add_chronogram_axis(
  p,
  tree,
  by = 5,
  digits = 0L,
  axis_offset_frac = 0.052,
  tick_height_frac = 0.01,
  label_offset_frac = 0.028,
  title_margin_top = 28,
  bar_size = 12,
  segment_size = 2.5
)
```

## Arguments

- p:

  A `ggtree` plot object.

- tree:

  Object of class `phylo`.

- by:

  Numeric. Interval in Ma between consecutive tick marks. Defaults to
  `5`.

- digits:

  Integer. Decimal precision for age labels. Defaults to `0L`.

- axis_offset_frac:

  Numeric. Vertical offset fraction for axis line. Defaults to `0.052`.

- tick_height_frac:

  Numeric. Vertical height fraction for tick marks. Defaults to `0.010`.

- label_offset_frac:

  Numeric. Vertical offset fraction for age text labels. Defaults to
  `0.028`.

- title_margin_top:

  Numeric. Top margin for axis title. Defaults to `28`.

- bar_size:

  Numeric. Axis title font size. Defaults to `12`.

- segment_size:

  Numeric. Text font size for age labels. Defaults to `2.5`.

## Value

Updated `ggtree` plot object containing the geological age axis layer.

## Examples

``` r
# A pruned copy of the distributed chronogram and the constraint table
ext <- function(f) system.file("extdata", f, package = "PhyloCactus")
tree <- ape::read.tree(ext("phylocactus_chronogram_ml.tree"))
genera <- c("Talinum", "Portulaca", "Leuenbergeria", "Pereskia", "Opuntia",
            "Cylindropuntia", "Tephrocactus", "Maihuenia")
tree <- ape::keep.tip(tree, tree$tip.label[sub("_.*", "", tree$tip.label) %in% genera])
constraints_tbl <- utils::read.csv(ext("cactus_constraints.csv"))
p <- add_chronogram_axis(ggtree::ggtree(tree), tree, by = 10)
```
