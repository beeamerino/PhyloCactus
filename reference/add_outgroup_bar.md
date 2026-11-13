# Mark the Outgroup of a ggtree Plot with Bars

The terminals outside the ingroup are usually not a clade of the tree,
so `geom_cladelabel()` cannot mark them. This function draws a bar over
each run of consecutive outgroup terminals in the plotting order and
writes `label` beside it, at the position the bars of
[`add_clade_labels_by_level()`](https://beeamerino.github.io/PhyloCactus/reference/add_clade_labels_by_level.md)
take with `align = TRUE` and the same `offset` and `offset_text`.

## Usage

``` r
add_outgroup_bar(
  p,
  ingroup_tips,
  label = "Outgroup",
  fontsize = 3,
  barsize = 0.45,
  offset = 0.03,
  offset_text = 0.004,
  fontface = "plain"
)
```

## Arguments

- p:

  A `ggtree` plot object.

- ingroup_tips:

  Character vector. Tip labels of the ingroup; every other tip is
  outgroup.

- label:

  Character. Text written beside each bar. Defaults to `"Outgroup"`.

- fontsize:

  Numeric. Font size of the label. Defaults to `3`.

- barsize:

  Numeric. Line width of the bar. Defaults to `0.45`.

- offset:

  Numeric. Horizontal offset of the bar, as a fraction of the width of
  the tree. Defaults to `0.03`.

- offset_text:

  Numeric. Offset of the label from the bar, as a fraction of the width
  of the tree. Defaults to `0.004`.

- fontface:

  Character. Font face of the label. Defaults to `"plain"`.

## Value

The `ggtree` plot with a segment layer and a text layer, or `p`
unchanged when no tip is outgroup.
