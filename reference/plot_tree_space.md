# Tree Space of a Comparison of Phylogenies

Draws the tree-space figure of
[`validate_phylogenies()`](https://beeamerino.github.io/PhyloCactus/reference/validate_phylogenies.md):
one panel per taxon set (A1 to A4) with the classical multidimensional
scaling (MDS) of the pairwise distances among the trees, and one legend
below.
[`validate_phylogenies()`](https://beeamerino.github.io/PhyloCactus/reference/validate_phylogenies.md)
calls it on the coordinates it computes; Tutorial 4 calls it on the
coordinates of the reference run distributed in `inst/extdata`
(`phylocactus_table_tree_space.csv`), so the figure is drawn without the
reference trees or the checklist.

## Usage

``` r
plot_tree_space(coordinates, metric = c("RF", "DPI"), tree_labels = NULL)
```

## Arguments

- coordinates:

  Data frame, or path to a CSV file, of tree-space coordinates as
  [`validate_phylogenies()`](https://beeamerino.github.io/PhyloCactus/reference/validate_phylogenies.md)
  writes them to `SUPP_tree_space_coordinates.csv`: columns
  `analysis_id`, `analysis_label`, `tree`, `Dimension_1`, `Dimension_2`
  and `metric`.

- metric:

  Character. Distance whose projection is drawn: `"RF"` (normalized
  Robinson-Foulds, the default) or `"DPI"` (Different Phylogenetic
  Information).

- tree_labels:

  Named character vector or `NULL`. Legend labels, as in
  [`validate_phylogenies()`](https://beeamerino.github.io/PhyloCactus/reference/validate_phylogenies.md).

## Value

A `ggplot` object (assembled with `cowplot`).

## Examples

``` r
if (FALSE) { # \dontrun{
coords <- system.file("extdata", "phylocactus_table_tree_space.csv", package = "PhyloCactus")
plot_tree_space(coords, metric = "RF")
} # }
```
