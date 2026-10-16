# Theme of the Figures of PhyloCactus

The one theme of the package figures, the theme written for the figures
of Tutorial 3 and moved into the package as it was:
[`theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
with a light grey panel and a white grid, italic axis titles, bold
legend and strip titles. Its colours come from
`phylocactus_colours("theme")`. The title of a figure is its caption,
not a title inside the plot. Trees use
[`ggtree::theme_tree()`](https://rdrr.io/pkg/ggtree/man/theme_tree.html)
instead.

## Usage

``` r
theme_phylocactus(base_size = 12, base_family = "")
```

## Arguments

- base_size:

  Numeric. Base font size in points. Defaults to `12`.

- base_family:

  Character. Base font family. Defaults to `""`.

## Value

A `ggplot2` theme.

## See also

[`phylocactus_colours()`](https://beeamerino.github.io/PhyloCactus/reference/phylocactus_colours.md).

## Examples

``` r
if (FALSE) { # \dontrun{
ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point() + theme_phylocactus()
} # }
```
