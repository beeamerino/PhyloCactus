# Panel Tags of the Composite Figures of PhyloCactus

Every composite figure of the package and of the tutorial scripts takes
its panel tags from this function, so that their style is changed in one
place when a journal asks for it: `A`, `B` by default;
`options(phylocactus.tag_suffix = ")")` gives `A)`, `B)`, and
`options(phylocactus.tag_levels = "a")` gives lower-case letters.

## Usage

``` r
phylocactus_panel_tags(...)
```

## Arguments

- ...:

  Further arguments of
  [`patchwork::plot_annotation()`](https://patchwork.data-imaginist.com/reference/plot_annotation.html)
  (a caption, a theme).

## Value

A `patchwork` annotation to add to a composite figure.

## See also

[`theme_phylocactus()`](https://beeamerino.github.io/PhyloCactus/reference/theme_phylocactus.md).

## Examples

``` r
if (FALSE) { # \dontrun{
patchwork::wrap_plots(p1, p2) + phylocactus_panel_tags()
options(phylocactus.tag_suffix = ")")
} # }
```
