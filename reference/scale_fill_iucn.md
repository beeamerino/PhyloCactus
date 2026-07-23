# Standardized IUCN Category Color Scale for ggplot2

Applies a standardized, publication-grade hex color palette representing
official IUCN Red List threat categories (EX, EW, CR, EN, VU, NT, LC,
DD, NE) onto \`ggplot2\` comparative phylogenetic graphics.

## Usage

``` r
scale_fill_iucn(...)
```

## Arguments

- ...:

  Additional arguments passed directly to
  \`ggplot2::scale_fill_manual\`.

## Value

A \`ggplot2\` manual scale fill layer object.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
ggplot(df, aes(x = species, fill = iucn_category)) +
  geom_bar() +
  scale_fill_iucn()
} # }
```
