# Colours of the Figures of PhyloCactus

Returns the colours of one role of the palette the package uses for
every figure. Each role does one job: `"categorical"` tells unordered
things apart (trees compared, methods, loci, nucleotides), in a fixed
order that is never cycled; `"outcome"` and `"state"` give the answers
of an identification on a diverging scale, blue for a right answer, grey
for none, red for a wrong one; `"sequential"` encodes a magnitude in
greys, from light to dark; `"support"` gives the classes of bootstrap
support of a node as Tutorial 3 draws them, with constrained nodes
apart; `"neutral"` holds the greys of text, lines, rules, backgrounds
and the marks of a single series (bars, histograms, segments, points);
`"theme"` holds the colours of
[`theme_phylocactus()`](https://beeamerino.github.io/PhyloCactus/reference/theme_phylocactus.md);
`"iucn"` keeps the colours of the IUCN Red List categories.

## Usage

``` r
phylocactus_colours(
  role = c("categorical", "outcome", "state", "sequential", "support", "neutral",
    "theme", "iucn")
)
```

## Arguments

- role:

  Character. One of `"categorical"`, `"outcome"`, `"state"`,
  `"sequential"`, `"support"`, `"neutral"`, `"theme"` or `"iucn"`.

## Value

A character vector of hexadecimal colours, named by category where the
role has categories.

## See also

[`theme_phylocactus()`](https://beeamerino.github.io/PhyloCactus/reference/theme_phylocactus.md),
[`scale_fill_iucn()`](https://beeamerino.github.io/PhyloCactus/reference/scale_fill_iucn.md).

## Examples

``` r
phylocactus_colours("outcome")
#> species correct   genus correct      abstention   wrong species     wrong genus 
#>       "#1c5cab"       "#86b6ef"       "#a8a79f"       "#ee8a85"       "#b2182b" 
phylocactus_colours("categorical")[1:3]
#>      blue    orange      aqua 
#> "#2a78d6" "#eb6834" "#1baf7a" 
```
