# Standardize Taxon Name Strings

Standardizes scientific taxon binomials by replacing whitespace, dashes,
and duplicate underscores with single underscores.

## Usage

``` r
standardize_taxon(x)
```

## Arguments

- x:

  Character vector of taxon binomial strings.

## Value

Character vector of formatted taxon names.

## Examples

``` r
standardize_taxon(c("Opuntia ficus-indica", " Pereskia  aculeata "))
#> [1] "Opuntia_ficus_indica" "Pereskia_aculeata"   
```
