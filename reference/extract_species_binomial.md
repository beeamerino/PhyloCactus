# Extract Standardized Species Binomial from Sequence Header

Internal helper to parse species binomials from FASTA headers formatted
as `Genus_species|accession|sid` or `Genus species`.

## Usage

``` r
extract_species_binomial(header)
```

## Arguments

- header:

  Character string of the sequence header.

## Value

Standardized species binomial formatted with underscores (e.g.,
`Copiapoa_cinerea`). A hyphenated epithet is kept whole and its hyphen
written as an underscore, the convention of the terminal labels and of
[`standardize_taxon()`](https://beeamerino.github.io/PhyloCactus/reference/standardize_taxon.md)
(0.5.5): `Opuntia_ficus-indica|KJ773783.1` returns
`Opuntia_ficus_indica`, not `Opuntia_ficus`.
