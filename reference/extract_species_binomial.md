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
`Copiapoa_cinerea`). A hyphenated epithet is kept whole, so
`Opuntia_ficus-indica|KJ773783.1` returns `Opuntia_ficus-indica` and not
`Opuntia_ficus`.
