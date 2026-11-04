# Fetch IUCN Red List Metadata via rredlist API

Queries the IUCN Red List API through `rredlist` and returns the latest
assessment of one species: the Red List category (e.g., CR, EN, VU) and
its criteria, the population trend, the countries of occurrence and of
endemism, habitats, threats, and the texts of the assessment. Enriching
phylogenetic datasets with IUCN conservation statuses enables
evolutionary distinctness and extinction vulnerability analyses.

## Usage

``` r
get_iucn_data(sp_name)
```

## Arguments

- sp_name:

  Character. Standardized binomial scientific name of the target species
  (e.g., `"Astrophytum asterias"`).

## Value

A single-row `tibble` containing extracted IUCN conservation metadata
fields.

## Details

The table returned is Red List Data under the IUCN Red List Terms of
Use: it is for the user's own analysis, queried with the user's own
token, and is not redistributed (sections 1 and 4 of the terms).
[`summarise_iucn_registry()`](https://beeamerino.github.io/PhyloCactus/reference/summarise_iucn_registry.md)
reduces it to the category of each species and to counts, which is what
the package distributes.

## Examples

``` r
if (FALSE) { # \dontrun{
get_iucn_data("Astrophytum asterias")
} # }
```
