# Fetch IUCN Red List Metadata via rredlist API

Queries the IUCN Red List API (\`rredlist\`) to retrieve verified
macroevolutionary conservation threat categories (e.g., CR, EN, VU),
population trends, habitat descriptions, and geographic endemism metrics
for a focal species. Enriching phylogenetic datasets with IUCN
conservation statuses enables evolutionary distinctness and extinction
vulnerability analyses.

## Usage

``` r
get_iucn_data(sp_name)
```

## Arguments

- sp_name:

  Character. Standardized binomial scientific name of the target species
  (e.g., \`"Astrophytum asterias"\`).

## Value

A single-row \`tibble\` containing extracted IUCN conservation metadata
fields.

## Examples

``` r
if (FALSE) { # \dontrun{
get_iucn_data("Astrophytum asterias")
} # }
```
