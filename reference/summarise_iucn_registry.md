# Reduce the IUCN Red List table of a run to the category of each species and to counts

Reads the species table that Module 11 writes with the assessments of
[`get_iucn_data()`](https://beeamerino.github.io/PhyloCactus/reference/get_iucn_data.md)
(`9_Visualization/tables/TABLE_species_all_data.csv`) and returns three
tables that hold no assessment content: the species registry with
`iucn_found`, `iucn_category` and `iucn_category_name` as its only IUCN
columns, and two count tables. Under the IUCN Red List Terms of Use
(version 3.1), IUCN places no restriction on the use of the category of
a named taxon (section 3), whereas the assessments themselves are not to
be redistributed (section 4); counts are published with the citation of
the Red List and its version (section 6). The figures of endemism,
threats and habitats of Tutorial 3 are drawn from the count tables.

## Usage

``` r
summarise_iucn_registry(
  species_table,
  redlist_version,
  access_date,
  output_dir = NULL
)
```

## Arguments

- species_table:

  A data frame, or the path to the CSV file, with the columns of the
  species registry and of
  [`get_iucn_data()`](https://beeamerino.github.io/PhyloCactus/reference/get_iucn_data.md).

- redlist_version:

  Character. Version of the Red List the assessments were taken from, as
  [`rredlist::rl_version()`](https://docs.ropensci.org/rredlist/reference/rl_version.html)
  reports it (e.g., `"2026-1"`).

- access_date:

  Character. Date the assessments were retrieved, `YYYY-MM-DD`.

- output_dir:

  Character or `NULL`. Folder where the three tables are written as
  `TABLE_species_registry.csv`, `TABLE_iucn_endemic_by_country.csv` and
  `TABLE_iucn_classes.csv`. `NULL` (the default) writes nothing.

## Value

A list of three tibbles:

- `registry`: the species table without its `iucn_` columns other than
  `iucn_found`, `iucn_category` and `iucn_category_name`, and without
  the checklist columns `checklist_fullName`, `checklist_author`,
  `checklist_rank`, `checklist_taxon` and `checklist_uuid` (the
  checklist of Korotkova et al. 2021 may be used but not redistributed);

- `endemic_by_country`: `country`, `iucn_category`, `n` (endemic species
  of each country in each category), `redlist_version`, `access_date`;

- `classes`: `field` (`threat`, `habitat` or `endemic_species_in_tree`),
  `class`, `n`, `redlist_version`, `access_date`. The row of
  `endemic_species_in_tree` (class `all`) is the number of endemic
  species, each counted once.

## Details

The counts cover the ingroup species in the maximum-likelihood tree
(`is_ingroup` and `in_ml_tree` true), as the figures of Tutorial 3 do.
Countries, threats and habitats are split on `"; "`, the separator
[`get_iucn_data()`](https://beeamerino.github.io/PhyloCactus/reference/get_iucn_data.md)
writes, so a species is counted once in every country, threat class or
habitat class it has.

## References

IUCN (2026). *The IUCN Red List of Threatened Species*. Version 2026-1.
<https://www.iucnredlist.org>

## Examples

``` r
if (FALSE) { # \dontrun{
summarise_iucn_registry("9_Visualization/tables/TABLE_species_all_data.csv",
                        redlist_version = "2026-1", access_date = "2026-10-15",
                        output_dir = "9_Visualization/tables")
} # }
```
