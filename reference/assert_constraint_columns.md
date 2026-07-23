# Assert Mandatory Constraint Table Columns

Verifies that the input taxonomic constraint data frame contains
required classification columns.

## Usage

``` r
assert_constraint_columns(constraints_tbl)
```

## Arguments

- constraints_tbl:

  Data frame containing taxonomic classification hierarchy.

## Value

Invisible \`TRUE\` if validation passes; throws error otherwise.
