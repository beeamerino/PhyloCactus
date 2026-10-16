# Build the Validation Folds of the Molecular Diagnostic Branch

Step 5 of the branch. Reads the library written by
[`finalize_barcoding_library()`](https://beeamerino.github.io/PhyloCactus/reference/finalize_barcoding_library.md),
builds the two partition schemes (E and G), checks the hard rule on
every fold and writes the fold tables and their summary. The training
set is not written: it is the rest of the library of the locus, and
writing it would multiply the size of the tables without adding
information.

## Usage

``` r
build_barcoding_folds(
  library_dir = file.path("11_barcoding", "4_library"),
  output_dir = file.path("11_barcoding", "5_folds"),
  schemes = c("species", "genus")
)
```

## Arguments

- library_dir:

  Character. Output directory of
  [`finalize_barcoding_library()`](https://beeamerino.github.io/PhyloCactus/reference/finalize_barcoding_library.md).

- output_dir:

  Character. Destination; refused inside a directory of the phylogeny.

- schemes:

  Character vector. Schemes to build; defaults to both.

## Value

Invisibly, a list with `folds` (by scheme) and `summary`. Writes
`TABLE_barcoding_folds_species.csv`, `TABLE_barcoding_folds_genus.csv`
and `TABLE_barcoding_folds_summary.csv`.

## Examples

``` r
if (FALSE) { # \dontrun{
build_barcoding_folds(
  library_dir = "11_barcoding/4_library",
  output_dir = "11_barcoding/5_folds"
)
} # }
```
