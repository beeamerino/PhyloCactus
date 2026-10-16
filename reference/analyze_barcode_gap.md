# Barcode Gap of the Molecular Diagnostic Branch

Step 6 of the branch. For every locus of the library, computes the
intraspecific distances by pairs and, for every sequence, the distance
to its nearest neighbour of another species, with the `raw` and K80
models. Derives the candidate threshold with this rule: the 95th
percentile of the intraspecific distribution, published next to the 5th
percentile of the interspecific one. When the first is above the second
the locus is declared to have no gap.

## Usage

``` r
analyze_barcode_gap(
  library_dir = file.path("11_barcoding", "4_library"),
  output_dir = file.path("11_barcoding", "6_gap"),
  models = c("raw", "K80"),
  min_comparable = 100L,
  paralog_loci = "pepC_like",
  figures = TRUE
)
```

## Arguments

- library_dir:

  Character. Output directory of
  [`finalize_barcoding_library()`](https://beeamerino.github.io/PhyloCactus/reference/finalize_barcoding_library.md).

- output_dir:

  Character. Destination; refused inside a directory of the phylogeny.

- models:

  Character vector. Distance models of
  [`ape::dist.dna()`](https://rdrr.io/pkg/ape/man/dist.dna.html); both
  are reported.

- min_comparable:

  Integer. Minimum number of positions where both sequences carry A, C,
  G or T for a pair to have a value. Pairs below it are counted apart
  and left out of every distribution and percentile; they are never
  imputed.

- paralog_loci:

  Character vector. Loci flagged in `possible_paralog`.

- figures:

  Logical. Write one figure per locus, with one panel per model.

## Value

Invisibly, a list with `summary` and `distances` (by locus). Writes
`TABLE_barcoding_gap_summary.csv`,
`TABLE_barcoding_distances_<locus>.csv` and
`FIG_barcoding_gap_<locus>.png`.

## Examples

``` r
if (FALSE) { # \dontrun{
analyze_barcode_gap(
  library_dir = "11_barcoding/4_library",
  output_dir = "11_barcoding/6_gap"
)
} # }
```
