# Summarise the Per-Class Metrics of the Molecular Diagnostic Branch

Step 10 of the branch. For each classifier at its operating point of
step 9, per locus and per scheme, it writes the metrics of the
validation: per-class recall and precision as distributions (median,
interquartile range, the classes with recall 0 and the classes never
assigned), the unassigned rate and the two kinds of misassignment as
separate figures, and the five context columns of the locus on every
row. It runs no classifier and computes no operating point: both come
from steps 7 and 9.

## Usage

``` r
summarise_barcoding_metrics(
  classifier_dir = file.path("11_barcoding", "7_classifier"),
  threshold_dir = file.path("11_barcoding", "9_threshold"),
  library_dir = file.path("11_barcoding", "4_library"),
  output_dir = file.path("11_barcoding", "11_metrics"),
  methods = c("nn_add", "idtaxa"),
  idtaxa_rule = "default_60",
  figures = TRUE
)
```

## Arguments

- classifier_dir:

  Character. Directory of step 7, with the prediction tables.

- threshold_dir:

  Character. Directory of step 9, with the operating and curve tables.

- library_dir:

  Character. Directory of step 4, with
  `TABLE_barcoding_library_summary.csv`.

- output_dir:

  Character. Where the tables and figures are written.

- methods:

  Character vector, any of `"nn_add"` (the nearest neighbour through the
  identification path, at `t*`) and `"idtaxa"`.

- idtaxa_rule:

  Character. The row of the IdTaxa operating table used as its operating
  point. Defaults to `"default_60"`, as decided on 2026-09-26.

- figures:

  Logical. Draw the curve figures. Defaults to `TRUE`.

## Value

Invisibly, a list with the class and summary tables of each method.

## Details

The class is the species in scheme E and the genus in scheme G. In
scheme E an assignment is a state 1; a state 2 with the right genus is
an honest genus, reported as its own rate, not as a correct species. In
scheme G an assignment is a state 2; a state 1 names a species that is
absent from the training set and is a misassignment even when its genus
is right. The unassigned, the wrong species and the wrong genus are
never summed into one error.

Precision is NA, not 0, for a class that is never assigned, and those
classes are counted. The interquartile range uses
[`stats::quantile()`](https://rdrr.io/r/stats/quantile.html) type 7.

With `figures = TRUE` it draws, for each locus, the curve of step 9 as
the share of correct assignments at the rank of the scheme against the
unassigned share, one panel per scheme, both classifiers, their
operating points marked.

## See also

[`sweep_barcoding_threshold()`](https://beeamerino.github.io/PhyloCactus/reference/sweep_barcoding_threshold.md),
[`classify_barcoding_folds()`](https://beeamerino.github.io/PhyloCactus/reference/classify_barcoding_folds.md).

## Examples

``` r
if (FALSE) { # \dontrun{
summarise_barcoding_metrics()
} # }
```
