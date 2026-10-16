# Per-genus discrimination of each locus in the validation folds

From the IdTaxa predictions of the validation folds (step 7), cut at the
operating threshold of step 9 (60 by default), counts per locus, scheme
and true genus the queries whose species was named correctly, the wrong
species, the right and wrong genera, and the unassigned, with their
rates. It says how well each locus tells the species of one genus apart
in the library, which the identification report shows beside the answer.

## Usage

``` r
summarise_barcoding_genus_discrimination(
  classifier_dir = file.path("11_barcoding", "7_classifier"),
  threshold_dir = file.path("11_barcoding", "9_threshold"),
  output_dir = file.path("11_barcoding", "11_metrics"),
  idtaxa_rule = "default_60"
)
```

## Arguments

- classifier_dir:

  Character. Directory of step 7, with
  `TABLE_barcoding_predictions_{species,genus}_idtaxa.csv`.

- threshold_dir:

  Character. Directory of step 9, with
  `TABLE_barcoding_threshold_operating_idtaxa.csv`.

- output_dir:

  Character. Where `TABLE_barcoding_genus_discrimination_idtaxa.csv` is
  written.

- idtaxa_rule:

  Character. Row of the operating table used, `"default_60"` by default.

## Value

Invisibly, the table.

## See also

[`summarise_barcoding_metrics()`](https://beeamerino.github.io/PhyloCactus/reference/summarise_barcoding_metrics.md),
[`report_barcoding_identification()`](https://beeamerino.github.io/PhyloCactus/reference/report_barcoding_identification.md).

## Examples

``` r
if (FALSE) { # \dontrun{
summarise_barcoding_genus_discrimination()
} # }
```
