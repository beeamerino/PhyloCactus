# Sweep the Threshold of Remoteness of the Molecular Diagnostic Branch

Step 9 of the branch. Without a rule of remoteness the nearest neighbour
names a cactus species for any query: in the outgroup control, 288 of
312 outgroup queries came out in state 1. This step applies the rule (a
query whose distance to its nearest neighbour exceeds a threshold is
state 3, reason `remoteness`) over every threshold observed, and fixes
an operating threshold per locus.

## Usage

``` r
sweep_barcoding_threshold(
  classifier_dir = file.path("11_barcoding", "7_classifier"),
  controls_dir = file.path("11_barcoding", "8_controls"),
  output_dir = file.path("11_barcoding", "9_threshold"),
  alignment = c("add", "library"),
  q = NULL,
  quantile_type = 1L,
  figures = TRUE,
  method = c("nn", "idtaxa"),
  species_q = 0.95
)
```

## Arguments

- classifier_dir:

  Character. Output directory of
  [`classify_barcoding_folds()`](https://beeamerino.github.io/PhyloCactus/reference/classify_barcoding_folds.md).

- controls_dir:

  Character. Output directory of
  [`run_barcoding_controls()`](https://beeamerino.github.io/PhyloCactus/reference/run_barcoding_controls.md).
  When it holds no `TABLE_barcoding_cn2_queries.csv` the rejection is
  reported as not measured everywhere.

- output_dir:

  Character. Destination; refused inside a directory of the phylogeny.

- alignment:

  Character. Which predictions of step 7 are swept: `"add"`, the
  default, reads the tables of `alignment = "add"`, measured the way a
  user query is measured; `"library"` reads those of the joint
  alignment.

- q:

  Numeric. Quantile of scheme G that sets the operating threshold.
  Defaults to `0.99` of the distances for `"nn"` and `0.01` of the genus
  confidences for `"idtaxa"`.

- quantile_type:

  Integer. Type of
  [`stats::quantile()`](https://rdrr.io/r/stats/quantile.html). Defaults
  to `1L`, whose result is always a distance observed in the data.

- figures:

  Logical. Write one figure per locus, PDF and PNG at 300 dpi, 180 mm
  wide.

- method:

  Character. `"nn"`, the nearest neighbour and its rule of remoteness,
  or `"idtaxa"`, IdTaxa and its confidence threshold. With `"idtaxa"`
  the tables read are `TABLE_barcoding_predictions_<scheme>_idtaxa.csv`
  and `TABLE_barcoding_cn2_queries_idtaxa.csv`, the curve runs over the
  confidence threshold, the operating threshold of a locus is the
  quantile `q` of `genus_confidence` of its scheme G queries (the rule
  of D1 turned around, since confidence rises where distance falls), and
  it is reported next to 60, DECIPHER's default, in the column `rule`.
  The tables written carry the suffix `_idtaxa`.

- species_q:

  Numeric. With `"idtaxa"` only: the quantile of the species confidence
  of the scheme G queries of each locus (species absent from the
  library) that sets its species threshold, never under 60. Written to
  `TABLE_barcoding_species_threshold_idtaxa.csv` and read by the
  identification functions through `threshold_dir`; the validation
  (step 10) stays at 60. Defaults to `0.95`.

## Value

Invisibly, a list with `curve` and `operating`. Writes
`TABLE_barcoding_threshold_curve.csv`,
`TABLE_barcoding_threshold_operating.csv` and, with `figures = TRUE`,
`FIG_barcoding_threshold_<locus>.pdf` and `.png`.

## Details

For every locus and scheme the curve reports the legitimate queries in
five disjoint outcomes (correct species, wrong species, correct genus,
wrong genus, not assigned; a wrong species and an honest genus are never
summed) and, when the outgroup control has comparable outgroup queries
for the locus, how many of them the rule rejects. Outgroup queries that
match the locus in neither direction are not comparable and stay out.

The operating threshold of a locus is the quantile `q` (type
`quantile_type` of
[`stats::quantile()`](https://rdrr.io/r/stats/quantile.html)) of the
distances of the legitimate queries of scheme G: a cactus whose species
is missing from the library while its genus is present, which is the
query the rule must not reject. It is fixed without the outgroup
queries, which then report the rejection with its exact binomial
interval; a locus with no comparable outgroup query reports it as not
measured.

## Examples

``` r
if (FALSE) { # \dontrun{
sweep_barcoding_threshold(
  classifier_dir = "11_barcoding/7_classifier",
  controls_dir = "11_barcoding/8_controls",
  output_dir = "11_barcoding/9_threshold"
)
} # }
```
