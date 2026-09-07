# Age of Every Calibrated Node Across a Range of Rate-Smoothing Values

Penalized likelihood requires a rate-smoothing parameter, and treePL
selects one by cross-validation. When the cross-validated minimum falls
on the edge of the tested grid, as it does for this dataset, the
selection is the boundary of the search rather than an optimum, and a
reader is entitled to ask whether the reported ages are an artefact of
that choice. This function answers the question directly: it dates the
same tree, under the same calibrations, at each smoothing value given,
and reports the age of every calibrated node in each run.

## Usage

``` r
report_smoothing_sensitivity(
  cfg_file,
  out_csv = NULL,
  smoothing_values = c(1e-04, 0.01, 1, 10, 100),
  treepl_bin = "treePL",
  timeout = 1800,
  work_dir = NULL
)
```

## Arguments

- cfg_file:

  Character. A treePL configuration carrying the tree, `numsites`, the
  optimisation parameters and the calibrations. The configuration
  written for the maximum-likelihood chronogram is the natural input.

- out_csv:

  Character or `NULL`. Where to write the table. Defaults to
  `TABLE_smoothing_sensitivity.csv` beside `cfg_file`.

- smoothing_values:

  Numeric vector of smoothing values to test. Defaults to five values
  spanning six orders of magnitude, which is wide enough that stability
  across them is informative.

- treepl_bin:

  Character. The treePL executable. Defaults to `"treePL"`.

- timeout:

  Numeric. Seconds allowed per run before it is abandoned and recorded
  as unconverged. Defaults to 1800.

- work_dir:

  Character or `NULL`. Directory for the intermediate configurations,
  logs and chronograms. Defaults to a `smoothing_sensitivity` directory
  beside `cfg_file`, so the evidence behind the table survives alongside
  it.

## Value

Invisibly, a data frame with one row per smoothing value and calibrated
node: `smoothing`, `smoothing_used`, `node`, `age_ma`, `seconds` and
`converged`.

## Details

The nodes come from the `mrca` lines of the configuration itself, so the
table covers exactly the nodes the analysis makes claims about and
cannot drift out of step with them.

Each run is checked with the same verification applied everywhere else:
the smoothing treePL reports in its log must match the smoothing
requested. Until 2026-09-02 the pipeline wrote the keyword `smoothing`,
which treePL does not recognise and discards without a message, so every
chronogram was produced at the built-in default of 10 and a table like
this one would have shown five identical rows. A run whose smoothing
cannot be confirmed is reported as such rather than tabulated as a
result.

## Examples

``` r
if (FALSE) { # \dontrun{
report_smoothing_sensitivity(
  cfg_file = "8_Dating/auto_results/ML_tree/configure_smooth_ML_tree"
)
} # }
```
