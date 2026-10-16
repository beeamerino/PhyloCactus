# Merge the Chunks of a Chunked Run of Step 7

Reads the tables that
`classify_barcoding_folds(chunk = k, n_chunks = n)` wrote in `chunks/`
of `output_dir`, checks them against the folds of step 5 and writes the
tables of the whole run, which are the tables a run without chunks would
have written, row for row and in the same order.

## Usage

``` r
merge_barcoding_chunks(
  library_dir = file.path("11_barcoding", "4_library"),
  folds_dir = file.path("11_barcoding", "5_folds"),
  output_dir = file.path("11_barcoding", "7_classifier"),
  method = c("nn", "idtaxa"),
  n_chunks,
  alignment = c("library", "add"),
  schemes = c("species", "genus"),
  loci = NULL,
  max_folds = NULL,
  seed = 1L
)
```

## Arguments

- library_dir, folds_dir, output_dir, method, alignment, schemes, loci,
  max_folds, seed:

  As in
  [`classify_barcoding_folds()`](https://beeamerino.github.io/PhyloCactus/reference/classify_barcoding_folds.md),
  with the values the chunks were run with.

- n_chunks:

  Integer. Number of chunks of the run.

## Value

Invisibly, a list of prediction tables by scheme. Writes
`TABLE_barcoding_predictions_<scheme>_<method>.csv` and
`TABLE_barcoding_timing_<method>.csv` in `output_dir`.

## Details

The merge refuses, naming what it found, when a chunk is missing, when a
fold appears that is not among the folds planned for the run (the folds
of step 5, with the declared sample of `max_folds`), when a query
appears more than once, and when a planned query is missing. The
arguments that choose the folds (`schemes`, `loci`, `max_folds`, `seed`)
must be those of the chunks. The timing table sums the seconds of the
chunks by locus and scheme; the chunk tables and their sessions are left
in `chunks/`.

## See also

[`classify_barcoding_folds()`](https://beeamerino.github.io/PhyloCactus/reference/classify_barcoding_folds.md),
[`generate_barcoding_job_scripts()`](https://beeamerino.github.io/PhyloCactus/reference/generate_barcoding_job_scripts.md).

## Examples

``` r
if (FALSE) { # \dontrun{
merge_barcoding_chunks(
  library_dir = "11_barcoding/4_library",
  folds_dir = "11_barcoding/5_folds",
  output_dir = "11_barcoding/7_classifier",
  method = "idtaxa", n_chunks = 60
)
} # }
```
