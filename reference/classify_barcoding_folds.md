# Classify the Validation Folds of the Molecular Diagnostic Branch

Step 7 of the branch. Reads the library of step 4 and the folds of step
5, classifies every query of every fold and writes the raw predictions,
one row per query. It aggregates nothing: no accuracy, no rate, no mean.
The negative controls (step 8) and the metrics (step 10) read these
tables afterwards.

## Usage

``` r
classify_barcoding_folds(
  library_dir = file.path("11_barcoding", "4_library"),
  folds_dir = file.path("11_barcoding", "5_folds"),
  output_dir = file.path("11_barcoding", "7_classifier"),
  method = c("nn", "idtaxa"),
  model = "raw",
  min_comparable = 100L,
  schemes = c("species", "genus"),
  loci = NULL,
  max_folds = NULL,
  seed = 1L,
  threshold = 60,
  chunk = NULL,
  n_chunks = NULL,
  alignment = c("library", "add"),
  mafft_exec = "mafft",
  mafft_opts = "--auto",
  workers = 1L,
  notify = FALSE,
  notify_to = NULL,
  notify_credentials = NULL
)
```

## Arguments

- library_dir:

  Character. Output directory of
  [`finalize_barcoding_library()`](https://beeamerino.github.io/PhyloCactus/reference/finalize_barcoding_library.md).

- folds_dir:

  Character. Output directory of
  [`build_barcoding_folds()`](https://beeamerino.github.io/PhyloCactus/reference/build_barcoding_folds.md).

- output_dir:

  Character. Destination; refused inside a directory of the phylogeny.

- method:

  Character. `"nn"`, the branch classifier, or `"idtaxa"`, the contrast.

- model:

  Character. Distance model of
  [`ape::dist.dna()`](https://rdrr.io/pkg/ape/man/dist.dna.html), for
  `"nn"`.

- min_comparable:

  Integer. Minimum number of comparable positions for a pair to have a
  value Pairs below it are left without value, never imputed.

- schemes:

  Character vector. `"species"` (scheme E) and `"genus"` (scheme G).

- loci:

  Character vector or `NULL` for every locus of the library.

- max_folds:

  Integer or `NULL`. Folds per locus and scheme; a declared random
  subset, for the contrast with IdTaxa, which retrains once per fold.
  `NULL` means all of them.

- seed:

  Integer. Seed of that subset, so it is reproducible, and, for
  `"idtaxa"`, the seed from which the seed of every query is derived
  with its locus, scheme, fold and `sid`, so that each answer is
  reproducible and independent of the rest of the run.

- threshold:

  Numeric. Confidence threshold applied to the output of
  [`DECIPHER::IdTaxa()`](https://rdrr.io/pkg/DECIPHER/man/IdTaxa.html),
  for `"idtaxa"`. IdTaxa itself always runs at threshold 0 and the
  confidence of each rank is written (`genus_confidence`,
  `species_confidence`), so any other threshold can be applied
  afterwards with the same result as a run at that threshold.

- chunk, n_chunks:

  Integers, or `NULL` for a run without chunks. With both, only chunk
  `chunk` of `n_chunks` is classified: inside each locus and scheme the
  folds, in the order of step 5, are dealt round-robin to the chunks.
  The chunk writes its tables in `chunks/` of `output_dir`, with the
  suffix `_chunk<k>of<n>`, together with its session; the tables of a
  whole run are then written by
  [`merge_barcoding_chunks()`](https://beeamerino.github.io/PhyloCactus/reference/merge_barcoding_chunks.md).
  For a job array on a cluster, see
  [`generate_barcoding_job_scripts()`](https://beeamerino.github.io/PhyloCactus/reference/generate_barcoding_job_scripts.md).

- alignment:

  Character. How the distance of a query is measured, for `"nn"`.
  `"library"`, the default, reads it from the joint alignment of the
  library, in which the query took part. `"add"` removes the query's
  gaps, orients it against the training set of its fold and adds it to
  that training alignment with `MAFFT --add --keeplength`, one call per
  query: the path of the outgroup control and of a query brought by a
  user. On 2026-09-25 the two paths agreed on 699 of 720 queries and 4
  answers changed. `"add"` writes its own tables, with the suffix
  `_add`, and never touches those of `"library"`.

- mafft_exec, mafft_opts:

  Command of the `MAFFT` binary and its options, for
  `alignment = "add"`.

- workers:

  Integer. `1`, the default, runs in this R session. More than 1 splits
  the folds into that many chunks (as `chunk` and `n_chunks` do) and
  runs them as parallel R processes of this machine, each writing its
  log in `chunks/` of `output_dir`; the progress is written to
  `PROGRESS_<method>.txt` of `output_dir` every minute and printed; when
  every worker has finished the chunks are merged with
  [`merge_barcoding_chunks()`](https://beeamerino.github.io/PhyloCactus/reference/merge_barcoding_chunks.md),
  into the same tables as a run of one process, and `notify` sends one
  email. A worker that fails stops the others and the run. On macOS the
  machine is kept awake with `caffeinate` while the run lasts. The
  workers load the installed package. Cannot be combined with `chunk`.

- notify:

  Logical. Send an email when the run ends, whether it finished or
  failed. Intended for `method = "idtaxa"`, which retrains once per fold
  and is measured in hours. A notification that cannot be sent is
  reported and ignored: it never fails the run. See
  [`send_run_notification()`](https://beeamerino.github.io/PhyloCactus/reference/send_run_notification.md)
  for what has to be configured. Defaults to `FALSE`.

- notify_to, notify_credentials:

  Passed to
  [`send_run_notification()`](https://beeamerino.github.io/PhyloCactus/reference/send_run_notification.md)
  as `to` and `credentials`. Both default to `NULL`, which reads the
  `MY_EMAIL` and `PHYLOCACTUS_SMTP_CREDS` environment variables.

## Value

Invisibly, a list of prediction tables by scheme. Writes
`TABLE_barcoding_predictions_<scheme>_<method>.csv` and
`TABLE_barcoding_timing_<method>.csv`, one row per locus and scheme with
the folds, the queries and the seconds they took; with
`alignment = "add"` both names carry the suffix `_add` and the
prediction tables add the columns `alignment` and `orientation`.

## Details

The output has three states: species, genus with the specific ambiguity
declared, and not assignable. The state comes from the structure of the
tie at the minimum distance, not from a threshold, so the threshold can
be swept later over `nn_distance`, `margin` and `confidence` without
classifying again.

## Examples

``` r
if (FALSE) { # \dontrun{
classify_barcoding_folds(
  library_dir = "11_barcoding/4_library",
  folds_dir = "11_barcoding/5_folds",
  output_dir = "11_barcoding/7_classifier"
)
} # }
```
