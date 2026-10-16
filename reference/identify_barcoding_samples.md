# Identify the samples of a sample sheet

Reads a sheet with one row per sample and routes each row by the matrix
of accepted inputs of the package
(`inst/extdata/barcoding_accepted_inputs.csv`): an admitted input goes
to
[`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md)
(`sanger`, `plastome`, `genbank`; the sequences of a row are one
sample),
[`identify_barcoding_assembly()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_assembly.md)
(`assembly`) or
[`identify_barcoding_reads()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_reads.md)
(`reads_paired`), and gets its report
([`report_barcoding_identification()`](https://beeamerino.github.io/PhyloCactus/reference/report_barcoding_identification.md))
with the declared species and the voucher; any other input is not run
and gets the message of its row. An error in one sample is recorded and
the others go on. A sample without a voucher is identified and flagged:
without it the leakage rule per specimen cannot be applied to it.

## Usage

``` r
identify_barcoding_samples(
  sheet,
  library_dir = file.path("11_barcoding", "4_library"),
  metrics_dir = file.path("11_barcoding", "11_metrics"),
  output_dir = file.path("11_barcoding", "10_identify", "samples"),
  threads = 4L
)
```

## Arguments

- sheet:

  Character or data frame. The sample sheet, a CSV path or a data frame.

- library_dir, metrics_dir:

  As in
  [`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md).

- output_dir:

  Character. Where the samples and the index are written.

- threads:

  Integer. Threads of the external tools (assemblies and reads).

## Value

Invisibly, the index table.

## Details

Columns of the sheet: `sample_id` (unique), `input_type` (a row of the
matrix), `file_1`, and optionally `file_2` (the second file of paired
reads), `declared_species` and `voucher`.

Writes one folder per identified sample in `output_dir` and
`TABLE_barcoding_samples_index.csv`: status (`identified`, `not_run`,
`error`), message, loci found, answers per state over the loci found,
species and genera named, the answer against the declared species
(`species agrees`, `genus agrees`, `disagrees`, `not assigned`,
`not declared`), whether the voucher is missing, and the report.

## See also

[`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md),
[`identify_barcoding_assembly()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_assembly.md),
[`identify_barcoding_reads()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_reads.md).

## Examples

``` r
if (FALSE) { # \dontrun{
identify_barcoding_samples("samples.csv")
} # }
```
