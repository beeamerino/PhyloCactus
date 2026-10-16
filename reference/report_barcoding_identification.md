# Build the identification report of one run

Reads `TABLE_barcoding_identify_<run_name>.csv` and `RUN_<run_name>.csv`
from `results_dir` and writes, for each query,
`REPORT_<run_name>_<query>.html`, a self-contained page written by the
package (no `pandoc`); an assembly, a set of reads, or several sequences
declared as one sample (argument `sample` of
[`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md))
get one `REPORT_<run_name>.html` that reads them together. For the run,
a set of CSV files: `REPORT_<run_name>_answer.csv`, `_alternatives.csv`,
`_agreement.csv`, `_sampling.csv` and `_provenance.csv`, and
`REPORT_<run_name>_checksums.csv` with the md5 of the table, the run
record and every file written.

## Usage

``` r
report_barcoding_identification(
  run_name,
  results_dir = file.path("11_barcoding", "10_identify"),
  library_dir = file.path("11_barcoding", "4_library"),
  metrics_dir = file.path("11_barcoding", "11_metrics"),
  output_dir = results_dir
)
```

## Arguments

- run_name:

  Character. The run.

- results_dir:

  Character. Where the table and the run record are.

- library_dir:

  Character. Directory of step 4, with the `LIB_<locus>.fasta` files.

- metrics_dir:

  Character or `NULL`. Directory of step 10, where
  `TABLE_barcoding_genus_discrimination_idtaxa.csv` is looked for.

- output_dir:

  Character. Where the report is written; `results_dir` by default.

## Value

Invisibly, the paths of the HTML reports.

## Details

The report has nine sections: (1) query and route; (2) answer per locus;
(3) reading across loci, descriptive only, since the loci are never
combined into one call plan); (4) alternatives under the threshold,
marked not assigned: IdTaxa gives one path per query, so the
alternatives are that path when it falls under the threshold and, in
state 2, the species of the named genus in the library; (5) how far to
trust each locus: the wrong-species rates of the validation and the
per-genus discrimination of
[`summarise_barcoding_genus_discrimination()`](https://beeamerino.github.io/PhyloCactus/reference/summarise_barcoding_genus_discrimination.md);
(6) sampling of the named genera in the library; (7) reads; (8) library,
models and software; (9) limits. The figures follow Zeng et al. (2026,
New Phytologist). The three identification functions call it at the end;
a report can be rebuilt from saved tables.

## See also

[`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md),
[`identify_barcoding_assembly()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_assembly.md),
[`identify_barcoding_reads()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_reads.md).

## Examples

``` r
if (FALSE) { # \dontrun{
report_barcoding_identification("query")
} # }
```
