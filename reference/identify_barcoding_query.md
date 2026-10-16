# Identify user sequences against the barcoding library (step 11 of Tutorial 5)

Classifies each query against the IdTaxa model of each locus of the
library, or of the one declared, and answers locus by locus with the
three states of the branch. It never combines loci into one answer.

## Usage

``` r
identify_barcoding_query(
  query,
  locus = NULL,
  library_dir = file.path("11_barcoding", "4_library"),
  metrics_dir = file.path("11_barcoding", "11_metrics"),
  threshold_dir = file.path("11_barcoding", "9_threshold"),
  genus_only_loci = "rbcL",
  output_dir = file.path("11_barcoding", "10_identify"),
  run_name = "query",
  threshold = 60,
  min_overlap = 100L,
  seed = 1L,
  processors = 1L,
  sample = NULL,
  report = TRUE
)
```

## Arguments

- query:

  A FASTA or GenBank flat-file path (a GenBank record is named by its
  accession), a
  [`Biostrings::DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  or a named character vector of sequences.

- locus:

  Character or `NULL`. Loci to test; `NULL` tests every locus of the
  library.

- library_dir:

  Character. Directory of step 4, with the `LIB_<locus>.fasta` files.

- metrics_dir:

  Character. Directory of step 10, with
  `TABLE_barcoding_metrics_summary_idtaxa.csv`.

- threshold_dir:

  Character or `NULL`. Directory of step 9, with
  `TABLE_barcoding_species_threshold_idtaxa.csv`, written by
  `sweep_barcoding_threshold(method = "idtaxa")`: the species threshold
  of each locus (step 9). A row that IdTaxa names at species rank keeps
  the species only when its species confidence reaches that threshold;
  otherwise it is named at genus rank, with reason
  `below_species_threshold`, and keeps in `species_idtaxa` and
  `species_confidence` what IdTaxa reached. The value used is in the
  column `species_threshold`. Without the table, or for a locus it
  lacks, the species threshold is `threshold`, with a message.

- genus_only_loci:

  Character. Loci that answer at genus rank only: a row of one of them
  that IdTaxa names at species rank, over its species threshold, is
  named at genus rank, with reason `genus_only_locus`, and keeps in
  `species_idtaxa` and `species_confidence` what IdTaxa reached.
  Defaults to `"rbcL"`: on the library of this version, built from
  GenBank, the species answers of `rbcL` were right in 11 of 34 (single
  records, sets of records and plastomes), and its species threshold of
  step 9 is the floor of 60. `character(0)` lets every locus name a
  species.

- output_dir:

  Character. Where the table and the model cache are written.

- run_name:

  Character. Suffix of the output table,
  `TABLE_barcoding_identify_<run_name>.csv`.

- threshold:

  Numeric. IdTaxa confidence threshold. Defaults to 60, the operating
  threshold of the validation folds.

- min_overlap:

  Integer. Minimum length of the region to classify. Defaults to 100.

- seed:

  Integer. Seed of the models and of the queries. Defaults to 1, as in
  Tutorial 5.

- processors:

  Integer. Processors for
  [`DECIPHER::IdTaxa()`](https://rdrr.io/pkg/DECIPHER/man/IdTaxa.html).

- sample:

  Character or `NULL`. A name for the specimen when the sequences of
  `query` are several loci of one specimen: they are then one sample,
  with one report that reads them together. `NULL` (the default): each
  sequence is a sample of its own.

- report:

  Logical. Write the identification report
  ([`report_barcoding_identification()`](https://beeamerino.github.io/PhyloCactus/reference/report_barcoding_identification.md)).
  The run record `RUN_<run_name>.csv` is written in any case.

## Value

Invisibly, the table of answers: one row per query and locus tested.

## Details

For each locus:

1.  The whole query is tested against the 20-mers of the locus on both
    strands, with the rule of the rest of the branch (share of 0.15,
    ratio of 3). If it matches, it is classified as it is
    (`path = "whole"`). This is the path of a Sanger sequence.

2.  If it does not match and it is longer than 1.5 times the longest
    library sequence of the locus, the region of the query with the most
    20-mer hits in a window of that width is cropped and tested again
    (`path = "crop"`). `other_windows` counts the other windows with at
    least half as many hits, since an inverted repeat or an rDNA array
    can hold more than one copy.

3.  The region is cut to the core of the locus, the span of alignment
    columns covered by half or more of its library sequences: beyond its
    first and last 20-mer of the core, the query keeps only as many
    bases as the core itself holds there (`core_trimmed` gives the bases
    removed). A segment held by one or two references only is left out,
    so it cannot name them.

4.  When IdTaxa assigns a genus (state 1 or 2) with 3 or more sequences
    in the locus, the region as it was before step 3 is cut to the core
    of that genus, with the same rule, and classified again if that is
    not the region already classified (`genus_core_trimmed` gives the
    bases removed from it). The core of a genus can be narrower than
    that of the locus (two long references of a genus of short sequences
    cannot then take the queries of its absent species) or wider (a
    genus whose sequences are all longer than the core of the locus
    keeps them).

5.  State 3 is returned with its reason: `no_overlap` (no match),
    `short_overlap` (a region shorter than `min_overlap`) or
    `low_confidence` (IdTaxa's genus confidence under `threshold`).

The model of each locus is IdTaxa trained on the whole library of the
locus with a fixed seed, the same model as the outgroup control of step
8. It is cached in `output_dir/models/` with the md5 of its library file
and trained again only when the library or the seed change. Each query
has its own seed, from `seed`, the locus and its identifier (the
accession of a `Genus_species|accession` name, the whole name
otherwise). The random state of the session is left as it was.

Each row carries, for its locus, the wrong-species rate of IdTaxa
measured in step 10 with the species present
(`validation_ws_rate_species_present`, scheme E) and absent
(`validation_ws_rate_species_absent`, scheme G), when that table exists
and was measured at the same threshold.

A locus whose library holds a single species cannot be trained: its rows
are state 3 with reason `single_species_library`.

Raw reads (FASTQ) are not accepted: paired Illumina reads go to
[`identify_barcoding_reads()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_reads.md),
a genome assembly to
[`identify_barcoding_assembly()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_assembly.md).
Aligned reads (BAM, CRAM, SAM) stop with the command that converts them
to FASTQ.

## See also

[`run_barcoding_controls()`](https://beeamerino.github.io/PhyloCactus/reference/run_barcoding_controls.md),
[`sweep_barcoding_threshold()`](https://beeamerino.github.io/PhyloCactus/reference/sweep_barcoding_threshold.md),
[`summarise_barcoding_metrics()`](https://beeamerino.github.io/PhyloCactus/reference/summarise_barcoding_metrics.md).

## Examples

``` r
if (FALSE) { # \dontrun{
identify_barcoding_query("my_sequences.fasta")
identify_barcoding_query("my_plastome.fasta", run_name = "plastome")
} # }
```
