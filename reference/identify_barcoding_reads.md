# Identify paired Illumina reads through route A

Step 11 for genome-skim or WGS reads. The reads are cleaned with
`fastp`; the plastome (`embplant_pt`) and the nrDNA (`embplant_nr`) are
assembled with GetOrganelle (`-R 15` and `-R 10`, `-k 21,45,65,85,105`);
the path sequences of each (the `graph1.1` path first) are identified
with
[`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md),
and for each locus the row with the longest region is kept.

## Usage

``` r
identify_barcoding_reads(
  reads1,
  reads2,
  library_dir = file.path("11_barcoding", "4_library"),
  metrics_dir = file.path("11_barcoding", "11_metrics"),
  threshold_dir = file.path("11_barcoding", "9_threshold"),
  genus_only_loci = "rbcL",
  output_dir = file.path("11_barcoding", "10_identify"),
  run_name = "reads",
  threads = 4L,
  keep_intermediate = FALSE,
  fastp = "fastp",
  getorganelle = "get_organelle_from_reads.py",
  getorg_path = Sys.getenv("GETORG_PATH", "~/.GetOrganelle"),
  runner = .bc_reads_runner,
  threshold = 60,
  min_overlap = 100L,
  seed = 1L,
  report = TRUE
)
```

## Arguments

- reads1, reads2:

  Character. The two FASTQ files of the run; `reads2 = NULL`
  (single-end) stops.

- library_dir, metrics_dir, threshold_dir, genus_only_loci, output_dir,
  threshold, min_overlap, seed:

  As in
  [`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md).

- run_name:

  Character. Name of the run folder and of the table.

- threads:

  Integer. Threads of `fastp` and GetOrganelle.

- keep_intermediate:

  Logical. Keep the cleaned reads and the GetOrganelle work folders.

- fastp, getorganelle:

  Character. The executables.

- getorg_path:

  Character. The GetOrganelle folder of databases (`GETORG_PATH`, by
  default `~/.GetOrganelle`).

- runner:

  Function `(cmd, args, stdout, stderr)` that runs an external step, as
  [`system2()`](https://rdrr.io/r/base/system2.html); the default runs
  the tools. The tests pass a fake one.

- report:

  Logical. Write the identification report
  ([`report_barcoding_identification()`](https://beeamerino.github.io/PhyloCactus/reference/report_barcoding_identification.md)).
  The run record `RUN_<run_name>.csv` is written in any case.

## Value

Invisibly, the table of
[`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md),
one row per locus, with the columns `contig` and `organelle` (`plastome`
or `nrdna`).

## Details

Accepted: paired Illumina reads, two FASTQ files, plain or gzipped.
Single-end reads stop (not tested in 0.5.0), and so do long reads
(median length above 600 bases over the first 10 000 reads; not admitted
in 0.5.0). Aligned reads (BAM, CRAM, SAM) stop with the command that
converts them to FASTQ. Target capture, RAD-seq and GBS, RNA-seq and
mixed amplicons are rejected in 0.5.0 (`barcoding_accepted_inputs.csv`)
but cannot be told from the files; they should not be given to this
function.

A GetOrganelle step that fails, or exits 0 without writing a path, is
reported with a message and its stderr is kept; the loci not found then
have reason `assembly_failed`, not `no_overlap`, since their absence was
not observed. Before any step runs, the executables and the GetOrganelle
seed and label databases (`embplant_pt`, `embplant_mt`, `embplant_nr`)
are checked.

Writes to `output_dir/<run_name>/`: `fastp.json`, `plastome_path.fasta`,
`nrdna_path.fasta`, the `STDERR_<step>.txt` of each step, `MANIFEST.csv`
(md5 of inputs and outputs, versions of the tools, commands) and
`TABLE_barcoding_identify_<run_name>.csv`. With
`keep_intermediate = FALSE` the cleaned reads, the fastp HTML report and
the GetOrganelle work folders are deleted; the input reads are never
touched. The identification models are cached in `output_dir/models`.

## See also

[`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md),
[`identify_barcoding_assembly()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_assembly.md).

## Examples

``` r
if (FALSE) { # \dontrun{
identify_barcoding_reads("SRR_1.fastq.gz", "SRR_2.fastq.gz", run_name = "SRR", threads = 8)
} # }
```
