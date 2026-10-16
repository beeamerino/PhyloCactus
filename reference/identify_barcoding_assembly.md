# Identify a genome assembly against the library

Step 11 for a genome assembly. For each locus of the library, its
sequences (without gaps) are mapped to the assembly with
`minimap2 -x asm20`; the hit with the most matching bases is kept, cut
with `samtools faidx` with `flank` bases on each side (within the
scaffold) and identified with the rules of
[`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md)
for that locus: crop, cut to the core of the locus and of the genus,
IdTaxa, three states. A locus with no hit is state 3, `no_overlap`.

## Usage

``` r
identify_barcoding_assembly(
  assembly,
  locus = NULL,
  library_dir = file.path("11_barcoding", "4_library"),
  metrics_dir = file.path("11_barcoding", "11_metrics"),
  threshold_dir = file.path("11_barcoding", "9_threshold"),
  genus_only_loci = "rbcL",
  output_dir = file.path("11_barcoding", "10_identify"),
  run_name = "assembly",
  flank = 500L,
  low_mapq = 20L,
  threads = 1L,
  minimap2 = "minimap2",
  samtools = "samtools",
  threshold = 60,
  min_overlap = 100L,
  seed = 1L,
  report = TRUE
)
```

## Arguments

- assembly:

  Character. Path to the assembly, FASTA, plain or gzipped.

- locus:

  Character or `NULL`. Loci to look for; `NULL`, every locus of the
  library.

- library_dir, metrics_dir, threshold_dir, genus_only_loci, output_dir,
  threshold, min_overlap, seed:

  As in
  [`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md).

- run_name:

  Character. The table is `TABLE_barcoding_identify_<run_name>.csv`.

- flank:

  Integer. Bases added on each side of the hit. Defaults to 500.

- low_mapq:

  Integer. Mapping quality under which a hit is flagged. Defaults to 20.

- threads:

  Integer. Threads of `minimap2`.

- minimap2, samtools:

  Character. The executables.

- report:

  Logical. Write the identification report
  ([`report_barcoding_identification()`](https://beeamerino.github.io/PhyloCactus/reference/report_barcoding_identification.md)).
  The run record `RUN_<run_name>.csv` is written in any case.

## Value

Invisibly, the table of
[`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md),
one row per locus, with the columns `scaffold`, `scaffold_length`,
`hit_start`, `hit_end`, `mapq`, `matching_bases`, `scaffolds_hit` and
`low_mapq`.

## Details

The best hit is kept even when the locus has more than one copy in the
assembly: the table gives the scaffold, its length, the mapping quality
and the number of scaffolds hit, and a mapping quality under `low_mapq`
is flagged and reported, since the region may be a second copy (for
instance a plastid insertion in the nucleus) and not the locus.

Requires `minimap2` and `samtools` on the PATH (for instance from
bioconda). The assembly may be gzipped; its folder is not written to
(the index goes to a temporary folder).

## See also

[`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md).

## Examples

``` r
if (FALSE) { # \dontrun{
identify_barcoding_assembly("GCA_024363205.1_genomic.fna.gz")
} # }
```
