# Concatenate Locus Alignments and Build Partition Coordinate Maps

Concatenates individual orthologous locus alignments end-to-end into a
unified multilocus supermatrix. Combining independent molecular loci
increases statistical power to resolve difficult ancestral nodes while
allowing partitioned substitution modeling to account for mutational
rate heterogeneity across molecular locus alignments. Exports one
coordinate map per partition format rather than one per downstream
program: a RAxML-style map (`PARTITION_raxml_style.txt`) read by
`RAxML-NG`, `ModelTest-NG` and `IQ-TREE`, and a NEXUS SETS block
(`PARTITION_nexus_charset.nex`). Configuration scripts for
`PartitionFinder2` and `MrBayes` are exported separately under their
program names.

## Usage

``` r
run_concatenation_pipeline(
  input_dir,
  output_dir,
  outgroup_pattern = NULL,
  min_coverage = 0.2,
  exclude_markers = NULL
)
```

## Arguments

- input_dir:

  Character. Path to directory containing curated, aligned locus FASTA
  files.

- output_dir:

  Character. Path to destination root directory for concatenated
  alignments and partition maps.

- outgroup_pattern:

  Character or `NULL`. Regular expression identifying outgroup
  terminals. When supplied,
  [`report_marker_group_coverage()`](https://beeamerino.github.io/PhyloCactus/reference/report_marker_group_coverage.md)
  runs before concatenation and its table is written to
  `logs_and_qc/SUPP_TABLE_marker_group_coverage.csv`. Defaults to
  `NULL`.

- min_coverage:

  Numeric. Fraction of non-gap, non-missing characters at which a
  terminal counts as covered by a marker, passed to
  [`report_marker_group_coverage()`](https://beeamerino.github.io/PhyloCactus/reference/report_marker_group_coverage.md).
  Defaults to `0.2`.

- exclude_markers:

  Character vector or `NULL`. Markers to leave out of the supermatrix,
  named as they appear in the alignment file names without the `Masked_`
  prefix and without the extension. The files are not modified or
  removed, so an excluded locus remains available for other analyses;
  only this supermatrix is built without it. The excluded markers are
  named in the run log and listed in
  `logs_and_qc/TABLE_markers_excluded.csv`. A name that matches no
  alignment raises a warning rather than failing, so a typo is visible
  instead of silent. Defaults to `NULL`.

## Value

A data frame containing supermatrix dimensions, taxon coverage, and
locus partition bounds.

## Examples

``` r
if (FALSE) { # \dontrun{
run_concatenation_pipeline(
  input_dir = "5_MAFFT_Cleaned/aligned_markers",
  output_dir = "6_Concatenated"
)

# A dating matrix without the nuclear locus that carries no outgroup coverage.
run_concatenation_pipeline(
  input_dir = "5_MAFFT_Cleaned/aligned_markers",
  output_dir = "6_Concatenated_dating",
  exclude_markers = "phyC"
)
} # }
```
