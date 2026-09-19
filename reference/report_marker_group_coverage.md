# Report How Each Marker Covers the Ingroup and the Outgroup

Counts, per locus, how many ingroup and how many outgroup terminals
carry real sequence, and warns when either side is empty.

## Usage

``` r
report_marker_group_coverage(
  input_dir,
  outgroup_pattern,
  min_coverage = 0.2,
  out_csv = NULL,
  exclude_markers = NULL
)
```

## Arguments

- input_dir:

  Character. Directory of aligned locus FASTA files, one per marker.

- outgroup_pattern:

  Character. Regular expression matched against terminal names, or a
  named character vector of regular expressions (e.g.
  `c(Anacampserotaceae = "...", Portulacaceae = "...", Talinaceae = "...")`)
  to additionally compute per-group coverage columns.

- min_coverage:

  Numeric. Fraction of non-gap, non-missing characters at which a
  terminal counts as covered by that marker. Defaults to `0.2`.

- out_csv:

  Character or `NULL`. Path to write the table to. Defaults to `NULL`.

- exclude_markers:

  Character vector or `NULL`. Markers to leave out of the report, named
  as they appear in `input_dir` without the file extension and without
  any `Masked_` prefix. Passed through by
  [`run_concatenation_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_concatenation_pipeline.md)
  so that the coverage table describes the partitions the supermatrix
  has, not every alignment the folder holds. Defaults to `NULL`.

## Value

Invisibly, a data frame with one row per marker: alignment length,
terminals and covered terminals on each side, and the median coverage of
each side.

## Details

Concatenation assumes that the loci being joined describe the same
terminals. A locus sampled almost entirely on one side of the root
breaks that assumption without breaking anything visible: it adds
columns the other side cannot share, and the branch lengths spanning the
bipartition are then estimated from the loci that remain. The branch
subtending the outgroup can then be badly underestimated, and the
deepest calibrated nodes can return their own bounds instead of an
estimate.

The function reports and does not block. Whether a locus of that kind
belongs in a given matrix depends on which analysis the matrix is for: a
nuclear locus with no outgroup coverage is unusable for dating and
valuable for species discrimination.

## Examples

``` r
if (FALSE) { # \dontrun{
report_marker_group_coverage(
  input_dir = "5_MAFFT_Cleaned/aligned_markers",
  outgroup_pattern = c(
    Anacampserotaceae = "^(Anacampseros|Grahamia|Talinopsis)_",
    Portulacaceae     = "^Portulaca_",
    Talinaceae        = "^(Talinum|Talinella)_"
  )
)
} # }
```
