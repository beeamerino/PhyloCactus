# Perform Joint Realignment Across Integrated Ingroup and Outgroup Sequences

Perform Positional Homology Realignment Across Locus Sequence Alignments

## Usage

``` r
run_joint_realignment(
  input_dir,
  output_fasta_dir,
  output_aln_dir,
  min_non_gap_fraction = 0.3,
  max_missing_fraction = 0.3,
  preserve_iupac = TRUE,
  protect_pattern = NULL,
  protect_markers = NULL,
  rooting_pattern = NULL
)
```

## Arguments

- input_dir:

  Character. Directory containing curated locus FASTA files (e.g.,
  `4_Cleaned/cleaned_markers_joint` or
  `4_Cleaned/cleaned_markers_ingroup`).

- output_fasta_dir:

  Character. Directory path to save output realigned FASTA sequence
  files and QC logs.

- output_aln_dir:

  Character. Directory path to store final masked aligned FASTA files.

- min_non_gap_fraction:

  Numeric. Minimum proportion of non-gap characters, relative to the
  masked alignment width, required to retain an individual sequence in a
  locus. Defaults to `0.30`.

- max_missing_fraction:

  Numeric. Maximum proportion of missing characters (`N`) tolerated per
  sequence. Defaults to `0.30`.

- preserve_iupac:

  Logical. Retain IUPAC ambiguity codes instead of collapsing them to
  `N` before realignment. Defaults to `TRUE`, matching
  [`run_alignment_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_alignment_pipeline.md);
  see that function for the rationale.

- protect_pattern:

  Character or `NULL`. Regular expression matched against sequence
  names; matching terminals that carry at least one non-gap character
  are retained regardless of the occupancy filters. Intended as a
  last-resort safeguard for rooting terminals whose sequences are
  legitimately short (for example
  `"^(Anacampseros|Grahamia|Talinopsis|Portulaca)_"`). Defaults to
  `NULL` (no exemption).

- protect_markers:

  Character vector or `NULL`. Names of the markers, as they appear in
  `input_dir` without the file extension, in which `protect_pattern` is
  honoured. `NULL`, the default, applies the exemption to every marker.
  Naming markers restricts it to the loci where a short outgroup
  sequence is worth its gap cost, instead of retaining every fragment of
  every terminal across the whole matrix: in the August 2026 dataset
  only `trnL_trnF` lost rooting terminals, and only there does the
  exemption buy anything.

- rooting_pattern:

  Character or `NULL`. Regular expression identifying the terminals the
  tree will be rooted on. It exempts nothing. Any matching terminal
  removed by the occupancy filters is named in a warning and flagged in
  `LOG_SEQ_FILTER_<marker>.csv`. A filter that deletes the rooting
  outgroup must say so at the moment it does it, not four modules
  downstream: the five *Portulaca* sequences of `trnL_trnF`, 297 bp
  against a threshold of about 353, were removed silently and the branch
  subtending the outgroup collapsed from roughly 600 expected
  substitutions to 0.03. Defaults to `NULL`.

## Value

A data frame containing compiled alignment summary statistics across all
processed markers.

## Details

Re-estimates positional homology alignments (`MAFFT`) across curated
locus FASTA files, performs alignment quality masking with `DECIPHER`,
filters low-occupancy sequences, and generates comprehensive alignment
statistics.

This is the step that removes individual taxa from a locus without
removing the locus itself. Sequences that entered Stage 4 as short
fragments, typically outgroup accessions that were masked separately
from the ingroup in
[`run_alignment_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_alignment_pipeline.md),
fail the occupancy filter here and disappear from the supermatrix.
Passing `mask_alignment_regions = FALSE` for the outgroup in Stage 2 so
that masking happens only once, jointly, at this stage, is preferable to
relaxing these thresholds or resorting to `protect_pattern`, because
retaining very short sequences inflates the gap fraction of the final
supermatrix and can destabilise the affected terminals.

## References

Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment
software version 7: Improvements in performance and usability.
*Molecular Biology and Evolution*, 30(4), 772–780.
[doi:10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)

## Examples

``` r
if (FALSE) { # \dontrun{
run_joint_realignment(
  input_dir = "4_Cleaned/cleaned_markers_joint",
  output_fasta_dir = "5_MAFFT_Cleaned",
  output_aln_dir = "5_MAFFT_Cleaned/aligned_markers"
)
} # }
```
