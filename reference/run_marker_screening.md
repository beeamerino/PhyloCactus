# Screen Locus Alignments for Substitution Saturation and Informativeness

Evaluates phylogenetic informativeness, sequence coverage, alignment
length, and substitution saturation across individual locus alignments.
Filtering out loci exhibiting high substitution saturation or severe
site length anomalies prevents systematic noise and long-branch
attraction (LBA) artifacts from distorting maximum-likelihood
supermatrix inference.

## Usage

``` r
run_marker_screening(
  fasta_folder,
  out_base,
  min_cols_to_evaluate = 50L,
  min_aln_len_to_retain = 200L,
  min_nseq_to_retain = 100L,
  outgroup_folder = NULL,
  max_marker_missing = 0.7,
  saturation_flag_cutoff = 0.3,
  saturation_keep_cutoff = 0.5,
  iqr_multiplier = 1.5,
  saturation_method = c("corrected", "legacy"),
  gap_handling = c("iupac", "legacy")
)
```

## Arguments

- fasta_folder:

  Character. Directory path containing aligned locus FASTA files.

- out_base:

  Character. Base destination directory for diagnostic plots and
  screened FASTA outputs (`filtered_markers/`).

- min_cols_to_evaluate:

  Integer. Minimum number of alignment columns required to compute
  saturation metrics. Defaults to `50L`.

- min_aln_len_to_retain:

  Integer. Minimum alignment length in base pairs required to retain a
  locus. Defaults to `200L`.

- min_nseq_to_retain:

  Integer. Minimum number of ingroup sequences required per locus,
  counted after outlier filtering. Defaults to `100L`.

- outgroup_folder:

  Character or `NULL`. Directory of aligned outgroup locus FASTA files,
  matched to the ingroup files by marker name. Used for **reporting
  only**: the summary table gains `n_outgroup` and `n_total` so a locus
  rejected here can be seen to carry outgroup data, but no retention
  decision depends on them. Outgroup markers with no ingroup counterpart
  simply leave those two columns empty; which ones they are, and what
  becomes of them, is reported by
  [`integrate_and_clean_markers()`](https://beeamerino.github.io/PhyloCactus/reference/integrate_and_clean_markers.md),
  which owns that decision and applies `marker_aliases` before comparing
  the two sets. Defaults to `NULL`.

  This module asks whether a locus resolves the ingroup radiation, and
  the answer cannot depend on how many outgroup accessions exist.
  `trnT-psbD` is the case that forced the distinction: 50 ingroup
  sequences against a threshold of 100, correctly rejected as an ingroup
  marker, while carrying 49 Portulaca accessions of 1347 bp that are the
  best outgroup coverage in the dataset. Bringing it back is a decision
  about connecting the two groups, which belongs to
  [`integrate_and_clean_markers()`](https://beeamerino.github.io/PhyloCactus/reference/integrate_and_clean_markers.md)
  and its `readmit_markers` argument, not to a sequence count here.

- max_marker_missing:

  Numeric. Maximum allowable missing data fraction per locus. Defaults
  to `0.7`.

- saturation_flag_cutoff:

  Numeric. Uncorrected p-distance vs. raw distance slope threshold to
  flag substitution saturation. Defaults to `0.3`.

- saturation_keep_cutoff:

  Numeric. Saturation slope cutoff threshold below which saturated loci
  are excluded. Defaults to `0.5`.

- iqr_multiplier:

  Numeric. Interquartile range (IQR) multiplier for identifying
  site-length outlier bounds. Defaults to `1.5`.

- saturation_method:

  Character. Method for saturation test distance calculation. Defaults
  to `"corrected"`.

- gap_handling:

  Character. Method for handling ambiguous and gap characters. Defaults
  to `"iupac"`.

## Value

A data frame summarizing saturation statistics, alignment dimensions,
and retention decisions across screened loci.

## Examples

``` r
if (FALSE) { # \dontrun{
run_marker_screening(
  fasta_folder = "2_MAFFT_Cactaceae/alignments",
  out_base = "3_Screening_Ingroup",
  min_aln_len_to_retain = 200L,
  min_nseq_to_retain = 50L
)
} # }
```
