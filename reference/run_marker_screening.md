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
  max_marker_missing = 0.7,
  saturation_flag_cutoff = 0.3,
  saturation_keep_cutoff = 0.5,
  iqr_multiplier = 1.5
)
```

## Arguments

- fasta_folder:

  Character. Directory path containing aligned locus FASTA files.

- out_base:

  Character. Base destination directory for diagnostic plots and
  screened FASTA outputs (\`filtered_markers/\`).

- min_cols_to_evaluate:

  Integer. Minimum number of alignment columns required to compute
  saturation metrics. Defaults to \`50L\`.

- min_aln_len_to_retain:

  Integer. Minimum alignment length in base pairs required to retain a
  locus. Defaults to \`200L\`.

- min_nseq_to_retain:

  Integer. Minimum number of sequences required per locus alignment.
  Defaults to \`100L\`.

- max_marker_missing:

  Numeric. Maximum allowable missing data fraction per locus. Defaults
  to \`0.7\`.

- saturation_flag_cutoff:

  Numeric. Uncorrected p-distance vs. raw distance slope threshold to
  flag substitution saturation. Defaults to \`0.3\`.

- saturation_keep_cutoff:

  Numeric. Saturation slope cutoff threshold below which saturated loci
  are excluded. Defaults to \`0.5\`.

- iqr_multiplier:

  Numeric. Interquartile range (IQR) multiplier for identifying
  site-length outlier bounds. Defaults to \`1.5\`.

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
