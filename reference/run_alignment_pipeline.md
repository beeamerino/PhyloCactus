# Execute Complete Alignment and Gap-Masking Pipeline

Orchestrates multiple sequence alignment (MSA) and automated
quality-control masking across orthologous sequence clusters. Primary
alignment hypotheses are inferred using \`MAFFT\` (Katoh & Standley,
2013). Subsequently, ambiguous sites, poorly aligned terminal fragments,
and non-homologous insertions are masked using the \`DECIPHER\`
framework (Wright, 2024), eliminating systematic noise while retaining
phylogenetically informative nucleotide positions for downstream
supermatrix assembly.

## Usage

``` r
run_alignment_pipeline(
  input_folder,
  output_dir,
  fasta_pattern = "\\.fasta$",
  mask_alignment_regions = TRUE,
  min_non_gap_fraction = 0.3,
  max_missing_fraction = 0.3,
  mafft_exec = "mafft",
  mafft_opts = "--auto"
)
```

## Arguments

- input_folder:

  Character. Directory path containing raw unaligned orthologous FASTA
  files.

- output_dir:

  Character. Root destination directory for output subfolders
  (\`alignments/\`, \`tables/\`, \`logs/\`).

- fasta_pattern:

  Character. Regular expression pattern matching target FASTA files.
  Defaults to \`"\\fasta\$"\`.

- mask_alignment_regions:

  Logical. Apply automated alignment masking via \`DECIPHER\`? Defaults
  to \`TRUE\`.

- min_non_gap_fraction:

  Numeric. Minimum allowable proportion of non-gap characters required
  to retain a site column. Defaults to \`0.30\`.

- max_missing_fraction:

  Numeric. Maximum allowable proportion of missing or ambiguous
  characters (\`N\`) allowed per sequence. Defaults to \`0.30\`.

- mafft_exec:

  Character. System command or full path to the executable \`MAFFT\`
  binary. Defaults to \`"mafft"\`.

- mafft_opts:

  Character. Command-line parameters passed directly to \`MAFFT\`.
  Defaults to \`"–auto"\`.

## Value

A data frame containing site length, missingness, and sequence retention
statistics across processed loci.

## References

Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment
software version 7: Improvements in performance and usability.
\*Molecular Biology and Evolution\*, 30(4), 772–780.
[doi:10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)

Wright, E. S. (2024). Fast and Flexible Search for Homologous Biological
Sequences with DECIPHER v3. \*The R Journal\*, 16(2), 191-200.
[doi:10.18129/B9.bioc.DECIPHER](https://doi.org/10.18129/B9.bioc.DECIPHER)

## Examples

``` r
if (FALSE) { # \dontrun{
run_alignment_pipeline(
  input_folder = "1_phylotaR_out_ingroup",
  output_dir = "2_MAFFT_Cactaceae",
  min_non_gap_fraction = 0.30,
  max_missing_fraction = 0.30
)
} # }
```
