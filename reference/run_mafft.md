# Infer Positional Homology via MAFFT Alignment

Establishes hypotheses of positional homology across unaligned
orthologous nucleotide sequence clusters. Positional homology alignment
is a crucial prerequisite for maximum-likelihood phylogenetic inference,
ensuring that corresponding nucleotide sites derived from common
evolutionary ancestry are aligned prior to substitution model
evaluation.

## Usage

``` r
run_mafft(
  input_fasta,
  output_fasta,
  mafft_exec = "mafft",
  mafft_opts = "--auto"
)
```

## Arguments

- input_fasta:

  Character. Path to unaligned input FASTA file.

- output_fasta:

  Character. Path to destination aligned FASTA output file.

- mafft_exec:

  Character. System command or full path to the executable \`MAFFT\`
  binary. Defaults to \`"mafft"\`.

- mafft_opts:

  Character. Command-line parameters passed directly to \`MAFFT\`.
  Defaults to \`"–auto"\`.

## Value

Invisible numeric exit status code (0 for successful alignment
completion).

## References

Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment
software version 7: Improvements in performance and usability.
\*Molecular Biology and Evolution\*, 30(4), 772–780.
[doi:10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)

## Examples

``` r
if (FALSE) { # \dontrun{
run_mafft(
  input_fasta = "raw_cluster.fasta",
  output_fasta = "aligned_cluster.fasta",
  mafft_opts = "--auto"
)
} # }
```
