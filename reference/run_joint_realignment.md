# Perform Joint Realignment Across Integrated Ingroup and Outgroup Sequences

Re-estimates positional homology alignments (\`MAFFT\`) after
aggregating ingroup and outgroup sequence markers. Joint realignment
accommodates positional variations introduced when combining outgroup
taxons with the focal ingroup clade.

## Usage

``` r
run_joint_realignment(input_dir, output_fasta_dir, output_aln_dir)
```

## Arguments

- input_dir:

  Character. Directory containing merged, taxonomically reconciled locus
  FASTA files.

- output_fasta_dir:

  Character. Directory path to save output realigned FASTA sequence
  files.

- output_aln_dir:

  Character. Directory path to store detailed realignment log reports.

## Value

Invisible NULL upon successful execution.

## References

Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment
software version 7: Improvements in performance and usability.
\*Molecular Biology and Evolution\*, 30(4), 772–780.
[doi:10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)

## Examples

``` r
if (FALSE) { # \dontrun{
run_joint_realignment(
  input_dir = "4_Cleaned/cleaned_markers",
  output_fasta_dir = "5_MAFFT_Cleaned",
  output_aln_dir = "5_MAFFT_Cleaned/aligned_markers"
)
} # }
```
