# Concatenate Locus Alignments and Build Partition Coordinate Maps

Concatenates individual orthologous locus alignments end-to-end into a
unified multilocus supermatrix. Combining independent genomic loci
increases statistical power to resolve difficult ancestral nodes while
allowing partitioned substitution modeling to account for mutational
rate heterogeneity across genomic regions. Exports partition files
compatible with \`RAxML-NG\`, \`IQ-TREE\`, \`PartitionFinder2\`, and
\`MrBayes\`.

## Usage

``` r
run_concatenation_pipeline(input_dir, output_dir)
```

## Arguments

- input_dir:

  Character. Path to directory containing curated, aligned locus FASTA
  files.

- output_dir:

  Character. Path to destination root directory for concatenated
  alignments and partition maps.

## Value

A data frame containing supermatrix dimensions, taxon coverage, and
genomic partition bounds.

## Examples

``` r
if (FALSE) { # \dontrun{
run_concatenation_pipeline(
  input_dir = "5_MAFFT_Cleaned/aligned_markers",
  output_dir = "6_Concatenated"
)
} # }
```
