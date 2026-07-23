# Preprocess Partitions and Validate Alignment Syntax

Validates PHYLIP alignment syntax and partition file coordinates using
\`RAxML-NG\` (Kozlov \*et al.\*, 2019). Verifies site ranges, formatting
compatibility, and data integrity prior to substitution model
evaluation.

## Usage

``` r
preprocess_partitions(
  phy_matrix,
  part_file,
  raxml_path,
  output_dir = dirname(phy_matrix),
  force_check = FALSE
)
```

## Arguments

- phy_matrix:

  Character. Path to input PHYLIP supermatrix file.

- part_file:

  Character. Path to input partition mapping text file.

- raxml_path:

  Character. System command or full path to executable \`RAxML-NG\`
  binary.

- output_dir:

  Character. Output directory for validated partition outputs. Defaults
  to \`dirname(phy_matrix)\`.

- force_check:

  Logical. Bypass validation check if cleaned partition output file
  already exists? Defaults to \`FALSE\`.

## Value

Character path to the cleaned partition map file ready for model
evaluation.

## References

Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A.
(2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum
likelihood phylogenetic inference. \*Bioinformatics\*, 35(21),
4453-4455.
[doi:10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)
