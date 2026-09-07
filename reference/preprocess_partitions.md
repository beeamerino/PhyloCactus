# Preprocess Partitions and Validate Alignment Syntax

Validates PHYLIP alignment syntax and partition file coordinates using
`RAxML-NG` (Kozlov *et al.*, 2019). Verifies site ranges, formatting
compatibility, and data integrity prior to substitution model
evaluation.

## Usage

``` r
preprocess_partitions(
  phy_matrix,
  part_file,
  raxml_path,
  output_dir = dirname(phy_matrix),
  prefix = "cactus",
  force_check = FALSE,
  model_handling = c("force_dna", "preserve")
)
```

## Arguments

- phy_matrix:

  Character. Path to input PHYLIP supermatrix file.

- part_file:

  Character. Path to input partition mapping text file in RAxML-style
  format (`PARTITION_raxml_style.txt`, exported by
  [`run_concatenation_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_concatenation_pipeline.md)).

- raxml_path:

  Character. System command or full path to executable `RAxML-NG`
  binary.

- output_dir:

  Character. Output directory for validated partition outputs. Defaults
  to `dirname(phy_matrix)`.

- prefix:

  Character. Basename prefix for every file written by this step.
  Defaults to `"cactus"`.

- force_check:

  Logical. Re-run validation even if the validated partition file
  already exists? Defaults to `FALSE`.

- model_handling:

  Character. Model field written to the validated partition map.
  `"force_dna"` (default) writes the datatype token `DNA`, leaving the
  substitution model to be selected by
  [`run_modeltest_ng()`](https://beeamerino.github.io/PhyloCactus/reference/run_modeltest_ng.md).
  `"preserve"` keeps whatever model string `RAxML-NG` emitted, which is
  only appropriate when the partition map is passed straight to
  `RAxML-NG` without model selection. `ModelTest-NG` cannot parse the
  `RAxML-NG` model syntax and aborts on it.

## Value

Character path to the validated partition map file ready for model
evaluation.

## References

Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A.
(2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum
likelihood phylogenetic inference. *Bioinformatics*, 35(21), 4453-4455.
[doi:10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

## See also

[`run_concatenation_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_concatenation_pipeline.md)
for the partition map this function consumes, and
[`run_modeltest_ng()`](https://beeamerino.github.io/PhyloCactus/reference/run_modeltest_ng.md)
for the model selection step that consumes its output.
