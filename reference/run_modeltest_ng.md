# Evaluate Nucleotide Substitution Models via ModelTest-NG

Evaluates nucleotide substitution model fit per predefined supermatrix
partition using `ModelTest-NG` (Darriba *et al.*, 2020). The model of
each partition is selected under the AICc, so substitution models can
differ among loci.

## Usage

``` r
run_modeltest_ng(
  modeltest_exec_path,
  aln_file,
  part_file,
  prefix = "MODELTEST_cactus_phylo",
  threads = 4,
  output_dir = NULL,
  seed = NULL
)
```

## Arguments

- modeltest_exec_path:

  Character. System command or full path to executable `ModelTest-NG`
  binary.

- aln_file:

  Character. Path to validated PHYLIP supermatrix file.

- part_file:

  Character. Path to the validated partition map written by
  [`preprocess_partitions()`](https://beeamerino.github.io/PhyloCactus/reference/preprocess_partitions.md).
  Its model field must carry the datatype token `DNA`; `ModelTest-NG`
  cannot parse `RAxML-NG` model syntax such as `GTR+FC+G4m+B` and aborts
  with a segmentation fault on it.

- prefix:

  Character. Output filename prefix. Defaults to
  `"MODELTEST_cactus_phylo"`.

- threads:

  Integer. Number of processing threads. Defaults to `4`.

- output_dir:

  Character or `NULL`. Directory for the outputs of `ModelTest-NG`,
  created if needed. `NULL` (default) writes them, as before 0.5.5,
  relative to the working directory.

- seed:

  Integer or `NULL`. Seed of the random number generator of
  `ModelTest-NG` (`--rngseed`). `NULL` (default) leaves it to the
  program, as before 0.5.5; give a value to make the run reproducible
  and record it.

## Value

Character path to the resulting partition file containing selected model
parameters (`.part.aicc`).

## References

Darriba, D., Posada, D., Kozlov, A. M., Stamatakis, A., Morel, B., &
Flouri, T. (2020). ModelTest-NG: a new and scalable tool for the
selection of DNA and protein evolutionary models. *Molecular Biology and
Evolution*, 37(1), 291-294.
[doi:10.1093/molbev/msz189](https://doi.org/10.1093/molbev/msz189)

## Examples

``` r
if (FALSE) { # \dontrun{
run_modeltest_ng(
  modeltest_exec_path = "modeltest-ng",
  aln_file = "ALIGNMENT_supermatrix.phy",
  part_file = "cactus_partitions_validated.txt"
)
} # }
```
