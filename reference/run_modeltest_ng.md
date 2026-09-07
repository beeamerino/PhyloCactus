# Evaluate Nucleotide Substitution Models via ModelTest-NG

Evaluates nucleotide substitution model fit per predefined supermatrix
partition using `ModelTest-NG` (Darriba *et al.*, 2020). Selecting
optimal substitution models under the AICc criterion controls for
mutational rate heterogeneity across molecular locus alignments,
mitigating systematic long-branch attraction (LBA) bias during
maximum-likelihood inference.

## Usage

``` r
run_modeltest_ng(
  modeltest_exec_path,
  aln_file,
  part_file,
  prefix = "MODELTEST_cactus_phylo",
  threads = 4
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
