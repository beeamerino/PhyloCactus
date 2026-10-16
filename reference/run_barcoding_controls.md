# Negative Controls of the Molecular Diagnostic Branch

Step 8 of the branch, and the gate in front of every real figure. Three
controls, none of which is a result of the method: the three are the
condition for a result of the method to be interpretable.

## Usage

``` r
run_barcoding_controls(
  library_dir = file.path("11_barcoding", "4_library"),
  folds_dir = file.path("11_barcoding", "5_folds"),
  output_dir = file.path("11_barcoding", "8_controls"),
  outgroup_dir = NULL,
  method = c("nn", "idtaxa"),
  permutations = 10L,
  seed = 1L,
  model = "raw",
  min_comparable = 100L,
  schemes = c("species", "genus"),
  loci = NULL,
  max_folds = NULL,
  controls = c("CN1", "CN2", "CN3"),
  mafft_exec = "mafft",
  mafft_opts = "--auto",
  threshold = 60
)
```

## Arguments

- library_dir:

  Character. Output directory of
  [`finalize_barcoding_library()`](https://beeamerino.github.io/PhyloCactus/reference/finalize_barcoding_library.md).

- folds_dir:

  Character. Output directory of
  [`build_barcoding_folds()`](https://beeamerino.github.io/PhyloCactus/reference/build_barcoding_folds.md).

- output_dir:

  Character. Destination; refused inside a directory of the phylogeny.

- outgroup_dir:

  Character or `NULL`. Directory with `<locus>.fasta` of the outgroup,
  as the assembly writes them. `NULL` skips the outgroup control and
  says so.

- method:

  Character. `"nn"` or `"idtaxa"`. With `"idtaxa"` only the outgroup
  control ("CN2") runs: the label permutation would retrain IdTaxa once
  per fold and permutation, and the resubstitution measures a bias of
  the nearest neighbour. The outgroup queries are oriented against the
  library with the rule of the nearest neighbour, classified by IdTaxa
  trained on the whole locus at threshold 0 with a seed per query, and
  written to `TABLE_barcoding_cn2_queries_idtaxa.csv`; the tables of the
  nearest neighbour are not touched. MAFFT is not needed.

- permutations:

  Integer. Permutations of the label permutation control; at least 10 in
  a real run.

- seed:

  Integer. Seed of the permutations, declared so they are reproducible.

- model, min_comparable:

  Distance model and minimum of comparable positions, as in step 7.

- schemes, loci:

  Schemes and loci; `NULL` for every locus of the library.

- max_folds:

  Integer or `NULL`. Declared subset of folds, for the expensive method.

- controls:

  Character vector. Which of `"CN1"`, `"CN2"` and `"CN3"` to run.

- mafft_exec:

  Character. Command or path of the `MAFFT` binary. The outgroup control
  adds every outgroup query to the alignment of the library with
  `--add --keeplength`, one call per query, so the columns of the
  library stay fixed and the distance of one query is comparable with
  the distance of another. Without `outgroup_dir` MAFFT is never called.

- mafft_opts:

  Character. Further options for `MAFFT`. `--add` and `--keeplength` are
  added by the function and are not optional.

- threshold:

  Numeric. Confidence threshold applied to IdTaxa's output, for
  `"idtaxa"`, as in
  [`classify_barcoding_folds()`](https://beeamerino.github.io/PhyloCactus/reference/classify_barcoding_folds.md).

## Value

Invisibly, a list with the tables. Writes
`TABLE_barcoding_cn1_permutations.csv`,
`TABLE_barcoding_cn1_verdict.csv`, `TABLE_barcoding_cn2_outgroup.csv`
and `TABLE_barcoding_cn3_resubstitution.csv`.

## Details

**Label permutation ("CN1").** The whole procedure is repeated with the
species labels permuted inside the training set. The accuracy must fall
to the base rate; if it does not, there is leakage. Written as one row
per permutation plus its verdict.

**Outgroup queries ("CN2"), alien to the clade.** Sequences of the
outgroup, which are not in the library, are passed as queries. The table
publishes how many outgroup **species** each locus had, and how many
sequences, zeros included. The size of the control is the species count:
GenBank holds population level studies, so a locus can carry dozens of
sequences of a single species and control almost nothing. A control over
one species is an anecdote and has to be read as one.

**Resubstitution ("CN3").** The same queries as the folds, in both
schemes, each one classified against its own training set with the query
put back into it. Nothing else changes between this figure and the
honest one, which is what makes the two comparable: until 2026-09-25
this control ran over the whole library instead, a different set of
queries, and the two figures could not be subtracted. Written with the
label `irreproducible_resubstitution` in the table itself, so the figure
cannot travel without it. **Its value is published only next to the
honest figure**, which does not exist yet.

## Examples

``` r
if (FALSE) { # \dontrun{
run_barcoding_controls(
  library_dir = "11_barcoding/4_library",
  folds_dir = "11_barcoding/5_folds",
  output_dir = "11_barcoding/8_controls",
  outgroup_dir = "11_barcoding/8_controls/outgroup/1_assembly"
)
} # }
```
