# Resolve Every File Path of a Phylogenetic Run from Its Naming Convention

Reconstructs the absolute path of every file a `PhyloCactus` inference
run produces, given only the output directory, the run prefix and the
supermatrix location. Because each path is derived rather than carried
in memory, any stage of the pipeline can be executed in a fresh R
session without first re-running the stages before it.

## Usage

``` r
resolve_run_paths(
  output_dir = "7_Phylogenetics",
  prefix = "cactus",
  supermatrix = file.path("6_Concatenated", "concatenated_alignments",
    "ALIGNMENT_supermatrix.phy"),
  require = NULL
)

# S3 method for class 'cactus_run_paths'
print(x, ...)
```

## Arguments

- output_dir:

  Character. Run directory holding the inference outputs. Defaults to
  `"7_Phylogenetics"`.

- prefix:

  Character. Run prefix shared by every file of the run. Defaults to
  `"cactus"`.

- supermatrix:

  Character. Path to the concatenated supermatrix exported by
  [`run_concatenation_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_concatenation_pipeline.md).
  Defaults to the Module 6 location.

- require:

  Character vector. Names of entries that must already exist on disk.
  Any that do not trigger an error naming the module responsible for
  writing them. Defaults to `NULL`.

- x:

  A `cactus_run_paths` object.

- ...:

  Ignored.

## Value

An object of class `cactus_run_paths`: a named list of absolute paths,
plus `was_reduced`, `ml_on_cluster`, and an `exists` list of logicals
reporting which files are present.

## Details

The pipeline is a sequence of stages whose outputs feed the next, and a
linear script holds those outputs in variables such as `analysed_phy` or
`ml_results`. Restarting R and resuming at a later stage clears those
variables, and the stage fails on a missing object rather than on a
missing file. Hard-coding the filenames instead removes the dependency
on session state but reintroduces the problem this naming convention
exists to solve: the same literal name is retyped at a dozen call sites,
drifts out of step with the run prefix, and continues to be read after
it has stopped describing its contents.

Deriving the paths from the convention keeps a single definition of
every filename while leaving each stage independently runnable. Two
entries are resolved by inspection rather than by convention alone,
because they record facts about the run rather than choices about
naming:

- `analysed_phy` is the reduced PHYLIP when `RAxML-NG` collapsed
  identical terminals during validation, and the supermatrix itself when
  it did not. `was_reduced` reports which occurred.

- `best_tree` and `ml_trees` resolve to the `ml_search/` subdirectory
  when the maximum likelihood search was submitted to a cluster, and to
  the run root when it was executed locally.

Passing `require` converts a missing input into an error naming the
stage that produces it, instead of letting the failure surface several
calls later as an unreadable file.

## See also

[`preprocess_partitions()`](https://beeamerino.github.io/PhyloCactus/reference/preprocess_partitions.md),
[`calculate_ml_tree()`](https://beeamerino.github.io/PhyloCactus/reference/calculate_ml_tree.md),
[`collect_bootstraps()`](https://beeamerino.github.io/PhyloCactus/reference/collect_bootstraps.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Start of any stage, in a fresh R session
paths <- resolve_run_paths(output_dir = "7_Phylogenetics", prefix = "cactus")
paths

# Refuse to start Module 9 unless Module 8 actually finished
paths <- resolve_run_paths(require = c("best_tree", "constraint_tree"))
} # }
```
