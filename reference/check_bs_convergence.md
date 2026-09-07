# Check Bootstrap Convergence Criterion in RAxML-NG

Evaluates whether the generated pool of non-parametric bootstrap trees
has achieved statistical convergence (`RAxML-NG`), using the autoMRE
bootstopping criterion at a permutation cutoff (typically 0.03).

## Usage

``` r
check_bs_convergence(
  raxml_bin_path,
  bs_trees_file,
  bs_cutoff = 0.03,
  bs_metric = c("fbp", "tbe"),
  seed = NULL,
  threads = 4,
  output_dir = dirname(bs_trees_file),
  prefix = "cactus_bs_convergence"
)
```

## Arguments

- raxml_bin_path:

  Character. System command or full path to executable `RAxML-NG`
  binary.

- bs_trees_file:

  Character. Path to concatenated bootstrap tree file (e.g.,
  `cactus_ALL_bootstraps.tree`).

- bs_cutoff:

  Numeric. Permutation cutoff threshold for convergence. Defaults to
  `0.03`.

- bs_metric:

  Character. Branch support metric for the bootstopping test: `"fbp"`
  (Felsenstein's Bootstrap Percentage, the default) or `"tbe"` (Transfer
  Bootstrap Expectation). Must match the metric being reported.

- seed:

  Integer. Random seed for reproducible convergence testing. Defaults to
  `NULL` (random).

- threads:

  Integer. Number of CPU threads. Defaults to `4`.

- output_dir:

  Character. Directory path to save convergence report logs. Defaults to
  `dirname(bs_trees_file)`.

- prefix:

  Character. Output file prefix. Defaults to `"cactus_bs_convergence"`.

## Value

A structured S3 object of class `cactus_bs_convergence` containing:

- `converged`: Logical indicating whether the autoMRE convergence
  criterion was satisfied.

- `trees_at_convergence`: Integer number of trees required to reach
  convergence (or `NA` if not converged).

- `trees_analyzed`: Total number of bootstrap trees evaluated.

- `cutoff`: Numeric convergence cutoff threshold.

- `metric`: Support metric used (`"fbp"` or `"tbe"`).

- `log_file`: Character path to the resulting RAxML-NG convergence log
  file (`.raxml.log`).

## Details

The criterion is metric-specific and the two metrics do not converge at
the same replicate count. TBE is the more stable summary and reaches the
cutoff earlier, so a pool declared converged under TBE is not thereby
converged under FBP. The default is FBP, matching
[`map_branch_supports()`](https://beeamerino.github.io/PhyloCactus/reference/map_branch_supports.md),
so that the convergence statement and the reported support values refer
to the same quantity. Run the test under whichever metric is being
reported, and under both when both are tabulated.
