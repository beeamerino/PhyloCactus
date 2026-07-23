# Check Bootstrap Convergence Criterion in RAxML-NG

Evaluates whether the generated pool of non-parametric bootstrap trees
has achieved statistical convergence (\`RAxML-NG\`). Convergence is
assessed using the MRE-based cutoff criterion (typically \<= 0.03),
ensuring that a sufficient number of bootstrap replicates were sampled.

## Usage

``` r
check_bs_convergence(
  raxml_bin_path,
  bs_trees_file,
  bs_cutoff = 0.03,
  seed = 111,
  threads = 4,
  output_dir = dirname(bs_trees_file),
  prefix = "cactus_bs_convergence"
)
```

## Arguments

- raxml_bin_path:

  Character. System command or full path to executable \`RAxML-NG\`
  binary.

- bs_trees_file:

  Character. Path to concatenated bootstrap tree file.

- bs_cutoff:

  Numeric. Permutation cutoff threshold for convergence. Defaults to
  \`0.03\`.

- seed:

  Integer. Random seed for reproducible convergence testing. Defaults to
  \`111\`.

- threads:

  Integer. Number of CPU threads. Defaults to \`4\`.

- output_dir:

  Character. Directory path to save convergence report logs. Defaults to
  \`dirname(bs_trees_file)\`.

- prefix:

  Character. Output file prefix. Defaults to
  \`"cactus_bs_convergence"\`.

## Value

Character path to the convergence log file (\`.raxml.log\`).
