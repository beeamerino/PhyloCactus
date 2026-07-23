# Collect and Concatenate Parallel Bootstrap Tree Outputs

Scans a target directory for chunked bootstrap output files
(\`.raxml.bootstraps\`) and concatenates them into a single Newick tree
file.

## Usage

``` r
collect_bootstraps(
  bs_dir,
  output_dir = bs_dir,
  prefix = "cactus_ALL_bootstraps"
)
```

## Arguments

- bs_dir:

  Character. Directory path containing chunked bootstrap output files.

- output_dir:

  Character. Output directory path to save concatenated bootstrap file.
  Defaults to \`bs_dir\`.

- prefix:

  Character. File output prefix. Defaults to
  \`"cactus_ALL_bootstraps"\`.

## Value

Character path to the concatenated bootstrap tree file (\`.tree\`).
