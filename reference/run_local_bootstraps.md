# Generate Non-Parametric Bootstrap Trees Locally

Performs non-parametric bootstrap resampling over supermatrix site
columns to infer a distribution of bootstrap tree topologies
(\`RAxML-NG\`). Evaluates topological variation under non-parametric
resampling to quantify node support via Transfer Bootstrap Expectation
(TBE).

## Usage

``` r
run_local_bootstraps(
  raxml_bin_path,
  aln_file,
  part_file,
  constraint_file,
  bs_trees = 500,
  outgroup = NULL,
  seed = 111,
  threads = 8,
  workers = 1,
  output_dir = dirname(aln_file),
  prefix = "cactus_bs"
)
```

## Arguments

- raxml_bin_path:

  Character. System command or full path to executable \`RAxML-NG\`
  binary.

- aln_file:

  Character. Path to input PHYLIP supermatrix alignment file.

- part_file:

  Character. Path to partition file specifying substitution models.

- constraint_file:

  Character. Path to Newick topological constraint scaffold file.

- bs_trees:

  Integer. Total number of non-parametric bootstrap trees to generate.
  Defaults to \`500\`.

- outgroup:

  Character. Optional outgroup taxon binomial to root bootstrap
  topologies. Defaults to \`NULL\`.

- seed:

  Integer. Random seed for reproducible bootstrap initialization.
  Defaults to \`111\`.

- threads:

  Integer. Number of CPU threads. Defaults to \`8\`.

- workers:

  Integer. Parallel worker process count. Defaults to \`1\`.

- output_dir:

  Character. Output directory for generated bootstrap trees. Defaults to
  \`dirname(aln_file)\`.

- prefix:

  Character. Output file prefix. Defaults to \`"cactus_bs"\`.

## Value

Character path to the output bootstrap trees file
(\`.raxml.bootstraps\`).

## References

Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A.
(2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum
likelihood phylogenetic inference. \*Bioinformatics\*, 35(21),
4453-4455.
[doi:10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)
