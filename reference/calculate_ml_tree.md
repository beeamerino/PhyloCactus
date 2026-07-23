# Infer Maximum-Likelihood Phylogeny under Constrained Search

Infers the maximum-likelihood evolutionary hypothesis explaining the
concatenated supermatrix under specified partition models and
topological constraint scaffolds using \`RAxML-NG\` (Kozlov \*et al.\*,
2019). Executes multiple independent tree searches starting from
randomized and parsimony starting trees to avoid local likelihood
Optima.

## Usage

``` r
calculate_ml_tree(
  raxml_bin_path,
  aln_file,
  part_file,
  constraint_file,
  outgroup = NULL,
  n_init_trees = "rand{25},pars{25}",
  seed = 111,
  n_workers = 1,
  threads = 4,
  output_dir = dirname(aln_file),
  prefix = "cactus_search"
)
```

## Arguments

- raxml_bin_path:

  Character. System command or full path to executable \`RAxML-NG\`
  binary.

- aln_file:

  Character. Path to input PHYLIP supermatrix alignment file.

- part_file:

  Character. Path to partition file specifying substitution models per
  partition.

- constraint_file:

  Character. Path to Newick topological constraint scaffold file.

- outgroup:

  Character. Optional outgroup taxon binomial to root the resulting tree
  topology. Defaults to \`NULL\`.

- n_init_trees:

  Character. Initial starting tree specifications. Defaults to
  \`"rand25,pars25"\` (25 random + 25 parsimony trees).

- seed:

  Integer. Random seed for reproducible tree search initialization.
  Defaults to \`111\`.

- n_workers:

  Integer. Parallel worker process count for RAxML-NG. Defaults to
  \`1\`.

- threads:

  Integer. Number of CPU threads per worker. Defaults to \`4\`.

- output_dir:

  Character. Directory path to save resulting maximum-likelihood tree
  files. Defaults to \`dirname(aln_file)\`.

- prefix:

  Character. Output file prefix. Defaults to \`"cactus_search"\`.

## Value

A named list containing paths to the best ML tree (\`.raxml.bestTree\`)
and all evaluated trees (\`.raxml.mlTrees\`).

## References

Kozlov, A. M., Darriba, D., Flouri, T., Morel, B., & Stamatakis, A.
(2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum
likelihood phylogenetic inference. \*Bioinformatics\*, 35(21),
4453-4455.
[doi:10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

## Examples

``` r
if (FALSE) { # \dontrun{
calculate_ml_tree(
  raxml_bin_path = "raxml-ng",
  aln_file = "ALIGNMENT_supermatrix.phy",
  part_file = "MODELTEST_cactus_phylo.part.aicc",
  constraint_file = "cactus_constraints.tree"
)
} # }
```
