# Synthesize Multifurcating Monophyly Constraint Scaffold

Constructs a Newick multifurcating constraint tree enforcing monophyly
of established higher taxonomic ranks (e.g., subfamilies, tribes).
Constrained maximum-likelihood searches restrict branch topology
exploration to scientifically verified monophyletic backbone clades,
preventing aberrant tree topologies when analyzing sparse supermatrices.

## Usage

``` r
build_constraint_scaffold(
  alignment_path,
  constraints_csv_path,
  output_dir = dirname(alignment_path)
)
```

## Arguments

- alignment_path:

  Character. Path to input PHYLIP supermatrix alignment file.

- constraints_csv_path:

  Character. Path to taxonomy CSV table mapping species binomials to
  taxonomic ranks.

- output_dir:

  Character. Directory path to save generated constraint scaffold file.
  Defaults to \`dirname(alignment_path)\`.

## Value

Character string path to the saved Newick constraint tree file
(\`cactus_constraints.tree\`).

## Examples

``` r
if (FALSE) { # \dontrun{
build_constraint_scaffold(
  alignment_path = "ALIGNMENT_supermatrix.phy",
  constraints_csv_path = "Cactaceae_taxonomy.csv"
)
} # }
```
