# Rescale Branch Lengths of a Phylogenetic Tree

Multiplies all edge lengths of a \`phylo\` tree object by a constant
scaling factor. Rescaling is used prior to penalized likelihood dating
(\`treePL\`) to prevent numerical underflow issues caused by small
substitution rates per site.

## Usage

``` r
rescale_tree(tree, factor = 100)
```

## Arguments

- tree:

  An object of class \`phylo\` representing a phylogenetic tree.

- factor:

  Numeric multiplier applied to all edge lengths. Defaults to \`100\`.

## Value

A rescaled \`phylo\` object with updated edge lengths.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ape)
tree <- rtree(10)
rescaled_tree <- rescale_tree(tree, factor = 100)
} # }
```
