# Rescale Branch Lengths of a Phylogenetic Tree

Multiplies all edge lengths of a `phylo` tree object by a constant
scaling factor.

## Usage

``` r
rescale_tree(tree, factor = 100)
```

## Arguments

- tree:

  An object of class `phylo` representing a phylogenetic tree.

- factor:

  Numeric multiplier applied to all edge lengths. Defaults to `100`.

## Value

A rescaled `phylo` object with updated edge lengths.

## Details

The factor does not change which branches `treePL` clamps. `treePL`
refuses to let a branch carry less than one expected substitution and
rewrites any shorter branch to `1/numsites`, so the clamp acts on the
substitution count, not on the raw length. Measured on the Cactaceae
supermatrix, 492 of 2088 branches are clamped, identically at factor 100
and at factor 1, once `numsites` is divided by the same factor as
[`automate_treePL()`](https://beeamerino.github.io/PhyloCactus/reference/automate_treePL.md)
does. What the factor does change is the scale of the rate parameters
and therefore the numerical conditioning of the optimisation.

Earlier versions of this documentation described the rescaling as a
guard against numerical underflow. That is not what it does.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ape)
tree <- rtree(10)
rescaled_tree <- rescale_tree(tree, factor = 100)
} # }
```
