# Resolve the Set of Rooting Terminals from a Tip Label Vector

Selects every terminal belonging to the outgroup lineage that carries
the root, returning them as a vector suitable for
[`ape::root()`](https://rdrr.io/pkg/ape/man/root.html) and for the
`outgroup` argument of the inference and dating functions.

## Usage

``` r
resolve_rooting_outgroup(tip_labels, pattern = "^Portulaca_")
```

## Arguments

- tip_labels:

  Character vector of terminal labels, typically `phylo$tip.label` or
  the row names of the supermatrix.

- pattern:

  Character. Regular expression matched against `tip_labels`. Defaults
  to `"^Portulaca_"`, the Portulacaceae sample of the reference
  Cactaceae dataset.

## Value

Character vector of matching terminals, sorted. Errors when no terminal
matches, since a silently empty rooting set would leave every downstream
tree unrooted.

## Details

Rooting is a property of the rooted tree, not of the search. `RAxML-NG`
returns an unrooted topology, and the root is imposed afterwards by
choosing the edge it sits on. Choosing a single terminal of a sampled
outgroup clade places that root *inside* the clade, which renders the
clade paraphyletic in the final tree and collapses its crown node onto
the root. Any calibration addressed by the MRCA of that clade then
silently lands on the root instead. Passing the whole clade places the
root on its stem edge and avoids both consequences.

The trailing underscore in the default pattern is required rather than
stylistic. Terminals use the underscore as binomial separator, so a bare
`"^Portulaca"` also matches *Portulacaria* (Didiereaceae), which is
neither Portulacaceae nor part of the intended rooting sample. The same
reasoning is applied in
[`run_marker_screening()`](https://beeamerino.github.io/PhyloCactus/reference/run_marker_screening.md)
and
[`run_concatenation_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_concatenation_pipeline.md).

## Examples

``` r
resolve_rooting_outgroup(c("Portulaca_oleracea", "Portulacaria_afra", "Opuntia_ficus-indica"))
#> [1] "Portulaca_oleracea"
```
