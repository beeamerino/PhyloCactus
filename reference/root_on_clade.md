# Root a Tree on a Clade of Terminals, Tolerating a Basal Polytomy

[`ape::root()`](https://rdrr.io/pkg/ape/man/root.html) decides monophyly
against the tree's *current* root. `RAxML-NG` writes unrooted Newick
with a basal trifurcation, and the terminals of the rooting clade
routinely fall on more than one branch of that trifurcation, so
[`ape::root()`](https://rdrr.io/pkg/ape/man/root.html) does not see them
as a clade and does not place the root where it was asked to. The
reference `cactus_support.raxml.support` is exactly this case:
`{Portulaca | rest}` is a valid bipartition of the unrooted topology,
yet three *Portulaca* terminals sit outside the largest basal branch.

## Usage

``` r
root_on_clade(phy, outgroup)
```

## Arguments

- phy:

  An object of class `phylo`.

- outgroup:

  Character vector of rooting terminals, typically from
  [`resolve_rooting_outgroup()`](https://beeamerino.github.io/PhyloCactus/reference/resolve_rooting_outgroup.md).

## Value

The rooted `phylo`. Errors, via
[`ape::root()`](https://rdrr.io/pkg/ape/man/root.html), when the
terminals are not a clade of the unrooted topology.

## Details

Rooting first on an ingroup terminal moves the current root into the
ingroup, which makes the outgroup a clade in the rooted sense whenever
`{outgroup | ingroup}` is a bipartition of the unrooted topology. The
second call then places the root on the intended stem edge.

The first step is unconditional.
[`ape::is.monophyletic()`](https://rdrr.io/pkg/ape/man/is.monophyletic.html)
cannot guard it, because on an unrooted tree it answers in the
bipartition sense and reports the terminals as monophyletic while
[`ape::root()`](https://rdrr.io/pkg/ape/man/root.html), which looks for
them as a clade of the tree as currently rooted, does not find them;
[`ape::root()`](https://rdrr.io/pkg/ape/man/root.html) then falls back
to their MRCA, which is the existing root, and the call silently does
nothing. That silent no-op is the failure this function exists to
remove.

Rooting cannot be delegated upstream to `RAxML-NG`. Under
time-reversible substitution models the likelihood is identical for
every rooting of the same unrooted topology (Felsenstein 1981), so
neither a constraint tree nor `--outgroup` can select a root: the data
carry no information about its position. The root is an outgroup
decision imposed after the search.

## Examples

``` r
tr <- ape::read.tree(text = "((A:1,B:1):1,Portulaca_fulgens:1,Portulaca_oleracea:1);")
rooted <- root_on_clade(tr, c("Portulaca_fulgens", "Portulaca_oleracea"))
ape::is.monophyletic(rooted, c("Portulaca_fulgens", "Portulaca_oleracea"))
#> [1] TRUE
```
