# Taxonomic IDs of the Outgroup of the Molecular Diagnostic Branch

The six genus-level NCBI Taxonomy IDs that the outgroup control queries
the library with: the three families that surround Cactaceae.

## Usage

``` r
barcoding_outgroup_taxids()
```

## Value

Character vector of six taxonomic IDs.

## Details

- **Anacampserotaceae:** *Talinopsis* (107598), *Grahamia* (107617),
  *Anacampseros* (107583).

- **Portulacaceae:** *Portulaca* (3582).

- **Talinaceae:** *Talinum* (107600), *Talinella* (108056).

*Talinopsis* is Anacampserotaceae, not Talinaceae, despite the name and
the adjacent id. *Amphipetalum* (1835425) is not requested: GenBank
holds no nucleotide records for it, so the request would only produce an
empty download.

**This is the branch's own copy**. The phylogeny keeps its list in
[`assemble_outgroup_phylotar()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_outgroup_phylotar.md)
and is not touched, so the branch runs without it. A test compares the
two and fails if they ever drift apart.

## Examples

``` r
barcoding_outgroup_taxids()
#> [1] "107598" "107617" "107583" "3582"   "107600" "108056"
```
