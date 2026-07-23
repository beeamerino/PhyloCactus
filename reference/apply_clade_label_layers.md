# Apply Multiple Taxonomic Clade Label Layers onto ggtree Plot

Sequentially overlays hierarchical taxonomic clade annotations (e.g.,
subfamilies, tribes, genera) onto a \`ggtree\` plot.

## Usage

``` r
apply_clade_label_layers(p, registry, layer_specs)
```

## Arguments

- p:

  A \`ggtree\` plot object.

- registry:

  Named list of node annotation data frames.

- layer_specs:

  Data frame defining visual style specifications per annotation layer.

## Value

Updated \`ggtree\` plot object.
