# Build Taxonomic Annotation Registry across Tree Levels

Constructs a comprehensive registry mapping internal nodes to
multi-level taxonomic annotations (subfamilies, tribes, subtribes,
genera).

## Usage

``` r
build_annotation_registry(
  tree,
  constraints_tbl,
  main_level4_values = NULL,
  supp_level4_min_tips = 10L,
  require_monophyly = TRUE
)
```

## Arguments

- tree:

  Object of class \`phylo\`.

- constraints_tbl:

  Data frame containing taxonomic classifications.

- main_level4_values:

  Character vector of main level 4 (genus) annotations. Defaults to
  \`NULL\`.

- supp_level4_min_tips:

  Integer. Minimum tip threshold for secondary level 4 annotations.
  Defaults to \`10L\`.

- require_monophyly:

  Logical. Enforce monophyly constraint? Defaults to \`TRUE\`.

## Value

A named list of annotation node data frames organized by taxonomic
level.
