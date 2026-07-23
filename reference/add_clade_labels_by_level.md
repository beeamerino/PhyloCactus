# Add Clade Labels to ggtree Plot by Taxonomic Level

Annotates clade labels and vertical bars onto a \`ggtree\` plot object
corresponding to monophyletic node annotations.

## Usage

``` r
add_clade_labels_by_level(
  p,
  label_tbl,
  fontsize = 3,
  barsize = 0.45,
  offset = 0.03,
  offset_text = 0.004,
  fontface = "plain",
  sort_desc = TRUE,
  angle = 0,
  align = TRUE
)
```

## Arguments

- p:

  A \`ggtree\` plot object.

- label_tbl:

  Data frame containing node label annotations.

- fontsize:

  Numeric. Font size for clade labels. Defaults to \`3\`.

- barsize:

  Numeric. Line width for clade annotation bars. Defaults to \`0.45\`.

- offset:

  Numeric. Horizontal offset fraction for clade bars. Defaults to
  \`0.03\`.

- offset_text:

  Numeric. Offset fraction for label text. Defaults to \`0.004\`.

- fontface:

  Character. Font face specification (e.g., \`"plain"\`, \`"bold"\`,
  \`"italic"\`). Defaults to \`"plain"\`.

- sort_desc:

  Logical. Sort clade labels in descending order by tip height? Defaults
  to \`TRUE\`.

- angle:

  Numeric. Rotation angle for text labels in degrees. Defaults to \`0\`.

- align:

  Logical. Align clade bars to common rightmost margin? Defaults to
  \`TRUE\`.

## Value

Updated \`ggtree\` plot object containing clade label layers.
