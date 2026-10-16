# ------------------------------------------------------------------------------
# Palette and theme of every figure of PhyloCactus.
#
# Every colour of the package is defined here, by the job it does, and nowhere else: the functions
# of R/, the scripts of the tutorials and the vignettes take their colours from phylocactus_colours()
# and their style from theme_phylocactus(). The colours were validated with the palette validator
# of the dataviz method on a light surface (#fcfcfb): the categorical order passes every adjacent
# check (worst colour-vision-deficiency Delta E 9.1, normal vision 19.6; the first three slots pass
# all pairs); colour is kept for categories, and magnitudes, presence and single series are
# drawn in the greys of the phylogeny figures of Tutorial 3, the reference: grey35
# for bars and points, grey95 against grey25 for absence and presence, gray60 for HPD bars,
# steelblue for constrained nodes; each arm of the outcome scale passes the ordinal check (one hue, monotone lightness,
# light end over 2:1); the neutral midpoint stands at 2.35:1.
# ------------------------------------------------------------------------------

#' Colours of the Figures of PhyloCactus
#'
#' Returns the colours of one role of the palette the package uses for every figure. Each role does
#' one job: `"categorical"` tells unordered things apart (trees compared, methods, loci, nucleotides),
#' in a fixed order that is never cycled; `"outcome"` and `"state"` give the answers of an
#' identification on a diverging scale, blue for a right answer, grey for none, red for a wrong
#' one; `"sequential"` encodes a magnitude in greys, from light to dark; `"support"` gives the classes
#' of bootstrap support of a node as Tutorial 3 draws them, with constrained nodes apart;
#' `"neutral"` holds the greys of text, lines, rules, backgrounds and the marks of a single series
#' (bars, histograms, segments, points); `"theme"` holds the colours of [theme_phylocactus()];
#' `"iucn"` keeps the colours of
#' the IUCN Red List categories.
#'
#' @param role Character. One of `"categorical"`, `"outcome"`, `"state"`, `"sequential"`,
#'   `"support"`, `"neutral"`, `"theme"` or `"iucn"`.
#' @return A character vector of hexadecimal colours, named by category where the role has
#'   categories.
#' @seealso [theme_phylocactus()], [scale_fill_iucn()].
#' @examples
#' phylocactus_colours("outcome")
#' phylocactus_colours("categorical")[1:3]
#' @export
phylocactus_colours <- function(role = c("categorical", "outcome", "state", "sequential", "support", "neutral", "theme",
                                          "iucn")) {
  role <- match.arg(role)
  switch(role,
    categorical = c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100",
                    magenta = "#e87ba4", green = "#008300", violet = "#6250d6", red = "#e34948"),
    outcome = c("species correct" = "#1c5cab", "genus correct" = "#86b6ef", "abstention" = "#a8a79f",
                "wrong species" = "#ee8a85", "wrong genus" = "#b2182b"),
    state = c("1" = "#1c5cab", "2" = "#86b6ef", "3" = "#a8a79f"),
    sequential = c("100" = "#F2F2F2", "200" = "#D9D9D9", "300" = "#BFBFBF", "400" = "#999999",
                   "500" = "#7F7F7F", "600" = "#595959", "700" = "#404040"),
    support = c(">= 90" = "#FFFFFF", "70-89" = "#BEBEBE", "< 70" = "#000000", "Constrained" = "#4682B4"),
    neutral = c(text = "#333333", line = "#999999", rule = "#E5E5E5", fill = "#F2F2F2", mark = "#595959"),
    theme = c(panel = "#F5F5F5", grid = "#FFFFFF", ticks = "#B3B3B3", axis_title = "#333333", axis_text = "#4D4D4D",
              title = "#262626", subtitle = "#666666", strip = "#E5E5E5"),
    iucn = c("EX" = "#000000", "EW" = "#5C155A", "CR" = "#D81E05", "EN" = "#FC7F3F", "VU" = "#F9E814",
             "NT" = "#CCE226", "LC" = "#60C659", "DD" = "#D1D1C6", "NE" = "#FFFFFF"))
}

#' Theme of the Figures of PhyloCactus
#'
#' The one theme of the package figures, the theme written for the figures of Tutorial 3 and
#' moved into the package as it was: `theme_minimal()` with a light
#' grey panel and a white grid, italic axis titles, bold legend and strip titles. Its colours come
#' from `phylocactus_colours("theme")`. The title of a figure is its caption, not a title inside the
#' plot. Trees use `ggtree::theme_tree()` instead.
#'
#' @param base_size Numeric. Base font size in points. Defaults to `12`.
#' @param base_family Character. Base font family. Defaults to `""`.
#' @return A `ggplot2` theme.
#' @seealso [phylocactus_colours()].
#' @examples
#' \dontrun{
#' ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point() + theme_phylocactus()
#' }
#' @export
theme_phylocactus <- function(base_size = 12, base_family = "") {
  th <- phylocactus_colours("theme")
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = th[["panel"]], color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.major = ggplot2::element_line(color = th[["grid"]], linewidth = 0.6),
      panel.grid.minor = ggplot2::element_line(color = th[["grid"]], linewidth = 0.3),
      axis.ticks = ggplot2::element_line(color = th[["ticks"]], linewidth = 0.4),
      axis.title = ggplot2::element_text(face = "italic", color = th[["axis_title"]]),
      axis.title.x = ggplot2::element_text(face = "italic", margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(face = "italic", margin = ggplot2::margin(r = 8)),
      axis.text = ggplot2::element_text(color = th[["axis_text"]]),
      plot.title = ggplot2::element_text(face = "bold", color = th[["title"]], hjust = 0),
      plot.subtitle = ggplot2::element_text(color = th[["subtitle"]], hjust = 0, margin = ggplot2::margin(b = 8)),
      legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.key = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.9)),
      strip.background = ggplot2::element_rect(fill = th[["strip"]], color = NA),
      strip.text = ggplot2::element_text(face = "bold", color = th[["axis_title"]], size = ggplot2::rel(0.95))
    )
}

#' Panel Tags of the Composite Figures of PhyloCactus
#'
#' Every composite figure of the package and of the tutorial scripts takes its panel tags from this
#' function, so that their style is changed in one place when a journal asks for it: `A`, `B` by
#' default; `options(phylocactus.tag_suffix = ")")` gives `A)`, `B)`, and
#' `options(phylocactus.tag_levels = "a")` gives lower-case letters.
#'
#' @param ... Further arguments of `patchwork::plot_annotation()` (a caption, a theme).
#' @return A `patchwork` annotation to add to a composite figure.
#' @seealso [theme_phylocactus()].
#' @examples
#' \dontrun{
#' patchwork::wrap_plots(p1, p2) + phylocactus_panel_tags()
#' options(phylocactus.tag_suffix = ")")
#' }
#' @export
phylocactus_panel_tags <- function(...) {
  patchwork::plot_annotation(tag_levels = getOption("phylocactus.tag_levels", "A"),
                             tag_prefix = getOption("phylocactus.tag_prefix", ""),
                             tag_suffix = getOption("phylocactus.tag_suffix", ""), ...)
}

#' The first n panel tags, in the style of phylocactus_panel_tags(), for the labels of cowplot
#' @noRd
.phylocactus_tag_labels <- function(n) {
  lev <- getOption("phylocactus.tag_levels", "A")
  base <- switch(lev, a = letters, A = LETTERS, "1" = as.character(seq_len(26L)), LETTERS)
  paste0(getOption("phylocactus.tag_prefix", ""), base[seq_len(n)], getOption("phylocactus.tag_suffix", ""))
}
