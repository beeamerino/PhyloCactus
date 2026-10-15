# ------------------------------------------------------------------------------
# Palette and theme of every figure of PhyloCactus (decisions E1 and E2 of BMM, 2026-10-15).
#
# Every colour of the package is defined here, by the job it does, and nowhere else: the functions
# of R/, the scripts of the tutorials and the vignettes take their colours from phylocactus_colours()
# and their style from theme_phylocactus(). The colours were validated with the palette validator
# of the dataviz method on a light surface (#fcfcfb): the categorical order passes every adjacent
# check (worst colour-vision-deficiency Delta E 9.1, normal vision 19.6; the first three slots pass
# all pairs); each arm of the outcome scale passes the ordinal check (one hue, monotone lightness,
# light end over 2:1); the neutral midpoint stands at 2.35:1.
# ------------------------------------------------------------------------------

#' Colours of the Figures of PhyloCactus
#'
#' Returns the colours of one role of the palette the package uses for every figure. Each role does
#' one job: `"categorical"` tells unordered things apart (trees compared, methods, loci, nucleotides),
#' in a fixed order that is never cycled; `"outcome"` and `"state"` give the answers of an
#' identification on a diverging scale, blue for a right answer, grey for none, red for a wrong
#' one; `"sequential"` encodes a magnitude from light to dark blue; `"support"` gives the classes
#' of bootstrap support of a node, filled for strong support, with constrained nodes apart;
#' `"neutral"` holds the greys of text, lines, rules and backgrounds; `"iucn"` keeps the colours of
#' the IUCN Red List categories.
#'
#' @param role Character. One of `"categorical"`, `"outcome"`, `"state"`, `"sequential"`,
#'   `"support"`, `"neutral"` or `"iucn"`.
#' @return A character vector of hexadecimal colours, named by category where the role has
#'   categories.
#' @seealso [theme_phylocactus()], [scale_fill_iucn()].
#' @examples
#' phylocactus_colours("outcome")
#' phylocactus_colours("categorical")[1:3]
#' @export
phylocactus_colours <- function(role = c("categorical", "outcome", "state", "sequential", "support", "neutral", "iucn")) {
  role <- match.arg(role)
  switch(role,
    categorical = c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100",
                    magenta = "#e87ba4", green = "#008300", violet = "#6250d6", red = "#e34948"),
    outcome = c("species correct" = "#1c5cab", "genus correct" = "#86b6ef", "abstention" = "#a8a79f",
                "wrong species" = "#ee8a85", "wrong genus" = "#b2182b"),
    state = c("1" = "#1c5cab", "2" = "#86b6ef", "3" = "#a8a79f"),
    sequential = c("100" = "#cde2fb", "200" = "#9ec5f4", "300" = "#6da7ec", "400" = "#3987e5",
                   "500" = "#256abf", "600" = "#184f95", "700" = "#0d366b"),
    support = c(">= 90" = "#0b0b0b", "70-89" = "#a8a79f", "< 70" = "#ffffff", "Constrained" = "#2a78d6"),
    neutral = c(text = "#52514e", line = "#a8a79f", rule = "#d6d5cf", fill = "#f0efec"),
    iucn = c("EX" = "#000000", "EW" = "#5C155A", "CR" = "#D81E05", "EN" = "#FC7F3F", "VU" = "#F9E814",
             "NT" = "#CCE226", "LC" = "#60C659", "DD" = "#D1D1C6", "NE" = "#FFFFFF"))
}

#' Theme of the Figures of PhyloCactus
#'
#' The one theme of the package figures: `theme_minimal()` with a recessive major grid, no minor
#' grid, black axis text, the legend at the bottom and no title inside the plot (the title of a
#' figure is its caption). Sizes are set for a double-column figure (180 mm), where the base size
#' lands between 8 and 10 points once printed. Trees use `ggtree::theme_tree()` instead.
#'
#' @param base_size Numeric. Base font size in points. Defaults to `9`.
#' @param base_family Character. Base font family. Defaults to `""`.
#' @return A `ggplot2` theme.
#' @seealso [phylocactus_colours()].
#' @examples
#' \dontrun{
#' ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point() + theme_phylocactus()
#' }
#' @export
theme_phylocactus <- function(base_size = 9, base_family = "") {
  n <- phylocactus_colours("neutral")
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size - 1, colour = n[["text"]]),
      plot.caption = ggplot2::element_text(size = base_size - 2, colour = n[["text"]], hjust = 0),
      plot.title.position = "plot",
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1, colour = "black"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = n[["rule"]]),
      legend.position = "bottom",
      legend.key.size = grid::unit(3.5, "mm"),
      legend.text = ggplot2::element_text(size = base_size - 1),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 2)
    )
}
