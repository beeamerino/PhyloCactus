# Decisions E1 to E6 of BMM (15-10, Phase 12): one visual and editorial standard for every figure and
# vignette. E1, every colour of the package comes from one palette module, by role, validated with
# the palette validator of the dataviz method (categorical order, a diverging outcome scale blue,
# grey and red, one blue sequential ramp, node support with filled symbols for strong support, the
# IUCN colours kept; neutral greys for lines, text and backgrounds). E2, one theme. E3, figure captions below the figure, "**Figure n.**", numbered
# in each vignette, no heading naming a figure. E4, table captions above, "Table n.". E5, Tutorials 1
# to 7 open with Abstract and close with References. Written before the code.

.es_tree <- function(...) {
  f <- testthat::test_path("..", "..", ...)
  skip_if(!file.exists(f), paste(file.path(...), "is not in this tree"))
  f
}

.es_r_files <- function() {
  d <- .es_tree("R")
  list.files(d, pattern = "\\.R$", full.names = TRUE)
}

.es_scripts <- function() {
  d <- system.file("scripts", package = "PhyloCactus")
  skip_if(!nzchar(d), "scripts not installed")
  list.files(d, pattern = "^tutorial-.*\\.R$", full.names = TRUE)
}

.es_vignettes <- function() {
  d <- .es_tree("vignettes")
  list.files(d, pattern = "\\.Rmd$", full.names = TRUE)
}

# Lines of R code only: roxygen and comments dropped
.es_code_lines <- function(f) {
  x <- readLines(f, warn = FALSE)
  x[!grepl("^\\s*#", x)]
}

# Chunks of a vignette: header options and code
.es_chunks <- function(f) {
  x <- readLines(f, warn = FALSE)
  out <- list(); inside <- FALSE
  for (l in x) {
    if (!inside && grepl("^```\\{r", l)) { inside <- TRUE; cur <- list(header = l, code = character(0)); next }
    if (inside && grepl("^```\\s*$", l)) { inside <- FALSE; out[[length(out) + 1L]] <- cur; next }
    if (inside) cur$code <- c(cur$code, l)
  }
  out
}

.es_evaluated <- function(ch) !grepl("eval\\s*=\\s*FALSE|include\\s*=\\s*FALSE", ch$header)

test_that("the palette gives the colours of E1 by role, in their fixed order", {
  pc <- PhyloCactus::phylocactus_colours
  expect_identical(unname(pc("categorical")),
                   c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"))
  expect_identical(pc("outcome"), c("species correct" = "#1c5cab", "genus correct" = "#86b6ef", "abstention" = "#a8a79f",
                                    "wrong species" = "#ee8a85", "wrong genus" = "#b2182b"))
  expect_identical(pc("state"), c("1" = "#1c5cab", "2" = "#86b6ef", "3" = "#a8a79f"))
  # Magnitudes in greys, light to dark, from the neutral fill to the text grey (amendment of BMM, 15-10:
  # colour only for categories; everything else in grey)
  # The greys of the phylogeny figures of Tutorial 3 before 15-10 are the reference (BMM, 15-10: "toda la parte
  # del vignette de visualizacion de la filogenia estaba bien en colores"): grey95 to grey25 for magnitudes
  expect_identical(unname(pc("sequential")), c("#F2F2F2", "#D9D9D9", "#BFBFBF", "#999999", "#7F7F7F", "#595959", "#404040"))
  # Node support exactly as Tutorial 3 drew it (white, grey, black, steelblue)
  expect_identical(pc("support"), c(">= 90" = "#FFFFFF", "70-89" = "#BEBEBE", "< 70" = "#000000", "Constrained" = "#4682B4"))
  # Lines, text and backgrounds of the figures take the neutral greys of the palette, not R colour names;
  # marks of a single series take the mark grey (amendment of E1 by BMM, 15-10: "gris como relleno, no azul")
  # Tutorial 3: text grey20, HPD bars gray60, strips grey90, absence grey95, bars and points grey35
  expect_identical(pc("neutral"), c(text = "#333333", line = "#999999", rule = "#E5E5E5", fill = "#F2F2F2",
                                    mark = "#595959"))
  iucn <- pc("iucn")
  expect_identical(names(iucn), c("EX", "EW", "CR", "EN", "VU", "NT", "LC", "DD", "NE"))
  expect_identical(unname(iucn[c("CR", "EN", "VU", "LC")]), c("#D81E05", "#FC7F3F", "#F9E814", "#60C659"))
  expect_error(pc("rainbow"))
})

test_that("the state colours of the report and of the plotting functions are the outcome scale (E1)", {
  expect_identical(PhyloCactus:::.bc_state_colours(), PhyloCactus::phylocactus_colours("state"))
  a <- data.frame(query = "q", locus = c("matK", "rbcL"), compartment = "plastid", state = c(1L, 2L),
                  predicted_species = c("Opuntia_a", NA), predicted_genus = "Opuntia", genus_confidence = c(90, 80),
                  species_confidence = c(80, 50), reason = c(NA, "genus_only_locus"), stringsAsFactors = FALSE)
  p <- PhyloCactus::plot_barcoding_identification_summary(list(x = a))
  sc <- Filter(function(s) "fill" %in% s$aesthetics, p$scales$scales)
  expect_length(sc, 1L)
  expect_identical(unname(sc[[1]]$palette(3)), unname(PhyloCactus::phylocactus_colours("state")))
})

test_that("presence and absence: grey25 against grey95 as in Tutorial 3, in every panel (E1, amendment of 15-10)", {
  n <- PhyloCactus::phylocactus_colours("neutral")
  s <- data.frame(genus = "Opuntia", species = c("Opuntia_a", "Opuntia_b"), locus = "matK", sequences = c(2L, 0L),
                  stringsAsFactors = FALSE)
  p <- PhyloCactus:::.bc_report_sampling_plot(s, "Opuntia_a")
  sc <- Filter(function(x) "fill" %in% x$aesthetics, p$scales$scales)[[1]]
  present <- PhyloCactus::phylocactus_colours("sequential")[["700"]]
  expect_identical(unname(sc$palette(2)[c("TRUE", "FALSE")]), unname(c(present, n[["fill"]])))
  t <- data.frame(path = "crop", core_trimmed = 0, genus_core_trimmed = 0, stringsAsFactors = FALSE)
  q <- PhyloCactus:::.bc_report_route_plot("plastome", t)
  sq <- Filter(function(x) "fill" %in% x$aesthetics, q$scales$scales)[[1]]
  expect_identical(unname(sq$palette(2)[c("run", "not run")]), unname(c(present, n[["fill"]])))
})

test_that("composite figures carry panel tags from one place, so their style changes in one line", {
  # A, B by default; options(phylocactus.tag_suffix = ")") gives A), B) in every composite figure
  a <- PhyloCactus::phylocactus_panel_tags()
  expect_identical(a$tag_levels, "A")
  expect_identical(a$tag_suffix, "")
  withr::local_options(phylocactus.tag_suffix = ")")
  expect_identical(PhyloCactus::phylocactus_panel_tags()$tag_suffix, ")")
  expect_identical(PhyloCactus:::.phylocactus_tag_labels(3), c("A)", "B)", "C)"))
})

test_that("every composite figure of R/ and of the scripts takes its tags from phylocactus_panel_tags()", {
  composite <- c("wrap_plots", "plot_grid")
  tagged <- c("phylocactus_panel_tags", ".phylocactus_tag_labels")
  units <- function(f) {
    ex <- as.list(parse(f, keep.source = FALSE))
    # In R/, a unit is a function definition; in a script, a top-level expression
    unlist(lapply(ex, function(e) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
                                       identical(e[[3]][[1]], as.name("function"))) list(e[[3]]) else list(e)),
           recursive = FALSE)
  }
  for (f in c(setdiff(.es_r_files(), .es_r_files()[basename(.es_r_files()) == "palette.R"]), .es_scripts())) {
    for (u in units(f)) {
      nm <- all.names(u)
      if (any(composite %in% nm)) expect_true(any(tagged %in% nm), info = paste(basename(f), paste(deparse(u)[1], collapse = "")))
    }
    expect_false(any(grepl("tag_levels\\s*=", .es_code_lines(f))), info = paste(basename(f), "tag_levels written by hand"))
  }
})

test_that("points are filled with the grey of the bars, and TB4 tells its classes apart with the outcome colours", {
  for (f in c(.es_scripts(), .es_vignettes())) {
    x <- .es_code_lines(f)
    pts <- grep("geom_point\\(", x, value = TRUE)
    pts <- pts[!grepl("aes\\(|data\\s*=", pts)]
    expect_true(all(grepl("mark", pts)), info = paste(basename(f), paste(pts[!grepl("mark", pts)], collapse = " | ")))
    # The same grey as the bars: no transparency on points (BMM, 15-10)
    expect_false(any(grepl("alpha", pts)), info = paste(basename(f), "points with alpha"))
  }
  s <- paste(readLines(system.file("scripts", "tutorial-5-cactus-phylogeny-barcoding.R", package = "PhyloCactus")), collapse = "\n")
  expect_true(grepl('phylocactus_colours("outcome")[c("species correct", "genus correct", "wrong species", "wrong genus")]', s, fixed = TRUE))
})

test_that("no colour literal outside the palette module (E1)", {
  hex <- "#[0-9A-Fa-f]{6}\\b"
  for (f in setdiff(.es_r_files(), .es_r_files()[basename(.es_r_files()) == "palette.R"])) {
    x <- .es_code_lines(f)
    # The CSS of the HTML report is not a figure
    if (basename(f) == "barcoding_report.R") x <- x[!grepl("(background|color|border[a-z-]*):|solid #", x)]
    expect_false(any(grepl(hex, x)), info = basename(f))
  }
  for (f in c(.es_scripts(), .es_vignettes())) expect_false(any(grepl(hex, readLines(f, warn = FALSE))), info = basename(f))
})

test_that("no colour scale or named colour outside the palette (E1)", {
  bad_scale <- "scale_(fill|colou?r)_(viridis|brewer|distiller|grey|gray|fermenter|hue)"
  # A colour given by its R name (grey35, steelblue): every name of grDevices::colors() but white and black
  named_cols <- setdiff(grDevices::colors(), c("white", "black"))
  colour_args <- function(x) {
    m <- regmatches(x, gregexpr("\\b(fill|colou?r)\\s*=\\s*\"[A-Za-z0-9]+\"", x))
    v <- sub("^.*\"([A-Za-z0-9]+)\"$", "\\1", unlist(m))
    v[v %in% named_cols]
  }
  files <- c(setdiff(.es_r_files(), .es_r_files()[basename(.es_r_files()) == "palette.R"]), .es_scripts(), .es_vignettes())
  for (f in files) {
    x <- .es_code_lines(f)
    expect_false(any(grepl(bad_scale, x)), info = paste(basename(f), "colour scale"))
    expect_identical(colour_args(x), character(0), info = basename(f))
  }
})

test_that("a single series is drawn in the mark grey, never in a categorical colour (E1, amendment of 15-10)", {
  # A categorical colour tells categories apart; given as the constant fill or colour of a geom it
  # names no category. Single bars, histograms, segments, points and HPD bars take neutral grey.
  constant <- "\\b(fill|colou?r)\\s*=\\s*phylocactus_colours\\(\"categorical\"\\)"
  for (f in c(.es_r_files(), .es_scripts(), .es_vignettes())) {
    expect_false(any(grepl(constant, .es_code_lines(f))), info = basename(f))
  }
})

test_that("one theme: theme_phylocactus() in the package, no other base theme (E2)", {
  # The theme of BMM written in the script of Tutorial 3 before 15-10, moved into the package as it was
  # (T1 to T3 of BMM, 15-10): theme_minimal at base size 12, a light grey panel with a white grid, italic
  # axis titles; its colours come from the theme role of the palette
  expect_identical(PhyloCactus::phylocactus_colours("theme"),
                   c(panel = "#F5F5F5", grid = "#FFFFFF", ticks = "#B3B3B3", axis_title = "#333333", axis_text = "#4D4D4D",
                     title = "#262626", subtitle = "#666666", strip = "#E5E5E5"))
  th <- PhyloCactus::theme_phylocactus()
  expect_s3_class(th, "theme")
  expect_identical(th$text$size, 12)
  expect_identical(th$panel.background$fill, "#F5F5F5")
  expect_identical(th$panel.grid.major$colour, "#FFFFFF")
  expect_identical(th$panel.grid.major$linewidth, 0.6)
  expect_identical(th$panel.grid.minor$linewidth, 0.3)
  expect_identical(th$axis.title$face, "italic")
  expect_identical(th$axis.title$colour, "#333333")
  expect_identical(th$axis.text$colour, "#4D4D4D")
  expect_identical(th$strip.background$fill, "#E5E5E5")
  expect_identical(th$strip.text$face, "bold")
  other <- "\\btheme_(minimal|bw|classic|light|void|gray|grey|linedraw|dark)\\("
  files <- c(setdiff(.es_r_files(), .es_r_files()[basename(.es_r_files()) == "palette.R"]), .es_scripts(), .es_vignettes())
  for (f in files) {
    x <- .es_code_lines(f)
    expect_false(any(grepl(other, x)), info = basename(f))
    expect_false(any(grepl("theme_phylocactus\\s*<-\\s*function", x)), info = paste(basename(f), "defines its own theme"))
  }
})

test_that("figures: a caption below, numbered, and no heading naming a figure or a table (E3)", {
  draws <- "include_graphics\\(|\\bggplot\\(|\\bggtree\\(|plot_barcoding_identification|wrap_plots\\(|print\\(p[0-9]*\\)"
  for (f in .es_vignettes()) {
    v <- readLines(f, warn = FALSE)
    expect_false(any(grepl("^#+ .*\\b(Figure|Table)\\b", v)), info = paste(basename(f), "heading"))
    # Tutorial 3 is the reference (BMM, 15-10): its captions keep the numbers of its text and of the
    # files its script writes (Figure_1 to Figure_17); Figures 7 and 9, the full-size trees, are not
    # rendered and have no caption. The other vignettes number their figures 1, 2, ...
    t3 <- grepl("^tutorial-3-", basename(f))
    expected <- if (t3) c(1:6, 8, 10:17) else NULL
    n <- 0L
    for (ch in Filter(.es_evaluated, .es_chunks(f))) {
      is_fig <- any(grepl(draws, ch$code)) || grepl("fig\\.(width|height)", ch$header)
      has_cap <- grepl("fig\\.cap\\s*=", ch$header)
      if (is_fig) expect_true(has_cap, info = paste(basename(f), ch$header))
      if (has_cap) {
        n <- n + 1L
        k <- if (t3) expected[n] else n
        expect_match(ch$header, paste0("fig\\.cap\\s*=\\s*\"\\*\\*Figure ", k, "\\.\\*\\* [A-Z]"), info = paste(basename(f), ch$header))
        lab <- regmatches(ch$header, regexpr("figure_[0-9]+_", ch$header))
        if (t3 && length(lab)) expect_identical(as.integer(gsub("\\D", "", lab)), as.integer(k), info = ch$header)
      }
    }
    if (t3) {
      expect_identical(n, length(expected), info = "Tutorial 3: 15 captions")
      for (k in c(7L, 9L)) expect_true(any(grepl(paste0("Figure ", k, " is not rendered on this page"), v, fixed = TRUE)), info = k)
    }
  }
})

test_that("tables: a caption above, numbered (E4)", {
  for (f in .es_vignettes()) {
    n <- 0L
    for (ch in Filter(.es_evaluated, .es_chunks(f))) {
      code <- paste(ch$code, collapse = " ")
      k <- gregexpr("kable\\(", code)[[1]]
      if (k[1] < 0) next
      caps <- regmatches(code, gregexpr("caption\\s*=\\s*\"Table [0-9]+\\. [A-Z]", code))[[1]]
      expect_length(caps, length(k))
      for (cp in caps) {
        n <- n + 1L
        expect_match(cp, paste0("Table ", n, "\\. "), info = basename(f))
      }
    }
  }
})

test_that("Tutorials 1 to 7 open with Abstract and close with References (E5)", {
  for (f in .es_vignettes()[grepl("^tutorial-[1-7]-", basename(.es_vignettes()))]) {
    h2 <- grep("^## ", readLines(f, warn = FALSE), value = TRUE)
    expect_identical(h2[1], "## Abstract", info = basename(f))
    expect_identical(h2[length(h2)], "## References", info = basename(f))
  }
})

# BMM (15-10): Tutorial 3 is the reference of the style, and running it must not change its figures.
# The titles and subtitles drawn inside its plots are kept as they were before 15-10.
test_that("the figures of Tutorial 3 keep their titles and subtitles (Tutorial 3 is the reference)", {
  titles <- c("Marker Coverage Across Species", "Distribution of Marker Completeness", "Marker Presence Heatmap",
              "Supermatrix Structure", "Marker Informativeness", "Species Composition",
              "Methodological status across checklist, supermatrix, and tree", "IUCN Categories",
              "Molecular Completeness vs IUCN Category",
              "Threatened Unsampled Species by Genus", "Priority Cactaceae taxa (CR, EN, VU) lacking molecular sequences",
              "Endemic Species by Country", "Endemic Species by Country and IUCN Category", "Most Common Threats",
              "Most Common Habitats")
  s <- paste(readLines(system.file("scripts", "tutorial-3-cactus-phylogeny-visualization.R", package = "PhyloCactus")),
             collapse = "\n")
  v <- paste(readLines(.es_tree("vignettes", "tutorial-3-cactus-phylogeny-visualization.Rmd")), collapse = "\n")
  for (t in titles) {
    expect_true(grepl(paste0("title = \"", t, "\""), s, fixed = TRUE), info = paste("script:", t))
    expect_true(grepl(paste0("title = \"", t, "\""), v, fixed = TRUE), info = paste("vignette:", t))
  }
  # The subtitle of Figure 12 keeps its words; its N is computed from the data (decision of BMM, 15-10,
  # D-b part 2: test-documentation-db-part2.R)
  n12 <- 'subtitle = paste0("Marker occupancy across IUCN threat categories for sampled ingroup taxa (N = ",'
  expect_true(grepl(n12, s, fixed = TRUE), info = "script: subtitle of Figure 12")
  expect_true(grepl(n12, v, fixed = TRUE), info = "vignette: subtitle of Figure 12")
})

# Decision TB2-a of BMM (15-10): TB2 of Tutorial 5 is drawn as the marker presence heatmap of Tutorial 3
# (Figure 3): the whole grid of species and loci, absent in grey95 and present in grey25, no ticks on
# the species axis, loci in italics, no grid, base size 10. The number of sequences stays in the table.
test_that("TB2 is drawn as the marker presence heatmap of Tutorial 3 (TB2-a)", {
  x <- readLines(system.file("scripts", "tutorial-5-cactus-phylogeny-barcoding.R", package = "PhyloCactus"), warn = FALSE)
  from <- grep("^# TB2\\.", x); to <- grep("save_fig\\(p_tb2", x)
  expect_length(from, 1L); expect_length(to, 1L)
  b <- paste(x[from:to], collapse = "\n")
  # The whole grid: table() gives every species and locus; no row dropped
  expect_false(grepl("Freq\\s*>\\s*0\\s*,\\s*\\]", b), info = "absent cells dropped")
  expect_true(grepl("present", b, fixed = TRUE))
  expect_true(grepl("fill = present", b, fixed = TRUE))
  expect_true(grepl(paste0('scale_fill_manual(values = c(phylocactus_colours("neutral")[["fill"]], ',
                           'phylocactus_colours("sequential")[["700"]]),'), b, fixed = TRUE))
  expect_true(grepl('name = "Present", labels = c("No", "Yes")', b, fixed = TRUE))
  expect_false(grepl("gradient|trans\\s*=", b), info = "a gradient of counts")
  expect_true(grepl("theme_phylocactus(base_size = 10)", b, fixed = TRUE))
  for (el in c("axis.text.y = element_blank()", "axis.ticks.y = element_blank()",
               "panel.grid.major = element_blank()", "panel.grid.minor = element_blank()",
               'axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")'))
    expect_true(grepl(el, b, fixed = TRUE), info = el)
})

# BMM (15-10): the figures of Tutorial 5 carry no caption inside the plot; the caption of a figure is
# the one of the vignette (E3).
test_that("no caption inside the plots of the tutorial scripts (E3, TB5, TB6, TB9)", {
  for (f in .es_scripts()) {
    x <- .es_code_lines(f)
    expect_false(any(grepl("caption\\s*=", x)), info = paste(basename(f), paste(grep("caption\\s*=", x, value = TRUE), collapse = " | ")))
  }
})
