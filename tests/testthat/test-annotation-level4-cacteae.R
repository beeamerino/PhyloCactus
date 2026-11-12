# BMM, 2026-11-12 (amendment 0.5.4, section 5.1): the stretches of Cacteae without a label in the
# chronogram take the names of the lineages of Guerrero et al. (2019, Figure 2). A name that ends in
# "clade" is an informal suprageneric name and is not set in italics; each level-4 group of Cacteae in
# the constraint table is a clade of the shipped chronogram, so that build_annotation_registry() draws it.

.label_faces <- function(labels) {
  tr <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  p <- ggtree::ggtree(tr)
  tbl <- data.frame(node = c(6L, 7L)[seq_along(labels)], label = labels, n_tips_target = 2L)
  p <- add_clade_labels_by_level(p, tbl, fontface = "plain")
  txt <- Filter(function(l) inherits(l$geom, "GeomText"), p$layers)
  unname(vapply(txt, function(l) l$aes_params$fontface, ""))
}

test_that("a level-4 name ending in 'clade' stays plain and a genus name is set in italics", {
  skip_if_not_installed("ggtree")
  expect_equal(.label_faces("Mammilloid clade"), "plain")
  expect_equal(.label_faces("Ferocactus"), "italic")
})

test_that("each level-4 group of Cacteae in the constraint table is a clade of the chronogram", {
  skip_if_not_installed("treeio")
  chrono <- system.file("extdata", "phylocactus_chronogram_hpd.tree", package = "PhyloCactus")
  path <- system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus")
  skip_if(!nzchar(chrono) || !nzchar(path), "chronogram or constraint table not installed")
  phy <- treeio::as.phylo(treeio::read.beast(chrono))
  cc <- utils::read.csv(path, stringsAsFactors = FALSE)
  cc <- cc[cc$Annot_level_3 == "Cacteae" & cc$Specie_name %in% phy$tip.label & nzchar(cc$Annot_level_4), ]
  for (g in unique(cc$Annot_level_4)) {
    tips <- cc$Specie_name[cc$Annot_level_4 == g]
    if (length(tips) < 2L) next
    expect_true(ape::is.monophyletic(phy, tips), info = g)
  }
})

# BMM, 2026-11-12: Figures 6 and 8 of Tutorial 3 mark the outgroup, which is not a clade of the tree, with a
# bar over each run of consecutive outgroup terminals, aligned with the bars of the level-4 groups.
test_that("add_outgroup_bar() draws one bar and one label per run of consecutive outgroup terminals", {
  skip_if_not_installed("ggtree")
  tr <- ape::read.tree(text = "((o1:1,(o2:1,(a:1,b:1):1):1):1,o3:1);")
  p <- ggtree::ggtree(tr)
  q <- add_outgroup_bar(p, ingroup_tips = c("a", "b"), label = "Outgroup")
  new <- q$layers[seq.int(length(p$layers) + 1L, length(q$layers))]
  seg <- Filter(function(l) inherits(l$geom, "GeomSegment"), new)
  txt <- Filter(function(l) inherits(l$geom, "GeomText"), new)
  expect_length(seg, 1L)
  expect_length(txt, 1L)
  ys <- sort(p$data$y[p$data$isTip & p$data$label %in% c("o1", "o2", "o3")])
  runs <- split(ys, cumsum(c(1, diff(ys) != 1)))
  expect_equal(nrow(seg[[1]]$data), length(runs))
  expect_equal(sort(seg[[1]]$data$y), sort(vapply(runs, min, 0)))
  expect_equal(unique(txt[[1]]$data$label), "Outgroup")
  expect_identical(add_outgroup_bar(p, ingroup_tips = tr$tip.label), p)
})
