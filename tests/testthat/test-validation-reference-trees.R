# Version 0.5.4 (decisions D18 and D19 of BMM, 11-11): the dated tree of Carruthers et al. (2026) is
# the fifth published reference tree of validate_phylogenies(), with its own colour, shape and default
# legend label; and under pruning_strategy = "common_set" the panel A3 of the tree-space figure, whose
# coordinates are those of A1, says so in its label. Written before the code.

test_that("the tree-space style has a sixth tree, Carruthers, with its own colour and shape", {
  st <- PhyloCactus:::.tree_space_style()
  trees <- c("FocalTree", "Zuntini", "Thompson", "Amaral", "deVos", "Carruthers")
  expect_identical(names(st$palette), trees)
  # The colours are the first six of the categorical palette, in its fixed order
  expect_identical(unname(st$palette), unname(phylocactus_colours("categorical")[1:6]))
  expect_setequal(names(st$shapes), trees)
  expect_false(anyDuplicated(unname(st$shapes)) > 0)
  expect_length(st$labels, 6L)
  # A label supplied by the user still replaces the default of the same name
  st2 <- PhyloCactus:::.tree_space_style(c(Carruthers = "Carruthers 2026"))
  expect_identical(as.character(st2$labels[[6]]), "Carruthers 2026")
})

test_that("plot_tree_space() draws a Carruthers tree without dropping its points", {
  co <- data.frame(
    analysis_id = "A1_species_all4", analysis_label = "A1",
    tree = c("FocalTree", "Zuntini", "Carruthers"),
    Dimension_1 = c(-0.1, 0.05, 0.05), Dimension_2 = c(0, 0.1, -0.1), metric = "RF"
  )
  expect_no_warning(p <- plot_tree_space(co, metric = "RF"))
  expect_s3_class(p, "ggplot")
})

# A synthetic comparison: three trees of ten species in four genera; the third tree lacks one
# species, so A3 (all species of the shared genera) holds more species than A1 in two trees.
.rt_inputs <- function() {
  d <- tempfile("rt_"); dir.create(d)
  sp <- c("Leuenbergeria_bleo", "Leuenbergeria_zinniiflora", "Opuntia_ficus-indica", "Opuntia_robusta",
          "Opuntia_stricta", "Cereus_jamacaru", "Cereus_hildmannianus", "Cereus_fernambucensis",
          "Parodia_graessneri", "Parodia_haselbergii")
  trees <- c(
    FocalTree = "((Leuenbergeria_bleo,Leuenbergeria_zinniiflora),((Opuntia_ficus-indica,(Opuntia_robusta,Opuntia_stricta)),((Cereus_jamacaru,(Cereus_hildmannianus,Cereus_fernambucensis)),(Parodia_graessneri,Parodia_haselbergii))));",
    Zuntini = "((Leuenbergeria_bleo,Leuenbergeria_zinniiflora),((((Opuntia_ficus-indica,Opuntia_robusta),Opuntia_stricta),(Parodia_graessneri,Parodia_haselbergii)),((Cereus_jamacaru,Cereus_hildmannianus),Cereus_fernambucensis)));",
    Carruthers = "((Leuenbergeria_bleo,Leuenbergeria_zinniiflora),((Opuntia_ficus-indica,Opuntia_robusta),((Cereus_jamacaru,Cereus_fernambucensis),(Parodia_graessneri,Parodia_haselbergii))));"
  )
  paths <- vapply(names(trees), function(nm) {
    f <- file.path(d, paste0(nm, ".tree")); writeLines(trees[[nm]], f); f
  }, character(1))
  genus <- sub("_.*", "", sp)
  ck <- file.path(d, "checklist.csv")
  utils::write.csv(data.frame(pureName = sub("_", " ", sp), Family = "Cactaceae"), ck, row.names = FALSE)
  cm <- file.path(d, "constraints.csv")
  utils::write.csv(data.frame(Genus = genus, Specie_name = sp, Clade = genus, Subfam = "Cactoideae",
                              Family = "Cactaceae", Annot_level_1 = NA, Annot_level_2 = NA,
                              Annot_level_3 = NA, Annot_level_4 = NA), cm, row.names = FALSE)
  list(dir = d, trees = as.list(paths), checklist = ck, constraints = cm)
}

.rt_run <- function(inp, strategy) {
  out <- file.path(inp$dir, paste0("out_", strategy))
  suppressMessages(suppressWarnings(validate_phylogenies(
    trees_mapping_list = inp$trees, checklist_csv = inp$checklist, constraints_map = inp$constraints,
    out_dir = out, min_tips = 4L, pruning_strategy = strategy
  )))
  list(coords = utils::read.csv(file.path(out, "tables", "SUPP_tree_space_coordinates.csv"), stringsAsFactors = FALSE),
       dist = utils::read.csv(file.path(out, "tables", "TABLE_pairwise_tree_distances_main.csv"), stringsAsFactors = FALSE))
}

test_that("under common_set the panel A3 has the coordinates of A1 and its label names A1", {
  skip_if_not_installed("TreeDist")
  inp <- .rt_inputs()
  res <- .rt_run(inp, "common_set")
  co <- res$coords[res$coords$metric == "RF", ]
  a1 <- co[co$analysis_id == "A1_species_all4", ]; a3 <- co[co$analysis_id == "A3_all_species_in_shared_genera", ]
  a3 <- a3[match(a1$tree, a3$tree), ]
  expect_equal(a3$Dimension_1, a1$Dimension_1)
  expect_equal(a3$Dimension_2, a1$Dimension_2)
  lab_a3 <- unique(a3$analysis_label)
  expect_length(lab_a3, 1L)
  expect_true(grepl("A1", sub("^A3", "", lab_a3), fixed = TRUE))
  # The pairwise table of A3 keeps all the species of the shared genera and its own label
  d3 <- res$dist[res$dist$analysis_id == "A3_all_species_in_shared_genera", ]
  expect_false(any(grepl("A1", sub("^A3", "", d3$analysis_label), fixed = TRUE)))
  expect_true(any(d3$n_common_tips == 10))
})

test_that("under pairwise the label of panel A3 does not name A1", {
  skip_if_not_installed("TreeDist")
  inp <- .rt_inputs()
  res <- .rt_run(inp, "pairwise")
  a3 <- res$coords[res$coords$analysis_id == "A3_all_species_in_shared_genera", ]
  expect_false(any(grepl("A1", sub("^A3", "", a3$analysis_label), fixed = TRUE)))
})
