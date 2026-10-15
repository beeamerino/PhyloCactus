# Decision of BMM (15-10, Phase 12, D-b part 2): the family filter of validate_phylogenies() keeps the
# species of Cactaceae in the checklist (Korotkova et al. 2021), not the species of the constraint
# table. Before, a species of Cactaceae in a published tree and absent from the focal sampling was
# removed with the outgroups (19 species of the tree of Zuntini et al. 2024 and 2 of the tree of
# de Vos et al. 2025 in the run of Phase 12). Written before the code.

.ff_checklist <- function() {
  f <- system.file("extdata", "CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx", package = "PhyloCactus")
  skip_if(!nzchar(f), "checklist not installed")
  f
}

test_that("the family names come from the sheet of the family in the checklist", {
  skip_if_not_installed("readxl")
  f <- .ff_checklist()
  cac <- PhyloCactus:::.checklist_family_names(f, family = "Cactaceae")
  # Species of Cactaceae absent from the constraint table are kept
  for (s in c("Acanthocereus_maculatus", "Chamaecereus_silvestrii", "Borzicactus_ventimigliae",
              "Parodia_graessneri", "Leuenbergeria_bleo"))
    expect_true(s %in% cac, info = s)
  # Species of Anacampserotaceae, in the checklist and in the constraint table, are removed
  for (s in c("Anacampseros_kurtzii", "Talinopsis_frutescens"))
    expect_false(s %in% cac, info = s)
  sheet <- readxl::read_excel(f, sheet = "Cactaceae")
  expect_setequal(cac, unique(PhyloCactus:::standardize_taxon(sheet$pureName)))
})

test_that("a CSV checklist gives the family from its Family column, and stops without one", {
  csv <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(pureName = c("Opuntia ficus-indica", "Anacampseros retusa"),
                              Family = c("Cactaceae", "Anacampserotaceae")), csv, row.names = FALSE)
  expect_identical(PhyloCactus:::.checklist_family_names(csv, family = "Cactaceae"), "Opuntia_ficus_indica")
  bad <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(pureName = "Opuntia ficus-indica"), bad, row.names = FALSE)
  expect_error(PhyloCactus:::.checklist_family_names(bad, family = "Cactaceae"), "family", ignore.case = TRUE)
})

test_that("the family filter of validate_phylogenies() uses the checklist, not the constraint table", {
  b <- paste(deparse(body(PhyloCactus::validate_phylogenies)), collapse = "\n")
  expect_true(grepl(".checklist_family_names(checklist_csv", b, fixed = TRUE))
  expect_false(grepl("setdiff(tr_before$tip.label, cactaceae_tips)", b, fixed = TRUE))
})
