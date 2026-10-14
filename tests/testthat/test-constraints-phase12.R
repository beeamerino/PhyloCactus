# Decision K1 of BMM (13-10, Phase 12): the supermatrix of the folder 0 of 29-09 holds 31 terminals
# that cactus_constraints.csv does not list (30 of Cactaceae and Talinum afrum), and
# build_constraint_scaffold() stops on any of them. Each is added with the assignment its genus has in
# the table; every one of the 15 genera has a single assignment there. Written before the data.

.k1_species <- c(
  "Airampoa_erectoclada", "Aylostera_albopectinata", "Aylostera_azurduyensis", "Aylostera_eos",
  "Aylostera_malochii", "Aylostera_mamillosa", "Aylostera_tuberosa", "Cochemiea_dioica", "Cochemiea_saboae",
  "Cochemiea_theresae", "Cylindropuntia_cedrosensis", "Cylindropuntia_parryi", "Ferocactus_pilosus",
  "Mammillaria_columbiana", "Mammillaria_varieaculeata", "Opuntia_matudae", "Oreocereus_doelzianus",
  "Pelecyphora_dasyacantha", "Pelecyphora_hesteri", "Pelecyphora_macromeris", "Pelecyphora_tuberculosa",
  "Pelecyphora_vivipara", "Rebutia_steinbachii", "Stenocactus_ochoterenianus", "Talinum_afrum",
  "Tephrocactus_paediophilus", "Weingartia_cardenasiana", "Weingartia_crispata", "Weingartia_purpurea",
  "Weingartia_tarijensis", "Xiquexique_bohlei")

.k1_table <- function() {
  path <- system.file("extdata", "cactus_constraints.csv", package = "PhyloCactus")
  skip_if(!nzchar(path), "cactus_constraints.csv not installed")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

test_that("the 31 terminals of the Phase 12 supermatrix are in the constraint table, once each (K1)", {
  cf <- .k1_table()
  expect_length(.k1_species, 31L)
  expect_true(all(.k1_species %in% trimws(cf$Specie_name)))
  expect_true(all(table(cf$Specie_name)[.k1_species] == 1L))
})

test_that("each added terminal carries the single assignment of its genus (K1)", {
  cf <- .k1_table()
  cols <- setdiff(names(cf), c("Genus", "Specie_name"))
  for (sp in .k1_species) {
    g <- sub("_.*$", "", sp)
    row <- cf[cf$Specie_name == sp, , drop = FALSE]
    expect_equal(nrow(row), 1L, info = sp)
    expect_equal(row$Genus, g, info = sp)
    others <- unique(cf[cf$Genus == g & cf$Specie_name != sp, cols, drop = FALSE])
    expect_equal(nrow(others), 1L, info = g)
    expect_equal(unname(unlist(row[, cols])), unname(unlist(others[1, cols])), info = sp)
  }
})
