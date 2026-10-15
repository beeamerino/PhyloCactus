# Decisions DB1 to DB5 of BMM (15-10, Phase 12, D-b, part 1): the figures of the documentation are
# taken from the clean run of Phase 12 (recount-db-figures.R, mode now, and recount-phase12-supports.R
# of the audit repository). (a) The figures of 0.4.5 and E1 that the recount replaced are absent from
# the vignettes, the scripts of the tutorials and the documentation, one expectation per figure, so
# that a stale figure cannot come back unnoticed; (b) Tutorial 5 shows the figures its script draws;
# the bootstrap call of Tutorial 2 states the resources measured in Phase 12; the Limitations of
# Tutorial 5 link the rates of Tutorial 6 (DB5). Written before the text.

.db_vignette <- function(name) {
  f <- system.file("doc", paste0(name, ".Rmd"), package = "PhyloCactus")
  if (!nzchar(f)) f <- testthat::test_path("..", "..", "vignettes", paste0(name, ".Rmd"))
  skip_if(!file.exists(f), paste(name, "not found"))
  paste(readLines(f, warn = FALSE), collapse = "\n")
}

.db_tree_file <- function(...) {
  f <- testthat::test_path("..", "..", ...)
  skip_if(!file.exists(f), paste(file.path(...), "is not in this tree"))
  f
}

.db_script <- function(name) {
  f <- system.file("scripts", paste0(name, ".R"), package = "PhyloCactus")
  skip_if(!nzchar(f), paste(name, "script not installed"))
  paste(readLines(f, warn = FALSE), collapse = "\n")
}

.db_absent <- function(text, stale, where) {
  for (s in stale) expect_false(grepl(s, text, fixed = TRUE), info = paste0(where, ": \"", s, "\""))
}

test_that("Tutorial 1 states no figure of 0.4.5 that the Phase 12 recount replaced (T1-a to T1-d)", {
  v <- .db_vignette("tutorial-1-cactus-phylogeny-prep")
  .db_absent(v, c(
    "99 unique accessions across 113 records", "16 ingroup clusters and 24 outgroup clusters",
    "`rbcL` and `matK` accessions of *Portulaca oleracea*", "163 ingroup sequences", "77 of `psbA-trnH` and 86 of `rpL16`",
    "368 ingroup sequences in place of 296", "677 in place of 600",
    "52 sequences of 52 *Portulaca* species, against 50 ingroup sequences", "0.040 of the outgroup 20-mers", "0.26 to 0.59",
    "605 bp in the ingroup", "16 outgroup terminals, 15 are Anacampserotaceae"), "Tutorial 1")
})

test_that("Tutorial 2 and map_branch_supports() state no figure of the 0.4.5 tree (T2-a, T2-b)", {
  v <- .db_vignette("tutorial-2-cactus-phylogeny-inference")
  .db_absent(v, c("34 of the 1021 internal bipartitions", "the other 987", "987 unconstrained nodes",
                  "median FBP of 0.220", "median TBE of 0.952"), "Tutorial 2 vignette")
  s <- .db_script("tutorial-2-cactus-phylogeny-inference")
  .db_absent(s, c("median FBP 0.470", "median TBE 0.783"), "Tutorial 2 script")
  for (f in c(.db_tree_file("R", "inference.R"), .db_tree_file("man", "map_branch_supports.Rd"))) {
    x <- paste(readLines(f, warn = FALSE), collapse = "\n")
    .db_absent(x, c("987 unconstrained nodes", "1024-terminal", "median TBE 0.783", "median FBP 0.470", "0.732",
                    "(median FBP 0.220", "703 nodes"), basename(f))
  }
})

test_that("the bootstrap call of Tutorial 2 states the resources measured in Phase 12 (T2-c)", {
  v <- strsplit(.db_vignette("tutorial-2-cactus-phylogeny-inference"), "\n", fixed = TRUE)[[1]]
  s <- strsplit(.db_script("tutorial-2-cactus-phylogeny-inference"), "\n", fixed = TRUE)[[1]]
  for (x in list(vignette = v, script = s)) {
    mem <- grep('cluster_mem = "16G"', x, fixed = TRUE, value = TRUE)
    tim <- grep('cluster_time = "24:00:00"', x, fixed = TRUE, value = TRUE)
    expect_length(mem, 1L); expect_length(tim, 1L)
    expect_match(mem, "used", fixed = TRUE)
    expect_match(tim, "used", fixed = TRUE)
  }
})

test_that("Tutorial 5 states the release library, not the figures of E1 (T5-a, T5-b, DB4, DB5)", {
  v <- .db_vignette("tutorial-5-cactus-phylogeny-barcoding")
  .db_absent(v, c("7738 sequences", "131 of them own sequences", "1209 species", "holds 1209",
                  "0.6 % of single records", "4 of 358 specimens", "2 of 93 plastomes"), "Tutorial 5")
  lim <- sub("^.*## Limitations", "", v)
  lim <- sub("\n## .*$", "", lim)
  expect_true(grepl("tutorial-6-cactus-barcoding-identification.html", lim, fixed = TRUE))
})

test_that("Tutorial 5 shows the figures its script draws, each present in vignettes/figures (T5-c)", {
  v <- .db_vignette("tutorial-5-cactus-phylogeny-barcoding")
  shown <- c("Figure_TB1_alignment_matK.png", "Figure_TB1_alignment_ITS.png", "Figure_TB2_species_by_locus.png",
             "Figure_TB3_funnel.png", "Figure_TB4_lengths_and_replicates.png", "Figure_TB5_accuracy_against_sampling.png",
             "Figure_TB6_accuracy_against_missing_data.png", "Figure_TB7_taxonomic_coverage.png",
             "Figure_TB9_chilean_genera.png")
  for (p in shown) expect_true(grepl(paste0("figures/", p), v, fixed = TRUE), info = p)
  d <- .db_tree_file("vignettes", "figures")
  for (p in shown) expect_true(file.exists(file.path(d, p)), info = p)
})
