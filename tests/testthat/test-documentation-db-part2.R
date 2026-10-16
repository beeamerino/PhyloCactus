# Decisions DB1 to DB5 of BMM (15-10, Phase 12, D-b, part 2): the figures of Tutorials 2 (temporal
# bootstraps and Module 10), 3 and 4 are taken from the recount after step 19
# (recount-db-figures.R, mode after19, of the audit repository, run after Tutorial 4 was run again
# with the family filter from the checklist). (a) The figures of 0.4.5 the recount replaced are absent
# from the vignettes, the scripts of the tutorials and the documentation; (b) the figures of the run
# are present; (c) the subtitle of Figure 12 of Tutorial 3 gives the species the figure plots,
# computed from the data (decision of BMM, 15-10). Written before the text.

.db2_vignette <- function(name) {
  f <- system.file("doc", paste0(name, ".Rmd"), package = "PhyloCactus")
  if (!nzchar(f)) f <- testthat::test_path("..", "..", "vignettes", paste0(name, ".Rmd"))
  skip_if(!file.exists(f), paste(name, "not found"))
  paste(readLines(f, warn = FALSE), collapse = "\n")
}

.db2_script <- function(name) {
  f <- system.file("scripts", paste0(name, ".R"), package = "PhyloCactus")
  skip_if(!nzchar(f), paste(name, "script not installed"))
  paste(readLines(f, warn = FALSE), collapse = "\n")
}

.db2_tree_file <- function(...) {
  f <- testthat::test_path("..", "..", ...)
  skip_if(!file.exists(f), paste(file.path(...), "is not in this tree"))
  paste(readLines(f, warn = FALSE), collapse = "\n")
}

.db2_absent <- function(text, stale, where) {
  for (s in stale) expect_false(grepl(s, text, fixed = TRUE), info = paste0(where, ": \"", s, "\""))
}

.db2_present <- function(text, fresh, where) {
  for (s in fresh) expect_true(grepl(s, text, fixed = TRUE), info = paste0(where, ": \"", s, "\""))
}

test_that("Tutorial 2 states the Module 10 of the Phase 12 run, not that of 0.4.5 (T2-d)", {
  v <- .db2_vignette("tutorial-2-cactus-phylogeny-inference")
  s <- .db2_script("tutorial-2-cactus-phylogeny-inference")
  .db2_absent(v, c("FBP 45;", "45.61", "10.3 units", "3.7 units", "221.12", "The reference dataset produces a plateau",
                   "The chronogram is dated at $10^{-12}$", "0.15 Ma", "0.25 Ma", "1.09 Ma", "0.89 Ma",
                   "the plateau of the reference run"), "Tutorial 2 vignette")
  .db2_absent(s, "the plateau of the reference run", "Tutorial 2 script")
  .db2_present(v, c("FBP 81", "45.53 Ma", "323.09", "interior minimum at $10^{-5}$", "0.62 Ma", "0.85 Ma", "0.93 Ma"),
               "Tutorial 2 vignette")
  stale_doc <- "Cactaceae dataset reaches a plateau below 1e-06"
  .db2_absent(.db2_tree_file("R", "treepl_cv.R"), stale_doc, "R/treepl_cv.R")
  .db2_absent(.db2_tree_file("man", "run_treePL_cv.Rd"), stale_doc, "man/run_treePL_cv.Rd")
})

test_that("the temporal bootstrap call of Tutorial 2 states the resources of Phase 12 (T2-c)", {
  v <- strsplit(.db2_vignette("tutorial-2-cactus-phylogeny-inference"), "\n", fixed = TRUE)[[1]]
  s <- strsplit(.db2_script("tutorial-2-cactus-phylogeny-inference"), "\n", fixed = TRUE)[[1]]
  for (x in list(vignette = v, script = s)) {
    from <- grep("generate_temporal_bootstrap_script(", x, fixed = TRUE)
    expect_length(from, 1L)
    blk <- x[from:min(length(x), from + 25L)]
    mem <- grep("cluster_mem =", blk, fixed = TRUE, value = TRUE)
    tim <- grep("cluster_time =", blk, fixed = TRUE, value = TRUE)
    expect_length(mem, 1L); expect_length(tim, 1L)
    expect_match(mem, "6.0 GiB used, 500 replicates", fixed = TRUE)
    expect_match(tim, "1 h 08 min used, 500 replicates", fixed = TRUE)
  }
})

test_that("Tutorial 3 states the registry, the tree and the gaps of the Phase 12 run (T3)", {
  v <- .db2_vignette("tutorial-3-cactus-phylogeny-visualization")
  .db2_absent(v, c("2,003", "1,060", "1,022", "1,024", "(986", "986)", "943", "134 of them", "*Parodia* (24 species)",
                   "*Mammillaria* (22)", "38 outgroup", "the 38 outgroups", "removed 36", "(36)", "N = 986"), "Tutorial 3 vignette")
  .db2_absent(.db2_script("tutorial-3-cactus-phylogeny-visualization"), "N = 986", "Tutorial 3 script")
  .db2_present(v, c("1,965 accepted species", "1,993 species", "1,218 of them (1,190 ingroup, 28 outgroup)",
                    "1,208 tips (1,180 ingroup, 28 outgroup)", "10 ingroup species", "775 accepted species",
                    "102 of them", "*Parodia* (18 species), *Melocactus* (14) and *Coleocephalocereus* (7)"),
               "Tutorial 3 vignette")
})

test_that("the subtitle of Figure 12 of Tutorial 3 counts the species the figure plots (decision of BMM, 15-10)", {
  for (x in list(script = .db2_script("tutorial-3-cactus-phylogeny-visualization"),
                 vignette = .db2_vignette("tutorial-3-cactus-phylogeny-visualization"))) {
    expect_true(grepl('subtitle = paste0("Marker occupancy across IUCN threat categories for sampled ingroup taxa (N = ",',
                      x, fixed = TRUE))
    expect_true(grepl("n_fig12 <- nrow(fig12_data)", x, fixed = TRUE))
  }
})

test_that("Tutorial 4 states the bipartitions of the Phase 12 run (T4)", {
  v <- .db2_vignette("tutorial-4-cactus-phylogeny-validation")
  .db2_absent(v, "1021 bipartitions", "Tutorial 4 vignette")
  .db2_present(v, "all 1205 bipartitions shared", "Tutorial 4 vignette")
})

# Decision W1 of BMM (15-10, D-b T5-e): the schema of the barcoding branch (TB8) is a figure of its own,
# drawn once by the script of the audit repository with the visual language of Fig1, and shown first in
# the Complete Pipeline Execution Workflow of Tutorial 5, before the Setup; Fig1 does not change.
test_that("Tutorial 5 opens its workflow with the schema of the branch, TB8 (W1, T5-e)", {
  f <- system.file("doc", "tutorial-5-cactus-phylogeny-barcoding.Rmd", package = "PhyloCactus")
  if (!nzchar(f)) f <- testthat::test_path("..", "..", "vignettes", "tutorial-5-cactus-phylogeny-barcoding.Rmd")
  skip_if(!file.exists(f), "Tutorial 5 vignette not found")
  v <- readLines(f, warn = FALSE)
  wf <- grep("^## Complete Pipeline Execution Workflow", v)
  setup <- grep("^### Setup", v)
  inc <- grep("figures/Figure_TB8_barcoding_workflow.png", v, fixed = TRUE)
  expect_length(wf, 1L); expect_length(setup, 1L); expect_length(inc, 1L)
  if (length(inc) == 1L && length(wf) == 1L && length(setup) == 1L) {
    expect_true(inc > wf && inc < setup, info = "TB8 between the workflow heading and the Setup")
    hdr <- v[max(grep("^```\\{r", v[seq_len(inc)]))]
    expect_match(hdr, 'fig.cap = "**Figure 1.** ', fixed = TRUE)
  }
  png <- testthat::test_path("..", "..", "vignettes", "figures", "Figure_TB8_barcoding_workflow.png")
  skip_if(!dir.exists(dirname(png)), "vignettes/figures is not in this tree")
  expect_true(file.exists(png))
})
