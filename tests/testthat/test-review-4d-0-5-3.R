# Closure of 0.5.2 into 0.5.3, phase 4d (review of the text of 0.5.3 by a session that did not write
# it, minutes of 27-10-2026 in the audit repository): the findings BMM accepted (M1 to M5, m1 to m20)
# are stated here where a test can state them. Written before the corrections.

.d4_root <- function(...) testthat::test_path("..", "..", ...)

.d4_flat <- function(...) {
  f <- .d4_root(...)
  skip_if(!file.exists(f), paste(file.path(...), "is not in this tree"))
  gsub("\\s+", " ", paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = " "))
}

.d4_extdata <- function(name) {
  f <- system.file("extdata", name, package = "PhyloCactus")
  skip_if(!nzchar(f), paste(name, "not installed"))
  f
}

.d4_t2 <- function() .d4_flat("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd")
.d4_t4 <- function() .d4_flat("vignettes", "tutorial-4-cactus-phylogeny-validation.Rmd")

test_that("Tutorial 2 does not say that the September result on the root was not repeated on this tree (M1)", {
  v <- .d4_t2()
  expect_false(grepl("has not been repeated on this tree", v, fixed = TRUE))
})

test_that("the ultrafast bootstrap of Amaral et al. (2022) is read at 95 in Tutorials 2 and 4 (M2)", {
  v2 <- .d4_t2()
  v4 <- .d4_t4()
  expect_false(grepl("the latter applied also", v2, fixed = TRUE))
  expect_false(grepl("contradicts with support three", v4, fixed = TRUE))
  expect_true(grepl("contradicts two relationships at that level with support", v4, fixed = TRUE))
  expect_true(grepl("ultrafast bootstrap 87, below the threshold of 95", v2, fixed = TRUE))
  expect_true(grepl("ultrafast bootstrap 87, below the", v4, fixed = TRUE))
  for (v in list(v2, v4)) expect_true(grepl("Hoang *et al*", v, fixed = TRUE))
})

test_that("Tutorial 6 states the wrong species answers per answer, from its own first row (M3)", {
  v <- .d4_flat("vignettes", "tutorial-6-cactus-barcoding-identification.Rmd")
  expect_false(grepl("rare but possible", v, fixed = TRUE))
  row <- regmatches(v, regexpr("own specimen set aside \\| locus \\| [0-9 /n.a]+\\([0-9]+\\) \\| [0-9 /n.a]+\\([0-9]+\\)", v))
  n <- as.integer(regmatches(row, gregexpr("[0-9]+", row))[[1]])
  # present: species / genus / none / wrong species / wrong genus (n); absent likewise
  present_species <- n[1]; present_wrong <- n[4]; absent_wrong <- n[10]
  s <- sprintf("wrong in %d of %d answers if the species is in the library", present_wrong, present_species + present_wrong)
  expect_true(grepl(s, v, fixed = TRUE), info = s)
  tot <- sprintf("%d of %d at the mix of this measurement", present_wrong + absent_wrong,
                 present_species + present_wrong + absent_wrong)
  expect_true(grepl(tot, v, fixed = TRUE), info = tot)
})

test_that("Tutorial 4 and its script say where each reference tree comes from, with the md5 the reference run read (M4)", {
  md5 <- c(amaral = "185926b62a0a01b3187750e664e0fcf7", thompson = "df479c182a14028aa710b508aebaa458",
           zuntini = "319e6985225858b836ab2ff37054b3e5", devos = "5ceb901a0555521da1a14f0815787362")
  where <- c("https://doi.org/10.6084/m9.figshare.17057816", "https://github.com/jamie-thompson/cactaceae",
             "treeoflife.all_support_values.4.0.tree", "Online Resource 4", "treepl_cacti.tre")
  for (src in list(c("vignettes", "tutorial-4-cactus-phylogeny-validation.Rmd"),
                   c("inst", "scripts", "tutorial-4-cactus-phylogeny-validation.R"))) {
    x <- do.call(.d4_flat, as.list(src))
    for (m in md5) expect_true(grepl(m, x, fixed = TRUE), info = paste(src[length(src)], m))
    for (w in where) expect_true(grepl(w, x, fixed = TRUE), info = paste(src[length(src)], w))
    expect_false(grepl("ultra_cacti_JT.tre", x, fixed = TRUE), info = src[length(src)])
  }
  v <- .d4_t4()
  # The pruning of the Kew tree, as the reference run's file was rebuilt (RF 0 on 217 tips)
  for (p in c("_Cactaceae_", "Trichocereus_macrogonus_var._pachanoi_ERR13173055", "keep.tip"))
    expect_true(grepl(p, v, fixed = TRUE), info = p)
  expect_false(grepl("The reference tree used here is the `ASTRAL-III` species tree estimated from their gene trees",
                     v, fixed = TRUE))
})

test_that("the distributed registry carries no checklist content beyond the names and the match (M5)", {
  sp <- utils::read.csv(.d4_extdata("phylocactus_table_species.csv"), nrows = 1, check.names = FALSE)
  dropped <- c("checklist_fullName", "checklist_author", "checklist_rank", "checklist_taxon", "checklist_uuid")
  expect_length(intersect(dropped, names(sp)), 0L)
  expect_true(all(c("species", "accepted_in_checklist", "status_pipeline") %in% names(sp)))
})

test_that("the table of schemes of Tutorial 2 does not say the schemes differ in the fixed age alone (m1)", {
  # The table of the five anchors, and the description of S3 (m2), were replaced by the three
  # calibration schemes in phase 4e.
  v <- .d4_t2()
  expect_false(grepl("differ in the anchor alone", v, fixed = TRUE))
})

test_that("the failed subsamples of the seed probe are stated where they fell (m3)", {
  for (src in list(c("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd"), c("R", "treepl_cv.R"),
                   c("man", "run_treePL_cv.Rd"))) {
    x <- do.call(.d4_flat, as.list(src))
    expect_false(grepl("next to 1e-05", x, fixed = TRUE), info = src[2])
  }
  h <- .d4_flat("man", "run_treePL_cv.Rd")
  expect_false(grepl("fail under most seeds", h, fixed = TRUE))
  expect_true(grepl("at 10 and above every one of five seeds failed in four to seven subsamples", h, fixed = TRUE))
})

test_that("FBP is written in percent in the text, as the distributed tree labels it (m5)", {
  skip_if_not_installed("ape")
  ml <- ape::read.tree(.d4_extdata("phylocactus_ml_tree.tree"))
  expect_identical(max(suppressWarnings(as.numeric(ml$node.label)), na.rm = TRUE), 100)
  for (src in list("README.md", c("vignettes", "PhyloCactus.Rmd"), c("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd"),
                   c("inst", "scripts", "tutorial-2-cactus-phylogeny-inference.R"))) {
    x <- do.call(.d4_flat, as.list(src))
    expect_false(grepl("FBP (< |of |)0\\.[0-9]", x), info = src[length(src)])
  }
  expect_false(grepl("support of 1.000", .d4_t2(), fixed = TRUE))
})

test_that("every work cited in a vignette is in its reference list (m6, m15)", {
  files <- list.files(.d4_root("vignettes"), pattern = "\\.Rmd$", full.names = TRUE)
  for (f in files) {
    x <- readLines(f, warn = FALSE, encoding = "UTF-8")
    i <- grep("^## References", x)
    if (!length(i)) next
    body <- paste(x[seq_len(i[1] - 1L)], collapse = " ")
    refs <- x[seq(i[1], length(x))]
    pat <- "(?:de )?([A-Z][A-Za-zÀ-ſ'-]+)(?: & [A-Z][A-Za-zÀ-ſ'-]+)? \\*et al\\.?\\*?\\.?,? \\(?((?:19|20)[0-9]{2})"
    m <- regmatches(body, gregexpr(pat, body, perl = TRUE))[[1]]
    keys <- unique(sub(pat, "\\1 \\2", m, perl = TRUE))
    for (k in keys) {
      parts <- strsplit(k, " ")[[1]]
      hit <- any(grepl(paste0("^- (de )?", parts[1], ".* ", parts[2], "\\."), refs))
      expect_true(hit, info = paste(basename(f), k))
    }
  }
})

test_that("the sweep script gives the full reference of the fossil pollen (m7)", {
  x <- .d4_flat("inst", "scripts", "analysis-calibration-schemes.R")
  expect_true(grepl("doi:10.1016/j.revpalbo.2025.105463", x, fixed = TRUE))
  expect_true(grepl("which this script does not evaluate", x, fixed = TRUE))
})

test_that("Tutorial 1 names the sheets of the checklist and their sources (m8), and the named groups of the constraint agree (m9)", {
  v <- .d4_flat("vignettes", "tutorial-1-cactus-phylogeny-prep.Rmd")
  expect_false(grepl("Cactaceae, Portulacaceae and Anacampserotaceae", v, fixed = TRUE))
  expect_true(grepl("Hernández-Ledesma", v, fixed = TRUE))
  cs <- utils::read.csv(.d4_extdata("cactus_constraints.csv"), stringsAsFactors = FALSE)
  cl <- unique(cs[, c("Clade", "Family")])
  n_cac <- length(unique(cl$Clade[cl$Family == "Cactaceae"]))
  n_out <- length(unique(cl$Clade[cl$Family != "Cactaceae"]))
  r <- .d4_flat("README.md")
  s <- sprintf("monophyly in %d nodes of Cactaceae and %d outgroup genera", n_cac, n_out)
  expect_true(grepl(s, r, fixed = TRUE), info = s)
  expect_false(grepl("21 named clades", .d4_t2(), fixed = TRUE))
})

test_that("the ground of the 55 exclusions is stated for the accessions it was read for (m12)", {
  for (src in list(c("vignettes", "tutorial-1-cactus-phylogeny-prep.Rmd"), c("man", "assemble_ingroup_phylotar.Rd"),
                   c("inst", "scripts", "tutorial-1-cactus-phylogeny-prep.R"))) {
    x <- do.call(.d4_flat, as.list(src))
    expect_false(grepl("were removed to resolve duplicate assignments", x, fixed = TRUE), info = src[2])
    # Tutorial 1 no longer gives the ground (BMM, 30-10); the help page and the script still do
    if (src[1] != "vignettes") expect_true(grepl("of the 44 found in the metadata of the mining, 43", x, fixed = TRUE), info = src[2])
  }
})

test_that("Tutorial 4 states the root of the Thompson tree and reads genus splits with the hedge they need (m13, m14)", {
  v <- .d4_t4()
  expect_true(grepl("its root is set at 37.24 Ma", v, fixed = TRUE))
  expect_false(grepl("a lack of signal that follows from", v, fixed = TRUE))
})

test_that("the script of Tutorial 2 carries the reading of the seed probe (m16)", {
  x <- .d4_flat("inst", "scripts", "tutorial-2-cactus-phylogeny-inference.R")
  expect_true(grepl("the cross-validation does not identify a smoothing value on this tree", x, fixed = TRUE))
})

test_that("the references are written as the publishers record them (m17)", {
  for (src in list("README.md", c("vignettes", "tutorial-1-cactus-phylogeny-prep.Rmd"),
                   c("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd"),
                   c("vignettes", "tutorial-3-cactus-phylogeny-visualization.Rmd"),
                   c("vignettes", "tutorial-4-cactus-phylogeny-validation.Rmd"),
                   c("vignettes", "tutorial-5-cactus-phylogeny-barcoding.Rmd"),
                   c("vignettes", "tutorial-6-cactus-barcoding-identification.Rmd"))) {
    x <- do.call(.d4_flat, as.list(src))
    expect_false(grepl("*Nature*, 629, 843", x, fixed = TRUE), info = src[length(src)])
    expect_false(grepl("Nature Ecology and Evolution", x, fixed = TRUE), info = src[length(src)])
    expect_false(grepl("Caryophyllales.org - A dynamic", x, fixed = TRUE), info = src[length(src)])
  }
})

test_that("the cross-validation warning says when the lowest value without failures is the selected one (m19)", {
  tmp <- withr::local_tempdir()
  grid <- 10^(1:-8)
  per <- matrix(300 / 10, nrow = length(grid), ncol = 10)
  per[grid %in% c(1e-5, 1e-4), 4] <- 1e36
  log <- unlist(lapply(seq_along(grid), function(i) {
    run <- cumsum(per[i, ])
    c(paste0("smoothing:", grid[i]), paste0("(chisq): ", signif(run, 8)),
      paste0("chisq: (", grid[i], ") ", signif(run[10], 8)))
  }))
  totals <- rowSums(per) + rev(seq_along(grid)) * 0.01
  writeLines(log, file.path(tmp, "cv.log"))
  writeLines(paste0("chisq: (", grid, ") ", signif(totals, 8)), file.path(tmp, "cv"))
  sel <- PhyloCactus:::.select_cv_smoothing(file.path(tmp, "cv"))
  w <- NULL
  withCallingHandlers(
    PhyloCactus:::.check_cv_subsamples(PhyloCactus:::.read_cv_subsamples(file.path(tmp, "cv.log")),
                                       sel$table, sel$smoothing, "run"),
    warning = function(cnd) { w <<- conditionMessage(cnd); invokeRestart("muffleWarning") })
  expect_false(is.null(w))
  expect_match(w, "), the selected value.", fixed = TRUE)
  h <- .d4_flat("man", "automate_treePL.Rd")
  expect_true(grepl("cv_failures", h, fixed = TRUE))
})

test_that("Tutorial 2 says once that the replicates share one topology (m20)", {
  v <- .d4_t2()
  expect_false(grepl("All dated replicates share one topology", v, fixed = TRUE))
  expect_true(grepl("does not change the ages", v, fixed = TRUE))
})

test_that("summarise_iucn_registry() leaves the checklist content out of the registry it writes (M5)", {
  sp <- data.frame(
    species = c("A_a", "A_b"), status_pipeline = c("sampled_and_dated", "unsampled_checklist_gap"),
    in_ml_tree = c(TRUE, FALSE), is_ingroup = c(TRUE, TRUE), in_checklist = c(TRUE, TRUE),
    accepted_in_checklist = c(TRUE, TRUE),
    checklist_pureName = c("A a", "A b"), checklist_fullName = c("A a L.", "A b Mill."),
    checklist_author = c("L.", "Mill."), checklist_rank = c("species", "species"),
    checklist_taxon = c("Taxon", "Taxon"), checklist_uuid = c("u1", "u2"),
    iucn_found = c(TRUE, TRUE), iucn_category = c("CR", "VU"),
    iucn_category_name = c("Critically Endangered", "Vulnerable"), iucn_is_endemic = c(TRUE, FALSE),
    iucn_endemic_locations = c("Mexico", NA), iucn_threat_names = c("Roads", NA),
    iucn_habitat_names = c("Desert", NA), stringsAsFactors = FALSE)
  out <- summarise_iucn_registry(sp, redlist_version = "2026-1", access_date = "2026-10-15")
  dropped <- c("checklist_fullName", "checklist_author", "checklist_rank", "checklist_taxon", "checklist_uuid")
  expect_length(intersect(dropped, names(out$registry)), 0L)
  expect_true(all(c("species", "in_checklist", "accepted_in_checklist", "checklist_pureName") %in% names(out$registry)))
})
