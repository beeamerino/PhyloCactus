# Closure of 0.5.2 into 0.5.3, phase 4c (audit of 26-10-2026, findings M7, M8, M9, M12): the results
# of the analyses of the closure are written where a reader meets them, and every number of them that
# inst/extdata can recompute is recomputed here: the genera that are not monophyletic in the distributed
# tree (the table of the anchors moved to test-redating-0-5-3.R in phase 4e). Written before the text.

.r53_root <- function(...) testthat::test_path("..", "..", ...)

.r53_lines <- function(...) {
  f <- .r53_root(...)
  skip_if(!file.exists(f), paste(file.path(...), "is not in this tree"))
  readLines(f, warn = FALSE, encoding = "UTF-8")
}

.r53_flat <- function(...) gsub("\\s+", " ", paste(.r53_lines(...), collapse = " "))

.r53_extdata <- function(name) {
  f <- system.file("extdata", name, package = "PhyloCactus")
  skip_if(!nzchar(f), paste(name, "not installed"))
  f
}

# The cells of the markdown table row that starts with `first`.
.r53_row <- function(lines, first) {
  r <- grep(paste0("^\\| ", first, " \\|"), lines, value = TRUE)
  expect_length(r, 1L)
  if (length(r) != 1L) return(rep(NA_character_, 8L))
  trimws(strsplit(sub("^\\|", "", sub("\\|\\s*$", "", r)), "\\|")[[1]])
}

# The table of the calibration schemes, its scheme 1 ages and the statements of 3e under scheme 1 are
# tested in test-redating-0-5-3.R and traced in the tracing table (phase 4e).

test_that("Tutorial 2 and its script give the TreeAnnotator command of the reference run (M8)", {
  for (x in list(.r53_flat("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd"),
                 .r53_flat("inst", "scripts", "tutorial-2-cactus-phylogeny-inference.R"))) {
    expect_match(x, "treeannotator -burnin 0 -height mean -topology MCC", fixed = TRUE)
    expect_match(x, "TreeAnnotator v2.7.7", fixed = TRUE)
  }
})

test_that("Tutorial 6 sets the baselines beside IdTaxa (M9, 3d)", {
  f <- .r53_flat("vignettes", "tutorial-6-cactus-barcoding-identification.Rmd")
  for (p in c("BLAST top hit", "nearest neighbour without threshold", "1091 / 345 / 88 / n.a. / 749 / 31 (2304)",
              "999 / 424 / 151 / n.a. / 659 / 71 (2304)", "0.593", "0.914")) {
    expect_match(f, p, fixed = TRUE, info = p)
  }
})

test_that("Tutorial 4 lists the genera that are not monophyletic in the distributed tree (M12)", {
  skip_if_not_installed("ape")
  ml <- ape::read.tree(.r53_extdata("phylocactus_ml_tree.tree"))
  cn <- ape::read.tree(.r53_extdata("phylocactus_constraint_tree.tree"))
  sp <- utils::read.csv(.r53_extdata("phylocactus_table_species.csv"), stringsAsFactors = FALSE)
  lab <- ml$tip.label
  tal <- grep("^(Talinum|Talinella)_", lab, value = TRUE)
  # The tree is rooted on Talinaceae: each bipartition is read as the clade on the side without it.
  clades <- function(tr) lapply(ape::prop.part(tr), function(x) {
    s <- tr$tip.label[x]
    if (any(tal %in% s)) setdiff(tr$tip.label, s) else s
  })[-1]
  key <- function(cl) vapply(cl, function(s) paste(sort(s), collapse = "|"), "")
  cl <- clades(ml)
  fbp <- suppressWarnings(as.numeric(ml$node.label))[-1]
  constrained <- key(cl) %in% key(clades(cn))
  ing <- setdiff(lab, sp$species[sp$species_class %in% "outgroup"])
  gen <- sub("_.*", "", ing)
  tg <- table(gen)
  g2 <- names(tg[tg >= 2])
  nonmono <- g2[!vapply(g2, function(g) paste(sort(ing[gen == g]), collapse = "|") %in% key(cl), TRUE)]
  supported <- nonmono[vapply(nonmono, function(g) {
    G <- ing[gen == g]
    any(vapply(seq_along(cl), function(i) {
      C <- cl[[i]]
      !constrained[i] && fbp[i] >= 70 && any(G %in% C) && !all(G %in% C) && !all(C %in% G)
    }, TRUE))
  }, TRUE)]
  v <- .r53_lines("vignettes", "tutorial-4-cactus-phylogeny-validation.Rmd")
  rows <- grep("^\\| \\*[A-Z][a-z]+\\* \\|", v, value = TRUE)
  cells <- lapply(rows, function(r) trimws(strsplit(sub("^\\|", "", sub("\\|\\s*$", "", r)), "\\|")[[1]]))
  genus <- gsub("\\*", "", vapply(cells, `[`, "", 1))
  expect_identical(sort(genus), sort(nonmono))
  expect_identical(as.integer(vapply(cells, `[`, "", 2)), as.integer(tg[genus]))
  expect_identical(sort(genus[grepl("^supported", vapply(cells, function(x) x[length(x)], ""))]), sort(supported))
  f <- .r53_flat("vignettes", "tutorial-4-cactus-phylogeny-validation.Rmd")
  n <- regmatches(f, regexec("Of the ([0-9]+) genera of Cactaceae with two species or more in the tree, ([0-9]+) are not monophyletic", f))[[1]]
  expect_identical(as.integer(n[-1]), c(length(g2), length(nonmono)))
  expect_match(f, "soft polytomy", fixed = TRUE)
})

test_that("Tutorial 4 states the data of the tree of Amaral et al. (2022), and the sweep script no longer calls the pollen the only Cactaceae fossil (M3, M7)", {
  f <- .r53_flat("vignettes", "tutorial-4-cactus-phylogeny-validation.Rmd")
  for (p in c("phylotaR", "IQ-TREE", "ultrafast bootstrap", "80 % missing data")) expect_match(f, p, fixed = TRUE, info = p)
  # Its supports are ultrafast bootstrap values (IQ-TREE), read as such in Tutorials 2 and 4.
  expect_match(.r53_flat("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd"), "ultrafast bootstrap 87", ignore.case = TRUE)
  expect_false(grepl("(bootstrap 87)", f, fixed = TRUE))
  s <- .r53_flat("inst", "scripts", "analysis-calibration-schemes.R")
  expect_false(grepl("only Cactaceae fossil", s, fixed = TRUE))
  expect_match(s, "Eopuntia", fixed = TRUE)
})

# Phase 3g (28-10) dated the five anchors of D14 without the three internal bounds; the re-dating of
# 3h and 4e replaced that section of Tutorial 2 (BMM, 30-10). The sources of the calibration schemes
# and of the reference ages stay in the reference list, with their DOI.
test_that("the sources of the calibration schemes and of the reference ages are in the reference list of Tutorial 2", {
  f <- .r53_flat("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd")
  expect_match(f, "Wu *et al*. 2026. Integrated analysis of fossils and molecular divergence time estimates a latest Jurassic origin of angiosperms. *Nature Plants*, 12(6), 1242–1251. <https://doi.org/10.1038/s41477-026-02311-x>", fixed = TRUE)
  expect_match(f, "<https://doi.org/10.1073/pnas.1100628108>", fixed = TRUE)
})
