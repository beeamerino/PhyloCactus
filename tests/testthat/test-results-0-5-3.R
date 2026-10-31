# Closure of 0.5.2 into 0.5.3, phase 4c (audit of 26-10-2026, findings M7, M8, M9, M12): the results
# of the analyses of the closure are written where a reader meets them, and every number of them that
# inst/extdata can recompute is recomputed here: the anchors and the S1 ages of the anchor table, and
# the genera that are not monophyletic in the distributed tree. Written before the text.

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

test_that("Tutorial 2 gives the anchor table with the anchors of calibrations_bounds.csv and the distributed S1 ages (M7, D14)", {
  skip_if_not_installed("ape")
  v <- .r53_lines("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd")
  cal <- utils::read.csv(.r53_extdata("calibrations_bounds.csv"), stringsAsFactors = FALSE)
  acp <- cal[cal$mrca == "ACP_root", ]
  scheme <- function(s) acp[grepl(paste0("scheme_", s, "($|_)"), acp$reference), ]
  head <- .r53_row(v, "Node")
  expect_identical(head[2], sprintf("S1 (%.2f, fixed)", scheme("S1")$min))
  expect_identical(head[3], sprintf("S2 (%.2f to %.2f)", scheme("S2")$min, scheme("S2")$max))
  for (k in 3:5) {
    expect_identical(head[k + 1], sprintf("S%d (%.2f, fixed)", k, scheme(paste0("S", k))$min), info = paste0("S", k))
  }
  t <- ape::read.tree(.r53_extdata("phylocactus_chronogram_ml.tree"))
  cc <- utils::read.csv(.r53_extdata("cactus_constraints.csv"), stringsAsFactors = FALSE)
  h <- ape::node.depth.edgelength(t); age <- max(h) - h
  mrca_age <- function(col, vals) age[ape::getMRCA(t, intersect(cc$Specie_name[cc[[col]] %in% vals], t$tip.label))]
  want <- c(`ACP crown` = mrca_age("Family", c("Cactaceae", "Anacampserotaceae", "Portulacaceae")),
            `Cactaceae crown` = mrca_age("Family", "Cactaceae"),
            `Cactoideae crown` = mrca_age("Subfam", "Cactoideae"),
            `Cacteae crown` = mrca_age("Clade", "Cacteae"),
            `Opuntioideae crown` = mrca_age("Subfam", "Opuntioideae"))
  for (n in names(want)) {
    s1 <- as.numeric(sub(" .*", "", .r53_row(v, n)[2]))
    expect_equal(s1, round(want[[n]], 2), info = n)
  }
  f <- .r53_flat("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd")
  for (p in c("conditional on the anchor", "20.98 to 23.62 Ma", "15 % of the replicates", "21.05 Ma")) {
    expect_match(f, p, fixed = TRUE, info = p)
  }
})

test_that("Tutorial 2 states what the cross-validation seeds showed (M7, 3e)", {
  f <- .r53_flat("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd")
  for (p in c("does not identify a smoothing value on this tree", "9.91 to 10.62 Ma", "16.76 to 18.35 Ma",
              "12.68 to 14.15 Ma", "mixes the effect of the smoothing with that of the optimizer")) {
    expect_match(f, p, fixed = TRUE, info = p)
  }
})

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
  s <- .r53_flat("inst", "scripts", "analysis-calibration-anchor-sweep.R")
  expect_false(grepl("only Cactaceae fossil", s, fixed = TRUE))
  expect_match(s, "Eopuntia", fixed = TRUE)
})

# Phase 3g (28-10): the five anchors dated without the three bounds of Hernandez-Hernandez et al.
# (2014). Its ages are in the probe results of the tutorial folder, not in inst/extdata, so the numbers
# are pinned as text here and traced to those files in the tracing table of the closure. The anchor
# section was shortened by BMM in 5b (28-10) to what the schemes give, without the arguments.
test_that("Tutorial 2 states what dating without the internal bounds showed (D14, 3g)", {
  f <- .r53_flat("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd")
  for (p in c("Dated without the three bounds of Hernández-Hernández *et al*., (2014), S1 changes the four Cactaceae nodes by 0.15 to 1.02 Ma",
              "it places nine of ten published crown ages",
              "inside their intervals, against eight under S4, two under S3 and one under S5",
              "the ACP crown is estimated at 46.08 Ma",
              "the Cactaceae nodes differ from S1 by 0.03 to 0.74 Ma")) {
    expect_match(f, p, fixed = TRUE, info = p)
  }
  # The section describes the schemes; the arguments for S1 stay in the minutes of 3b and 3g.
  for (p in c("the effect of the anchor exceeds that of the smoothing", "their ages are the bounds, not estimates",
              "the agreement does not show that S1 is correct")) {
    expect_false(grepl(p, f, fixed = TRUE), info = p)
  }
  # Every source the statement cites is in the reference list, with its DOI.
  expect_match(f, "Wu *et al*. 2026. Integrated analysis of fossils and molecular divergence time estimates a latest Jurassic origin of angiosperms. *Nature Plants*, 12(6), 1242–1251. <https://doi.org/10.1038/s41477-026-02311-x>", fixed = TRUE)
  expect_match(f, "<https://doi.org/10.1073/pnas.1100628108>", fixed = TRUE)
})
