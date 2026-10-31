# Closure of 0.5.2 into 0.5.3, phase 4b (audit of 26-10-2026, findings M1, M2, M3, M4, M6, M9, M10):
# the limits the method and the data keep are stated where a reader meets them, and every number of
# those statements that inst/extdata can recompute is recomputed here. Written before the text.

.d53_root <- function(...) testthat::test_path("..", "..", ...)

.d53_read <- function(...) {
  f <- .d53_root(...)
  skip_if(!file.exists(f), paste(file.path(...), "is not in this tree"))
  gsub("\\s+", " ", paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
}

.d53_vignette <- function(n) .d53_read("vignettes", paste0("tutorial-", n, ".Rmd"))

.d53_extdata <- function(name) {
  f <- system.file("extdata", name, package = "PhyloCactus")
  skip_if(!nzchar(f), paste(name, "not installed"))
  f
}

.d53_capture <- function(text, pattern) {
  m <- regmatches(text, regexec(pattern, text))[[1]]
  expect_true(length(m) > 1L, info = pattern)
  if (length(m) < 2L) return(rep(NA_integer_, 6L))
  as.integer(gsub(",", "", m[-1]))
}

.d53_registry <- function() {
  sp <- utils::read.csv(.d53_extdata("phylocactus_table_species.csv"), stringsAsFactors = FALSE)
  skip_if(!all(c("species_class", "retained_markers") %in% names(sp)), "registry columns absent")
  sp
}

.d53_loci <- function(sp) {
  m <- utils::read.csv(.d53_extdata("phylocactus_table_markers.csv"), stringsAsFactors = FALSE)
  intersect(gsub("-", "_", m[[1]]), names(sp))
}

t1 <- function() .d53_vignette("1-cactus-phylogeny-prep")
t2 <- function() .d53_vignette("2-cactus-phylogeny-inference")

test_that("Tutorial 1 and the README say that names are matched to accepted names only (M1)", {
  v <- t1()
  expect_false(grepl("synonyms, obsolete names, infraspecific designations and orthographic variants are resolved", v, fixed = TRUE))
  expect_false(grepl("resolve botanical synonymies", v, fixed = TRUE))
  expect_false(grepl("resolve botanical synonymies", .d53_read("inst", "scripts", "tutorial-1-cactus-phylogeny-prep.R"), fixed = TRUE))
  n <- .d53_capture(v, "rejected ([0-9,]+) of the ([0-9,]+) GenBank names of the ingroup, which carry ([0-9,]+) of its ([0-9,]+) records")
  expect_identical(n, c(329L, 1529L, 962L, 5527L))
  expect_match(v, "are rejected, not resolved", fixed = TRUE)
  r <- .d53_read("README.md")
  expect_match(r, "accepted names only", fixed = TRUE)
})

test_that("Tutorials 1 and 2 state the sampling of the outgroup and what the checks of the constraint assert (M2)", {
  sp <- .d53_registry()
  og <- sp[sp$species_class %in% "outgroup", ]
  v <- t1()
  n <- .d53_capture(v, "([0-9]+) of the ([0-9]+) outgroup terminals carry a single locus")
  expect_identical(n, c(sum(og$retained_markers == 1L), nrow(og)))
  tal <- og[grepl("^(Talinum|Talinella)_", og$species), ]
  n <- .d53_capture(v, "of the ([0-9]+) Talinaceae, ([0-9]+) carry a single locus")
  expect_identical(n, c(nrow(tal), sum(tal$retained_markers == 1L)))
  w <- t2()
  expect_match(w, "assert relationships that the constraint imposes", fixed = TRUE)
  expect_match(.d53_read("inst", "scripts", "tutorial-2-cactus-phylogeny-inference.R"),
               "assert relationships that the constraint imposes", fixed = TRUE)
})

test_that("Tutorials 2 and 4 state the constraint, its conflicts and the result without the bipartition of core Cactoideae II (M3)", {
  w <- t2()
  for (p in c("constrained backbone", "core Cactoideae II", "Amaral", "bootstrap 87",
              "Rhipsalideae (56 terminals)")) {
    expect_match(w, p, fixed = TRUE, info = p)
  }
  v <- .d53_vignette("4-cactus-phylogeny-validation")
  for (p in c("imposed by the constraint", "RF without the constrained bipartitions", "comparison tree only",
              "DPI cannot be split")) {
    expect_match(v, p, fixed = TRUE, info = p)
  }
  # The RF of the reference run quoted beside the RF without the constraint is the distributed one.
  d <- utils::read.csv(.d53_extdata("phylocactus_table_tree_distances.csv"), stringsAsFactors = FALSE)
  quoted <- as.numeric(regmatches(v, regexec("RF of A1 against Zuntini \\*et al\\*\\. \\(2024\\) rises from ([0-9.]+)", v))[[1]][2])
  num <- unlist(d[vapply(d, is.numeric, logical(1))])
  expect_true(any(abs(num - quoted) < 5e-4, na.rm = TRUE), info = "the published RF quoted is in the distributed table")
})

test_that("Tutorial 1 states what the concatenated tree estimates (M4)", {
  sp <- .d53_registry()
  loc <- .d53_loci(sp)
  plastid <- setdiff(loc, c("ITS", "pepC_like"))
  v <- t1()
  expect_match(v, "total-evidence estimate", fixed = TRUE)
  expect_identical(.d53_capture(v, "([0-9]+) of its ([0-9]+) partitions are plastid"), c(length(plastid), length(loc)))
})

test_that("Tutorial 1 states the rule of the distance filter and what it removed (M6)", {
  v <- t1()
  expect_false(grepl("detects sequence length outliers using the interquartile range", v, fixed = TRUE))
  expect_match(v, "Q3 + 1.5 IQR", fixed = TRUE)
  expect_identical(.d53_capture(v, "removed ([0-9]+) sequences from the ([0-9]+) loci retained"), c(290L, 14L))
  expect_match(v, "235 cells", fixed = TRUE)
})

test_that("Tutorial 6 leads with the independent estimate and gives the precision per answer (M9)", {
  v <- .d53_vignette("6-cactus-barcoding-identification")
  # Row 1 of the rates table: species correct / genus correct / no answer / conflict (n.a.) / wrong
  # species / wrong genus, for the present and the absent strata.
  strata <- regmatches(v, gregexpr("[0-9]+ / [0-9]+ / [0-9]+ / n\\.a\\. / [0-9]+ / [0-9]+", v))[[1]][1:2]
  cnt <- lapply(strata, function(s) as.integer(regmatches(s, gregexpr("[0-9]+", s))[[1]]))
  species_answers <- cnt[[1]][1] + cnt[[1]][4] + cnt[[2]][1] + cnt[[2]][4]
  wrong <- cnt[[1]][4] + cnt[[2]][4]
  expect_identical(.d53_capture(v, "([0-9]+) of the ([0-9]+) species answers of the first row are wrong"),
                   c(wrong, species_answers))
  expect_match(v, "set on the library they are evaluated on", fixed = TRUE)
  expect_true(regexpr("Of 94 Sanger records from outside the library", v) < regexpr("\\| Input \\| Unit \\|", v))
})

test_that("Tutorial 2 states the variation among the maximum-likelihood searches (M10)", {
  w <- t2()
  for (p in c("50 distinct topologies", "0.2246", "487 of 911", "4 of 11")) {
    expect_match(w, p, fixed = TRUE, info = p)
  }
})
