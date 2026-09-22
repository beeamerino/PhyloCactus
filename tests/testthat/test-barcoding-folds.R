# Tests of the validation partitions of the molecular diagnostic branch (11_barcoding/5_folds/).
# Written before the functions they test. Phase 3 of PhyloCactus 0.5.0.
#
# Hard rule of the validation plan (sec. 2): no accession may be in the training set and in the
# evaluation set of the same fold. Scheme E leaves one sequence out and asks for the species; scheme G
# leaves one whole species out and asks for the genus. The folds are built on the library of
# 4_library/ (representatives), are deterministic, and are written to disk by step 5.

.folds_library <- function() {
  data.frame(
    sid     = c("A1.1", "A2.1", "B1.1", "C1.1", "C2.1", "D1.1", "E1.1", "E2.1"),
    species = c("Opuntia_robusta", "Opuntia_robusta", "Opuntia_stricta", "Cereus_jamacaru",
                "Cereus_jamacaru", "Mammillaria_elongata", "Opuntia_robusta", "Opuntia_stricta"),
    genus   = c("Opuntia", "Opuntia", "Opuntia", "Cereus", "Cereus", "Mammillaria", "Opuntia", "Opuntia"),
    locus   = c("matK", "matK", "matK", "matK", "matK", "matK", "ITS", "ITS"),
    stringsAsFactors = FALSE
  )
}

test_that("scheme E leaves one sequence out, keeps the hard rule and never evaluates a species without replica", {
  lib <- .folds_library()

  folds <- .barcoding_folds(lib, scheme = "species")

  for (f in folds) {
    expect_length(intersect(f$train_ids, f$test_ids), 0L)
    expect_length(f$test_ids, 1L)
    expect_identical(f$target, "species")
    expect_setequal(c(f$train_ids, f$test_ids), lib$sid[lib$locus == f$locus])
  }
  # Queries: every sequence of a species with two or more sequences in the locus, once each
  queries <- vapply(folds, function(f) paste(f$locus, f$test_ids), character(1))
  expect_setequal(queries, c("matK A1.1", "matK A2.1", "matK C1.1", "matK C2.1"))
  expect_equal(anyDuplicated(queries), 0L)
  # A species with one sequence in the locus trains and is never a query
  expect_false(any(vapply(folds, function(f) f$stratum, character(1)) %in%
                     c("Opuntia_stricta", "Mammillaria_elongata")))
  expect_true(all(vapply(folds[vapply(folds, function(f) f$locus, character(1)) == "matK"],
                         function(f) all(c("B1.1", "D1.1") %in% f$train_ids), logical(1))))
})

test_that("scheme G leaves a whole species out, only in genera with two or more species in the locus", {
  lib <- .folds_library()

  folds <- .barcoding_folds(lib, scheme = "genus")

  strata <- vapply(folds, function(f) paste(f$locus, f$stratum), character(1))
  expect_setequal(strata, c("matK Opuntia_robusta", "matK Opuntia_stricta",
                            "ITS Opuntia_robusta", "ITS Opuntia_stricta"))
  for (f in folds) {
    d <- lib[lib$locus == f$locus, , drop = FALSE]
    expect_identical(f$target, "genus")
    expect_setequal(f$test_ids, d$sid[d$species == f$stratum])
    # No sequence and no species of the evaluation set is in the training set
    expect_length(intersect(f$train_ids, f$test_ids), 0L)
    expect_false(f$stratum %in% d$species[d$sid %in% f$train_ids])
    # The genus is still represented in the training set by another species
    g <- d$genus[d$species == f$stratum][1]
    expect_true(g %in% d$genus[d$sid %in% f$train_ids])
  }
  # Cereus and Mammillaria have one species in matK: they train and are never evaluated
  expect_false(any(vapply(folds, function(f) f$stratum, character(1)) %in%
                     c("Cereus_jamacaru", "Mammillaria_elongata")))
})

test_that("the folds are deterministic and do not depend on the random seed", {
  lib <- .folds_library()
  set.seed(1L)
  a <- .barcoding_folds(lib, scheme = "species")
  set.seed(99L)
  b <- .barcoding_folds(lib, scheme = "species")
  expect_identical(a, b)
  set.seed(7L)
  expect_identical(.barcoding_folds(lib, scheme = "genus"), .barcoding_folds(lib, scheme = "genus"))
})

test_that("the checker reports a leak introduced on purpose, and nothing on the real folds", {
  lib <- .folds_library()
  clean_e <- .barcoding_folds(lib, scheme = "species")
  clean_g <- .barcoding_folds(lib, scheme = "genus")

  expect_equal(nrow(.barcoding_check_folds(clean_e, lib)), 0L)
  expect_equal(nrow(.barcoding_check_folds(clean_g, lib)), 0L)

  # Scheme E with the query kept in the training set
  leaky_e <- clean_e
  leaky_e[[1]]$train_ids <- c(leaky_e[[1]]$train_ids, leaky_e[[1]]$test_ids)
  v <- .barcoding_check_folds(leaky_e, lib)
  expect_equal(nrow(v), 1L)
  expect_equal(v$motivo, "accesion_en_entrenamiento_y_evaluacion")

  # Scheme G with another sequence of the evaluated species kept in the training set
  leaky_g <- clean_g
  i <- which(vapply(clean_g, function(f) f$stratum == "Opuntia_robusta" && f$locus == "matK", logical(1)))[1]
  leaky_g[[i]]$test_ids <- setdiff(leaky_g[[i]]$test_ids, "A2.1")
  leaky_g[[i]]$train_ids <- c(leaky_g[[i]]$train_ids, "A2.1")
  v <- .barcoding_check_folds(leaky_g, lib)
  expect_equal(nrow(v), 1L)
  expect_equal(v$motivo, "especie_evaluada_en_entrenamiento")
})

test_that("label permutation touches the training set only, keeps the labels and is reproducible", {
  lib <- .folds_library()
  fold <- .barcoding_folds(lib, scheme = "species")[[1]]
  labels <- stats::setNames(lib$species, lib$sid)

  p1 <- .barcoding_permute_labels(fold, labels, seed = 42L)
  p2 <- .barcoding_permute_labels(fold, labels, seed = 42L)

  expect_identical(p1, p2)
  expect_identical(p1[fold$test_ids], labels[fold$test_ids])
  expect_setequal(names(p1), c(fold$train_ids, fold$test_ids))
  expect_equal(sort(unname(p1[fold$train_ids])), sort(unname(labels[fold$train_ids])))
  # It does not leave the session's random state altered
  set.seed(3L)
  before <- .Random.seed
  invisible(.barcoding_permute_labels(fold, labels, seed = 42L))
  expect_identical(.Random.seed, before)
})

test_that("CN1 on a synthetic fixture: permuted labels fall to the base rate, and a leak inflates accuracy to 1", {
  set.seed(11L)
  base <- function() paste(sample(c("A", "C", "G", "T"), 120, replace = TRUE), collapse = "")
  vary <- function(s, k) { for (p in sample(seq_len(120), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1); s }

  # With species signal: three species, four sequences each, close within the species
  spp <- c("Opuntia_robusta", "Cereus_jamacaru", "Mammillaria_elongata")
  seqs <- c(); lib <- NULL
  for (i in seq_along(spp)) {
    b <- base()
    for (j in 1:4) {
      sid <- paste0(substr(spp[i], 1, 1), i, j, ".1")
      seqs[sid] <- vary(b, 2)
      lib <- rbind(lib, data.frame(sid = sid, species = spp[i], genus = sub("_.*$", "", spp[i]),
                                   locus = "matK", stringsAsFactors = FALSE))
    }
  }
  labels <- stats::setNames(lib$species, lib$sid)
  folds <- .barcoding_folds(lib, scheme = "species")
  accuracy <- function(folds, seed = NULL) {
    ok <- vapply(folds, function(f) {
      labs <- if (is.null(seed)) labels else .barcoding_permute_labels(f, labels, seed = seed)
      pred <- .bc_nn_predict(seqs[f$train_ids], labs[f$train_ids], seqs[[f$test_ids]])
      identical(pred, unname(labels[f$test_ids]))
    }, logical(1))
    mean(ok)
  }
  base_rate <- max(table(labels[vapply(folds, function(f) f$test_ids, character(1))])) / length(folds)

  expect_gte(accuracy(folds), 0.8)
  perm <- vapply(1:10, function(k) accuracy(folds, seed = k), numeric(1))
  expect_lte(mean(perm), base_rate + 3 * stats::sd(perm))

  # Without signal: leak-free folds give chance, folds with the query kept in training give 1
  set.seed(12L)
  noise <- stats::setNames(vapply(lib$sid, function(s) base(), character(1)), lib$sid)
  seqs <- noise
  clean <- accuracy(folds)
  leaky_folds <- lapply(folds, function(f) { f$train_ids <- c(f$train_ids, f$test_ids); f })
  expect_equal(accuracy(leaky_folds), 1)
  expect_lt(clean, 1)
  expect_gt(nrow(.barcoding_check_folds(leaky_folds, lib)), 0L)
})

test_that("step 5 writes the fold tables and the summary, and names the step that is missing", {
  skip_if_not_installed("Biostrings")
  tmp <- withr::local_tempdir()
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir)
  d <- .folds_library()
  for (l in unique(d$locus)) {
    x <- d[d$locus == l, ]
    Biostrings::writeXStringSet(
      Biostrings::DNAStringSet(stats::setNames(rep(paste(rep("ACGT", 30), collapse = ""), nrow(x)),
                                               paste0(x$species, "|", x$sid))),
      file.path(lib_dir, paste0("LIB_", l, ".fasta")))
  }

  out <- suppressMessages(build_barcoding_folds(library_dir = lib_dir, output_dir = file.path(tmp, "5_folds")))

  e <- utils::read.csv(file.path(tmp, "5_folds", "TABLE_barcoding_folds_species.csv"), stringsAsFactors = FALSE)
  g <- utils::read.csv(file.path(tmp, "5_folds", "TABLE_barcoding_folds_genus.csv"), stringsAsFactors = FALSE)
  s <- utils::read.csv(file.path(tmp, "5_folds", "TABLE_barcoding_folds_summary.csv"), stringsAsFactors = FALSE)

  expect_setequal(paste(e$locus, e$sid), c("matK A1.1", "matK A2.1", "matK C1.1", "matK C2.1"))
  expect_equal(length(unique(paste(g$locus, g$pliegue))), 4L)
  expect_setequal(g$sid[g$locus == "matK" & g$estrato == "Opuntia_robusta"], c("A1.1", "A2.1"))
  m <- s[s$locus == "matK" & s$esquema == "species", ]
  expect_equal(c(m$pliegues, m$especies_evaluables, m$especies_solo_entrenamiento), c(4L, 2L, 2L))
  mg <- s[s$locus == "matK" & s$esquema == "genus", ]
  expect_equal(c(mg$pliegues, mg$generos_evaluables, mg$generos_solo_entrenamiento), c(2L, 1L, 2L))
  expect_equal(nrow(out$summary), nrow(s))

  expect_error(build_barcoding_folds(library_dir = file.path(tmp, "none"), output_dir = file.path(tmp, "5_folds")),
               "finalize_barcoding_library")
  expect_error(build_barcoding_folds(library_dir = lib_dir, output_dir = file.path("4_Cleaned", "folds")),
               "phylogeny")
})
