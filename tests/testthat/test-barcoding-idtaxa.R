# Tests of IdTaxa in steps 7, 8 and 9 of the molecular diagnostic branch. Written before the code.
# Phase 6B of PhyloCactus 0.5.0, decisions K1 to K3, K6 and K7 of BMM (26-09).
#
# Measured on 26-09 with DECIPHER 3.9.4, before this file was written: LearnTaxa() is deterministic;
# IdTaxa() gives other confidences at every call unless a seed is set just before it; several
# queries in one IdTaxa() call do not give what one call per query gives, even with the same seed;
# and a run at threshold 0 cut afterwards at t gives exactly what a run at threshold t gives. So:
# one seed per query, set before its own IdTaxa() call; one training per fold; IdTaxa always at
# threshold 0, with the confidence of each rank written, and the threshold applied afterwards.
#
# IdTaxa does not run under devtools::test() (5A minutes, section 4). The tests that call it carry
# the probe of test-barcoding-classifier.R and run under R CMD check and with the installed package.
# The tests of step 9 build their tables by hand and run everywhere.

.idt_runs <- function() {
  tryCatch({
    tr <- Biostrings::DNAStringSet(c(a = "ACGTACGTAACCGGTTAACC", b = "ACGTACGTACCCGGTTAACC",
                                     c = "TTTTGGGGAATTCCGGAATT", d = "TTTTGGGGAATTCCGGAATG"))
    lt <- DECIPHER::LearnTaxa(tr, c("Root;Opuntia;Opuntia_robusta", "Root;Opuntia;Opuntia_stricta",
                                    "Root;Cereus;Cereus_jamacaru", "Root;Cereus;Cereus_horrida"),
                              verbose = FALSE)
    DECIPHER::IdTaxa(Biostrings::DNAStringSet("ACGTACGTAACCGGTTAACC"), lt, strand = "top",
                     threshold = 60, processors = 1L, verbose = FALSE)
    TRUE
  }, error = function(e) FALSE)
}

.idt_skip <- function() {
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  if (!suppressWarnings(suppressMessages(.idt_runs()))) {
    skip("DECIPHER::IdTaxa() does not run in this loading context: match() over XStringSet fails to dispatch inside IdTaxa. The same call works in a plain R session.")
  }
}

.idt_vary <- function(s, k) {
  for (p in sample(seq_len(nchar(s)), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1)
  s
}

# Three genera of two species, three sequences each, in two loci: IdTaxa has a genus and a species
# to find in both schemes, and scheme G has two species per genus to leave out.
.idt_fixture <- function(tmp) {
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir, showWarnings = FALSE, recursive = TRUE)
  set.seed(31L)
  root <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")
  genera <- c("Opuntia", "Cereus", "Mammillaria")
  epithets <- c("alpha", "beta")
  seqs <- list()
  for (l in c("matK", "rbcL")) {
    x <- character(0)
    for (g in seq_along(genera)) {
      gb <- .idt_vary(root, 40)
      for (e in seq_along(epithets)) {
        sb <- .idt_vary(gb, 10)
        for (r in 1:3) {
          x[paste0(genera[g], "_", epithets[e], "|", l, "_", substr(genera[g], 1, 1), e, r, ".1")] <- .idt_vary(sb, 2)
        }
      }
    }
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), file.path(lib_dir, paste0("LIB_", l, ".fasta")))
    seqs[[l]] <- x
  }
  folds_dir <- file.path(tmp, "5_folds")
  suppressMessages(build_barcoding_folds(library_dir = lib_dir, output_dir = folds_dir))
  list(library_dir = lib_dir, folds_dir = folds_dir, seqs = seqs, root = root)
}

.idt_run <- function(f, out_dir, ...) {
  suppressWarnings(suppressMessages(utils::capture.output(
    classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir,
                             output_dir = out_dir, method = "idtaxa", ...))))
  invisible(out_dir)
}

.idt_md5 <- function(dir) {
  unname(tools::md5sum(file.path(dir, paste0("TABLE_barcoding_predictions_", c("species", "genus"),
                                             "_idtaxa.csv"))))
}

# ---- K1: reproducible ---------------------------------------------------------------------------

test_that("two IdTaxa runs of the same folds write the same tables, byte for byte", {
  .idt_skip()
  tmp <- withr::local_tempdir()
  f <- .idt_fixture(tmp)

  set.seed(1L)
  .idt_run(f, file.path(tmp, "a"))
  set.seed(2L); stats::runif(10)
  .idt_run(f, file.path(tmp, "b"))

  expect_identical(.idt_md5(file.path(tmp, "a")), .idt_md5(file.path(tmp, "b")))
})

test_that("an IdTaxa run leaves the session's random state as it found it", {
  .idt_skip()
  tmp <- withr::local_tempdir()
  f <- .idt_fixture(tmp)
  set.seed(3L)
  before <- get(".Random.seed", envir = globalenv())
  .idt_run(f, file.path(tmp, "a"))
  expect_true(identical(get(".Random.seed", envir = globalenv()), before))
})

test_that("the answer of a query does not depend on the other queries of the run", {
  .idt_skip()
  tmp <- withr::local_tempdir()
  f <- .idt_fixture(tmp)
  .idt_run(f, file.path(tmp, "all"))
  g <- utils::read.csv(file.path(tmp, "all", "TABLE_barcoding_predictions_genus_idtaxa.csv"),
                       stringsAsFactors = FALSE)
  # The last query of a fold with several, classified on its own with the seed of that query
  r0 <- g[g$fold == g$fold[1], , drop = FALSE]
  row <- r0[nrow(r0), ]
  x <- f$seqs[[row$locus]]
  sid <- sub("^[^|]*\\|", "", names(x))
  sp <- sub("\\|.*$", "", names(x))
  names(x) <- sid
  folds <- utils::read.csv(file.path(f$folds_dir, "TABLE_barcoding_folds_genus.csv"), stringsAsFactors = FALSE)
  test_ids <- folds$sid[folds$locus == row$locus & folds$fold == row$fold]
  train <- setdiff(sid, test_ids)

  alone <- suppressWarnings(.bc_classify_idtaxa(
    x[train], stats::setNames(sp, sid)[train], x[[row$sid]],
    query_seed = .bc_query_seed(1L, row$locus, "genus", row$fold, row$sid)))

  # The table went through a CSV, which keeps 15 significant digits
  expect_equal(alone$genus_confidence, row$genus_confidence, tolerance = 1e-12)
  expect_equal(alone$species_confidence, row$species_confidence, tolerance = 1e-12)
  expect_identical(alone$state, row$state)
})

test_that("IdTaxa in three chunks, merged, gives the tables of one run without chunks", {
  .idt_skip()
  tmp <- withr::local_tempdir()
  f <- .idt_fixture(tmp)
  .idt_run(f, file.path(tmp, "one"))
  parts <- file.path(tmp, "parts")
  for (k in 1:3) .idt_run(f, parts, chunk = k, n_chunks = 3L)
  suppressMessages(utils::capture.output(
    merge_barcoding_chunks(library_dir = f$library_dir, folds_dir = f$folds_dir, output_dir = parts,
                           method = "idtaxa", n_chunks = 3L)))
  expect_identical(.idt_md5(parts), .idt_md5(file.path(tmp, "one")))
})

# ---- The training draws from the generator too (pilot of 26-09) ----------------------------------

# Found by the pilot of 26-09 on real data: two identical runs gave different tables although each
# IdTaxa() call had its seed. LearnTaxa() classifies its own training sequences in rounds, with
# sample(), and where the labels conflict, as GenBank's do, the training it returns depends on the
# generator. The fixture of this file had clean labels, so it never showed. Added before the fix.
.idt_noisy <- function() {
  set.seed(31L)
  root <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")
  x <- character(0); lab <- character(0)
  for (g in 1:4) {
    gb <- .idt_vary(root, 30)
    for (e in 1:6) {
      k <- sample(0:3, 1)
      sb <- if (k == 0) gb else .idt_vary(gb, k)
      for (r in 1:3) { x <- c(x, .idt_vary(sb, 2)); lab <- c(lab, paste0(c("Opuntia", "Cereus", "Mammillaria", "Echinopsis")[g], "_sp", e)) }
    }
  }
  sw <- c(3, 20, 41); lab[sw] <- lab[sw + 4]
  names(x) <- paste0("s", seq_along(x), ".1")
  list(x = x, lab = stats::setNames(lab, names(x)))
}

test_that("the training of a fold is reproducible, although LearnTaxa draws from the generator", {
  .idt_skip()
  n <- .idt_noisy()
  prints <- vapply(1:6, function(s) {
    set.seed(s)
    paste(round(.bc_idtaxa_train(n$x, n$lab)$fraction, 6), collapse = ",")
  }, character(1))
  # The fact the fix answers: without a seed, the training depends on the session's generator
  expect_gt(length(unique(prints)), 1L)

  set.seed(1L); a <- .bc_idtaxa_train(n$x, n$lab, train_seed = 11L)
  set.seed(2L); b <- .bc_idtaxa_train(n$x, n$lab, train_seed = 11L)
  expect_identical(a, b)
})

test_that("two runs on a library with conflicting labels write the same tables", {
  .idt_skip()
  tmp <- withr::local_tempdir()
  n <- .idt_noisy()
  lib_dir <- file.path(tmp, "4_library"); dir.create(lib_dir)
  y <- n$x
  names(y) <- paste0(n$lab, "|matK_", names(n$x))
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(y), file.path(lib_dir, "LIB_matK.fasta"))
  folds_dir <- file.path(tmp, "5_folds")
  suppressMessages(build_barcoding_folds(library_dir = lib_dir, output_dir = folds_dir))
  f <- list(library_dir = lib_dir, folds_dir = folds_dir)
  set.seed(1L); .idt_run(f, file.path(tmp, "a"))
  set.seed(2L); .idt_run(f, file.path(tmp, "b"))
  expect_identical(.idt_md5(file.path(tmp, "a")), .idt_md5(file.path(tmp, "b")))
})

# ---- K2: one training per fold -------------------------------------------------------------------

test_that("a training learnt once per fold gives the same answer as a training per query", {
  .idt_skip()
  tmp <- withr::local_tempdir()
  f <- .idt_fixture(tmp)
  x <- f$seqs$matK
  sid <- sub("^[^|]*\\|", "", names(x)); sp <- stats::setNames(sub("\\|.*$", "", names(x)), sid)
  names(x) <- sid
  train <- sid[-c(1, 4)]
  trained <- .bc_idtaxa_train(x[train], sp[train], train_seed = 3L)
  for (q in sid[c(1, 4)]) {
    each <- suppressWarnings(.bc_classify_idtaxa(x[train], sp[train], x[[q]], query_seed = 7L, train_seed = 3L))
    once <- suppressWarnings(.bc_classify_idtaxa(x[train], sp[train], x[[q]], query_seed = 7L, trained = trained))
    expect_identical(once, each)
  }
})

# ---- K3: threshold 0, every rank kept -------------------------------------------------------------

test_that("IdTaxa runs at threshold 0 and the three states come from the confidences written", {
  .idt_skip()
  tmp <- withr::local_tempdir()
  f <- .idt_fixture(tmp)
  .idt_run(f, file.path(tmp, "t60"))
  .idt_run(f, file.path(tmp, "t80"), threshold = 80)

  for (sc in c("species", "genus")) {
    a <- utils::read.csv(file.path(tmp, "t60", paste0("TABLE_barcoding_predictions_", sc, "_idtaxa.csv")),
                         stringsAsFactors = FALSE)
    b <- utils::read.csv(file.path(tmp, "t80", paste0("TABLE_barcoding_predictions_", sc, "_idtaxa.csv")),
                         stringsAsFactors = FALSE)
    expect_true(all(c("genus_idtaxa", "species_idtaxa", "genus_confidence", "species_confidence") %in% names(a)))
    # IdTaxa ran at 0 both times: the raw answer is the same whatever threshold was asked for
    expect_identical(a[, c("genus_idtaxa", "species_idtaxa", "genus_confidence", "species_confidence")],
                     b[, c("genus_idtaxa", "species_idtaxa", "genus_confidence", "species_confidence")])
    expect_false(anyNA(a$genus_confidence))
    for (d in list(list(tab = a, t = 60), list(tab = b, t = 80))) {
      p <- d$tab; t <- d$t
      want <- ifelse(p$species_confidence >= t, 1L, ifelse(p$genus_confidence >= t, 2L, 3L))
      expect_identical(p$state, want)
      expect_identical(as.character(p$predicted_species[p$state == 1L]),
                       as.character(p$species_idtaxa[p$state == 1L]))
      expect_true(all(is.na(p$predicted_species[p$state != 1L])))
      expect_identical(as.character(p$predicted_genus[p$state %in% 1:2]),
                       as.character(p$genus_idtaxa[p$state %in% 1:2]))
      expect_true(all(is.na(p$predicted_genus[p$state == 3L])))
    }
  }
})

# ---- K6: CN2 with IdTaxa --------------------------------------------------------------------------

.idt_outgroup <- function(tmp, root) {
  # root is a sequence of the library: an outgroup query that shares no 20-mer with the locus would
  # be noise and not an outgroup (same reason as the fixture of test-barcoding-controls.R)
  og <- file.path(tmp, "outgroup")
  dir.create(og, showWarnings = FALSE, recursive = TRUE)
  set.seed(12L)
  s <- .idt_vary(root, 15)
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(s)))
  alien <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")
  x <- stats::setNames(c(s, rc, alien),
                       paste0(c("Portulaca_amilis", "Portulaca_amilis", "Talinum_paniculatum"),
                              "|matK_h", 1:3, ".1"))
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), file.path(og, "matK.fasta"))
  og
}

test_that("CN2 with IdTaxa classifies the comparable outgroup queries and leaves the no-match ones out", {
  .idt_skip()
  tmp <- withr::local_tempdir()
  f <- .idt_fixture(tmp)
  og <- .idt_outgroup(tmp, unname(f$seqs$matK[1]))
  out_dir <- file.path(tmp, "8_controls")
  dir.create(out_dir)
  # A table of the nearest neighbour already in place must come out untouched
  nn_file <- file.path(out_dir, "TABLE_barcoding_cn2_queries.csv")
  writeLines("\"locus\",\"sid\"\n\"matK\",\"x\"", nn_file)
  before <- tools::md5sum(nn_file)

  suppressWarnings(suppressMessages(utils::capture.output(
    run_barcoding_controls(library_dir = f$library_dir, folds_dir = f$folds_dir, output_dir = out_dir,
                           outgroup_dir = og, method = "idtaxa", controls = "CN2", loci = "matK"))))

  q <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_cn2_queries_idtaxa.csv"), stringsAsFactors = FALSE)
  expect_equal(nrow(q), 3L)
  expect_true(all(c("orientation", "state", "genus_idtaxa", "species_idtaxa", "genus_confidence",
                    "species_confidence") %in% names(q)))
  expect_identical(q$orientation, c("forward", "reverse", "no_match"))
  nm <- q[q$orientation == "no_match", ]
  expect_equal(nm$state, 3L)
  expect_equal(nm$reason, "no_match")
  expect_true(is.na(nm$genus_confidence))
  # The same sequence on either strand gets the same taxa. Not the same confidences: the seed of a
  # query is derived from its sid, and the two copies have different sids (changed on 2026-09-26,
  # before the code, when the seed rule of K1 was written into the function)
  expect_identical(q$genus_idtaxa[1], q$genus_idtaxa[2])
  expect_identical(q$species_idtaxa[1], q$species_idtaxa[2])
  expect_false(is.na(q$genus_confidence[1]))
  expect_false(is.na(q$genus_confidence[2]))
  expect_identical(tools::md5sum(nn_file), before)
})

test_that("CN1 and CN3 are not run with IdTaxa", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  f <- .idt_fixture(tmp)
  for (cn in c("CN1", "CN3")) {
    expect_error(run_barcoding_controls(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                        output_dir = file.path(tmp, "8_controls"), method = "idtaxa",
                                        controls = cn),
                 "nearest neighbour")
  }
})

# ---- K7: the IdTaxa threshold in step 9 -----------------------------------------------------------

# Tables of step 7 and CN2 written by hand, consistent with threshold 60, so step 9 runs everywhere
.idt_hand_tables <- function(tmp) {
  cls <- file.path(tmp, "7_classifier"); ctl <- file.path(tmp, "8_controls")
  dir.create(cls); dir.create(ctl)
  mk <- function(scheme, gconf, sconf) {
    n <- length(gconf)
    true_sp <- rep(c("Opuntia_alpha", "Cereus_beta"), length.out = n)
    g_id <- sub("_.*", "", true_sp)
    s_id <- if (scheme == "species") true_sp else sub("alpha|beta", "gamma", true_sp)
    p <- data.frame(locus = "matK", scheme = scheme, fold = seq_len(n), stratum = NA, sid = paste0("s", seq_len(n)),
                    true_species = true_sp, true_genus = g_id, method = "idtaxa",
                    genus_idtaxa = g_id, species_idtaxa = s_id,
                    genus_confidence = gconf, species_confidence = sconf, stringsAsFactors = FALSE)
    p$state <- ifelse(sconf >= 60, 1L, ifelse(gconf >= 60, 2L, 3L))
    p$predicted_species <- ifelse(p$state == 1L, s_id, NA)
    p$predicted_genus <- ifelse(p$state %in% 1:2, g_id, NA)
    p$candidates <- NA; p$nn_distance <- NA_real_; p$margin <- NA_real_
    p$confidence <- ifelse(p$state == 1L, sconf, ifelse(p$state == 2L, gconf, NA))
    p$reason <- ifelse(p$state == 3L, "low_confidence", NA)
    p
  }
  e <- mk("species", c(99, 95, 90, 85, 70, 65, 55, 40), c(98, 90, 80, 62, 50, 45, 30, 20))
  g <- mk("genus", c(97, 92, 88, 81, 77, 64, 58, 35, 30, 12), c(70, 40, 30, 20, 20, 10, 5, 5, 3, 2))
  utils::write.csv(e, file.path(cls, "TABLE_barcoding_predictions_species_idtaxa.csv"), row.names = FALSE)
  utils::write.csv(g, file.path(cls, "TABLE_barcoding_predictions_genus_idtaxa.csv"), row.names = FALSE)
  og <- data.frame(locus = "matK", sid = paste0("o", 1:5), query_species = "Portulaca_amilis",
                   orientation = c("forward", "forward", "reverse", "forward", "no_match"),
                   genus_idtaxa = "Opuntia", species_idtaxa = "Opuntia_alpha",
                   genus_confidence = c(80, 50, 20, 10, NA), species_confidence = c(60, 30, 10, 5, NA),
                   stringsAsFactors = FALSE)
  og$state <- ifelse(!is.na(og$species_confidence) & og$species_confidence >= 60, 1L,
                     ifelse(!is.na(og$genus_confidence) & og$genus_confidence >= 60, 2L, 3L))
  og$predicted_species <- ifelse(og$state == 1L, og$species_idtaxa, NA)
  og$predicted_genus <- ifelse(og$state %in% 1:2, og$genus_idtaxa, NA)
  og$reason <- ifelse(og$orientation == "no_match", "no_match", ifelse(og$state == 3L, "low_confidence", NA))
  utils::write.csv(og, file.path(ctl, "TABLE_barcoding_cn2_queries_idtaxa.csv"), row.names = FALSE)
  list(cls = cls, ctl = ctl, e = e, g = g, og = og)
}

test_that("the confidence rule gives back the table at 60 and sends everything to state 3 above every confidence", {
  tmp <- withr::local_tempdir()
  h <- .idt_hand_tables(tmp)
  for (p in list(h$e, h$g)) {
    at60 <- .bc_apply_confidence(p, 60)
    expect_identical(at60$state, p$state)
    expect_identical(at60$predicted_species, p$predicted_species)
    expect_identical(at60$predicted_genus, p$predicted_genus)
    top <- .bc_apply_confidence(p, 101)
    expect_true(all(top$state == 3L))
    expect_true(all(is.na(top$predicted_genus)))
    zero <- .bc_apply_confidence(p, 0)
    expect_true(all(zero$state == 1L))
    # The confidences themselves are never rewritten
    expect_identical(top$genus_confidence, p$genus_confidence)
  }
})

test_that("step 9 with IdTaxa fixes t* at the 0.01 quantile of scheme G genus confidence and reports 60 beside it", {
  tmp <- withr::local_tempdir()
  h <- .idt_hand_tables(tmp)
  out_dir <- file.path(tmp, "9_threshold")
  res <- suppressMessages(utils::capture.output(
    r <- sweep_barcoding_threshold(classifier_dir = h$cls, controls_dir = h$ctl, output_dir = out_dir,
                                   method = "idtaxa", figures = FALSE)))

  expect_true(file.exists(file.path(out_dir, "TABLE_barcoding_threshold_curve_idtaxa.csv")))
  expect_true(file.exists(file.path(out_dir, "TABLE_barcoding_threshold_operating_idtaxa.csv")))
  # The tables of the nearest neighbour are not written by an IdTaxa sweep
  expect_false(file.exists(file.path(out_dir, "TABLE_barcoding_threshold_operating.csv")))

  op <- r$operating
  expect_setequal(unique(op$rule), c("quantile", "default_60"))
  expect_equal(nrow(op), 4L)
  t_star <- unname(stats::quantile(h$g$genus_confidence, 0.01, type = 1))
  expect_equal(unique(op$threshold[op$rule == "quantile"]), t_star)
  expect_equal(unique(op$threshold[op$rule == "default_60"]), 60)

  # Outgroup: four comparable queries; rejected are those that end in state 3
  q60 <- op[op$rule == "default_60" & op$scheme == "genus", ]
  expect_equal(q60$outgroup_comparable, 4L)
  expect_equal(q60$outgroup_rejected, 3L)
  qs <- op[op$rule == "quantile" & op$scheme == "genus", ]
  expect_equal(qs$outgroup_rejected, sum(c(80, 50, 20, 10) < t_star))

  # The curve at 60 counts the table as it is
  cv <- r$curve
  e60 <- cv[cv$scheme == "species" & cv$threshold == 60, ]
  expect_equal(nrow(e60), 1L)
  cat5 <- table(factor(.bc_outcome_category(h$e, "species"), levels = .bc_outcome_levels()))
  expect_equal(unname(unlist(e60[, .bc_outcome_levels()])), as.integer(cat5))
})
