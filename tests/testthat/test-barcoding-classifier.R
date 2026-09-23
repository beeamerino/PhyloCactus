# Tests of the classifier of the molecular diagnostic branch (11_barcoding/7_classifier/).
# Written before the functions they test. Phase 5A of PhyloCactus 0.5.0.
#
# The three states of the validation plan (sec. 1) do not come from a threshold: they come from the
# structure of the tie at the minimum distance. One species there is state 1; several species of one
# genus is state 2, with the candidates declared; a tie across genera is state 3. The threshold is
# swept later over the two scores, in Phase 6, without running the classifier again.
#
# Nothing here measures accuracy on real data. Phase 5A implements the classifier and does not look
# at a single real figure; that is the whole point of splitting 5A from 5B.

.cls_species <- function() {
  stats::setNames(
    c("Opuntia_robusta", "Opuntia_robusta", "Opuntia_stricta", "Cereus_jamacaru", "Opuntia_robusta"),
    c("A1.1", "A2.1", "B1.1", "C1.1", "Q1.1"))
}

# Distance matrix where only the distances from the query Q1.1 matter; the rest are far away
.cls_matrix <- function(d) {
  ids <- c("A1.1", "A2.1", "B1.1", "C1.1", "Q1.1")
  m <- matrix(0.5, length(ids), length(ids), dimnames = list(ids, ids))
  diag(m) <- 0
  for (k in names(d)) {
    m["Q1.1", k] <- d[[k]]
    m[k, "Q1.1"] <- d[[k]]
  }
  m
}

.cls_train <- c("A1.1", "A2.1", "B1.1", "C1.1")

# p-distance by hand, to compare the two methods on the same fixture without pulling in ape
.cls_pdist <- function(seqs) {
  ch <- lapply(seqs, function(s) strsplit(s, "", fixed = TRUE)[[1]])
  n <- length(seqs)
  m <- matrix(0, n, n, dimnames = list(names(seqs), names(seqs)))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i < j) {
        d <- mean(ch[[i]] != ch[[j]])
        m[i, j] <- d
        m[j, i] <- d
      }
    }
  }
  m
}

test_that("one species at the minimum distance gives state 1, with a positive margin", {
  m <- .cls_matrix(c(A1.1 = 0.01, A2.1 = 0.02, B1.1 = 0.10, C1.1 = 0.20))

  r <- .bc_classify_nn(m, .cls_train, "Q1.1", .cls_species())

  expect_equal(nrow(r), 1L)
  expect_equal(r$estado, 1L)
  expect_equal(r$especie_predicha, "Opuntia_robusta")
  expect_equal(r$genero_predicho, "Opuntia")
  expect_equal(r$candidatas, "Opuntia_robusta")
  expect_equal(r$distancia_vecino, 0.01)
  # Margin: nearest sequence of another species (B1.1, 0.10) minus the minimum
  expect_equal(r$margen, 0.09)
  expect_true(is.na(r$motivo))
})

test_that("a tie between species of one genus gives state 2, never state 1, and the candidates are declared", {
  m <- .cls_matrix(c(A1.1 = 0, A2.1 = 0.30, B1.1 = 0, C1.1 = 0.20))

  r <- .bc_classify_nn(m, .cls_train, "Q1.1", .cls_species())

  expect_equal(r$estado, 2L)
  expect_true(is.na(r$especie_predicha))
  expect_equal(r$genero_predicho, "Opuntia")
  expect_equal(r$candidatas, "Opuntia_robusta|Opuntia_stricta")
  expect_equal(r$distancia_vecino, 0)
  # Another species sits at the same distance: the margin is exactly 0
  expect_equal(r$margen, 0)
})

test_that("a tie across genera gives state 3, with the candidates and no genus", {
  m <- .cls_matrix(c(A1.1 = 0, A2.1 = 0.30, B1.1 = 0.40, C1.1 = 0))

  r <- .bc_classify_nn(m, .cls_train, "Q1.1", .cls_species())

  expect_equal(r$estado, 3L)
  expect_true(is.na(r$especie_predicha))
  expect_true(is.na(r$genero_predicho))
  expect_equal(r$candidatas, "Cereus_jamacaru|Opuntia_robusta")
  expect_equal(r$margen, 0)
})

test_that("a query with no distance that has a value is state 3 with the reason declared", {
  m <- .cls_matrix(c(A1.1 = NA, A2.1 = NA, B1.1 = NA, C1.1 = NA))

  r <- .bc_classify_nn(m, .cls_train, "Q1.1", .cls_species())

  expect_equal(r$estado, 3L)
  expect_equal(r$motivo, "sin_posiciones_comparables")
  expect_true(is.na(r$distancia_vecino))
  expect_true(is.na(r$margen))
  expect_true(is.na(r$candidatas) || r$candidatas == "")
})

test_that("the classifier refuses a query that is in its own training set", {
  m <- .cls_matrix(c(A1.1 = 0.01, A2.1 = 0.02, B1.1 = 0.10, C1.1 = 0.20))

  expect_error(.bc_classify_nn(m, c(.cls_train, "Q1.1"), "Q1.1", .cls_species()), "training set")
  # And the folds of Phase 3 never put it there
  lib <- data.frame(sid = c("A1.1", "A2.1", "B1.1", "B2.1"),
                    species = c("Opuntia_robusta", "Opuntia_robusta", "Opuntia_stricta", "Opuntia_stricta"),
                    genus = "Opuntia", locus = "matK", stringsAsFactors = FALSE)
  for (f in .barcoding_folds(lib, scheme = "species")) {
    expect_false(f$test_ids %in% f$train_ids)
  }
})

test_that("the prediction does not depend on the order of the training sequences", {
  m <- .cls_matrix(c(A1.1 = 0.01, A2.1 = 0.02, B1.1 = 0.10, C1.1 = 0.20))
  sp <- .cls_species()

  a <- .bc_classify_nn(m, .cls_train, "Q1.1", sp)
  b <- .bc_classify_nn(m, rev(.cls_train), "Q1.1", sp)
  d <- .bc_classify_nn(m[rev(rownames(m)), rev(colnames(m))], .cls_train, "Q1.1", sp)

  expect_identical(a, b)
  expect_identical(a, d)
})

# Probe of 2026-09-22. The same IdTaxa call works in a plain R session and fails under
# devtools::test() with pkgload, with a dispatch error of match() over XStringSet raised inside
# IdTaxa itself. Checked by hand: BiocGenerics::match() on a DNAStringSet works in that session, and
# the body of .bc_classify_idtaxa() run line by line outside the package returns the right species.
# The mechanism of the loading context is NOT identified, so this probe states the fact and no more:
# it runs the test wherever IdTaxa can run, R CMD check included, and skips where it cannot.
.cls_idtaxa_runs <- function() {
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

test_that("IdTaxa answers with the same columns and the same three states as the nearest neighbour", {
  skip_if_not_installed("DECIPHER")
  if (!suppressWarnings(suppressMessages(.cls_idtaxa_runs()))) {
    skip("DECIPHER::IdTaxa() does not run in this loading context: match() over XStringSet fails to dispatch inside IdTaxa. The same call works in a plain R session.")
  }
  set.seed(5L)
  base <- function() paste(sample(c("A", "C", "G", "T"), 400, replace = TRUE), collapse = "")
  vary <- function(s, k) { for (p in sample(seq_len(400), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1); s }
  spp <- c("Opuntia_robusta", "Cereus_jamacaru", "Mammillaria_elongata")
  seqs <- character(0); sp <- character(0)
  for (i in seq_along(spp)) {
    b <- base()
    for (j in 1:5) {
      sid <- paste0(substr(spp[i], 1, 1), i, j, ".1")
      seqs[sid] <- vary(b, 3)
      sp[sid] <- spp[i]
    }
  }
  train <- setdiff(names(seqs), "O11.1")

  r <- suppressWarnings(suppressMessages(
    .bc_classify_idtaxa(seqs[train], sp[train], seqs[["O11.1"]])))
  nn <- .bc_classify_nn(.cls_pdist(seqs), train, "O11.1", sp)

  expect_setequal(names(r), names(nn))
  expect_true(r$estado %in% c(1L, 2L, 3L))
  expect_equal(r$estado, 1L)
  expect_equal(r$especie_predicha, "Opuntia_robusta")
})

test_that("step 7 writes one row per query, nothing aggregated, and names the step that is missing", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir)
  set.seed(9L)
  base <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")
  vary <- function(s, k) { for (p in sample(seq_len(300), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1); s }
  for (l in c("matK", "rbcL")) {
    x <- c(`Opuntia_robusta|A1.1` = vary(base, 2), `Opuntia_robusta|A2.1` = vary(base, 3),
           `Opuntia_stricta|B1.1` = vary(base, 25), `Opuntia_stricta|B2.1` = vary(base, 28),
           `Cereus_jamacaru|C1.1` = vary(base, 60), `Cereus_horrida|C2.1` = vary(base, 62))
    names(x) <- sub("\\|", paste0("|", l, "_"), names(x))
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), file.path(lib_dir, paste0("LIB_", l, ".fasta")))
  }
  folds_dir <- file.path(tmp, "5_folds")
  suppressMessages(build_barcoding_folds(library_dir = lib_dir, output_dir = folds_dir))
  out_dir <- file.path(tmp, "7_classifier")

  suppressMessages(classify_barcoding_folds(library_dir = lib_dir, folds_dir = folds_dir,
                                            output_dir = out_dir, method = "nn"))

  e <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_predictions_species_nn.csv"),
                       stringsAsFactors = FALSE)
  g <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_predictions_genus_nn.csv"),
                       stringsAsFactors = FALSE)
  fe <- utils::read.csv(file.path(folds_dir, "TABLE_barcoding_folds_species.csv"), stringsAsFactors = FALSE)
  fg <- utils::read.csv(file.path(folds_dir, "TABLE_barcoding_folds_genus.csv"), stringsAsFactors = FALSE)

  # One row per query of the folds, no more and no fewer
  expect_equal(nrow(e), nrow(fe))
  expect_equal(nrow(g), nrow(fg))
  expect_setequal(paste(e$locus, e$sid), paste(fe$locus, fe$sid))
  expect_true(all(c("locus", "esquema", "pliegue", "sid", "especie_verdadera", "genero_verdadero",
                    "metodo", "estado", "especie_predicha", "genero_predicho", "candidatas",
                    "distancia_vecino", "margen", "motivo") %in% names(e)))
  expect_true(all(e$estado %in% c(1L, 2L, 3L)))
  # Nothing aggregated: this phase does not look at any accuracy
  expect_false(any(grepl("exactitud|acierto|precision|exhaustividad|tasa|media|porcentaje",
                         names(e), ignore.case = TRUE)))

  expect_error(classify_barcoding_folds(library_dir = lib_dir, folds_dir = file.path(tmp, "none"),
                                        output_dir = out_dir),
               "build_barcoding_folds")
  expect_error(classify_barcoding_folds(library_dir = lib_dir, folds_dir = folds_dir,
                                        output_dir = file.path("4_Cleaned", "cls")),
               "phylogeny")
})

# Decision of BMM, 2026-09-23: the branch keeps the two conventions the phylogeny already has for a
# run measured in hours. A closing banner with the cactus, so the end of the analysis is visible in
# the console, and the email notification of R/notify.R, because the IdTaxa contrast is measured in
# hours and nobody is going to sit and watch it. Both follow automate_treePL(): the notification
# never fails the run, and a run that fails raises its own error, unchanged.

.cls_fixture <- function(tmp) {
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir, showWarnings = FALSE)
  set.seed(9L)
  base <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")
  vary <- function(s, k) { for (p in sample(seq_len(300), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1); s }
  for (l in c("matK", "rbcL")) {
    x <- c(`Opuntia_robusta|A1.1` = vary(base, 2), `Opuntia_robusta|A2.1` = vary(base, 3),
           `Opuntia_stricta|B1.1` = vary(base, 25), `Opuntia_stricta|B2.1` = vary(base, 28),
           `Cereus_jamacaru|C1.1` = vary(base, 60), `Cereus_horrida|C2.1` = vary(base, 62))
    names(x) <- sub("\\|", paste0("|", l, "_"), names(x))
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), file.path(lib_dir, paste0("LIB_", l, ".fasta")))
  }
  folds_dir <- file.path(tmp, "5_folds")
  suppressMessages(build_barcoding_folds(library_dir = lib_dir, output_dir = folds_dir))
  list(library_dir = lib_dir, folds_dir = folds_dir)
}

test_that("step 7 closes with the banner of the package, the cactus included", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  f <- .cls_fixture(tmp)
  out_dir <- file.path(tmp, "7_classifier")

  salida <- utils::capture.output(
    suppressMessages(classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                              output_dir = out_dir, method = "nn")))
  texto <- paste(salida, collapse = "\n")

  expect_match(texto, "Barcoding Classification Complete", fixed = TRUE)
  expect_match(texto, "\U0001f335", fixed = TRUE)
  expect_match(texto, out_dir, fixed = TRUE)
  expect_match(texto, "nn", fixed = TRUE)
  # The banner counts rows written, which is not a measure of anything: this phase has no accuracy
  expect_false(grepl("accuracy|exactitud|correct", texto, ignore.case = TRUE))
})

test_that("the notification never fails the run of step 7", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  f <- .cls_fixture(tmp)
  out_dir <- file.path(tmp, "7_classifier_notify")
  withr::local_envvar(MY_EMAIL = "")

  expect_message(
    suppressWarnings(utils::capture.output(
      classify_barcoding_folds(library_dir = f$library_dir, folds_dir = f$folds_dir,
                               output_dir = out_dir, method = "nn", notify = TRUE))),
    "No recipient address")

  # The run finished and wrote its tables even though nothing could be sent
  expect_true(file.exists(file.path(out_dir, "TABLE_barcoding_predictions_species_nn.csv")))
  expect_true(file.exists(file.path(out_dir, "TABLE_barcoding_predictions_genus_nn.csv")))
})

test_that("step 7 still fails with its own error when notify is on", {
  # A notification that swallowed a failure would be worse than no notification: an eight hour run
  # would look finished. Same rule as automate_treePL().
  tmp <- withr::local_tempdir()
  withr::local_envvar(MY_EMAIL = "")

  expect_error(
    suppressMessages(classify_barcoding_folds(library_dir = file.path(tmp, "none"),
                                              folds_dir = file.path(tmp, "none"),
                                              output_dir = file.path(tmp, "out"),
                                              notify = TRUE)),
    "finalize_barcoding_library")
})
