# Tests of the negative controls of the molecular diagnostic branch (11_barcoding/8_controls/).
# Written before the functions they test. Phase 5B of PhyloCactus 0.5.0.
#
# The three controls run before any real figure (validation plan, sec. 4). CN1 says whether there is
# leakage, CN2 whether the classifier knows how to keep quiet in front of something that is not a
# cactus, and CN3 measures the bias that retired the section in v0.4.2. None of the three is a
# result of the method: the three are the condition for a result of the method to be interpretable.
#
# The metric functions are written here because CN1 needs them. In this phase they are applied only
# to permuted labels, to outgroup queries and to resubstitution. Not once to the true labels.

.ctl_predictions <- function(esquema = "species") {
  # A prediction table by hand, with the three states and the cases the plan separates
  data.frame(
    locus = "matK", esquema = esquema, pliegue = 1:6, estrato = NA_character_,
    sid = paste0("q", 1:6),
    especie_verdadera = c("Opuntia_robusta", "Opuntia_robusta", "Opuntia_robusta",
                          "Opuntia_stricta", "Cereus_jamacaru", "Cereus_jamacaru"),
    genero_verdadero = c("Opuntia", "Opuntia", "Opuntia", "Opuntia", "Cereus", "Cereus"),
    metodo = "nn",
    estado = c(1L, 1L, 2L, 3L, 1L, 2L),
    especie_predicha = c("Opuntia_robusta", "Opuntia_stricta", NA, NA, "Cereus_jamacaru", NA),
    genero_predicho = c("Opuntia", "Opuntia", "Opuntia", NA, "Cereus", "Opuntia"),
    candidatas = c("Opuntia_robusta", "Opuntia_stricta", "Opuntia_robusta|Opuntia_stricta",
                   "Cereus_jamacaru|Opuntia_robusta", "Cereus_jamacaru",
                   "Opuntia_robusta|Opuntia_stricta"),
    distancia_vecino = 0, margen = 0, confianza = NA_real_, motivo = NA_character_,
    stringsAsFactors = FALSE
  )
}

.ctl_library <- function(tmp) {
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir, showWarnings = FALSE, recursive = TRUE)
  set.seed(4L)
  vary <- function(s, k) { for (p in sample(seq_len(300), k)) substr(s, p, p) <- sample(c("A","C","G","T"), 1); s }
  base <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")
  for (l in c("matK", "rbcL")) {
    a <- vary(base, 5); b <- vary(base, 60); cc <- vary(base, 120)
    x <- c(vary(a, 1), vary(a, 2), vary(b, 1), vary(b, 2), vary(cc, 1), vary(cc, 2))
    names(x) <- paste0(c("Opuntia_robusta", "Opuntia_robusta", "Opuntia_stricta", "Opuntia_stricta",
                         "Cereus_jamacaru", "Cereus_horrida"), "|", l, "_s", 1:6, ".1")
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), file.path(lib_dir, paste0("LIB_", l, ".fasta")))
  }
  folds_dir <- file.path(tmp, "5_folds")
  suppressMessages(build_barcoding_folds(library_dir = lib_dir, output_dir = folds_dir))
  list(library_dir = lib_dir, folds_dir = folds_dir, base = base)
}

.ctl_mafft_ok <- function() {
  # A probe, not an assumption. Since 2026-09-25 CN2 aligns every query against the library, so the
  # controls need MAFFT at run time, which they did not need before.
  isTRUE(tryCatch(.bc_assert_mafft("mafft"), error = function(e) FALSE))
}

.ctl_outgroup <- function(tmp, base, largos = c(240, 300, 355)) {
  # The outgroup as it comes out of step 1: unaligned, and every sequence of a different length.
  # The first version of this fixture wrote three sequences of exactly the 300 bases of the library,
  # so the truncation to min(ncol) never truncated anything and the defect stayed hidden until the
  # real data hit it on 2026-09-25.
  og <- file.path(tmp, "outgroup")
  dir.create(og, showWarnings = FALSE, recursive = TRUE)
  set.seed(8L)
  ajeno <- function(n) {
    s <- base
    # Fifteen substitutions, not a hundred: measured on the real data of 2026-09-25, the outgroup
    # shares between 0.23 and 0.82 of its 20-mers with the library of its locus, so a fixture that
    # shares none of them would not be an outgroup, it would be noise
    for (p in sample(seq_len(300), 15)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1)
    if (n > 300) paste0(s, paste(sample(c("A", "C", "G", "T"), n - 300, replace = TRUE), collapse = ""))
    else substr(s, 1, n)
  }
  x <- stats::setNames(
    vapply(largos, ajeno, character(1)),
    paste0(c("Portulaca_amilis", "Talinum_paniculatum", "Anacampseros_filamentosa"),
           "|matK_o", seq_along(largos), ".1"))
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), file.path(og, "matK.fasta"))
  og
}

test_that("the metric counts a hit as the plan defines it, and never a candidate set as a species hit", {
  e <- .bc_metric_accuracy(.ctl_predictions("species"))
  g <- .bc_metric_accuracy(.ctl_predictions("genus"))

  # Species: only state 1 with the right species. Rows 1 and 5 hit; row 2 is a confident error,
  # rows 3 and 6 are genus with ambiguity, row 4 is not assignable
  expect_equal(e$consultas, 6L)
  expect_equal(e$aciertos, 2L)
  expect_equal(e$exactitud, 2 / 6)
  # Genus: rows 1, 2, 3 and 5 hit; row 6 gives the wrong genus, row 4 gives none
  expect_equal(g$aciertos, 4L)
  expect_equal(g$exactitud, 4 / 6)
})

test_that("the base rate is the majority class among the queries evaluated, not of the library", {
  e <- .bc_metric_accuracy(.ctl_predictions("species"))
  g <- .bc_metric_accuracy(.ctl_predictions("genus"))

  # Three of the six queries are Opuntia_robusta; four of the six are Opuntia
  expect_equal(e$tasa_base, 3 / 6)
  expect_equal(g$tasa_base, 4 / 6)
})

test_that("CN1 declares leakage when the permuted accuracy does not fall, and not when it falls", {
  # Permuted accuracies around the base rate: no leakage
  limpio <- .bc_cn1_verdict(c(0.34, 0.30, 0.36, 0.33, 0.31, 0.35, 0.32, 0.34, 0.33, 0.30),
                            tasa_base = 1 / 3)
  # Permuted accuracies that stay high: the labels are reaching the classifier through some other way
  con_fuga <- .bc_cn1_verdict(c(0.98, 0.99, 1.00, 0.97, 0.99, 1.00, 0.98, 0.99, 0.98, 1.00),
                              tasa_base = 1 / 3)

  expect_false(limpio$hay_fuga)
  expect_true(con_fuga$hay_fuga)
  expect_equal(limpio$permutaciones, 10L)
  # The threshold is the rule of the plan, written out: base rate plus three standard deviations
  expect_equal(limpio$umbral, 1 / 3 + 3 * stats::sd(c(0.34, 0.30, 0.36, 0.33, 0.31, 0.35, 0.32, 0.34, 0.33, 0.30)))
})

test_that("CN1 writes one row per permutation and its verdict, with the base rate beside it", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  f <- .ctl_library(tmp)
  out_dir <- file.path(tmp, "8_controls")

  suppressMessages(run_barcoding_controls(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                          output_dir = out_dir, permutations = 5L, seed = 1L))

  p <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_cn1_permutations.csv"), stringsAsFactors = FALSE)
  v <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_cn1_verdict.csv"), stringsAsFactors = FALSE)

  expect_setequal(unique(p$permutacion), 1:5)
  expect_true(all(c("locus", "esquema", "metodo", "permutacion", "exactitud", "tasa_base") %in% names(p)))
  expect_true(all(c("media", "sd", "umbral", "hay_fuga", "tasa_base", "permutaciones") %in% names(v)))
  expect_equal(nrow(p), nrow(v) * 5L)
  # The folds of this fixture have no leakage
  expect_false(any(v$hay_fuga))
  # And not a single column carries the accuracy with the true labels: that is Phase 6
  expect_false(any(grepl("verdadera|real|honesta", names(v), ignore.case = TRUE)))
})

.ctl_dnabin <- function(x, nombre = "consulta") {
  ape::as.DNAbin(matrix(strsplit(tolower(x), "")[[1]], nrow = 1, dimnames = list(nombre, NULL)))
}

.ctl_rc <- function(x) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
}

.ctl_outgroup_hebras <- function(tmp, base) {
  # Three queries: one forward, the same one reverse complemented, and one that belongs to no
  # locus at all. The second is the one that matters: a control that gives a different answer to
  # the same sequence depending on how GenBank stored it is not measuring the classifier.
  og <- file.path(tmp, "outgroup_hebras")
  dir.create(og, showWarnings = FALSE, recursive = TRUE)
  set.seed(12L)
  s <- base
  for (p in sample(seq_len(300), 15)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1)
  ajena <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")
  x <- stats::setNames(c(s, .ctl_rc(s), ajena),
                       paste0(c("Portulaca_amilis", "Portulaca_amilis", "Talinum_paniculatum"),
                              "|matK_h", 1:3, ".1"))
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), file.path(og, "matK.fasta"))
  og
}

# Added on 2026-09-25, after the probe over rbcL. Four of the thirteen outgroup sequences,
# MH767707.1 to MH767710.1, are deposited on the opposite strand: counted over the files, they
# share 0.000 of their 12-mers with the library forward and 0.747 reverse complemented. They came
# out at distance 0.48 while their conspecifics sat at 0.011, which is orientation and not biology.
# MAFFT --add does not turn anything, and the branch already knew: .normalise_strand() exists for
# exactly this, and says in writing why --adjustdirection is not used. CN2 skipped it because the
# outgroup is read from step 1, which is raw. The library itself is clean, 0 reversed out of 5332.
#
# The rule is now shared instead of copied, so these tests come in two parts: one that pins what
# .normalise_strand() decides, so the extraction cannot change it, and the ones that orient a query
# against the library. The reference pool is built from the library alone, never from other
# queries, which is why .bc_strand_pool() takes the library and .bc_orient_to_library() takes the
# pool already built.

test_that("the strand rule keeps deciding what it decided before it was shared", {
  skip_if_not_installed("Biostrings")
  set.seed(11L)
  base <- paste(sample(c("A", "C", "G", "T"), 400, replace = TRUE), collapse = "")
  vary <- function(s, k) { for (p in sample(seq_len(400), k)) substr(s, p, p) <- sample(c("A","C","G","T"), 1); s }
  d_directa <- vary(base, 10)
  x <- c(a = base, b = vary(base, 8), c = vary(base, 12),
         d = .ctl_rc(d_directa),
         e = paste(sample(c("A", "C", "G", "T"), 400, replace = TRUE), collapse = ""))

  r <- .normalise_strand(Biostrings::DNAStringSet(x))

  expect_equal(r$log$Action,
               c("kept", "kept", "kept", "reverse_complemented", "no_match_either_direction"))
  expect_equal(r$log$Seq, names(x))
  # What it turns comes back in the direction of the marker, base by base
  expect_equal(as.character(r$dna[["d"]]), d_directa)
  # What matches neither direction is returned untouched, and only reported
  expect_equal(as.character(r$dna[["e"]]), unname(x["e"]))
})

test_that("a query is turned against the library, and one that matches no locus is only marked", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  f <- .ctl_library(tmp)
  lib <- ape::read.dna(file.path(f$library_dir, "LIB_matK.fasta"), format = "fasta", as.matrix = TRUE)
  pool <- .bc_strand_pool(lib)

  set.seed(13L)
  directa <- f$base
  ajena <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")

  d <- .bc_orient_to_library(.ctl_dnabin(directa), pool)
  r <- .bc_orient_to_library(.ctl_dnabin(.ctl_rc(directa)), pool)
  a <- .bc_orient_to_library(.ctl_dnabin(ajena), pool)

  expect_equal(d$orientacion, "directa")
  expect_equal(r$orientacion, "reversa")
  expect_equal(a$orientacion, "sin_coincidencia")
  # The turned one comes back as the forward one, base by base
  expect_identical(as.character(r$query), as.character(d$query))
  # The one that matches neither is returned as it arrived
  expect_identical(as.character(a$query), as.character(.ctl_dnabin(ajena)))
  # And the name it travels under does not change
  expect_equal(rownames(r$query), "consulta")
})

test_that("CN2 answers the same for a sequence and for its reverse complement, and says which it was", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  skip_if_not(.ctl_mafft_ok(), "MAFFT is not available")
  tmp <- withr::local_tempdir()
  f <- .ctl_library(tmp)
  og <- .ctl_outgroup_hebras(tmp, f$base)
  out_dir <- file.path(tmp, "8_controls")

  suppressMessages(run_barcoding_controls(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                          output_dir = out_dir, outgroup_dir = og,
                                          controls = "CN2", loci = "matK"))

  q <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_cn2_queries.csv"), stringsAsFactors = FALSE)
  expect_true("orientacion" %in% names(q))
  expect_equal(q$orientacion, c("directa", "reversa", "sin_coincidencia"))
  # The same sequence, stored either way, has to give the same answer. This is the test that the
  # real run of 2026-09-25 would have failed
  expect_equal(q$distancia_vecino[2], q$distancia_vecino[1])
  expect_equal(q$especie_predicha[2], q$especie_predicha[1])
  expect_equal(q$estado[2], q$estado[1])
})

test_that("the CN2 summary counts apart what was turned and what matches no locus", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  skip_if_not(.ctl_mafft_ok(), "MAFFT is not available")
  tmp <- withr::local_tempdir()
  f <- .ctl_library(tmp)
  og <- .ctl_outgroup_hebras(tmp, f$base)
  out_dir <- file.path(tmp, "8_controls")

  suppressMessages(run_barcoding_controls(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                          output_dir = out_dir, outgroup_dir = og,
                                          controls = "CN2", loci = "matK"))

  cn2 <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_cn2_outgroup.csv"), stringsAsFactors = FALSE)
  expect_true(all(c("comparables", "invertidas", "sin_coincidencia") %in% names(cn2)))
  expect_equal(cn2$consultas, 3L)
  expect_equal(cn2$invertidas, 1L)
  expect_equal(cn2$sin_coincidencia, 1L)
  # A sequence that matches the locus in neither direction is not a query of this control: it is
  # counted, it is published, and it stays out of the proportions and out of the species
  expect_equal(cn2$comparables, 2L)
  expect_equal(cn2$estado_1 + cn2$estado_2 + cn2$estado_3, 2L)
  expect_equal(cn2$especies, 1L)
})

# Added on 2026-09-25, after CN2 crashed on the real data with "DNA sequences in list not of the
# same length". The crash was the small half. .bc_cn2_queries() truncated both alignments to
# min(ncol) and compared them position by position, which assumes that column i of one alignment is
# homologous to column i of another. It is not. In the real run the library of matK had 573
# sequences of 2464 columns and the outgroup sequences ran from 326 to 1533 bases, so the comparison
# would have been between column 1 of the library and column 1 of a sequence that starts somewhere
# else. The fixture hid it by writing the outgroup with the same 300 bases as the library.
#
# The query is now added to the alignment of the library with MAFFT --add --keeplength, through the
# run_mafft() that the branch already uses in step 2. The columns of the library stay fixed, which is
# what makes the distance of one query comparable with the distance of another, and the insertions
# of the query are discarded. Decision of BMM, 2026-09-25: one call per query, so that no query is
# ever aligned in the company of another, which is how a real user query will arrive.

test_that("a query of a different length enters the library alignment without moving one column", {
  skip_if_not_installed("ape")
  skip_if_not(.ctl_mafft_ok(), "MAFFT is not available")
  tmp <- withr::local_tempdir()
  f <- .ctl_library(tmp)
  lib <- ape::read.dna(file.path(f$library_dir, "LIB_matK.fasta"), format = "fasta", as.matrix = TRUE)
  rownames(lib) <- .bc_parse_header(labels(lib))$sid

  # The query is the first sequence of the library with its first 40 bases removed. Its homology is
  # known, and any comparison by position displaces it by exactly 40 columns
  s <- as.character(lib)[1, ]
  query <- ape::as.DNAbin(matrix(s[s != "-"][-seq_len(40)], nrow = 1,
                                 dimnames = list("consulta", NULL)))

  al <- .bc_align_to_library(lib, query)

  expect_equal(ncol(al), ncol(lib))
  expect_equal(nrow(al), nrow(lib) + 1L)
  expect_true("consulta" %in% rownames(al))
  # The library comes back exactly as it went in, row by row and column by column
  expect_identical(as.character(al[rownames(lib), , drop = FALSE]), as.character(lib))
  # And the query lands on the columns it came from, not forty positions to the left
  d <- .bc_classifier_matrix(al, "raw", 100L)
  expect_lt(d["consulta", rownames(lib)[1]], 0.01)
})

test_that("the alignment of a query survives a directory whose name carries a space", {
  # run_mafft() splits its options by whitespace before pasting the input file, so a path with a
  # space in it is a real way to break this call and not a hypothetical one
  skip_if_not_installed("ape")
  skip_if_not(.ctl_mafft_ok(), "MAFFT is not available")
  tmp <- file.path(withr::local_tempdir(), "con espacio")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  f <- .ctl_library(tmp)
  lib <- ape::read.dna(file.path(f$library_dir, "LIB_matK.fasta"), format = "fasta", as.matrix = TRUE)
  rownames(lib) <- .bc_parse_header(labels(lib))$sid
  s <- as.character(lib)[2, ]
  query <- ape::as.DNAbin(matrix(s[s != "-"][-seq_len(25)], nrow = 1,
                                 dimnames = list("consulta", NULL)))

  al <- .bc_align_to_library(lib, query)
  expect_equal(ncol(al), ncol(lib))
  d <- .bc_classifier_matrix(al, "raw", 100L)
  expect_lt(d["consulta", rownames(lib)[2]], 0.01)
})

test_that("CN2 names MAFFT when it is missing, before CN1 spends a single permutation", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  f <- .ctl_library(tmp)
  og <- .ctl_outgroup(tmp, f$base)
  out_dir <- file.path(tmp, "8_controls")

  expect_error(
    suppressMessages(run_barcoding_controls(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                            output_dir = out_dir, outgroup_dir = og,
                                            permutations = 2L, mafft_exec = "mafft_que_no_existe")),
    "MAFFT")
  # The check belongs at the start: otherwise the run dies after an hour of permutations
  expect_false(file.exists(file.path(out_dir, "TABLE_barcoding_cn1_permutations.csv")))
})

# Rewritten on 2026-09-23, after the first version failed. The expectation was wrong, not the code.
# .bc_classify_nn() has no rule of remoteness by design (decision of Phase 5A): the three states come
# from the structure of the tie, and the threshold is swept later over the scores. So a query that
# resembles nothing still has a nearest neighbour, and when that neighbour belongs to a single
# species the answer is state 1. CN2 over the nearest neighbour therefore cannot return "not
# assignable" by distance, and asserting that it does would be adapting the test to a result nobody
# implemented.
#
# What CN2 can prove, and what the threshold of Phase 6 needs, is the other half of sec. 4 of the
# plan: whether the alien queries fall far enough from the library for any threshold to separate
# them. That is what this test asserts, and the table has to publish the distance to make it
# checkable.

test_that("CN2 publishes the distance of every alien query, and they fall beyond the library", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  skip_if_not(.ctl_mafft_ok(), "MAFFT is not available")
  tmp <- withr::local_tempdir()
  f <- .ctl_library(tmp)
  og <- .ctl_outgroup(tmp, f$base)
  out_dir <- file.path(tmp, "8_controls")

  # The fixture is unaligned and of three different lengths, which is the state the outgroup of
  # step 1 is really in. If this stops being true the test stops testing anything
  largos <- nchar(as.character(Biostrings::readDNAStringSet(file.path(og, "matK.fasta"))))
  expect_equal(length(unique(largos)), 3L)

  suppressMessages(run_barcoding_controls(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                          output_dir = out_dir, outgroup_dir = og, permutations = 2L))

  cn2 <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_cn2_outgroup.csv"), stringsAsFactors = FALSE)
  q <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_cn2_queries.csv"), stringsAsFactors = FALSE)

  # Every locus of the library appears, including the one with no outgroup sequence at all: a
  # control over zero sequences is not a control, and it has to say so rather than disappear
  expect_setequal(cn2$locus, c("matK", "rbcL"))
  expect_equal(cn2$consultas[cn2$locus == "matK"], 3L)
  expect_equal(cn2$consultas[cn2$locus == "rbcL"], 0L)
  expect_true(all(c("estado_1", "estado_2", "estado_3", "mediana_distancia_vecino",
                    "especies") %in% names(cn2)))

  # Finding of the real run, 2026-09-25: the outgroup of GenBank carries population level studies.
  # trnS-trnG came back with 83 sequences, all of Talinopsis frutescens, and 117 of the 235 ITS
  # sequences are that same species. Counting queries would make trnS-trnG look like the best
  # controlled locus when it is the worst: one species tested 83 times. So the table publishes the
  # species beside the queries, and the acta cites the species as the size of the control.
  expect_equal(cn2$especies[cn2$locus == "matK"], 3L)
  expect_equal(cn2$especies[cn2$locus == "rbcL"], 0L)

  # One row per alien query, with the distance that Phase 6 will sweep the threshold over
  expect_equal(nrow(q), 3L)
  expect_true(all(c("locus", "sid", "especie_consulta", "estado", "distancia_vecino") %in% names(q)))
  expect_setequal(q$especie_consulta,
                  c("Portulaca_amilis", "Talinum_paniculatum", "Anacampseros_filamentosa"))
  # All three arrive in the direction of the library, and the table says so
  expect_equal(q$orientacion, rep("directa", 3))

  # And the point of the control: the alien queries fall beyond anything the library holds inside a
  # species. If they did not, no threshold could ever separate them
  dna <- ape::read.dna(file.path(f$library_dir, "LIB_matK.fasta"), format = "fasta")
  h <- .bc_parse_header(labels(dna))
  rownames(dna) <- h$sid
  g <- .compute_barcode_gap(.bc_distance_matrix(dna, "raw"), stats::setNames(h$species, h$sid),
                            locus = "matK", model = "raw")
  expect_true(all(q$distancia_vecino > max(g$intra$distancia)))
})

# Added on 2026-09-25, after the real run of the three controls. CN3 measured resubstitution over
# the whole library, 693 sequences in ITS, while the honest accuracy of Phase 5A was measured over
# the queries of the folds, 488 in ITS. Two different query sets, so putting one figure beside the
# other is not the comparison CN3 promises to Phase 6. And the scheme was written as "species" in
# the code, so the genus scheme had no resubstitution counterpart at all.
#
# Decision of BMM, 2026-09-25: CN3 runs over the same queries as the folds, in both schemes, and
# the training set of each query is the one of its own fold plus the query itself. That way the
# only difference between the inflated figure and the honest one is the query being inside its
# training set, which is the whole point of the control.
#
# What is lost is worth writing down: the whole library version had an external check, its accuracy
# equalled one minus the proportion of zero interspecific distances of the gap table of Phase 4, in
# the twelve loci to machine precision. That check was run on 2026-09-25 and is recorded in the
# acta; the paired version cannot be checked that way.

test_that("CN3 evaluates the queries of the folds, in both schemes, and nothing else", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  f <- .ctl_library(tmp)
  out_dir <- file.path(tmp, "8_controls")

  suppressMessages(run_barcoding_controls(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                          output_dir = out_dir, controls = "CN3"))
  cn3 <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_cn3_resubstitution.csv"),
                         stringsAsFactors = FALSE)

  expect_setequal(cn3$esquema, c("species", "genus"))
  expect_true(all(cn3$etiqueta == "irreproducible_resustitucion"))

  # The number of queries is counted from the folds file, not from the function under test
  for (sc in c("species", "genus")) {
    tab <- utils::read.csv(file.path(f$folds_dir, paste0("TABLE_barcoding_folds_", sc, ".csv")),
                           stringsAsFactors = FALSE)
    for (l in unique(cn3$locus)) {
      expect_equal(cn3$consultas[cn3$locus == l & cn3$esquema == sc],
                   sum(tab$locus == l))
    }
  }
  # An earlier version of this test also asserted that the count never equals the size of the
  # library of the locus. That assertion was wrong and it was removed on 2026-09-25 after it failed:
  # in the genus scheme of this fixture the three species have a congener, so every one of the six
  # sequences is a query and coinciding with six is legitimate. What matters is already pinned
  # above, where the number of queries is counted from the folds file.
})

test_that("resubstitution hands the query the answer the honest run cannot give it", {
  # A species with a single sequence in the locus. With the query out of the training set nobody of
  # its species is left and the honest run cannot name it; with the query in, the nearest neighbour
  # is the query itself at distance zero and the answer is its own species. That gap is the bias
  # CN3 exists to measure, and it is why its figure never travels without its label.
  ids <- c("A1.1", "B1.1", "B2.1")
  m <- matrix(c(0, 0.20, 0.21,
                0.20, 0, 0.01,
                0.21, 0.01, 0), 3, 3, dimnames = list(ids, ids))
  sp <- stats::setNames(c("Opuntia_robusta", "Opuntia_stricta", "Opuntia_stricta"), ids)

  honesta <- .bc_classify_nn(m, c("B1.1", "B2.1"), "A1.1", sp)
  resust <- .bc_classify_nn(m, c("B1.1", "B2.1", "A1.1"), "A1.1", sp, allow_resubstitution = TRUE)

  expect_false(identical(honesta$especie_predicha, "Opuntia_robusta"))
  expect_equal(resust$especie_predicha, "Opuntia_robusta")
  expect_equal(resust$distancia_vecino, 0)
})

test_that("CN3 has to ask for the resubstitution out loud, and the guard holds without it", {
  ids <- c("A1.1", "A2.1", "B1.1")
  m <- matrix(c(0, 0.01, 0.20,
                0.01, 0, 0.21,
                0.20, 0.21, 0), 3, 3, dimnames = list(ids, ids))
  sp <- stats::setNames(c("Opuntia_robusta", "Opuntia_robusta", "Opuntia_stricta"), ids)

  expect_error(.bc_classify_nn(m, ids, "A1.1", sp), "training set")

  r <- .bc_classify_nn(m, ids, "A1.1", sp, allow_resubstitution = TRUE)
  expect_equal(r$estado, 1L)
  expect_equal(r$especie_predicha, "Opuntia_robusta")
  expect_equal(r$distancia_vecino, 0)
})

test_that("the assembly keeps the GenBank names when no checklist is given, and does not change with one", {
  nombres <- c("Portulaca amilis", "Talinum paniculatum", "Opuntia robusta")
  lista <- c("Opuntia robusta", "Opuntia stricta")

  sin_lista <- .bc_resolve_names(nombres, checklist_path = NULL)
  con_lista <- .bc_resolve_names(nombres, checklist_path = lista)

  # Without a checklist nothing is dropped: the names come back cleaned, with underscores
  expect_equal(sin_lista, c("Portulaca_amilis", "Talinum_paniculatum", "Opuntia_robusta"))
  # With a checklist the behaviour is the one of the closed phase: what does not match is NA
  expect_equal(con_lista, c(NA, NA, "Opuntia_robusta"))
  # And the default of the assembly is still the Cactaceae checklist, not NULL
  expect_null(formals(assemble_barcoding_dataset)$checklist_path)
  expect_match(deparse(body(assemble_barcoding_dataset)), "CactaceaeFullList", all = FALSE)
})

test_that("step 8 writes the three control tables, labels CN3, and names the step that is missing", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")
  tmp <- withr::local_tempdir()
  f <- .ctl_library(tmp)
  out_dir <- file.path(tmp, "8_controls")

  suppressMessages(run_barcoding_controls(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                          output_dir = out_dir, permutations = 2L))

  expect_true(file.exists(file.path(out_dir, "TABLE_barcoding_cn1_permutations.csv")))
  expect_true(file.exists(file.path(out_dir, "TABLE_barcoding_cn1_verdict.csv")))
  expect_true(file.exists(file.path(out_dir, "TABLE_barcoding_cn3_resubstitution.csv")))
  cn3 <- utils::read.csv(file.path(out_dir, "TABLE_barcoding_cn3_resubstitution.csv"), stringsAsFactors = FALSE)
  # CN3 carries its label in the table itself, so it cannot travel without it
  expect_true("etiqueta" %in% names(cn3))
  expect_true(all(cn3$etiqueta == "irreproducible_resustitucion"))

  expect_error(run_barcoding_controls(library_dir = f$library_dir, folds_dir = file.path(tmp, "none"),
                                      output_dir = out_dir),
               "build_barcoding_folds")
  expect_error(run_barcoding_controls(library_dir = f$library_dir, folds_dir = f$folds_dir,
                                      output_dir = file.path("4_Cleaned", "ctl")),
               "phylogeny")
})
