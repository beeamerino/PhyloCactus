# Tests of the remoteness threshold of the molecular diagnostic branch (11_barcoding/9_threshold/).
# Written before the code. Phase 6A of PhyloCactus 0.5.0, proposal of 2026-09-25, decisions D1 to D4.
#
# Without a rule of remoteness the nearest neighbour names a cactus species for anything: in CN2, 288
# of 312 alien queries came out in state 1. The rule sends a query whose distance to its nearest
# neighbour exceeds a threshold to state 3. The threshold of each locus is the quantile 0.99, type 1,
# of the distances of the legitimate queries of scheme G (a cactus whose species is missing from the
# library and whose genus is not), and it is fixed without looking at the alien queries, which only
# report afterwards how many it rejects.

.thr_pred <- function(esquema = "species") {
  # Legitimate queries of one locus by hand, with every category the plan separates
  data.frame(
    locus = "matK", esquema = esquema, pliegue = 1:6, estrato = NA_character_,
    sid = paste0("q", 1:6),
    especie_verdadera = c("Opuntia_robusta", "Opuntia_robusta", "Opuntia_robusta",
                          "Opuntia_stricta", "Cereus_jamacaru", "Cereus_jamacaru"),
    genero_verdadero = c("Opuntia", "Opuntia", "Opuntia", "Opuntia", "Cereus", "Cereus"),
    metodo = "nn", alineamiento = "add", orientacion = "directa",
    estado = c(1L, 1L, 2L, 3L, 2L, 3L),
    especie_predicha = c("Opuntia_robusta", "Opuntia_stricta", NA, NA, NA, NA),
    genero_predicho = c("Opuntia", "Opuntia", "Opuntia", NA, "Opuntia", NA),
    candidatas = NA_character_,
    distancia_vecino = c(0.01, 0.02, 0.03, 0.04, 0.05, NA),
    margen = 0, confianza = NA_real_,
    motivo = c(NA, NA, NA, "empate_entre_generos", NA, "sin_posiciones_comparables"),
    stringsAsFactors = FALSE
  )
}

.thr_alien <- function() {
  # Alien queries as CN2 writes them, one of them with no homology to the locus
  data.frame(
    locus = "matK", sid = paste0("o", 1:5),
    especie_consulta = c("Portulaca_amilis", "Portulaca_amilis", "Talinum_paniculatum",
                         "Anacampseros_filamentosa", "Talinopsis_frutescens"),
    orientacion = c("directa", "directa", "reversa", "directa", "sin_coincidencia"),
    estado = c(1L, 1L, 1L, 3L, 3L),
    especie_predicha = c("Opuntia_robusta", "Opuntia_stricta", "Cereus_jamacaru", NA, NA),
    genero_predicho = c("Opuntia", "Opuntia", "Cereus", NA, NA),
    candidatas = NA_character_,
    distancia_vecino = c(0.015, 0.06, 0.08, 0.09, NA),
    margen = 0, confianza = NA_real_,
    motivo = c(NA, NA, NA, "empate_entre_generos", "sin_coincidencia"),
    stringsAsFactors = FALSE
  )
}

.thr_dirs <- function(tmp, with_cn2 = TRUE) {
  cls <- file.path(tmp, "7_classifier")
  ctl <- file.path(tmp, "8_controls")
  dir.create(cls, recursive = TRUE, showWarnings = FALSE)
  dir.create(ctl, recursive = TRUE, showWarnings = FALSE)
  e <- .thr_pred("species")
  g <- .thr_pred("genus")
  g$distancia_vecino <- c(0.02, 0.03, 0.06, 0.07, 0.08, 0.10)
  g$motivo[6] <- NA
  g$estado[6] <- 2L
  g$genero_predicho[6] <- "Cereus"
  # A second locus with no alien query at all
  e2 <- e; e2$locus <- "rbcL"
  g2 <- g; g2$locus <- "rbcL"
  utils::write.csv(rbind(e, e2), file.path(cls, "TABLE_barcoding_predictions_species_nn_add.csv"), row.names = FALSE)
  utils::write.csv(rbind(g, g2), file.path(cls, "TABLE_barcoding_predictions_genus_nn_add.csv"), row.names = FALSE)
  if (with_cn2) {
    utils::write.csv(.thr_alien(), file.path(ctl, "TABLE_barcoding_cn2_queries.csv"), row.names = FALSE)
  }
  list(classifier_dir = cls, controls_dir = ctl, out = file.path(tmp, "9_threshold"))
}

test_that("the rule of remoteness sends only the distant queries to state 3 and keeps every earlier reason", {
  p <- .thr_pred()
  r <- .bc_apply_remoteness(p, 0.025)
  # At or below the threshold nothing changes
  expect_identical(r[1:2, ], p[1:2, ])
  # Beyond it: state 3, reason lejania, no species and no genus
  expect_equal(r$estado[3], 3L)
  expect_equal(r$motivo[3], "lejania")
  expect_true(is.na(r$genero_predicho[3]))
  expect_equal(r$estado[5], 3L)
  expect_equal(r$motivo[5], "lejania")
  # A query already in state 3 keeps its own reason, and one with no distance is left as it was
  expect_equal(r$motivo[4], "empate_entre_generos")
  expect_equal(r$motivo[6], "sin_posiciones_comparables")
  # The distance is never altered: the rule reads it, it does not rewrite it
  expect_identical(r$distancia_vecino, p$distancia_vecino)
})

test_that("the five categories are disjoint and add up to the queries, in both schemes", {
  for (sc in c("species", "genus")) {
    p <- .thr_pred(sc)
    cat5 <- .bc_outcome_category(p, sc)
    expect_length(cat5, nrow(p))
    expect_true(all(cat5 %in% c("especie_correcta", "especie_equivocada", "genero_correcto",
                                "genero_equivocado", "no_asignada")))
  }
  e <- .bc_outcome_category(.thr_pred("species"), "species")
  expect_equal(e, c("especie_correcta", "especie_equivocada", "genero_correcto", "no_asignada",
                    "genero_equivocado", "no_asignada"))
})

test_that("in scheme G no query is ever a correct species, even when it is forced to state 1", {
  p <- .thr_pred("genus")
  # The fixture names the true species in state 1; in scheme G that species is not in the training
  # set, so a species name there is always a false one
  g <- .bc_outcome_category(p, "genus")
  expect_false(any(g == "especie_correcta"))
  expect_equal(g[1:2], c("especie_equivocada", "especie_equivocada"))
})

test_that("above every distance the curve reproduces step 7 row by row, and below every distance nothing is assigned", {
  p <- .thr_pred()
  expect_identical(.bc_apply_remoteness(p, 1), p)
  bajo <- .bc_apply_remoteness(p, -1)
  expect_true(all(bajo$estado == 3L))
  expect_true(all(is.na(bajo$especie_predicha)))
})

test_that("the share of unassigned queries never falls as the threshold goes down", {
  cv <- .bc_threshold_curve(.thr_pred(), "species", ajenas = .thr_alien())
  cv <- cv[order(cv$umbral, decreasing = TRUE), ]
  expect_true(all(diff(cv$no_asignada) >= 0))
  expect_true(all(diff(cv$ajenas_rechazadas) >= 0))
  # One row per distance observed, legitimate or alien, plus zero
  obs <- c(.thr_pred()$distancia_vecino, .thr_alien()$distancia_vecino[.thr_alien()$orientacion != "sin_coincidencia"])
  expect_setequal(cv$umbral, sort(unique(c(0, obs[!is.na(obs)]))))
  # The five categories add up to the queries at every threshold
  expect_true(all(cv$especie_correcta + cv$especie_equivocada + cv$genero_correcto +
                    cv$genero_equivocado + cv$no_asignada == cv$consultas))
})

test_that("alien queries with no homology stay out of the curve, and a locus with only those is not measured", {
  a <- .thr_alien()
  cv <- .bc_threshold_curve(.thr_pred(), "species", ajenas = a)
  expect_true(all(cv$ajenas_comparables == 4L))
  # An alien query already in state 3 counts as rejected at any threshold
  expect_true(all(cv$ajenas_rechazadas >= 1L))

  solo <- a[a$orientacion == "sin_coincidencia", , drop = FALSE]
  cv0 <- .bc_threshold_curve(.thr_pred(), "species", ajenas = solo)
  expect_true(all(cv0$ajenas_comparables == 0L))
  expect_true(all(is.na(cv0$ajenas_rechazadas)))
})

test_that("the operating threshold is fixed without the alien queries", {
  skip_if_not_installed("withr")
  tmp1 <- withr::local_tempdir()
  tmp2 <- withr::local_tempdir()
  d1 <- .thr_dirs(tmp1, with_cn2 = TRUE)
  d2 <- .thr_dirs(tmp2, with_cn2 = FALSE)
  suppressMessages(sweep_barcoding_threshold(classifier_dir = d1$classifier_dir, controls_dir = d1$controls_dir,
                                             output_dir = d1$out, figures = FALSE))
  suppressMessages(sweep_barcoding_threshold(classifier_dir = d2$classifier_dir, controls_dir = d2$controls_dir,
                                             output_dir = d2$out, figures = FALSE))
  o1 <- utils::read.csv(file.path(d1$out, "TABLE_barcoding_threshold_operating.csv"), stringsAsFactors = FALSE)
  o2 <- utils::read.csv(file.path(d2$out, "TABLE_barcoding_threshold_operating.csv"), stringsAsFactors = FALSE)
  expect_equal(o1$umbral, o2$umbral)
  # Without CN2 the rejection is not measured anywhere, and it says so
  expect_true(all(o2$rechazo_estado == "no medido"))
  # With CN2, matK is measured and rbcL, which has no alien query, is not
  expect_equal(unique(o1$rechazo_estado[o1$locus == "matK"]), "medido")
  expect_equal(unique(o1$rechazo_estado[o1$locus == "rbcL"]), "no medido")
  expect_true(all(is.na(o1$rechazo[o1$locus == "rbcL"])))
})

test_that("the operating threshold is the declared quantile of scheme G, type 1, and the table says which", {
  g <- .thr_pred("genus")
  g$distancia_vecino <- c(0.02, 0.03, 0.06, 0.07, 0.08, 0.10)
  op <- .bc_operating_point(g, q = 0.5, type = 1L)
  # Type 1 returns a distance that exists in the data
  expect_equal(op, unname(stats::quantile(g$distancia_vecino, 0.5, type = 1)))
  expect_true(op %in% g$distancia_vecino)
  expect_error(.bc_operating_point(.thr_pred("species"), q = 0.99, type = 1L), "genus")

  skip_if_not_installed("withr")
  tmp <- withr::local_tempdir()
  d <- .thr_dirs(tmp)
  suppressMessages(sweep_barcoding_threshold(classifier_dir = d$classifier_dir, controls_dir = d$controls_dir,
                                             output_dir = d$out, q = 0.99, quantile_type = 1L, figures = FALSE))
  o <- utils::read.csv(file.path(d$out, "TABLE_barcoding_threshold_operating.csv"), stringsAsFactors = FALSE)
  expect_true(all(c("locus", "esquema", "umbral", "q", "tipo_cuantil", "consultas", "especie_correcta",
                    "especie_equivocada", "genero_correcto", "genero_equivocado", "no_asignada",
                    "ajenas_comparables", "ajenas_rechazadas", "rechazo", "rechazo_ic_inf",
                    "rechazo_ic_sup", "rechazo_estado") %in% names(o)))
  expect_true(all(o$q == 0.99))
  expect_true(all(o$tipo_cuantil == 1L))
  # One row per locus and scheme, and the threshold is the same in both schemes of a locus
  expect_equal(nrow(o), 4L)
  expect_equal(length(unique(o$umbral[o$locus == "matK"])), 1L)
})

test_that("the interval of the rejection is the exact binomial one", {
  for (xn in list(c(0L, 7L), c(7L, 7L), c(12L, 13L))) {
    ic <- .bc_clopper_pearson(xn[1], xn[2])
    ref <- stats::binom.test(xn[1], xn[2])$conf.int
    expect_equal(unname(ic), as.numeric(ref), tolerance = 1e-12)
  }
  expect_true(all(is.na(.bc_clopper_pearson(0L, 0L))))
})

test_that("the sweep writes its tables and figures, and names the step that is missing", {
  skip_if_not_installed("withr")
  tmp <- withr::local_tempdir()
  d <- .thr_dirs(tmp)
  suppressMessages(sweep_barcoding_threshold(classifier_dir = d$classifier_dir, controls_dir = d$controls_dir,
                                             output_dir = d$out, figures = TRUE))
  expect_true(file.exists(file.path(d$out, "TABLE_barcoding_threshold_curve.csv")))
  expect_true(file.exists(file.path(d$out, "TABLE_barcoding_threshold_operating.csv")))
  for (l in c("matK", "rbcL")) {
    expect_true(file.exists(file.path(d$out, paste0("FIG_barcoding_threshold_", l, ".pdf"))))
    expect_true(file.exists(file.path(d$out, paste0("FIG_barcoding_threshold_", l, ".png"))))
  }
  cv <- utils::read.csv(file.path(d$out, "TABLE_barcoding_threshold_curve.csv"), stringsAsFactors = FALSE)
  expect_setequal(unique(paste(cv$locus, cv$esquema)),
                  c("matK species", "matK genus", "rbcL species", "rbcL genus"))

  expect_error(sweep_barcoding_threshold(classifier_dir = file.path(tmp, "none"),
                                         controls_dir = d$controls_dir, output_dir = d$out),
               "classify_barcoding_folds")
  expect_error(sweep_barcoding_threshold(classifier_dir = d$classifier_dir, controls_dir = d$controls_dir,
                                         output_dir = file.path("4_Cleaned", "thr")),
               "phylogeny")
})
