# Tests of Tutorial 5 as Phase 11 runs it (decisions V1 and V11 of BMM, 29-09). Written before the
# changes. IdTaxa is run once, on Leftraru, on the frozen library; the tutorial reads its merged
# tables through the switch `idtaxa_tables` and never trains IdTaxa itself. The metrics of 6D are a
# step of the tutorial.

.t5_calls <- function() {
  f <- system.file("scripts", "tutorial-5-cactus-phylogeny-barcoding.R", package = "PhyloCactus")
  skip_if(!nzchar(f), "Tutorial 5 not installed")
  list(file = f, text = readLines(f, warn = FALSE), exprs = as.list(parse(f)))
}
.t5_find <- function(exprs, fun) {
  out <- list()
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(e[[1]], as.name(fun))) out[[length(out) + 1L]] <<- e
      for (a in as.list(e)[-1]) walk(a)
    }
  }
  for (e in exprs) walk(e)
  out
}

test_that("Tutorial 5 declares idtaxa_tables as NULL and no longer samples IdTaxa folds (V1, V11)", {
  t5 <- .t5_calls()
  assigns <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")), t5$exprs)
  lhs <- vapply(assigns, function(e) deparse(e[[2]]), "")
  expect_true("idtaxa_tables" %in% lhs)
  expect_null(eval(assigns[[which(lhs == "idtaxa_tables")]][[3]]))
  expect_false("idtaxa_max_folds" %in% lhs)
})

test_that("Tutorial 5 never trains IdTaxa: no classifier nor control is called with method idtaxa (V11)", {
  t5 <- .t5_calls()
  for (fun in c("classify_barcoding_folds", "run_barcoding_controls")) {
    for (e in .t5_find(t5$exprs, fun)) expect_false(identical(e$method, "idtaxa"))
  }
})

test_that("with idtaxa_tables set, Tutorial 5 copies the merged tables and sweeps and measures IdTaxa (V11)", {
  t5 <- .t5_calls()
  txt <- paste(t5$text, collapse = "\n")
  expect_match(txt, "if (!is.null(idtaxa_tables))", fixed = TRUE)
  expect_match(txt, "TABLE_barcoding_cn2_queries_idtaxa.csv", fixed = TRUE)
  sweeps <- .t5_find(t5$exprs, "sweep_barcoding_threshold")
  expect_true(any(vapply(sweeps, function(e) identical(e$method, "idtaxa"), logical(1))))
})

test_that("Tutorial 5 runs the metrics of 6D (V1)", {
  t5 <- .t5_calls()
  expect_gte(length(.t5_find(t5$exprs, "summarise_barcoding_metrics")), 1L)
})

test_that("Tutorial 5 no longer quotes the running times of the library of the 0.4.5 mining (V1)", {
  txt <- paste(.t5_calls()$text, collapse = "\n")
  for (s in c("9 h 11 min", "34 h", "29 475 s", "5928 folds")) expect_false(grepl(s, txt, fixed = TRUE))
})

test_that("the comment of step 9 says that the IdTaxa sweep writes the species threshold of each locus (R-c, C6)", {
  tut <- system.file("scripts", "tutorial-5-cactus-phylogeny-barcoding.R", package = "PhyloCactus")
  skip_if(!nzchar(tut), "Tutorial 5 is not installed")
  txt <- paste(readLines(tut, warn = FALSE), collapse = "\n")
  expect_true(grepl("TABLE_barcoding_species_threshold_idtaxa.csv", txt, fixed = TRUE))
})
