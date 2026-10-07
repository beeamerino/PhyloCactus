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
      # An empty argument, as in x[i, ], is the missing symbol: only calls are walked into
      args <- as.list(e)[-1]
      for (i in seq_along(args)) if (is.call(args[[i]])) walk(args[[i]])
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

# Phase D (decisions D3 and D4 of BMM, 07-10): step 11 identifies the example plastome of extdata,
# the routes that need the user's own data are shown but not run, and the script draws the figures
# TB1 to TB7 and TB9 into 11_barcoding/figures/ (TB8 is drawn once, outside the script).

test_that("step 11 identifies the example plastome distributed in extdata (D3, D5)", {
  t5 <- .t5_calls()
  calls <- .t5_find(t5$exprs, "identify_barcoding_query")
  expect_gte(length(calls), 1L)
  txt <- paste(vapply(calls, function(e) paste(deparse(e), collapse = ""), ""), collapse = "\n")
  expect_match(txt, "barcoding_example_plastome_MN517611.1.gb", fixed = TRUE)
  expect_match(txt, "system.file", fixed = TRUE)
})

test_that("the routes that need the user's data are shown as comments, not run (D3)", {
  t5 <- .t5_calls()
  for (fun in c("identify_barcoding_samples", "identify_barcoding_assembly", "identify_barcoding_reads")) {
    expect_length(.t5_find(t5$exprs, fun), 0L)
    expect_true(any(grepl(paste0("^#.*", fun, "\\("), t5$text)))
  }
})

test_that("the script draws TB1 to TB7 and TB9 into 11_barcoding/figures (D4)", {
  t5 <- .t5_calls()
  txt <- paste(t5$text, collapse = "\n")
  expect_match(txt, "11_barcoding/figures", fixed = TRUE)
  for (k in c(1:7, 9)) expect_match(txt, sprintf("Figure_TB%d_", k), fixed = TRUE)
  expect_false(grepl("Figure_TB8_", txt, fixed = TRUE))
})
