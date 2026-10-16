# Decision P11 of BMM (13-10, Phase 12; backlog request of 29-09): the length limits of the mining
# (100 and 5000 bases, decision L1 of 28-09) are arguments of the three mining functions, passed to
# the shared phylotaR call, and step 1 of Tutorial 1 shows them, with the email notification, as
# explicit arguments. Defaults unchanged. Written before the code.

test_that("the three mining functions take mnsql and mxsql, 100 and 5000 by default (P11)", {
  for (fun in c("assemble_ingroup_phylotar", "assemble_outgroup_phylotar", "assemble_barcoding_dataset")) {
    f <- formals(get(fun, envir = asNamespace("PhyloCactus")))
    expect_true(all(c("mnsql", "mxsql") %in% names(f)), info = fun)
    expect_equal(eval(f$mnsql), 100L, info = fun)
    expect_equal(eval(f$mxsql), 5000L, info = fun)
  }
})

test_that("each mining function passes its limits to the shared phylotaR call (P11)", {
  find_call <- function(e) {
    if (!is.call(e)) return(NULL)
    if (identical(e[[1]], as.name(".phylotar_load_or_mine"))) return(e)
    args <- as.list(e)[-1]
    for (i in seq_along(args)) if (is.call(args[[i]])) { r <- find_call(args[[i]]); if (!is.null(r)) return(r) }
    NULL
  }
  for (fun in c("assemble_ingroup_phylotar", "assemble_outgroup_phylotar", "assemble_barcoding_dataset")) {
    cl <- find_call(body(get(fun, envir = asNamespace("PhyloCactus"))))
    expect_false(is.null(cl), info = fun)
    expect_identical(cl$mnsql, as.name("mnsql"), info = fun)
    expect_identical(cl$mxsql, as.name("mxsql"), info = fun)
  }
})

test_that("step 1 of Tutorial 1 shows the notification and the length limits as explicit arguments (P11)", {
  f <- system.file("scripts", "tutorial-1-cactus-phylogeny-prep.R", package = "PhyloCactus")
  skip_if(!nzchar(f), "Tutorial 1 not installed")
  exprs <- as.list(parse(f))
  found <- list()
  walk <- function(e) {
    if (is.call(e)) {
      if (as.character(e[[1]])[1] %in% c("assemble_ingroup_phylotar", "assemble_outgroup_phylotar")) {
        found[[length(found) + 1L]] <<- e
      }
      args <- as.list(e)[-1]
      for (i in seq_along(args)) if (is.call(args[[i]])) walk(args[[i]])
    }
  }
  for (e in exprs) walk(e)
  expect_length(found, 2L)
  for (e in found) expect_true(all(c("notify", "mnsql", "mxsql") %in% names(e)))
})
