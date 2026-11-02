# Decision P11 of BMM (13-10, Phase 12; backlog request of 29-09): the length limits of the mining
# (100 and 5000 bases, decision L1 of 28-09) are arguments of the three mining functions, passed to
# the shared phylotaR call. Defaults unchanged. Written before the code.

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
