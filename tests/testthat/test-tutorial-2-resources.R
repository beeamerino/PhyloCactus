# Decisions T1 to T3 of BMM (14-10, Phase 12): the cluster resources of Tutorial 2 fitted to the matrix
# of the folder 0 of 29-09 (1208 analysed terminals, 13 700 sites). The ML search of Phase 12 used
# 15.5 GiB of 16G and 1 h 15 min; at that scale a bootstrap task of 500 replicates would pass 24 h.
# T1: ML search with 24G. T2: bootstraps in 4 tasks of 250 replicates, 16G, 24 h. T3: temporal
# bootstraps with 24G. The script and the vignette carry the same arguments. Written before the change.

.t2_sources <- function() {
  script <- system.file("scripts", "tutorial-2-cactus-phylogeny-inference.R", package = "PhyloCactus")
  vig <- system.file("doc", "tutorial-2-cactus-phylogeny-inference.Rmd", package = "PhyloCactus")
  if (!nzchar(vig)) vig <- testthat::test_path("..", "..", "vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd")
  out <- list()
  if (nzchar(script)) out$script <- as.list(parse(script))
  if (file.exists(vig)) {
    lines <- readLines(vig, warn = FALSE)
    inside <- FALSE; code <- character(0)
    for (l in lines) {
      if (!inside && grepl("^```\\{r", l)) { inside <- TRUE; next }
      if (inside && grepl("^```\\s*$", l)) { inside <- FALSE; next }
      if (inside) code <- c(code, l)
    }
    out$vignette <- as.list(parse(text = code))
  }
  skip_if(length(out) < 2L, "Tutorial 2 script or vignette not found")
  out
}

.t2_find <- function(exprs, fun) {
  found <- list()
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(e[[1]], as.name(fun))) found[[length(found) + 1L]] <<- e
      args <- as.list(e)[-1]
      for (i in seq_along(args)) if (is.call(args[[i]])) walk(args[[i]])
    }
  }
  for (e in exprs) walk(e)
  found
}

test_that("the ML search script asks for 24G and 2 h (T1)", {
  for (src in names(s <- .t2_sources())) {
    calls <- .t2_find(s[[src]], "generate_ml_search_script")
    expect_length(calls, 1L)
    expect_identical(calls[[1]]$cluster_mem, "24G", info = src)
    expect_identical(calls[[1]]$cluster_time, "02:00:00", info = src)
  }
})

test_that("the bootstraps run as 4 tasks of 250 replicates, 16G and 24 h each (T2)", {
  for (src in names(s <- .t2_sources())) {
    calls <- .t2_find(s[[src]], "generate_bootstrap_script")
    expect_length(calls, 1L)
    expect_equal(calls[[1]]$bs_per_rep, 250, info = src)
    expect_equal(calls[[1]]$max_reps, 4, info = src)
    expect_identical(calls[[1]]$cluster_mem, "16G", info = src)
    expect_identical(calls[[1]]$cluster_time, "24:00:00", info = src)
  }
})

test_that("the temporal bootstraps ask for 24G and 4 h (T3)", {
  for (src in names(s <- .t2_sources())) {
    calls <- .t2_find(s[[src]], "generate_temporal_bootstrap_script")
    expect_length(calls, 1L)
    expect_identical(calls[[1]]$cluster_mem, "24G", info = src)
    expect_identical(calls[[1]]$cluster_time, "04:00:00", info = src)
  }
})

test_that("the resource comments no longer quote the 0.4.5 measurements (T1, T2)", {
  script <- system.file("scripts", "tutorial-2-cactus-phylogeny-inference.R", package = "PhyloCactus")
  skip_if(!nzchar(script), "Tutorial 2 not installed")
  txt <- paste(readLines(script, warn = FALSE), collapse = "\n")
  for (s in c("~11 GB max used", "~41 min real runtime", "~5.3 GB max used", "~17 h real runtime", "2 x 500 = 1000")) {
    expect_false(grepl(s, txt, fixed = TRUE), info = s)
  }
})
