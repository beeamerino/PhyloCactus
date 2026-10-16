# Decisions DB6 (V5-a) and DB7 of BMM (15-10, Phase 12, D-b part 1b): the vignette of Tutorial 5 is
# a step-by-step tutorial like Tutorials 1 and 2: a Setup with the switches of the script, one section
# per step (1 to 11) with the call the script runs, the same arguments and values, the blocks that
# read the IdTaxa tables of the cluster, and a section on running IdTaxa on a computing cluster. The
# tests of D-a read only the script; these read the chunks of the vignette. Written before the text.

.t5v_exprs_script <- function() {
  f <- system.file("scripts", "tutorial-5-cactus-phylogeny-barcoding.R", package = "PhyloCactus")
  skip_if(!nzchar(f), "Tutorial 5 script not installed")
  txt <- readLines(f, warn = FALSE)
  from <- grep("^# Step 1:", txt)
  to <- grep("^# Results", txt)
  expect_length(from, 1L); expect_length(to, 1L)
  ex <- parse(f, keep.source = TRUE)
  first <- vapply(attr(ex, "srcref"), function(s) s[1], integer(1))
  as.list(ex[first > from & first < to])
}

.t5v_exprs_vignette <- function() {
  f <- system.file("doc", "tutorial-5-cactus-phylogeny-barcoding.Rmd", package = "PhyloCactus")
  if (!nzchar(f)) f <- testthat::test_path("..", "..", "vignettes", "tutorial-5-cactus-phylogeny-barcoding.Rmd")
  skip_if(!file.exists(f), "Tutorial 5 vignette not found")
  lines <- readLines(f, warn = FALSE)
  inside <- FALSE; code <- character(0)
  for (l in lines) {
    if (!inside && grepl("^```\\{r", l)) { inside <- TRUE; next }
    if (inside && grepl("^```\\s*$", l)) { inside <- FALSE; next }
    if (inside) code <- c(code, l)
  }
  list(lines = lines, exprs = as.list(parse(text = code, keep.source = FALSE)))
}

.t5v_calls <- function(exprs, funs) {
  found <- list()
  walk <- function(e) {
    if (is.call(e)) {
      if (is.name(e[[1]]) && as.character(e[[1]]) %in% funs) found[[length(found) + 1L]] <<- e
      args <- as.list(e)[-1]
      for (i in seq_along(args)) if (is.call(args[[i]])) walk(args[[i]])
    }
  }
  for (e in exprs) walk(e)
  found
}

.t5v_exported <- function() getNamespaceExports("PhyloCactus")

test_that("the vignette of Tutorial 5 has a Setup and one section per step, 1 to 11 (DB6)", {
  v <- .t5v_exprs_vignette()$lines
  expect_true(any(grepl("^### Setup", v)))
  for (k in 1:11) expect_true(any(grepl(paste0("^### Step ", k, ": "), v)), info = paste("Step", k))
})

test_that("every call of the steps of the script is in the vignette with the same arguments (DB6, DB7)", {
  s <- .t5v_calls(.t5v_exprs_script(), .t5v_exported())
  expect_gte(length(s), 13L)
  v <- .t5v_calls(.t5v_exprs_vignette()$exprs, .t5v_exported())
  used <- rep(FALSE, length(v))
  for (e in s) {
    hit <- which(!used & vapply(v, function(x) identical(x, e), logical(1)))
    expect_true(length(hit) >= 1L, info = paste(deparse(e), collapse = " "))
    if (length(hit)) used[hit[1]] <- TRUE
  }
})

test_that("the Setup of the vignette sets the three switches of the script, and the IdTaxa blocks are shown (DB6)", {
  f <- system.file("scripts", "tutorial-5-cactus-phylogeny-barcoding.R", package = "PhyloCactus")
  skip_if(!nzchar(f), "Tutorial 5 script not installed")
  sw <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
                 as.character(e[[2]]) %in% c("notify_email", "idtaxa_tables", "n_workers"), as.list(parse(f)))
  expect_length(sw, 3L)
  vx <- .t5v_exprs_vignette()$exprs
  for (e in sw) expect_true(any(vapply(vx, function(x) identical(x, e), logical(1))), info = as.character(e[[2]]))
  is_idtaxa_if <- function(x) is.call(x) && identical(x[[1]], as.name("if")) &&
    identical(x[[2]], quote(!is.null(idtaxa_tables)))
  n_script <- sum(vapply(.t5v_exprs_script(), is_idtaxa_if, logical(1)))
  expect_gte(n_script, 3L)
  expect_gte(sum(vapply(vx, is_idtaxa_if, logical(1))), n_script)
})

test_that("the vignette shows how to run IdTaxa on a computing cluster (DB6)", {
  vg <- .t5v_exprs_vignette()
  g <- .t5v_calls(vg$exprs, "generate_barcoding_job_scripts")
  expect_gte(length(g), 1L)
  ok_names <- names(formals(PhyloCactus::generate_barcoding_job_scripts))
  for (e in g) {
    a <- names(as.list(e))[-1]
    expect_true(all(nzchar(a)), info = "every argument named")
    expect_true(all(a %in% ok_names), info = paste(setdiff(a, ok_names), collapse = ", "))
    expect_identical(e$method, "idtaxa")
  }
  expect_true(any(grepl("merge_barcoding_chunks", vg$lines, fixed = TRUE)))
})
