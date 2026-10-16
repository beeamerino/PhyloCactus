# Decisions V1 to V4 and T6-1 to T6-7 of BMM (14-10 and 15-10, Phase 12): a new Tutorial 6 identifies
# public samples of each kind of input with the release library, each example downloaded in a block of
# its own and identified with a path to its file, as a user writes it; its figures are drawn from tables
# of extdata by the two plotting functions; the function reference becomes Tutorial 7, and the old
# address is redirected; the rates of identification move from Tutorial 5 to Tutorial 6. Examples (V3,
# 15-10): LT745702.1, PZ714536.1, the specimen of Coryphantha elephantidens subsp. bumamma (OK340410.1,
# OK340349.1, OK340287.1, OK340224.1), the specimen of Xiquexique gounellei (OP189074.1 and its loci),
# MN517611.1, MT387452.1, a sample sheet, GCA_024363205.1 and SRR26398994. Written before the tutorial.

.t6_name <- "tutorial-6-cactus-barcoding-identification"
.t6_examples <- c("LT745702.1", "PZ714536.1", "OK340410.1", "OP189074.1", "MN517611.1", "MT387452.1",
                  "GCA_024363205.1", "SRR26398994")

.t6_script <- function() {
  f <- system.file("scripts", paste0(.t6_name, ".R"), package = "PhyloCactus")
  skip_if(!nzchar(f), "Tutorial 6 script not installed")
  list(text = readLines(f, warn = FALSE), exprs = as.list(parse(f)))
}

.t6_vignette_file <- function(name) {
  f <- system.file("doc", paste0(name, ".Rmd"), package = "PhyloCactus")
  if (!nzchar(f)) f <- testthat::test_path("..", "..", "vignettes", paste0(name, ".Rmd"))
  f
}

.t6_vignette <- function() {
  f <- .t6_vignette_file(.t6_name)
  skip_if(!file.exists(f), "Tutorial 6 vignette not found")
  readLines(f, warn = FALSE)
}

.t6_find <- function(exprs, fun) {
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

.t6_root <- function(file) {
  f <- testthat::test_path("..", "..", file)
  skip_if(!file.exists(f), paste(file, "is not in this tree"))
  f
}

test_that("Tutorial 6 has a script and a vignette (V1)", {
  expect_true(nzchar(system.file("scripts", paste0(.t6_name, ".R"), package = "PhyloCactus")))
  expect_true(file.exists(.t6_vignette_file(.t6_name)))
})

test_that("the script identifies each kind of input with the function of the package (T6-1)", {
  s <- .t6_script()
  expect_gte(length(.t6_find(s$exprs, "identify_barcoding_query")), 4L)
  expect_gte(length(.t6_find(s$exprs, "identify_barcoding_samples")), 1L)
  expect_gte(length(.t6_find(s$exprs, "identify_barcoding_assembly")), 1L)
  expect_gte(length(.t6_find(s$exprs, "identify_barcoding_reads")), 1L)
  # Specimens of several loci are one sample each
  q <- .t6_find(s$exprs, "identify_barcoding_query")
  expect_gte(sum(vapply(q, function(e) !is.null(e$sample), logical(1))), 2L)
})

test_that("no identification of the script reads its file from the package: a user's path is shown (T6-1)", {
  s <- .t6_script()
  for (fun in c("identify_barcoding_query", "identify_barcoding_samples", "identify_barcoding_assembly",
                "identify_barcoding_reads")) {
    for (e in .t6_find(s$exprs, fun)) {
      expect_false(any(grepl("system.file", deparse(e), fixed = TRUE)), info = paste(deparse(e), collapse = " "))
    }
  }
})

test_that("every example is downloaded from its public source in the script (T6-1, C-a, V3)", {
  s <- .t6_script()
  txt <- paste(s$text, collapse = "\n")
  for (a in .t6_examples) expect_true(grepl(a, txt, fixed = TRUE), info = a)
  # rentrez is in Suggests, so the script calls it as rentrez::entrez_fetch()
  expect_true(grepl("rentrez::entrez_fetch(", txt, fixed = TRUE))
  expect_true(grepl("ebi.ac.uk/ena", txt, fixed = TRUE))
  expect_true(grepl("ftp.ncbi.nlm.nih.gov/genomes", txt, fixed = TRUE))
})

test_that("the script and the vignette hold no path of a local machine", {
  s <- .t6_script()
  expect_false(any(grepl("/Users/|/home/", s$text)))
  expect_false(any(grepl("/Users/|/home/", .t6_vignette())))
})

test_that("the vignette draws its figures with the two plotting functions from tables of extdata (F2, F4, T6-5)", {
  v <- paste(.t6_vignette(), collapse = "\n")
  expect_true(grepl("plot_barcoding_identification(", v, fixed = TRUE))
  expect_true(grepl("plot_barcoding_identification_summary(", v, fixed = TRUE))
  used <- unique(regmatches(v, gregexpr("barcoding_example_[A-Za-z0-9_.]+_(answer|sampling)\\.csv", v))[[1]])
  expect_gte(length(used), 2L)
  m <- utils::read.csv(system.file("extdata", "MANIFEST.csv", package = "PhyloCactus"), stringsAsFactors = FALSE)
  expect_true(all(used %in% m$canonical), info = paste(setdiff(used, m$canonical), collapse = ", "))
})

test_that("the vignette links every example to its full report on the website (F3)", {
  v <- paste(.t6_vignette(), collapse = "\n")
  for (a in .t6_examples) expect_true(grepl(paste0("reports/", a), v, fixed = TRUE), info = a)
})

test_that("the vignette states that no Sanger record from outside the library was named at species rank (option A)", {
  v <- paste(.t6_vignette(), collapse = "\n")
  expect_true(grepl("94", v, fixed = TRUE))
  expect_true(grepl("below_species_threshold", v, fixed = TRUE))
})

test_that("the rates of identification are in Tutorial 6, and Tutorial 5 links to it (T6-7)", {
  v6 <- paste(.t6_vignette(), collapse = "\n")
  expect_true(grepl("## How Well This Version Identifies", v6, fixed = TRUE))
  f5 <- .t6_vignette_file("tutorial-5-cactus-phylogeny-barcoding")
  skip_if(!file.exists(f5), "Tutorial 5 vignette not found")
  v5 <- paste(readLines(f5, warn = FALSE), collapse = "\n")
  expect_false(grepl("## How Well This Version Identifies", v5, fixed = TRUE))
  expect_true(grepl(paste0(.t6_name, ".html"), v5, fixed = TRUE))
})

test_that("the function reference is Tutorial 7 and its old address is redirected (V1, T6-7)", {
  expect_true(file.exists(.t6_vignette_file("tutorial-7-cactus-phylogeny-functions")))
  expect_false(file.exists(testthat::test_path("..", "..", "vignettes", "tutorial-6-cactus-phylogeny-functions.Rmd")))
  y <- paste(readLines(.t6_root("_pkgdown.yml"), warn = FALSE), collapse = "\n")
  expect_true(grepl("redirects:", y, fixed = TRUE))
  expect_true(grepl("articles/tutorial-6-cactus-phylogeny-functions.html", y, fixed = TRUE))
  expect_true(grepl("articles/tutorial-7-cactus-phylogeny-functions.html", y, fixed = TRUE))
  expect_true(grepl(paste0("articles/", .t6_name, ".html"), y, fixed = TRUE))
})

test_that("README and Get Started link the two tutorials at their new addresses (T6-7)", {
  for (f in c("README.md", file.path("vignettes", "PhyloCactus.Rmd"))) {
    x <- paste(readLines(.t6_root(f), warn = FALSE), collapse = "\n")
    expect_true(grepl("tutorial-7-cactus-phylogeny-functions.html", x, fixed = TRUE), info = f)
    expect_true(grepl(paste0(.t6_name, ".html"), x, fixed = TRUE), info = f)
    expect_false(grepl("tutorial-6-cactus-phylogeny-functions.html", x, fixed = TRUE), info = f)
  }
})

# Decision E7-a of BMM (15-10): Example 7 had no figure and no result; the index of its sample sheet
# is distributed and shown as a table, and the text says that the figures of its three samples are
# those of Examples 1, 2 and 5. Written before the vignette is changed.
test_that("Example 7 shows the index of its sample sheet from extdata, as a table (E7-a)", {
  v <- .t6_vignette()
  i <- grep("^## Example 7", v); j <- grep("^## Example 8", v)
  expect_length(i, 1L); expect_length(j, 1L)
  e7 <- paste(v[i:(j - 1L)], collapse = "\n")
  expect_true(grepl("barcoding_example_samples_index.csv", e7, fixed = TRUE))
  expect_true(grepl("kable(", e7, fixed = TRUE))
  expect_true(grepl("Examples 1, 2 and 5", e7, fixed = TRUE))
  m <- utils::read.csv(system.file("extdata", "MANIFEST.csv", package = "PhyloCactus"), stringsAsFactors = FALSE)
  expect_true("barcoding_example_samples_index.csv" %in% m$canonical)
})
