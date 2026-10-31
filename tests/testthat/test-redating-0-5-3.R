# Closure of 0.5.2 into 0.5.3, phase 4e (re-dating of 3h, BMM 30-10 and 31-10): Module 10 dates
# under three calibration schemes, each one fixed age on ACP_root and no calibration inside
# Cactaceae; scheme 1 (Ramirez-Barahona et al. 2020, 41.82 Ma) is the distributed chronogram. What
# the distributed files or the code can recompute is recomputed here: the fixed ages of the table of
# Tutorial 2 from calibrations_bounds.csv, the scheme 1 column and the root age from the distributed
# chronogram, the TreeAnnotator paths from MANIFEST.csv, the smoothing and the seed from
# PROVENANCE.csv, and the calls of Module 10 in the vignette and its script. The ages of schemes 2
# and 3 are not distributed; they are traced to the run of 3h in the tracing table of the closure.
# Written before the change.

.rd_root <- function(...) testthat::test_path("..", "..", ...)

.rd_lines <- function(...) {
  f <- .rd_root(...)
  skip_if(!file.exists(f), paste(file.path(...), "is not in this tree"))
  readLines(f, warn = FALSE, encoding = "UTF-8")
}

.rd_extdata <- function(name) {
  f <- system.file("extdata", name, package = "PhyloCactus")
  skip_if(!nzchar(f), paste(name, "not installed"))
  f
}

.rd_row <- function(lines, first) {
  r <- grep(paste0("^\\| ", first, " \\|"), lines, value = TRUE)
  expect_length(r, 1L)
  if (length(r) != 1L) return(rep(NA_character_, 4L))
  trimws(strsplit(sub("^\\|", "", sub("\\|\\s*$", "", r)), "\\|")[[1]])
}

.rd_t2 <- function() .rd_lines("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd")

# Code of the R chunks of the vignette, and of the script, as parsed expressions.
.rd_code <- function() {
  lines <- .rd_t2()
  inside <- FALSE; code <- character(0)
  for (l in lines) {
    if (!inside && grepl("^```\\{r", l)) { inside <- TRUE; next }
    if (inside && grepl("^```\\s*$", l)) { inside <- FALSE; next }
    if (inside) code <- c(code, l)
  }
  list(vignette = as.list(parse(text = code)),
       script = as.list(parse(file = .rd_root("inst", "scripts", "tutorial-2-cactus-phylogeny-inference.R"))))
}

.rd_calls <- function(exprs, fun) {
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

.rd_assign <- function(exprs, name) {
  for (e in exprs) if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(name))) return(e[[3]])
  NULL
}

test_that("the table of Tutorial 2 gives the fixed ages of the three schemes and the scheme 1 ages of the distributed chronogram", {
  skip_if_not_installed("ape")
  cal <- utils::read.csv(.rd_extdata("calibrations_bounds.csv"), stringsAsFactors = FALSE)
  fixed <- vapply(1:3, function(k) cal$min[grepl(paste0("_calibration_scheme_", k, "$"), cal$reference)], 0)
  v <- .rd_t2()
  head <- .rd_row(v, "Node")
  expect_identical(head[2:4], sprintf("Scheme %d (%s Ma)", 1:3, c(sprintf("%.2f", fixed[1:2]), sprintf("%.1f", fixed[3]))))

  t <- ape::read.tree(.rd_extdata("phylocactus_chronogram_ml.tree"))
  cc <- utils::read.csv(.rd_extdata("cactus_constraints.csv"), stringsAsFactors = FALSE)
  h <- ape::node.depth.edgelength(t); age <- max(h) - h
  crown <- function(col, val) age[ape::getMRCA(t, intersect(cc$Specie_name[cc[[col]] %in% val], t$tip.label))]
  want <- c(`Root (ACPT crown)` = max(h),
            `Cactaceae crown` = crown("Family", "Cactaceae"),
            `Cactoideae crown` = crown("Subfam", "Cactoideae"),
            `Cacteae crown` = crown("Clade", "Cacteae"),
            `Opuntioideae crown` = crown("Subfam", "Opuntioideae"))
  for (n in names(want)) {
    s1 <- as.numeric(sub(" .*", "", .rd_row(v, gsub("([()])", "\\\\\\1", n))[2]))
    expect_equal(s1, round(want[[n]], 2), info = n)
  }
  # The fixed node of the distributed chronogram carries the age of scheme 1.
  acp <- crown("Family", c("Cactaceae", "Anacampserotaceae", "Portulacaceae"))
  expect_equal(acp, fixed[1], tolerance = 1e-5)
  # The root age given beside ACPT_root is that of the distributed chronogram.
  m <- regmatches(paste(v, collapse = " "), regexec("its age, ([0-9.]+) Ma under scheme 1", paste(v, collapse = " ")))[[1]]
  expect_equal(as.numeric(m[2]), round(max(h), 2))
})

test_that("Module 10 dates the scheme chosen in the code, scheme 1 by default, the active row of calibrations_bounds.csv", {
  cal <- utils::read.csv(.rd_extdata("calibrations_bounds.csv"), stringsAsFactors = FALSE)
  active_scheme <- as.integer(sub(".*_calibration_scheme_([0-9]+)$", "\\1", cal$reference[cal$used_in_analysis]))
  for (src in .rd_code()) {
    expect_equal(eval(.rd_assign(src, "calibration_scheme")), active_scheme)
    dirs <- eval(.rd_assign(src, "scheme_dirs"))
    expect_length(dirs, 3L)
    # The folder of scheme 1 is the one the distributed chronogram was published from.
    man <- utils::read.csv(.rd_extdata("MANIFEST.csv"), stringsAsFactors = FALSE)
    src_ml <- man$source_path[man$canonical == "phylocactus_chronogram_ml.tree"]
    expect_identical(src_ml, file.path("8_Dating", dirs[active_scheme], "BestTree_treePL.tree"))
  }
})

test_that("the TreeAnnotator command of Tutorial 2 writes the summary MANIFEST.csv says was distributed", {
  man <- utils::read.csv(.rd_extdata("MANIFEST.csv"), stringsAsFactors = FALSE)
  hpd <- man$source_path[man$canonical == "phylocactus_chronogram_hpd.tree"]
  cmd <- grep("^treeannotator ", .rd_t2(), value = TRUE)
  expect_length(cmd, 1L)
  words <- strsplit(cmd, "\\s+")[[1]]
  expect_identical(words[length(words)], hpd)
  expect_identical(words[length(words) - 1L], file.path(dirname(hpd), "bsTree_treePL.tree"))
  # The node count the text gives is that of the distributed summary.
  skip_if_not_installed("ape")
  tr <- ape::read.nexus(.rd_extdata("phylocactus_chronogram_hpd.tree"))
  n <- as.integer(regmatches(paste(.rd_t2(), collapse = " "),
                             regexec("phylocactus_chronogram_hpd.tree`, ([0-9]+) nodes", paste(.rd_t2(), collapse = " ")))[[1]][2])
  expect_identical(n, ape::Ntip(tr) + ape::Nnode(tr))
})

test_that("PROVENANCE.csv records the scheme, smoothing and treePL seed of the chronogram, the ones Module 10 uses", {
  pv <- utils::read.csv(.rd_extdata("PROVENANCE.csv"), stringsAsFactors = FALSE)
  row <- pv[pv$canonical == "phylocactus_chronogram_ml.tree", ]
  expect_equal(nrow(row), 1L)
  for (k in c("calibration_scheme", "smoothing_treepl", "seed_treepl")) expect_true(k %in% names(pv), info = k)
  cal <- utils::read.csv(.rd_extdata("calibrations_bounds.csv"), stringsAsFactors = FALSE)
  expect_identical(row$calibration_scheme, cal$reference[cal$used_in_analysis])
  for (src in .rd_code()) {
    call <- .rd_calls(src, "automate_treePL")
    expect_length(call, 1L)
    expect_identical(as.character(eval(call[[1]]$seed)), as.character(row$seed_treepl))
  }
  expect_equal(as.numeric(row$smoothing_treepl), 1e-05)
})

test_that("the sensitivity step of Module 10 reads the uncalibrated crowns, in the vignette and the script", {
  for (src in .rd_code()) {
    call <- .rd_calls(src, "report_smoothing_sensitivity")
    expect_length(call, 1L)
    nodes <- call[[1]]$nodes
    expect_false(is.null(nodes))
    expect_identical(names(as.list(nodes))[-1],
                     c("root", "Cactaceae_crown", "Cactoideae_crown", "Cacteae_crown", "Opuntioideae_crown"))
  }
  code <- .rd_code()
  expect_identical(deparse(.rd_calls(code$vignette, "report_smoothing_sensitivity")[[1]]),
                   deparse(.rd_calls(code$script, "report_smoothing_sensitivity")[[1]]))
  expect_identical(deparse(.rd_calls(code$vignette, "automate_treePL")[[1]]),
                   deparse(.rd_calls(code$script, "automate_treePL")[[1]]))
})
