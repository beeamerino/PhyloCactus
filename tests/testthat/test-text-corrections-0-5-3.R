# Closure of 0.5.2 into 0.5.3, phase 4a (audit of 26-10-2026, findings M5, M7, M8, M10, M13, M14,
# m1 to m9, n1): every number the text states that inst/extdata can recompute is recomputed here,
# the numbers that Tutorial 1 and its script both state are the same, help pages carry no number of
# an earlier run, the intervals of the dated replicates are called HPD only where TreeAnnotator names
# its column, the references are the ones the publishers record, and Tutorial 2 passes the base
# seeds of the reference run. Written before the text.

.t53_root <- function(...) testthat::test_path("..", "..", ...)

.t53_read <- function(...) {
  f <- .t53_root(...)
  skip_if(!file.exists(f), paste(file.path(...), "is not in this tree"))
  paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

.t53_flat <- function(x) gsub("\\s+", " ", x)

# A script with the comment marks removed, so that a sentence broken over comment lines reads whole.
.t53_script <- function(name) .t53_flat(gsub("(^|\n)\\s*#+ ?", "\\1", .t53_read("inst", "scripts", name)))

.t53_extdata <- function(name) {
  f <- system.file("extdata", name, package = "PhyloCactus")
  skip_if(!nzchar(f), paste(name, "not installed"))
  f
}

# The groups a pattern captures in a text, as character, failing when the pattern is absent.
.t53_capture <- function(text, pattern, info = pattern) {
  m <- regmatches(text, regexec(pattern, text))[[1]]
  expect_true(length(m) > 1L, info = info)
  if (length(m) < 2L) return(rep(NA_character_, 8L))
  gsub(",", "", m[-1])
}

# The text of a call, from its name to the parenthesis that closes it; of the occurrences of the
# name, the longest, so that a mention in prose ("automate_treePL()") is not taken for the call.
.t53_call <- function(text, fun) {
  starts <- gregexpr(paste0(fun, "("), text, fixed = TRUE)[[1]]
  expect_true(starts[1] > 0L, info = fun)
  if (starts[1] < 0L) return("")
  calls <- vapply(starts, function(start) {
    chars <- strsplit(substring(text, start + nchar(fun)), "")[[1]]
    depth <- cumsum((chars == "(") - (chars == ")"))
    substring(text, start, start + nchar(fun) + which(depth == 0L)[1] - 1L)
  }, "")
  calls[which.max(nchar(calls))]
}

# Free bipartitions of the distributed ML tree (not in the constraint tree) and their FBP.
.t53_free_support <- function() {
  skip_if_not_installed("ape")
  ml <- ape::unroot(ape::read.tree(.t53_extdata("phylocactus_ml_tree.tree")))
  cn <- ape::unroot(ape::read.tree(.t53_extdata("phylocactus_constraint_tree.tree")))
  ref <- sort(ml$tip.label)[1]
  key <- function(tr) {
    vapply(ape::prop.part(tr), function(x) {
      s <- tr$tip.label[x]
      if (ref %in% s) s <- setdiff(tr$tip.label, s)
      paste(sort(s), collapse = "|")
    }, "")[-1]
  }
  km <- key(ml)
  sup <- suppressWarnings(as.numeric(ml$node.label))[-1]
  list(all = length(km), constrained = length(key(cn)), free = sup[!km %in% key(cn)])
}

test_that("the README states the bipartitions of the distributed ML and constraint trees (M5)", {
  s <- .t53_free_support()
  r <- .t53_flat(.t53_read("README.md"))
  n <- as.integer(.t53_capture(r, "fixes ([0-9,]+) of the ([0-9,]+) internal bipartitions"))
  expect_identical(n[1], s$constrained)
  expect_identical(n[2], s$all)
  expect_identical(as.integer(.t53_capture(r, "the remaining ([0-9,]+) bipartitions are estimated")),
                   length(s$free))
})

test_that("the help page of integrate_publication_tree() states the free nodes under FBP 70 of the distributed tree (M5)", {
  s <- .t53_free_support()
  h <- .t53_flat(.t53_read("man", "integrate_publication_tree.Rd"))
  n <- as.integer(.t53_capture(h, "collapsing ([0-9,]+) of ([0-9,]+) free nodes"))
  expect_identical(n[1], sum(s$free < 70))
  expect_identical(n[2], length(s$free))
})

test_that("Tutorial 2 states the median FBP of the parents of single-locus terminals (M5)", {
  skip_if_not_installed("ape")
  ml <- ape::read.tree(.t53_extdata("phylocactus_ml_tree.tree"))
  sp <- utils::read.csv(.t53_extdata("phylocactus_table_species.csv"), stringsAsFactors = FALSE)
  skip_if(!"retained_markers" %in% names(sp), "the distributed species table has no retained_markers column")
  loci <- sp$retained_markers[match(ml$tip.label, sp$species)]
  parent <- ml$edge[match(seq_len(ape::Ntip(ml)), ml$edge[, 2]), 1]
  fbp <- suppressWarnings(as.numeric(ml$node.label[parent - ape::Ntip(ml)]))  # percent, as the tree labels it
  v <- .t53_flat(.t53_read("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd"))
  x <- .t53_capture(v, "single locus, whose parent nodes have a median FBP of ([0-9.]+), against ([0-9.]+) for terminals with two loci or more")
  digits <- nchar(sub("^[0-9]*\\.", "", x))
  expect_equal(as.numeric(x[1]), round(stats::median(fbp[loci == 1L], na.rm = TRUE), digits[1]))
  expect_equal(as.numeric(x[2]), round(stats::median(fbp[loci >= 2L], na.rm = TRUE), digits[2]))
  n1 <- as.integer(.t53_capture(v, "([0-9,]+) of the ([0-9,]+) terminals are represented by a single locus"))
  expect_identical(n1, c(sum(loci == 1L, na.rm = TRUE), ape::Ntip(ml)))
})

test_that("Tutorial 1 and its script state the same numbers and one reason for excluding phyC (M5)", {
  v <- .t53_flat(.t53_read("vignettes", "tutorial-1-cactus-phylogeny-prep.Rmd"))
  s <- .t53_script("tutorial-1-cactus-phylogeny-prep.R")
  same <- function(pattern) {
    a <- .t53_capture(v, pattern, info = paste("vignette:", pattern))
    b <- .t53_capture(s, pattern, info = paste("script:", pattern))
    expect_identical(b, a, info = pattern)
  }
  same("and so are ([0-9]+) ingroup sequences")
  same("([0-9]+) of ([0-9]+) in `?rpL16`?")
  same("([0-9]+) of ([0-9]+) in `?psbA-trnH`?")
  same("held ([0-9]+) ingroup sequences and was filtered out")
  same("\\(([0-9]+) sequences, of ([0-9]+) \\*?Portulaca\\*? species")
  same("([0-9.]+) of the outgroup 20-mers occur in the ingroup sequences, all of them on the opposite strand")
  same("holds ([0-9]+) outgroup terminals \\(([0-9]+) Anacampserotaceae, ([0-9]+) Portulacaceae")
  same("only ([0-9.]+) of the outgroup 20-mers occur in the ingroup sequences")
  expect_false(grepl("carries no Portulacaceae terminal", s, fixed = TRUE))
  conc <- .t53_flat(.t53_read("R", "concatenation.R"))
  expect_false(grepl("carries no Portulacaceae", conc, fixed = TRUE))
  expect_false(grepl("nuclear locus that carries no outgroup coverage", conc, fixed = TRUE))
  # The help page of run_marker_screening() gives the ingroup count of the vignette.
  h <- .t53_flat(.t53_read("man", "run_marker_screening.Rd"))
  expect_identical(.t53_capture(h, "trnT-psbD\\}? is the case that forced the distinction: ([0-9]+) ingroup sequences"),
                   .t53_capture(v, "held ([0-9]+) ingroup sequences and was filtered out"))
})

test_that("no help page carries the numbers of an earlier run (M5)", {
  d <- .t53_root("man")
  skip_if(!dir.exists(d), "man/ is not in this tree")
  rd <- vapply(list.files(d, "\\.Rd$", full.names = TRUE),
               function(f) .t53_flat(paste(readLines(f, warn = FALSE), collapse = " ")), "")
  for (p in c("477 of 2044", "703 of 987", "50 ingroup sequences", "the honest figure, which does not exist yet")) {
    expect_false(any(grepl(p, rd, fixed = TRUE)), info = p)
  }
  expect_false(any(grepl("\\b(1021|987)\\b", rd)), info = "1021 or 987 bipartitions of the 0.4.5 tree")
})

test_that("Tutorial 3 and the manifest call the intervals HPD only in the name of the TreeAnnotator column (M7)", {
  for (x in list(.t53_read("vignettes", "tutorial-3-cactus-phylogeny-visualization.Rmd"),
                 .t53_read("inst", "scripts", "tutorial-3-cactus-phylogeny-visualization.R"))) {
    expect_false(grepl("HPD", gsub("height_0.95_HPD", "", x, fixed = TRUE), fixed = TRUE))
  }
  m <- utils::read.csv(.t53_extdata("MANIFEST.csv"), stringsAsFactors = FALSE)
  expect_false(any(grepl("HPD", m$description, fixed = TRUE)))
})

test_that("Tutorial 2 states the replicates dated, and the branches treePL lengthened (M7)", {
  v <- .t53_flat(.t53_read("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd"))
  n <- as.integer(.t53_capture(v, "([0-9]+) of the ([0-9]+) temporal replicates, drawn at random"))
  expect_match(.t53_call(v, "automate_treePL"), paste0("num_bs = ", n[1]), fixed = TRUE)
  expect_match(.t53_call(v, "generate_temporal_bootstrap_script"), paste0("bs_trees = ", n[2]), fixed = TRUE)
  expect_match(v, "563 branches \\(300 of them terminal\\)")
  expect_false(grepl("are floors", v, fixed = TRUE))
  # "On one random partition of the data" described the curve of the 15-10 run; under scheme 1 the
  # selection rests on six of the ten random subsamples (phase 4e, test-redating-0-5-3.R).
})

test_that("Tutorial 3 states the threatened species without sequences and their genera from the registry (M5)", {
  sp <- utils::read.csv(.t53_extdata("phylocactus_table_species.csv"), stringsAsFactors = FALSE)
  skip_if(!all(c("status_pipeline", "iucn_category") %in% names(sp)), "registry columns absent")
  gap <- sp[sp$status_pipeline == "unsampled_checklist_gap", ]
  thr <- gap[gap$iucn_category %in% c("CR", "EN", "VU"), ]
  v <- .t53_flat(.t53_read("vignettes", "tutorial-3-cactus-phylogeny-visualization.Rmd"))
  n <- as.integer(.t53_capture(v, "The other ([0-9,]+) accepted species have no sequence data, and ([0-9,]+) of them are listed as threatened"))
  expect_identical(n, c(nrow(gap), nrow(thr)))
  s <- .t53_capture(v, "among the [0-9,]+ accepted species without sequence data, by genus\\. (.*?) have the most\\.")[1]
  named <- regmatches(s, gregexpr("\\*[A-Z][a-z]+\\* \\([0-9]+", s))[[1]]
  got <- stats::setNames(as.integer(sub(".*\\(", "", named)), gsub("\\*| \\(.*", "", named))
  tab <- sort(table(sub("[_ ].*", "", thr$species)), decreasing = TRUE)
  # Every genus at or above the smallest count named is named, ties included.
  want <- tab[tab >= min(got)]
  expect_identical(sort(names(got)), sort(names(want)))
  expect_identical(unname(got[names(want)]), as.integer(want))
})

test_that("Figure 2 of Tutorial 3 has an axis break for every partition of the supermatrix (m9)", {
  m <- utils::read.csv(.t53_extdata("phylocactus_table_markers.csv"), stringsAsFactors = FALSE)
  for (x in list(.t53_read("vignettes", "tutorial-3-cactus-phylogeny-visualization.Rmd"),
                 .t53_read("inst", "scripts", "tutorial-3-cactus-phylogeny-visualization.R"))) {
    b <- as.integer(.t53_capture(x, "breaks = seq\\(1, ([0-9]+), 1\\)"))
    expect_identical(b, nrow(m))
  }
})

test_that("the references are the ones the publishers record (M8, m1 to m6)", {
  docs <- c(README = .t53_read("README.md"),
            vapply(list.files(.t53_root("vignettes"), "\\.Rmd$", full.names = TRUE),
                   function(f) paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n"), ""))
  expect_false(any(grepl("Helfrich", docs, fixed = TRUE)))
  expect_false(any(grepl("7114", docs, fixed = TRUE)))
  t2 <- .t53_read("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd")
  expect_match(t2, "Bouckaert *et al*. 2019. BEAST 2.5", fixed = TRUE)
  expect_match(t2, "https://doi.org/10.1371/journal.pcbi.1006650", fixed = TRUE)
  expect_match(t2, "Hurvich & Tsai 1989.", fixed = TRUE)
  expect_match(t2, "Flouri *et al*. 2015. The Phylogenetic Likelihood Library", fixed = TRUE)
  expect_match(t2, "Flouri *et al*., 2015)", fixed = TRUE)
  t1 <- .t53_read("vignettes", "tutorial-1-cactus-phylogeny-prep.Rmd")
  expect_false(grepl("Cactineae", t1, fixed = TRUE))
  for (ref in c("Nyffeler & Eggli 2010.", "Ramírez-Barahona *et al*. 2020.", "Zuntini *et al*. 2024.",
                "https://doi.org/10.32614/RJ-2024-022")) {
    expect_match(t1, ref, fixed = TRUE, info = ref)
  }
  expect_false(grepl("B9.bioc.DECIPHER", t1, fixed = TRUE))
  expect_false(grepl("B9.bioc.DECIPHER", .t53_read("R", "alignment.R"), fixed = TRUE))
  for (f in c("README.md", file.path("vignettes", "tutorial-4-cactus-phylogeny-validation.Rmd"))) {
    expect_match(.t53_read(f), "*Nature Communications*, 15(1), 7282.", fixed = TRUE, info = f)
  }
  expect_match(.t53_read("R", "barcoding_report.R"), "10.1111/nph.71494", fixed = TRUE)
})

test_that("the manifest names no path outside the two repositories and no forbidden term (M10, n1)", {
  m <- utils::read.csv(.t53_extdata("MANIFEST.csv"), stringsAsFactors = FALSE)
  txt <- unlist(m, use.names = FALSE)
  expect_false(any(grepl("plantlab-core", txt, fixed = TRUE)))
  expect_false(any(grepl("phylogenom", txt, ignore.case = TRUE)))
  expect_false(any(grepl("^(/|~)", txt)))
})

test_that("Tutorial 2 and its script pass the base seeds of the reference run (M10)", {
  for (x in list(.t53_read("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd"),
                 .t53_read("inst", "scripts", "tutorial-2-cactus-phylogeny-inference.R"))) {
    expect_match(.t53_call(x, "generate_bootstrap_script"), "base_seed = 2150", fixed = TRUE)
    expect_match(.t53_call(x, "generate_temporal_bootstrap_script"), "base_seed = 87786", fixed = TRUE)
  }
})

test_that("the generality claim of Tutorial 1 is restricted to the steps the code allows (M13)", {
  t1 <- .t53_flat(.t53_read("vignettes", "tutorial-1-cactus-phylogeny-prep.Rmd"))
  expect_false(grepl("without modifying the package", t1, fixed = TRUE))
})

test_that("the README lists the external programs with the versions of the reference run (M14)", {
  lines <- strsplit(.t53_read("README.md"), "\n")[[1]]
  versions <- c(MAFFT = "7.526", `ModelTest-NG` = "0.1.7", `RAxML-NG` = "1.2.2", treePL = "2022.04.06",
                `BLAST+` = "2.17.0", TreeAnnotator = "2.7.7", minimap2 = "2.31",
                samtools = "1.24", fastp = "1.3.7", GetOrganelle = "1.7.7.1")
  for (p in names(versions)) {
    has <- grepl(p, lines, fixed = TRUE) & grepl(versions[[p]], lines, fixed = TRUE)
    expect_true(any(has), info = paste(p, versions[[p]]))
  }
  d <- read.dcf(.t53_root("DESCRIPTION"))
  expect_match(.t53_flat(d[1, "SystemRequirements"]), "TreeAnnotator", fixed = TRUE)
})

test_that("the outgroup records without a FASTA file are reported beside those cut (m7)", {
  t1 <- .t53_flat(.t53_read("vignettes", "tutorial-1-cactus-phylogeny-prep.Rmd"))
  expect_match(t1, "51 outgroup records were cut", fixed = TRUE)
  expect_match(t1, "15 outgroup records", fixed = TRUE)
})

test_that("the polytomies are said to belong to integrate_publication_tree(), which Tutorial 3 does not call (m8)", {
  for (x in list(.t53_read("README.md"), .t53_read("vignettes", "PhyloCactus.Rmd"))) {
    expect_match(.t53_flat(x), "integrate_publication_tree()", fixed = TRUE)
  }
  h <- .t53_flat(.t53_read("man", "integrate_publication_tree.Rd"))
  expect_match(h, "constraints_path}{Not used", fixed = TRUE)
})

test_that("the exclusions without a stated ground are described as such (M11)", {
  t1 <- .t53_flat(.t53_read("vignettes", "tutorial-1-cactus-phylogeny-prep.Rmd"))
  expect_false(grepl("on taxonomic or sequence-quality grounds that cannot be recovered", t1, fixed = TRUE))
  ex <- rbind(utils::read.csv(.t53_extdata("manual_exclusions_ingroup.csv"), stringsAsFactors = FALSE)[, c("locus", "reason")],
              utils::read.csv(.t53_extdata("manual_exclusions_outgroup.csv"), stringsAsFactors = FALSE)[, c("locus", "reason")])
  n <- as.integer(.t53_capture(t1, "([0-9]+) were removed by manual conflict resolution"))
  expect_identical(n, sum(grepl("manual conflict resolution", ex$reason, fixed = TRUE)))
})

test_that("the script of the calibration schemes and Tutorial 2 describe the fixed age of scheme 1 (D14, 4e)", {
  # The sweep of the five anchors of D14 became analysis-calibration-schemes.R in phase 4e; the
  # September reason for a fixed anchor (0.4.x tree) is no longer given.
  s <- .t53_script("analysis-calibration-schemes.R")
  expect_false(grepl("A fixed anchor removes that artifact", s, fixed = TRUE))
  v <- .t53_flat(.t53_read("vignettes", "tutorial-2-cactus-phylogeny-inference.Rmd"))
  expect_match(v, "median height", fixed = TRUE)
  expect_match(v, "maximum clade credibility tree", fixed = TRUE)
})
