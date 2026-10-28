# Version 0.5.3 closes the audit of 26-10-2026 (closure plan, phase 5a). The package carries the
# version and NEWS opens with it and says what left inst/extdata and how a user obtains each item;
# PROVENANCE.csv records the seeds of the reference run read from its logs, the versions of the
# programs, the version of the package that published the run and the checklist it read; the counts
# of the exclusion lists in the help pages of assemble_*_phylotar() and the table of Tutorial 4 with
# and without the constrained bipartitions are recomputed from inst/extdata. Written before the change.

.r53_root <- function(...) testthat::test_path("..", "..", ...)

.r53_read <- function(...) {
  f <- .r53_root(...)
  skip_if(!file.exists(f), paste(file.path(...), "is not in this tree"))
  readLines(f, warn = FALSE, encoding = "UTF-8")
}

.r53_flat <- function(x) gsub("\\s+", " ", paste(x, collapse = "\n"))

.r53_extdata <- function(name) {
  f <- system.file("extdata", name, package = "PhyloCactus")
  skip_if(!nzchar(f), paste(name, "not installed"))
  f
}

.r53_capture <- function(text, pattern) {
  m <- regmatches(text, regexec(pattern, text))[[1]]
  expect_true(length(m) > 1L, info = pattern)
  if (length(m) < 2L) return(rep(NA_character_, 8L))
  m[-1]
}

# The entry of 0.5.3 in NEWS, from its heading to the heading of 0.5.2.
.r53_news_entry <- function() {
  news <- .r53_read("NEWS.md")
  end <- grep("^# PhyloCactus 0\\.5\\.2$", news)
  expect_length(end, 1L)
  if (length(end) != 1L) return("")
  .r53_flat(news[seq_len(end - 1L)])
}

test_that("the package is version 0.5.3 and NEWS opens with it", {
  expect_identical(as.character(utils::packageVersion("PhyloCactus")), "0.5.3")
  news <- .r53_read("NEWS.md")
  expect_identical(news[1], "# PhyloCactus 0.5.3")
  expect_length(grep("^# PhyloCactus 0\\.5\\.2$", news), 1L)
})

test_that("NEWS of 0.5.3 says what left inst/extdata and how a user obtains each item", {
  e <- .r53_news_entry()
  for (p in c("phylocactus_checklist()", "PHYLOCACTUS_CHECKLIST", "Korotkova et al. (2021)",
              "Caryophyllales.org", "0ab7855cef6fe7db0e0105646d5fee08",
              "PHYLOCACTUS_REFERENCE_TREES", "get_iucn_data()", "summarise_iucn_registry()",
              "iucn_found", "iucn_category", "iucn_category_name",
              "phylocactus_table_iucn_endemic_by_country.csv", "phylocactus_table_iucn_classes.csv",
              "phylocactus_table_tree_space.csv", "plot_tree_space()", "run_treePL_cv()",
              "PROVENANCE.csv")) {
    expect_true(grepl(p, e, fixed = TRUE), info = p)
  }
  for (f in c("reference_amaral_2022.tree", "reference_thompson_2024.tree", "reference_zuntini_2024.tree",
              "reference_devos_2025_astral.tree", "CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx")) {
    expect_true(grepl(f, e, fixed = TRUE), info = f)
  }
})

test_that("no text plans work for a version 0.6", {
  for (f in c(list(c("NEWS.md"), c("README.md")),
              lapply(list.files(.r53_root("vignettes"), pattern = "\\.Rmd$"), function(x) c("vignettes", x)))) {
    txt <- .r53_flat(do.call(.r53_read, as.list(f)))
    expect_false(grepl("planned for (version )?0\\.6", txt), info = paste(f, collapse = "/"))
  }
})

test_that("PROVENANCE.csv records the seeds of the reference run read from its logs (M10)", {
  p <- utils::read.csv(.r53_extdata("PROVENANCE.csv"), stringsAsFactors = FALSE, colClasses = "character")
  seeds <- c(seed_ml_search = "1111", seed_bootstrap = "2150; 4300; 6450; 8600",
             seed_temporal_bootstrap = "87786", seed_treepl = "79992967")
  for (s in names(seeds)) {
    expect_true(s %in% names(p), info = s)
    if (s %in% names(p)) expect_identical(unique(p[[s]]), unname(seeds[s]), info = s)
  }
})

test_that("PROVENANCE.csv records the versions of the programs and of the package that published the run (M10)", {
  p <- utils::read.csv(.r53_extdata("PROVENANCE.csv"), stringsAsFactors = FALSE, colClasses = "character")
  versions <- c(mafft = "^7\\.526( \\(|$)", modeltest_ng = "^0\\.1\\.7( \\(|$)", phylotar = "^1\\.3\\.0( \\(|$)",
                blast = "^2\\.17\\.0( \\(|$)", decipher = "3\\.6\\.0.*3\\.8\\.1|3\\.8\\.1.*3\\.6\\.0",
                treeannotator = "^2\\.7\\.7$", run_package_version = "^0\\.5\\.0$")
  for (v in names(versions)) {
    expect_true(v %in% names(p), info = v)
    if (v %in% names(p)) {
      expect_length(unique(p[[v]]), 1L)
      expect_match(unique(p[[v]])[1], versions[[v]], info = v)
    }
  }
  # The versions already recorded before 0.5.3 stay.
  expect_true(all(c("raxml_ng", "treepl", "R_version", "package_version") %in% names(p)))
})

test_that("PROVENANCE.csv names the checklist the run read, with its md5 (Korotkova et al. 2021)", {
  p <- utils::read.csv(.r53_extdata("PROVENANCE.csv"), stringsAsFactors = FALSE, colClasses = "character")
  expect_true(all(c("checklist_file", "checklist_md5") %in% names(p)))
  skip_if(!all(c("checklist_file", "checklist_md5") %in% names(p)))
  expect_identical(unique(p$checklist_file), "CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx")
  expect_identical(unique(p$checklist_md5), "0ab7855cef6fe7db0e0105646d5fee08")
  expect_identical(unique(p$checklist_md5), PhyloCactus:::.phylocactus_checklist_md5)
})

test_that("the help pages of assemble_*_phylotar() state the counts of the distributed exclusion lists", {
  ing <- utils::read.csv(.r53_extdata("manual_exclusions_ingroup.csv"), stringsAsFactors = FALSE)
  out <- utils::read.csv(.r53_extdata("manual_exclusions_outgroup.csv"), stringsAsFactors = FALSE)
  all <- rbind(ing[, c("locus", "sid", "reason")], out[, c("locus", "sid", "reason")])
  conflict <- all$reason == "manual conflict resolution"
  expected <- c(nrow(all), sum(!conflict), sum(conflict),
                nrow(all), nrow(ing), length(unique(ing$locus)), nrow(out), length(unique(out$locus)))
  expect_identical(expected[1:3], c(102L, 47L, 55L))
  for (src in list(c("man", "assemble_ingroup_phylotar.Rd"), c("R", "assembly.R"))) {
    txt <- gsub("\\\\(code|strong)\\{([^}]*)\\}", "\\2", .r53_flat(do.call(.r53_read, as.list(src))))
    a <- .r53_capture(txt, "Of their ([0-9]+) accessions, ([0-9]+) state their ground")
    b <- .r53_capture(txt, "and ([0-9]+) carry only \"manual conflict resolution\"")
    c <- .r53_capture(txt, paste0("they hold ([0-9]+) accessions, ([0-9]+) in ([0-9]+) ingroup loci ",
                                  "and ([0-9]+) in ([0-9]+) outgroup loci"))
    expect_identical(as.integer(c(a[1:2], b[1], c[1:5])), expected, info = paste(src, collapse = "/"))
  }
  # The help page of the outgroup sends the reader to that of the ingroup, where the counts are.
  o <- .r53_flat(.r53_read("man", "assemble_outgroup_phylotar.Rd"))
  expect_true(grepl("assemble_ingroup_phylotar", o, fixed = TRUE))
})

# The table of Tutorial 4, "RF with and without the constrained bipartitions" (phase 4b, M3). The RF
# and the shared terminals are those of the distributed distances. The imposed bipartitions and the
# RF without them need the reference trees, which are not distributed; on every pair of the reference
# run both pruned trees are fully resolved (n - 3 bipartitions each; recount m3_tutorial4_distances.R
# of the audit, columns F and R), so with n shared terminals, RF the published distance, k the imposed
# bipartitions also in the reference and C all the imposed ones on the set (k <= C <= the 34 of the
# constraint tree), RF without them is (2(n - 3) RF - (C - k)) / (2(n - 3) - C - k). Each row must
# be that value, at three decimals, for some such C; and the shared bipartitions are (n - 3)(1 - RF).
.r53_rf_table <- function() {
  v <- .r53_read("vignettes", "tutorial-4-cactus-phylogeny-validation.Rmd")
  rows <- grep("^\\| A[1-4] \\|", v, value = TRUE)
  cells <- lapply(strsplit(rows, "\\|"), function(x) trimws(x[-1]))
  data.frame(set = vapply(cells, `[`, "", 1), reference = gsub(" ", "", vapply(cells, `[`, "", 2)),
             n = as.integer(vapply(cells, `[`, "", 3)), imposed = as.integer(vapply(cells, `[`, "", 4)),
             rf = as.numeric(vapply(cells, `[`, "", 5)), rf_free = as.numeric(vapply(cells, `[`, "", 6)),
             stringsAsFactors = FALSE)
}

test_that("the RF table of Tutorial 4 gives the distributed distances for its 14 pairs (M3)", {
  t4 <- .r53_rf_table()
  expect_identical(nrow(t4), 14L)
  d <- utils::read.csv(.r53_extdata("phylocactus_table_tree_distances.csv"), stringsAsFactors = FALSE)
  d <- d[d$metric == "RF" & (d$tree_1 == "FocalTree" | d$tree_2 == "FocalTree"), ]
  d$reference <- ifelse(d$tree_1 == "FocalTree", d$tree_2, d$tree_1)
  d$set <- sub("_.*$", "", d$analysis_id)
  expect_setequal(paste(t4$set, t4$reference), paste(d$set, d$reference))
  i <- match(paste(t4$set, t4$reference), paste(d$set, d$reference))
  expect_false(anyNA(i))
  expect_identical(t4$n, as.integer(d$n_common_tips[i]))
  expect_equal(t4$rf, round(d$distance[i], 3), tolerance = 1e-9)
})

test_that("the imposed bipartitions and the RF without them in Tutorial 4 agree with the distributed distances (M3)", {
  skip_if_not_installed("ape")
  t4 <- .r53_rf_table()
  d <- utils::read.csv(.r53_extdata("phylocactus_table_tree_distances.csv"), stringsAsFactors = FALSE)
  d <- d[d$metric == "RF" & (d$tree_1 == "FocalTree" | d$tree_2 == "FocalTree"), ]
  rf <- d$distance[match(paste(t4$set, t4$reference),
                         paste(sub("_.*$", "", d$analysis_id), ifelse(d$tree_1 == "FocalTree", d$tree_2, d$tree_1)))]
  cn <- ape::unroot(ape::read.tree(.r53_extdata("phylocactus_constraint_tree.tree")))
  n_constraint <- length(ape::prop.part(cn)) - 1L
  expect_identical(n_constraint, 34L)
  m <- 2 * (t4$n - 3)
  for (j in seq_len(nrow(t4))) {
    C <- t4$imposed[j]:n_constraint
    pred <- round((m[j] * rf[j] - (C - t4$imposed[j])) / (m[j] - C - t4$imposed[j]), 3)
    expect_true(any(abs(pred - t4$rf_free[j]) < 1e-9), info = paste(t4$set[j], t4$reference[j]))
  }
  shared <- round((t4$n - 3) * (1 - rf))
  share <- round(100 * t4$imposed / shared)
  v <- .r53_flat(.r53_read("vignettes", "tutorial-4-cactus-phylogeny-validation.Rmd"))
  r <- as.integer(.r53_capture(v, "([0-9]+) to ([0-9]+) of the bipartitions the focal tree shares"))
  expect_identical(r[1:2], range(t4$imposed))
  small <- t4$set %in% c("A1", "A4")
  large <- t4$set %in% c("A2", "A3") & t4$reference %in% c("Thompson", "Amaral")
  r <- as.integer(.r53_capture(v, "([0-9]+) to ([0-9]+) % of the shared bipartitions in the small sets"))
  expect_identical(r[1:2], as.integer(range(share[small])))
  r <- as.integer(.r53_capture(v, "and ([0-9]+) to ([0-9]+) % in the large sets"))
  expect_identical(r[1:2], as.integer(range(share[large])))
  r <- as.numeric(.r53_capture(v, "Every distance rises, by ([0-9.]+) to ([0-9.]+);"))
  expect_equal(r[1:2], range(round(t4$rf_free - t4$rf, 3)), tolerance = 1e-9)
  # The order of the references by distance holds in A1 to A3; in A4 Zuntini and Thompson exchange.
  for (s in c("A1", "A2", "A3")) {
    k <- t4$set == s
    expect_identical(order(t4$rf[k]), order(t4$rf_free[k]), info = s)
  }
  a4 <- t4[t4$set == "A4", ]
  expect_true(a4$rf[a4$reference == "Zuntini"] < a4$rf[a4$reference == "Thompson"])
  expect_true(a4$rf_free[a4$reference == "Zuntini"] > a4$rf_free[a4$reference == "Thompson"])
})
