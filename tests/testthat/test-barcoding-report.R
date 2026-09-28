# Tests of the identification report of the branch (Phase 9A; decisions IR1 to IR6 of BMM, 28-09).
# Written before the code. The report reads the identification table and the run record; the tests
# of the report write both by hand, so they run no classifier. The tests that go through the three
# identification functions run IdTaxa and skip where it does not run (as in test-barcoding-identify.R).

.rep_write <- function(d, path) utils::write.csv(d, path, row.names = FALSE)

# A library of three loci and three genera; Cereus has no trnL-trnF of C. beta, so a cell is empty
.rep_library <- function(tmp) {
  lib <- file.path(tmp, "4_library")
  dir.create(lib, recursive = TRUE)
  put <- function(locus, heads) {
    s <- stats::setNames(rep("ACGTACGTACGTACGTACGT", length(heads)), heads)
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(s), file.path(lib, paste0("LIB_", locus, ".fasta")))
  }
  put("matK", c("Opuntia_alpha|m1.1", "Opuntia_alpha|m2.1", "Opuntia_beta|m3.1", "Cereus_alpha|m4.1",
                "Cereus_beta|m5.1", "Mammillaria_alpha|m6.1"))
  put("rbcL", c("Opuntia_alpha|r1.1", "Opuntia_beta|r2.1", "Cereus_alpha|r3.1", "Mammillaria_alpha|r4.1"))
  put("trnL-trnF", c("Opuntia_alpha|t1.1", "Cereus_alpha|t2.1", "Cereus_alpha|t3.1"))
  put("ITS", c("Opuntia_alpha|i1.1", "Cereus_alpha|i2.1"))
  lib
}

.rep_row <- function(query, locus, state, species, genus, candidates, genus_idtaxa, species_idtaxa, gc, sc, reason,
                     region_length = 300L) {
  data.frame(query = query, locus = locus, query_length = 300L, path = if (identical(reason, "no_overlap")) "none" else "whole",
             region_start = 1L, region_end = region_length, region_length = if (identical(reason, "no_overlap")) NA_integer_ else region_length,
             other_windows = 0L, orientation = "forward", state = state, predicted_species = species,
             predicted_genus = genus, candidates = candidates, genus_idtaxa = genus_idtaxa, species_idtaxa = species_idtaxa,
             genus_confidence = gc, species_confidence = sc, reason = reason, core_trimmed = 0L, genus_core_trimmed = 0L,
             threshold = 60, validation_ws_rate_species_present = 0.02, validation_ws_rate_species_absent = 0.05,
             stringsAsFactors = FALSE)
}

# One query answered in four ways: a species (matK), a genus (rbcL), nothing over the threshold
# (trnL-trnF) and no overlap (ITS)
.rep_table <- function(query = "q1") {
  rbind(.rep_row(query, "ITS", 3L, NA, NA, NA, NA, NA, NA, NA, "no_overlap"),
        .rep_row(query, "matK", 1L, "Opuntia_alpha", "Opuntia", "Opuntia_alpha", "Opuntia", "Opuntia_alpha", 95, 80, NA),
        .rep_row(query, "rbcL", 2L, NA, "Opuntia", "Opuntia_alpha|Opuntia_beta", "Opuntia", "Opuntia_beta", 85, 40, NA),
        .rep_row(query, "trnL-trnF", 3L, NA, NA, NA, "Cereus", "Cereus_alpha", 45, 30, "low_confidence"))
}

.rep_record <- function(input_type = "sequences") {
  data.frame(kind = c("input", "input", "input", "setting", "setting", "software", "software"),
             name = c("path", "md5", "input_type", "threshold", "seed", "R", "PhyloCactus"),
             value = c("query.fasta", "0123456789abcdef0123456789abcdef", input_type, "60", "1",
                       "R version 4.5.3", "0.4.5.9000"), stringsAsFactors = FALSE)
}

.rep_fixture <- function(tmp, table = .rep_table(), record = .rep_record(), run = "run1") {
  res <- file.path(tmp, "10_identify")
  dir.create(res, recursive = TRUE)
  .rep_write(table, file.path(res, paste0("TABLE_barcoding_identify_", run, ".csv")))
  .rep_write(record, file.path(res, paste0("RUN_", run, ".csv")))
  list(results_dir = res, library_dir = .rep_library(tmp), metrics_dir = file.path(tmp, "11_metrics"), run = run)
}

.rep_run <- function(f) {
  suppressMessages(invisible(utils::capture.output(
    r <- report_barcoding_identification(f$run, results_dir = f$results_dir, library_dir = f$library_dir,
                                         metrics_dir = f$metrics_dir))))
  r
}

.rep_html <- function(f, query = "q1") {
  paste(readLines(file.path(f$results_dir, paste0("REPORT_", f$run, "_", query, ".html")), warn = FALSE), collapse = "\n")
}

.rep_csv <- function(f, part) {
  utils::read.csv(file.path(f$results_dir, paste0("REPORT_", f$run, "_", part, ".csv")), stringsAsFactors = FALSE)
}

test_that("the internal base64 encoder gives the values of RFC 4648", {
  enc <- PhyloCactus:::.bc_base64
  expect_identical(enc(charToRaw("Man")), "TWFu")
  expect_identical(enc(charToRaw("Ma")), "TWE=")
  expect_identical(enc(charToRaw("M")), "TQ==")
  expect_identical(enc(charToRaw("foobar")), "Zm9vYmFy")
  expect_identical(enc(raw(0)), "")
})

test_that("the report is one self-contained HTML per query with the nine sections", {
  skip_if_not_installed("Biostrings")
  f <- .rep_fixture(withr::local_tempdir())
  .rep_run(f)
  h <- .rep_html(f)
  for (i in 1:9) expect_true(grepl(sprintf("id=\"section-%d\"", i), h, fixed = TRUE), info = i)
  # Self-contained: every image is embedded, nothing is loaded from a file or the network
  src <- regmatches(h, gregexpr("(src|href)=\"[^\"]*\"", h))[[1]]
  expect_true(length(src) > 0L)
  expect_true(all(grepl("^(src|href)=\"(data:|#)", src)))
  expect_true(grepl("data:image/png;base64,", h, fixed = TRUE))
})

test_that("the report never states a combined call across loci (E8)", {
  skip_if_not_installed("Biostrings")
  f <- .rep_fixture(withr::local_tempdir())
  .rep_run(f)
  h <- tolower(.rep_html(f))
  for (w in c("consensus", "combined identification", "overall identification", "final identification")) {
    expect_false(grepl(w, h, fixed = TRUE), info = w)
  }
  ag <- .rep_csv(f, "agreement")
  expect_false(any(grepl("consensus|combined", names(ag), ignore.case = TRUE)))
})

test_that("the answer table of the report repeats the identification table, locus by locus", {
  skip_if_not_installed("Biostrings")
  f <- .rep_fixture(withr::local_tempdir())
  .rep_run(f)
  a <- .rep_csv(f, "answer")
  t <- .rep_table()
  a <- a[match(t$locus, a$locus), ]
  expect_identical(a$state, t$state)
  expect_identical(a$predicted_species, t$predicted_species)
  expect_identical(a$predicted_genus, t$predicted_genus)
  expect_equal(a$genus_confidence, t$genus_confidence)
  expect_equal(a$species_confidence, t$species_confidence)
  expect_identical(a$compartment[a$locus == "ITS"], "nuclear")
  expect_identical(a$compartment[a$locus == "matK"], "plastid")
})

test_that("the alternatives are the best path under the threshold and the candidates of state 2, never assigned", {
  skip_if_not_installed("Biostrings")
  f <- .rep_fixture(withr::local_tempdir())
  .rep_run(f)
  al <- .rep_csv(f, "alternatives")
  expect_true(all(!al$assigned))
  b <- al[al$kind == "best_path_under_threshold", ]
  expect_setequal(b$locus, c("rbcL", "trnL-trnF"))
  expect_identical(b$species[b$locus == "trnL-trnF"], "Cereus_alpha")
  expect_equal(b$genus_confidence[b$locus == "trnL-trnF"], 45)
  expect_equal(b$species_confidence[b$locus == "rbcL"], 40)
  c2 <- al[al$kind == "candidate_of_named_genus", ]
  expect_setequal(c2$species, c("Opuntia_alpha", "Opuntia_beta"))
  expect_true(all(c2$locus == "rbcL"))
  # The section is headed with the threshold
  h <- .rep_html(f)
  s4 <- sub(".*id=\"section-4\"(.*?)id=\"section-5\".*", "\\1", h)
  expect_true(grepl("60", s4, fixed = TRUE))
  expect_true(grepl("not assigned", s4, fixed = TRUE))
})

test_that("the reading across loci counts the loci that name each taxon, assigned or not", {
  skip_if_not_installed("Biostrings")
  f <- .rep_fixture(withr::local_tempdir())
  .rep_run(f)
  ag <- .rep_csv(f, "agreement")
  op <- ag[ag$rank == "genus" & ag$taxon == "Opuntia", ]
  expect_equal(op$n_loci, 2L)
  expect_equal(op$n_assigned, 2L)
  ce <- ag[ag$rank == "genus" & ag$taxon == "Cereus", ]
  expect_equal(ce$n_loci, 1L)
  expect_equal(ce$n_assigned, 0L)
  sp <- ag[ag$rank == "species" & ag$taxon == "Opuntia_alpha", ]
  expect_equal(sp$n_assigned, 1L)
  none <- ag[ag$rank == "none", ]
  expect_identical(none$loci, "ITS")
})

test_that("the sampling of the named genera in the library shows every species and the empty cells", {
  skip_if_not_installed("Biostrings")
  f <- .rep_fixture(withr::local_tempdir())
  .rep_run(f)
  s <- .rep_csv(f, "sampling")
  expect_setequal(unique(s$genus), c("Opuntia", "Cereus"))
  expect_false(any(grepl("Mammillaria", s$species)))
  expect_equal(s$sequences[s$species == "Opuntia_alpha" & s$locus == "matK"], 2L)
  expect_equal(s$sequences[s$species == "Cereus_alpha" & s$locus == "trnL-trnF"], 2L)
  expect_equal(s$sequences[s$species == "Cereus_beta" & s$locus == "trnL-trnF"], 0L)
})

test_that("the provenance repeats the run record, and the reads section applies only to reads", {
  skip_if_not_installed("Biostrings")
  f <- .rep_fixture(withr::local_tempdir())
  .rep_run(f)
  p <- .rep_csv(f, "provenance")
  expect_identical(p$value[p$name == "md5"], "0123456789abcdef0123456789abcdef")
  h <- .rep_html(f)
  s7 <- sub(".*id=\"section-7\"(.*?)id=\"section-8\".*", "\\1", h)
  expect_true(grepl("Not applicable", s7, fixed = TRUE))
})

test_that("a report rebuilt from the saved tables is identical, byte for byte", {
  skip_if_not_installed("Biostrings")
  f <- .rep_fixture(withr::local_tempdir())
  .rep_run(f)
  files <- list.files(f$results_dir, pattern = "^REPORT_", full.names = TRUE)
  first <- tools::md5sum(files)
  .rep_run(f)
  expect_identical(tools::md5sum(files), first)
})

test_that("a run of two queries gives two HTML files and one set of tables with the query", {
  skip_if_not_installed("Biostrings")
  f <- .rep_fixture(withr::local_tempdir(), table = rbind(.rep_table("q1"), .rep_table("q2")))
  .rep_run(f)
  expect_true(all(file.exists(file.path(f$results_dir, paste0("REPORT_run1_", c("q1", "q2"), ".html")))))
  a <- .rep_csv(f, "answer")
  expect_setequal(unique(a$query), c("q1", "q2"))
  expect_equal(nrow(a), 8L)
})

test_that("the per-genus discrimination of each locus is counted from the 6B predictions at 60", {
  tmp <- withr::local_tempdir()
  cls <- file.path(tmp, "7_classifier"); thr <- file.path(tmp, "9_threshold")
  dir.create(cls); dir.create(thr)
  row <- function(scheme, sid, true_species, gi, si, gc, sc) {
    data.frame(locus = "matK", scheme = scheme, fold = sid, stratum = true_species, sid = paste0("s", sid),
               true_species = true_species, true_genus = sub("_.*", "", true_species), method = "idtaxa",
               state = 3L, predicted_species = NA_character_, predicted_genus = NA_character_, candidates = NA_character_,
               confidence = NA_real_, reason = "low_confidence", genus_idtaxa = gi, species_idtaxa = si,
               genus_confidence = gc, species_confidence = sc, stringsAsFactors = FALSE)
  }
  .rep_write(rbind(row("species", 1, "A_a", "A", "A_a", 95, 90), row("species", 2, "A_b", "A", "A_a", 99, 70),
                   row("species", 3, "B_c", "B", "B_d", 80, 40), row("species", 4, "B_c", "B", "B_c", 30, 20)),
             file.path(cls, "TABLE_barcoding_predictions_species_idtaxa.csv"))
  .rep_write(rbind(row("genus", 1, "A_a", "A", "A_b", 90, 50), row("genus", 2, "B_c", "A", "A_a", 70, 20)),
             file.path(cls, "TABLE_barcoding_predictions_genus_idtaxa.csv"))
  .rep_write(data.frame(locus = "matK", scheme = rep(c("species", "genus"), each = 2), rule = rep(c("quantile", "default_60"), 2),
                        threshold = rep(c(10, 60), 2)), file.path(thr, "TABLE_barcoding_threshold_operating_idtaxa.csv"))
  out <- file.path(tmp, "11_metrics")
  suppressMessages(invisible(utils::capture.output(
    summarise_barcoding_genus_discrimination(classifier_dir = cls, threshold_dir = thr, output_dir = out))))
  g <- utils::read.csv(file.path(out, "TABLE_barcoding_genus_discrimination_idtaxa.csv"), stringsAsFactors = FALSE)
  a <- g[g$scheme == "species" & g$genus == "A", ]
  expect_equal(c(a$queries, a$correct_species, a$wrong_species, a$unassigned), c(2L, 1L, 1L, 0L))
  expect_equal(a$wrong_species_rate, 0.5)
  b <- g[g$scheme == "species" & g$genus == "B", ]
  expect_equal(c(b$queries, b$correct_genus, b$unassigned), c(2L, 1L, 1L))
  bg <- g[g$scheme == "genus" & g$genus == "B", ]
  expect_equal(c(bg$queries, bg$wrong_genus), c(1L, 1L))
  expect_true(all(g$threshold == 60))
})

test_that("the report shows the per-genus discrimination of the genus of the answer, or says it is absent", {
  skip_if_not_installed("Biostrings")
  f <- .rep_fixture(withr::local_tempdir())
  .rep_run(f)
  s5 <- sub(".*id=\"section-5\"(.*?)id=\"section-6\".*", "\\1", .rep_html(f))
  expect_true(grepl("not available", s5, fixed = TRUE))

  g <- .rep_fixture(withr::local_tempdir())
  dir.create(g$metrics_dir)
  .rep_write(data.frame(locus = c("matK", "rbcL"), scheme = "species", genus = "Opuntia", threshold = 60, queries = c(12L, 9L),
                        correct_species = c(5L, 2L), wrong_species = c(1L, 0L), correct_genus = c(4L, 5L), wrong_genus = 0L,
                        unassigned = c(2L, 2L), correct_species_rate = c(5 / 12, 2 / 9), wrong_species_rate = c(1 / 12, 0),
                        correct_genus_rate = c(4 / 12, 5 / 9), unassigned_rate = c(2 / 12, 2 / 9)),
             file.path(g$metrics_dir, "TABLE_barcoding_genus_discrimination_idtaxa.csv"))
  .rep_run(g)
  s5 <- sub(".*id=\"section-5\"(.*?)id=\"section-6\".*", "\\1", .rep_html(g))
  expect_false(grepl("not available", s5, fixed = TRUE))
  expect_true(grepl("Opuntia", s5, fixed = TRUE))
})

# ---- Through the identification functions (IdTaxa) -----------------------------------------------

.rep_idtaxa_runs <- function() {
  tryCatch({
    tr <- Biostrings::DNAStringSet(c(a = "ACGTACGTAACCGGTTAACC", b = "ACGTACGTACCCGGTTAACC",
                                     c = "TTTTGGGGAATTCCGGAATT", d = "TTTTGGGGAATTCCGGAATG"))
    lt <- DECIPHER::LearnTaxa(tr, c("Root;Opuntia;Opuntia_robusta", "Root;Opuntia;Opuntia_stricta",
                                    "Root;Cereus;Cereus_jamacaru", "Root;Cereus;Cereus_horrida"), verbose = FALSE)
    DECIPHER::IdTaxa(Biostrings::DNAStringSet("ACGTACGTAACCGGTTAACC"), lt, strand = "top", threshold = 60,
                     processors = 1L, verbose = FALSE)
    TRUE
  }, error = function(e) FALSE)
}

.rep_skip_idtaxa <- function() {
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("Biostrings")
  if (!suppressWarnings(suppressMessages(.rep_idtaxa_runs()))) {
    skip("DECIPHER::IdTaxa() does not run in this loading context: match() over XStringSet fails to dispatch inside IdTaxa. The same call works in a plain R session.")
  }
}

.rep_random <- function(n) paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
.rep_vary <- function(s, k) {
  for (p in sample(seq_len(nchar(s)), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1)
  s
}
.rep_idn_library <- function(tmp) {
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir, recursive = TRUE)
  set.seed(41L)
  seqs <- list()
  for (l in c("matK", "rbcL")) {
    root <- .rep_random(300); x <- character(0)
    for (g in c("Opuntia", "Cereus")) {
      gb <- .rep_vary(root, 6)
      for (e in c("alpha", "beta")) {
        sb <- .rep_vary(gb, 3)
        for (r in 1:3) x[paste0(g, "_", e, "|", l, "_", substr(g, 1, 1), e, r, ".1")] <- .rep_vary(sb, 1)
      }
    }
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), file.path(lib_dir, paste0("LIB_", l, ".fasta")))
    seqs[[l]] <- x
  }
  list(library_dir = lib_dir, seqs = seqs, out = file.path(tmp, "10_identify"), tmp = tmp)
}

test_that("identify_barcoding_query() writes the run record and the report, unless report = FALSE", {
  .rep_skip_idtaxa()
  f <- .rep_idn_library(withr::local_tempdir())
  fa <- file.path(f$tmp, "query.fasta")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(q1 = unname(f$seqs$matK[1]))), fa)
  suppressMessages(invisible(utils::capture.output(
    identify_barcoding_query(fa, library_dir = f$library_dir, metrics_dir = NULL, output_dir = f$out, run_name = "r"))))
  rec <- utils::read.csv(file.path(f$out, "RUN_r.csv"), stringsAsFactors = FALSE)
  expect_identical(rec$value[rec$name == "md5"], unname(tools::md5sum(fa)))
  expect_identical(rec$value[rec$name == "input_type"], "sequences")
  for (l in c("matK", "rbcL")) {
    expect_identical(rec$value[rec$kind == "library" & rec$name == paste0("LIB_", l, ".fasta")],
                     unname(tools::md5sum(file.path(f$library_dir, paste0("LIB_", l, ".fasta")))))
  }
  expect_true(all(c("R", "PhyloCactus", "DECIPHER", "Biostrings") %in% rec$name[rec$kind == "software"]))
  expect_identical(rec$value[rec$name == "threshold"], "60")
  expect_true(file.exists(file.path(f$out, "REPORT_r_q1.html")))

  suppressMessages(invisible(utils::capture.output(
    identify_barcoding_query(fa, library_dir = f$library_dir, metrics_dir = NULL, output_dir = f$out, run_name = "s",
                             report = FALSE))))
  expect_true(file.exists(file.path(f$out, "RUN_s.csv")))
  expect_false(any(grepl("^REPORT_s", list.files(f$out))))
})

test_that("a GenBank flat file is recorded as such in the run record", {
  .rep_skip_idtaxa()
  f <- .rep_idn_library(withr::local_tempdir())
  s <- tolower(unname(f$seqs$rbcL[4]))
  gb <- file.path(f$tmp, "q.gb")
  writeLines(c("LOCUS       AB000001 300 bp    DNA     linear   PLN 01-JAN-2026", "VERSION     AB000001.1",
               "FEATURES             Location/Qualifiers", "     source          1..300",
               "                     /organism=\"Opuntia beta\"", "ORIGIN",
               vapply(seq(1, 300, 60), function(i) sprintf("%9d %s", i, substr(s, i, i + 59)), ""), "//"), gb)
  suppressMessages(invisible(utils::capture.output(
    identify_barcoding_query(gb, library_dir = f$library_dir, metrics_dir = NULL, output_dir = f$out, run_name = "g"))))
  rec <- utils::read.csv(file.path(f$out, "RUN_g.csv"), stringsAsFactors = FALSE)
  expect_identical(rec$value[rec$name == "input_type"], "genbank")
  expect_true(file.exists(file.path(f$out, "REPORT_g_AB000001.1.html")))
})

test_that("an assembly or a set of reads gives one report for the sample, reading its regions together", {
  skip_if_not_installed("Biostrings")
  t <- .rep_table("acc")
  t$query <- paste0("acc__", t$locus)
  for (type in c("assembly", "reads")) {
    f <- .rep_fixture(withr::local_tempdir(), table = t, record = .rep_record(type))
    .rep_run(f)
    html <- list.files(f$results_dir, pattern = "\\.html$")
    expect_identical(html, "REPORT_run1.html")
    ag <- .rep_csv(f, "agreement")
    op <- ag[ag$rank == "genus" & ag$taxon == "Opuntia", ]
    expect_equal(op$n_loci, 2L)
    expect_identical(unique(ag$query), "run1")
    a <- .rep_csv(f, "answer")
    expect_setequal(a$query, t$query)
  }
})

test_that("the report names the kind of data at its start and the library it used", {
  skip_if_not_installed("Biostrings")
  kinds <- list(sequences = "Sanger-type sequence", assembly = "Nuclear genome assembly",
                reads = "Genome skimming or whole-genome sequencing reads")
  for (type in names(kinds)) {
    f <- .rep_fixture(withr::local_tempdir(), record = .rep_record(type))
    .rep_run(f)
    h <- paste(readLines(list.files(f$results_dir, pattern = "\\.html$", full.names = TRUE)[1], warn = FALSE), collapse = "\n")
    top <- sub("id=\"section-2\".*", "", h)
    expect_true(grepl(kinds[[type]], top, fixed = TRUE), info = type)
    expect_true(grepl("Sanger records of GenBank (phylotaR)", top, fixed = TRUE), info = type)
  }
  # A query of plastome length is named as such
  t <- .rep_table()
  t$query_length <- 150000L
  f <- .rep_fixture(withr::local_tempdir(), table = t, record = .rep_record("genbank"))
  .rep_run(f)
  top <- sub("id=\"section-2\".*", "", .rep_html(f))
  expect_true(grepl("Complete plastome", top, fixed = TRUE))
})

test_that("the report carries the cactus of PhyloCactus in its title", {
  skip_if_not_installed("Biostrings")
  f <- .rep_fixture(withr::local_tempdir())
  .rep_run(f)
  h <- .rep_html(f)
  expect_true(grepl("\U0001F335", h, fixed = TRUE))
})

test_that("a genus confidence under the first stripe still gives a report", {
  skip_if_not_installed("Biostrings")
  t <- .rep_table()
  t$genus_confidence[t$locus == "trnL-trnF"] <- 1.9
  f <- .rep_fixture(withr::local_tempdir(), table = t)
  expect_error(.rep_run(f), NA)
  expect_true(file.exists(file.path(f$results_dir, "REPORT_run1_q1.html")))
})

test_that("the report says which loci were found in the query and whether they were declared", {
  skip_if_not_installed("Biostrings")
  f <- .rep_fixture(withr::local_tempdir())
  .rep_run(f)
  top <- sub("id=\"section-2\".*", "", .rep_html(f))
  expect_true(grepl("Loci found", top, fixed = TRUE))
  expect_true(grepl("matK, rbcL, trnL-trnF", top, fixed = TRUE))
  expect_true(grepl("detected by overlap with the library", top, fixed = TRUE))
  rec <- rbind(.rep_record(), data.frame(kind = "setting", name = "loci_declared", value = "matK"))
  g <- .rep_fixture(withr::local_tempdir(), record = rec)
  .rep_run(g)
  top <- sub("id=\"section-2\".*", "", .rep_html(g))
  expect_true(grepl("declared by the user: matK", top, fixed = TRUE))
})

test_that("the report ends with the version of PhyloCactus and how to cite it", {
  skip_if_not_installed("Biostrings")
  f <- .rep_fixture(withr::local_tempdir())
  .rep_run(f)
  h <- .rep_html(f)
  foot <- sub(".*<footer", "<footer", h)
  expect_true(grepl("<footer", h, fixed = TRUE))
  expect_true(grepl(as.character(utils::packageVersion("PhyloCactus")), foot, fixed = TRUE))
  expect_true(grepl("How to cite", foot, fixed = TRUE))
  expect_true(grepl("https://github.com/beeamerino/PhyloCactus", foot, fixed = TRUE))
  expect_false(grepl("????", foot, fixed = TRUE))
})

test_that("identify_barcoding_query() records the loci declared by the user, or none", {
  .rep_skip_idtaxa()
  f <- .rep_idn_library(withr::local_tempdir())
  q <- c(q1 = unname(f$seqs$matK[1]))
  suppressMessages(invisible(utils::capture.output(
    identify_barcoding_query(q, library_dir = f$library_dir, metrics_dir = NULL, output_dir = f$out, run_name = "a"))))
  suppressMessages(invisible(utils::capture.output(
    identify_barcoding_query(q, locus = "matK", library_dir = f$library_dir, metrics_dir = NULL, output_dir = f$out, run_name = "b"))))
  a <- utils::read.csv(file.path(f$out, "RUN_a.csv"), stringsAsFactors = FALSE)
  b <- utils::read.csv(file.path(f$out, "RUN_b.csv"), stringsAsFactors = FALSE)
  expect_identical(a$value[a$name == "loci_declared"], "none")
  expect_identical(b$value[b$name == "loci_declared"], "matK")
})
