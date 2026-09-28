# Tests of the identification of paired Illumina reads through route A (Phase 8; decisions RA1 to RA6
# of BMM, 28-09). Written before the code. The external steps (fastp, GetOrganelle) go through the
# runner argument: here a fake runner writes the files the real tools would write, so the tests run
# offline without the tools. IdTaxa runs only under R CMD check and with the installed package (as in
# test-barcoding-identify.R).

.rds_runs <- function() {
  tryCatch({
    tr <- Biostrings::DNAStringSet(c(a = "ACGTACGTAACCGGTTAACC", b = "ACGTACGTACCCGGTTAACC",
                                     c = "TTTTGGGGAATTCCGGAATT", d = "TTTTGGGGAATTCCGGAATG"))
    lt <- DECIPHER::LearnTaxa(tr, c("Root;Opuntia;Opuntia_robusta", "Root;Opuntia;Opuntia_stricta",
                                    "Root;Cereus;Cereus_jamacaru", "Root;Cereus;Cereus_horrida"),
                              verbose = FALSE)
    DECIPHER::IdTaxa(Biostrings::DNAStringSet("ACGTACGTAACCGGTTAACC"), lt, strand = "top",
                     threshold = 60, processors = 1L, verbose = FALSE)
    TRUE
  }, error = function(e) FALSE)
}

.rds_skip <- function() {
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("Biostrings")
  if (!suppressWarnings(suppressMessages(.rds_runs()))) {
    skip("DECIPHER::IdTaxa() does not run in this loading context: match() over XStringSet fails to dispatch inside IdTaxa. The same call works in a plain R session.")
  }
}

.rds_random <- function(n) paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
.rds_vary <- function(s, k) {
  for (p in sample(seq_len(nchar(s)), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1)
  s
}

# The library of test-barcoding-identify.R, the two read files of a run and the seed databases
.rds_fixture <- function(tmp, read_length = 150L) {
  lib_dir <- file.path(tmp, "4_library")
  dir.create(lib_dir, showWarnings = FALSE, recursive = TRUE)
  set.seed(41L)
  genera <- c("Opuntia", "Cereus", "Mammillaria")
  epithets <- c("alpha", "beta")
  seqs <- list()
  for (l in c("matK", "rbcL", "trnL-trnF")) {
    root <- .rds_random(300)
    x <- character(0)
    for (g in seq_along(genera)) {
      gb <- .rds_vary(root, 6)
      for (e in seq_along(epithets)) {
        sb <- .rds_vary(gb, 3)
        for (r in 1:3) {
          x[paste0(genera[g], "_", epithets[e], "|", l, "_", substr(genera[g], 1, 1), e, r, ".1")] <-
            .rds_vary(sb, 1)
        }
      }
    }
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(x), file.path(lib_dir, paste0("LIB_", l, ".fasta")))
    seqs[[l]] <- x
  }
  reads <- file.path(tmp, c("run_1.fastq", "run_2.fastq"))
  for (f in reads) {
    writeLines(unlist(lapply(1:20, function(i) c(paste0("@r", i), .rds_random(read_length), "+",
                                                 strrep("I", read_length)))), f)
  }
  go <- file.path(tmp, "GetOrganelle")
  for (d in c("SeedDatabase", "LabelDatabase")) {
    dir.create(file.path(go, d), recursive = TRUE)
    for (db in c("embplant_pt", "embplant_mt", "embplant_nr")) writeLines(">x\nACGT", file.path(go, d, paste0(db, ".fasta")))
  }
  list(library_dir = lib_dir, seqs = seqs, reads = reads, getorg_path = go, out = file.path(tmp, "10_identify"), tmp = tmp)
}

# A fake runner: fastp copies the reads and writes its report; GetOrganelle writes a path holding
# matK and rbcL (plastome) or trnL-trnF (nrDNA) of Opuntia alpha between random flanks. With
# fail = "embplant_nr", that step exits 0 and writes no path, as GetOrganelle did on 27-09.
.rds_runner <- function(f, fail = character(0)) {
  calls <- character(0)
  runner <- function(cmd, args, stdout = "", stderr = "") {
    calls <<- c(calls, basename(cmd))
    arg <- function(flag) args[which(args == flag)[1] + 1L]
    if ("--version" %in% args) return(if (grepl("fastp", cmd)) "fastp 0.0.0" else "GetOrganelle v0.0.0")
    if (is.character(stderr) && nzchar(stderr)) writeLines("log of the fake step", stderr)
    if (grepl("fastp", cmd)) {
      file.copy(arg("-i"), arg("-o")); file.copy(arg("-I"), arg("-O"))
      writeLines("{}", arg("-j")); writeLines("<html></html>", arg("-h"))
      return(0L)
    }
    db <- arg("-F"); o <- arg("-o")
    dir.create(o, recursive = TRUE, showWarnings = FALSE)
    writeLines("work", file.path(o, "extended_1_paired.fq"))
    if (db %in% fail) return(0L)
    set.seed(43L)
    s <- if (db == "embplant_pt") {
      c(paste0(.rds_random(3000), unname(f$seqs$matK[1]), .rds_random(3000), unname(f$seqs$rbcL[1]), .rds_random(3000)))
    } else c(paste0(.rds_random(2000), unname(f$seqs[["trnL-trnF"]][1]), .rds_random(2000)))
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(contig = s)),
                                file.path(o, paste0(db, ".K105.complete.graph1.1.path_sequence.fasta")))
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(contig = .rds_random(500))),
                                file.path(o, paste0(db, ".K105.complete.graph1.2.path_sequence.fasta")))
    0L
  }
  list(runner = runner, calls = function() calls)
}

.rds_run <- function(f, runner, ...) {
  res <- NULL
  msg <- character(0)
  txt <- suppressWarnings(utils::capture.output(withCallingHandlers(
    res <- identify_barcoding_reads(f$reads[1], f$reads[2], library_dir = f$library_dir, metrics_dir = NULL,
                                    output_dir = f$out, run_name = "run", getorg_path = f$getorg_path,
                                    runner = runner, ...),
    message = function(m) {
      msg <<- c(msg, conditionMessage(m))
      invokeRestart("muffleMessage")
    })))
  list(tab = res, text = c(txt, msg), dir = file.path(f$out, "run"))
}

test_that("the paths of both assemblies are identified, one row per locus, with contig and organelle", {
  .rds_skip()
  f <- .rds_fixture(withr::local_tempdir())
  x <- .rds_run(f, .rds_runner(f)$runner)
  r <- x$tab
  expect_setequal(r$locus, c("matK", "rbcL", "trnL-trnF"))
  expect_equal(nrow(r), 3L)
  expect_true(all(c("contig", "organelle") %in% names(r)))
  expect_identical(r$organelle[r$locus == "matK"], "plastome")
  expect_identical(r$organelle[r$locus == "trnL-trnF"], "nrdna")
  expect_true(all(r$predicted_genus == "Opuntia"))
  expect_true(all(r$state %in% 1:2))
  for (fl in c("fastp.json", "plastome_path.fasta", "nrdna_path.fasta", "MANIFEST.csv",
               "TABLE_barcoding_identify_run.csv")) {
    expect_true(file.exists(file.path(x$dir, fl)), info = fl)
  }
  # The graph1.1 path is the one kept
  pp <- Biostrings::readDNAStringSet(file.path(x$dir, "plastome_path.fasta"))
  expect_equal(length(pp), 1L)
  expect_gt(Biostrings::width(pp), 9000L)
})

test_that("a GetOrganelle step that exits 0 without a path is reported, not read as absence", {
  .rds_skip()
  f <- .rds_fixture(withr::local_tempdir())
  x <- .rds_run(f, .rds_runner(f, fail = "embplant_nr")$runner)
  t <- x$tab[x$tab$locus == "trnL-trnF", ]
  expect_equal(t$state, 3L)
  expect_identical(t$reason, "assembly_failed")
  expect_true(any(grepl("embplant_nr", x$text) & grepl("no path", x$text)))
  expect_true(file.exists(file.path(x$dir, "STDERR_getorganelle_embplant_nr.txt")))
  expect_true(all(x$tab$state[x$tab$locus != "trnL-trnF"] %in% 1:2))
})

test_that("a missing seed database stops with the command that installs it, before any step runs", {
  skip_if_not_installed("Biostrings")
  f <- .rds_fixture(withr::local_tempdir())
  unlink(file.path(f$getorg_path, "SeedDatabase", "embplant_mt.fasta"))
  fr <- .rds_runner(f)
  expect_error(.rds_run(f, fr$runner), "get_organelle_config.py -a embplant_pt,embplant_mt,embplant_nr", fixed = TRUE)
  expect_length(fr$calls(), 0L)
})

test_that("single-end reads and long reads stop with their message", {
  skip_if_not_installed("Biostrings")
  f <- .rds_fixture(withr::local_tempdir())
  expect_error(identify_barcoding_reads(f$reads[1], NULL, library_dir = f$library_dir, output_dir = f$out,
                                        getorg_path = f$getorg_path, runner = .rds_runner(f)$runner),
               "Single-end")
  g <- .rds_fixture(withr::local_tempdir(), read_length = 5000L)
  expect_error(identify_barcoding_reads(g$reads[1], g$reads[2], library_dir = g$library_dir, output_dir = g$out,
                                        getorg_path = g$getorg_path, runner = .rds_runner(g)$runner),
               "long reads")
})

test_that("without keep_intermediate only the files of M3 stay, and the input reads are untouched", {
  .rds_skip()
  f <- .rds_fixture(withr::local_tempdir())
  md5_in <- unname(tools::md5sum(f$reads))
  x <- .rds_run(f, .rds_runner(f)$runner)
  kept <- list.files(x$dir, recursive = TRUE)
  expect_false(any(grepl("^clean_", kept)))
  expect_false(any(grepl("^getorganelle_", kept)))
  expect_false("fastp.html" %in% kept)
  expect_identical(unname(tools::md5sum(f$reads)), md5_in)

  g <- .rds_fixture(withr::local_tempdir())
  y <- .rds_run(g, .rds_runner(g)$runner, keep_intermediate = TRUE)
  kept <- list.files(y$dir, recursive = TRUE)
  expect_true(any(grepl("^clean_", kept)))
  expect_true(any(grepl("^getorganelle_embplant_pt/", kept)))
})

test_that("the manifest gives the md5 of the inputs and outputs and the versions of the tools", {
  .rds_skip()
  f <- .rds_fixture(withr::local_tempdir())
  x <- .rds_run(f, .rds_runner(f)$runner)
  m <- utils::read.csv(file.path(x$dir, "MANIFEST.csv"), stringsAsFactors = FALSE)
  expect_true(all(c("kind", "name", "value") %in% names(m)))
  expect_identical(m$value[m$kind == "input" & m$name == basename(f$reads[1])], unname(tools::md5sum(f$reads[1])))
  expect_identical(m$value[m$kind == "output" & m$name == "plastome_path.fasta"],
                   unname(tools::md5sum(file.path(x$dir, "plastome_path.fasta"))))
  expect_true(any(m$kind == "tool" & grepl("fastp", m$value)))
  expect_true(any(m$kind == "tool" & grepl("GetOrganelle", m$value)))
  expect_true(any(m$kind == "command" & grepl("-F embplant_pt", m$value)))
})

test_that("a missing external tool stops with its name before anything runs", {
  skip_if_not_installed("Biostrings")
  f <- .rds_fixture(withr::local_tempdir())
  expect_error(identify_barcoding_reads(f$reads[1], f$reads[2], library_dir = f$library_dir, output_dir = f$out,
                                        getorg_path = f$getorg_path, fastp = "no_such_fastp"),
               "no_such_fastp")
  expect_false(dir.exists(f$out))
})
