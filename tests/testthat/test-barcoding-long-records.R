# Tests of the cut of long records in step 1 (Phase 10, decision L2 of BMM, 28-09). Written before
# the code. The mining keeps records of up to 5000 bases (L1); in every locus, a record longer than
# 2000 bases is cut, before alignment, to the extent that the records of 2000 bases or fewer of that
# locus cover, so that a locus keeps one extent (the trnK intron around matK, for instance).

.lr_random <- function(n) paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
.lr_vary <- function(s, k) {
  for (p in sample(seq_len(nchar(s)), k)) substr(s, p, p) <- sample(c("A", "C", "G", "T"), 1)
  s
}
.lr_rc <- function(s) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(s)))

# One locus: a core of 800 bases; ten short records covering parts of it; the long records hold the
# whole core between flanks of 700 bases that no short record holds.
.lr_fixture <- function() {
  set.seed(71L)
  core <- .lr_random(800)
  starts <- c(1, 1, 1, 101, 201, 1, 51, 1, 301, 1)
  ends <- c(800, 800, 700, 800, 800, 600, 800, 800, 800, 750)
  short <- vapply(seq_along(starts), function(i) .lr_vary(substr(core, starts[i], ends[i]), 8), "")
  names(short) <- paste0("S", seq_along(short), ".1")
  long_core <- .lr_vary(core, 8)
  long <- c(L1.1 = paste0(.lr_random(700), long_core, .lr_random(700)),
            L2.1 = .lr_rc(paste0(.lr_random(700), .lr_vary(core, 8), .lr_random(700))),
            L3.1 = .lr_random(2300))
  seqs <- c(short, long)
  registry <- data.frame(sid = names(seqs), species = paste0("Opuntia_sp", seq_along(seqs)),
                         genus = "Opuntia", locus = "matK", cluster_id = 1L,
                         genbank_name = paste0("Opuntia_sp", seq_along(seqs)), stringsAsFactors = FALSE)
  list(registry = registry, sequences = seqs, core = core, long_core = long_core)
}

test_that("a record over 2000 bases is cut to the extent of the short records of its locus (L2)", {
  f <- .lr_fixture()
  out <- .bc_cut_long_records(f$registry, f$sequences)
  s <- out$sequences[["L1.1"]]
  # The flanks are gone and the core is kept: the cut lies within a few bases of the core
  expect_lt(abs(nchar(s) - 800), 30)
  expect_true(grepl(substr(f$long_core, 30, 770), s, fixed = TRUE))
  tab <- out$table
  expect_true(all(c("sid", "locus", "length", "region_start", "region_end", "region_length", "action") %in% names(tab)))
  r <- tab[tab$sid == "L1.1", ]
  expect_identical(r$action, "cut")
  expect_equal(r$length, 2200L)
  expect_lt(abs(r$region_start - 701), 30)
  expect_identical(substr(f$sequences[["L1.1"]], r$region_start, r$region_end), s)
})

test_that("a long record on the other strand is cut on its own strand (L2)", {
  f <- .lr_fixture()
  out <- .bc_cut_long_records(f$registry, f$sequences)
  s <- out$sequences[["L2.1"]]
  expect_lt(abs(nchar(s) - 800), 30)
  r <- out$table[out$table$sid == "L2.1", ]
  expect_identical(substr(f$sequences[["L2.1"]], r$region_start, r$region_end), s)
})

test_that("a long record that shares nothing with the short records is left out and listed (L2)", {
  f <- .lr_fixture()
  out <- .bc_cut_long_records(f$registry, f$sequences)
  expect_false("L3.1" %in% out$registry$sid)
  expect_false("L3.1" %in% names(out$sequences))
  expect_identical(out$table$action[out$table$sid == "L3.1"], "no_overlap")
})

test_that("records of 2000 bases or fewer and loci without long records are left untouched (L2)", {
  f <- .lr_fixture()
  its <- c(I1.1 = .lr_random(600), I2.1 = .lr_random(650))
  reg <- rbind(f$registry, data.frame(sid = names(its), species = c("Cereus_a", "Cereus_b"), genus = "Cereus",
                                      locus = "ITS", cluster_id = 2L, genbank_name = c("Cereus_a", "Cereus_b"),
                                      stringsAsFactors = FALSE))
  out <- .bc_cut_long_records(reg, c(f$sequences, its))
  short <- paste0("S", 1:10, ".1")
  expect_identical(out$sequences[short], f$sequences[short])
  expect_identical(out$sequences[names(its)], its)
  expect_false(any(c(short, names(its)) %in% out$table$sid))
  expect_identical(nrow(out$registry), nrow(reg) - 1L)
})

test_that("a stray 20-mer of the core far in a flank does not stretch the cut (L2)", {
  f <- .lr_fixture()
  s <- f$sequences[["L1.1"]]
  substr(s, 100, 119) <- substr(f$long_core, 400, 419)
  f$sequences[["L1.1"]] <- s
  out <- .bc_cut_long_records(f$registry, f$sequences)
  expect_lt(abs(nchar(out$sequences[["L1.1"]]) - 800), 30)
})

test_that("step 1 cuts the long records before the extras and writes the table of the cuts (L2)", {
  b <- paste(deparse(body(assemble_barcoding_dataset)), collapse = "\n")
  expect_match(b, ".bc_cut_long_records(", fixed = TRUE)
  expect_match(b, "TABLE_barcoding_long_records_cut.csv", fixed = TRUE)
  expect_lt(regexpr(".bc_cut_long_records(", b, fixed = TRUE), regexpr(".bc_add_extra_records(", b, fixed = TRUE))
})
