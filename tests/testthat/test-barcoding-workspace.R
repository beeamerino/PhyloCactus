# Tests of the shared phylotaR workspace step and of the cluster-to-locus assignment used by both
# branches (decisions A1 and B1 of BMM, 2026-09-21). Written before the functions they test.
#
# A1: the phylotaR download is one internal function, called by assemble_ingroup_phylotar() and by
# assemble_barcoding_dataset(), so the two branches make literally the same call and share
# 0_phylotaR_raw_Ingroup/ whichever runs first. B1: the barcoding branch assigns clusters to loci
# itself, with the same enrichment code as the phylogeny and its own metadata cache.

test_that("the phylotaR search term is the one the phylogeny has always used", {
  expected <- paste0(
    "NOT predicted[TI] ",
    "NOT \"whole genome shotgun\"[TI] ",
    "NOT unverified[TI] ",
    "NOT \"synthetic construct\"[Organism] ",
    "NOT refseq[filter] ",
    "NOT TSA[Keyword] ",
    "NOT \"sp.\"[TI] ",
    "NOT \"sp.\"[Organism] ",
    "NOT \"sp\"[Organism] ",
    "NOT \"sp\"[Organism] ",
    "NOT \"aff.\"[TI] ",
    "NOT \"aff\"[Organism] ",
    "NOT \"cf.\"[TI] ",
    "NOT \"cf\"[Organism] ",
    "NOT \"var.\"[TI] ",
    "NOT \"var\"[TI] ",
    "NOT \"var\"[Organism] ",
    "NOT \"var.\"[Organism] ",
    "NOT \"variety\"[TI] ",
    "NOT \"subsp.\"[TI] ",
    "NOT \"subsp\"[TI] ",
    "NOT \"subsp.\"[Organism] ",
    "NOT \"subsp\"[Organism] ",
    "NOT \"subspecies\"[Organism] ",
    "NOT \"x\"[Organism] ",
    "NOT \" x \"[Organism]"
  )
  expect_identical(.phylotar_search_terms(), expected)
})

test_that("an existing workspace is read and nothing is downloaded", {
  calls <- character(0)
  obj <- .phylotar_load_or_mine(
    "wd", preferred_parent = "3593",
    reader   = function(wd) { calls <<- c(calls, "read"); "PHYLOTA" },
    setup_fn = function(...) { calls <<- c(calls, "setup"); invisible(NULL) },
    run_fn   = function(...) { calls <<- c(calls, "run"); invisible(NULL) }
  )
  expect_identical(obj, "PHYLOTA")
  expect_identical(calls, "read")
})

test_that("a missing workspace is mined with the shared parameters, then read", {
  calls <- character(0)
  setup_args <- NULL
  n_read <- 0L
  obj <- .phylotar_load_or_mine(
    "wd", preferred_parent = "3593", ncbi_dr = "/opt/blast/bin",
    reader   = function(wd) {
      n_read <<- n_read + 1L
      calls <<- c(calls, "read")
      if (n_read == 1L) stop("no workspace") else "PHYLOTA"
    },
    setup_fn = function(...) { setup_args <<- list(...); calls <<- c(calls, "setup") },
    run_fn   = function(...) { calls <<- c(calls, "run") }
  )
  expect_identical(obj, "PHYLOTA")
  expect_identical(calls, c("read", "setup", "run", "read"))
  expect_identical(setup_args$wd, "wd")
  expect_identical(setup_args$txid, "3593")
  expect_identical(setup_args$ncbi_dr, "/opt/blast/bin")
  expect_identical(setup_args$mncvrg, 80)
  expect_identical(setup_args$srch_trm, .phylotar_search_terms())
})

test_that("force_download mines even when a workspace exists", {
  calls <- character(0)
  .phylotar_load_or_mine(
    "wd", force_download = TRUE, ncbi_dr = "/opt/blast/bin",
    reader   = function(wd) { calls <<- c(calls, "read"); "PHYLOTA" },
    setup_fn = function(...) calls <<- c(calls, "setup"),
    run_fn   = function(...) calls <<- c(calls, "run")
  )
  expect_identical(calls, c("setup", "run", "read"))
})

test_that("both branches obtain the workspace and the locus assignment through the same functions", {
  body_of <- function(f) paste(deparse(body(f)), collapse = "\n")
  for (f in list(assemble_ingroup_phylotar, assemble_barcoding_dataset)) {
    b <- body_of(f)
    expect_match(b, ".phylotar_load_or_mine(", fixed = TRUE)
    expect_match(b, ".cp_enrich_cluster_summary(", fixed = TRUE)
  }
})

test_that("the cluster enrichment maps the dominant marker through genes_map and keeps unmapped names", {
  smmry <- tibble::tibble(ID = c(0L, 1L), Seed = c("S0.1", "S1.1"),
                          top_marker = c("trnL-trnF", "rpS3"))
  seed_meta <- tibble::tibble(Seed = c("S0.1", "S1.1"),
                              Species = c("Opuntia robusta", "Cereus jamacaru"),
                              Description = c("Opuntia robusta trnL-trnF intergenic spacer", "Cereus jamacaru rpS3 gene"))
  markers <- c("trnL-trnF", "rpS3")
  lookup <- data.frame(search = c("trnL-trnF"), replace = c("trnL-trnF"), stringsAsFactors = FALSE)

  out <- .cp_enrich_cluster_summary(
    smmry, seed_meta, pattern = cp_build_pattern(markers),
    gene_lookup = dplyr::select(lookup, search_gene = search, Gene_std = replace),
    marker_lookup = dplyr::select(lookup, search_marker = search, Marker_std = replace)
  )

  expect_equal(out$Marker_std, c("trnL-trnF", "rpS3"))
  expect_equal(out$Species, c("Opuntia robusta", "Cereus jamacaru"))
  expect_equal(out$Genes_text, c("trnL-trnF", "rpS3"))
})

test_that("seed metadata are cached in the branch and only missing seeds are fetched", {
  tmp <- withr::local_tempdir()
  cache <- file.path(tmp, "cache", "CACHE_SEED_METADATA_BARCODING.csv")
  fetched <- list()
  fetcher <- function(sids, ...) {
    fetched[[length(fetched) + 1L]] <<- sids
    tibble::tibble(Seed = sids, Species = paste("sp", sids), Description = paste("desc", sids))
  }

  a <- .bc_seed_metadata_cached(c("S1.1", "S2.1"), cache, fetcher = fetcher)
  b <- .bc_seed_metadata_cached(c("S2.1", "S3.1", "S1.1"), cache, fetcher = fetcher)

  expect_equal(fetched, list(c("S1.1", "S2.1"), "S3.1"))
  expect_equal(b$Seed, c("S2.1", "S3.1", "S1.1"))
  expect_equal(b$Description, c("desc S2.1", "desc S3.1", "desc S1.1"))
  expect_true(file.exists(cache))
})

test_that("the branch assignment is compared with the phylogeny map, and every difference is reported", {
  branch    <- data.frame(ID = c("0", "11", "45", "7"), Marker_std = c("trnL-trnF", "rbcL", "rbcL", "matK"), stringsAsFactors = FALSE)
  phylogeny <- data.frame(ID = c("0", "7", "9"), Marker_std = c("trnL-trnF", "ITS", "rpL16"), stringsAsFactors = FALSE)

  cmp <- .bc_compare_marker_maps(branch, phylogeny)

  get <- function(id) cmp$state[cmp$ID == id]
  expect_equal(get("0"), "same")
  expect_equal(get("11"), "branch_only")
  expect_equal(get("45"), "branch_only")
  expect_equal(get("7"), "different")
  expect_equal(get("9"), "phylogeny_only")
  expect_equal(nrow(cmp), 5L)
})

# Decision of BMM, 2026-09-23: the branch has to run without the phylogeny branch. Its outgroup, the
# one CN2 queries with, therefore has to be obtainable the same way the ingroup is, by mining GenBank
# or reading the cache, and not by assuming that somebody already ran the phylogeny.
#
# The phylogeny mines its outgroup with six genus-level taxids across three families, with the
# taxonomic reasoning written in Tutorial 1. phylotaR accepts a vector and turns it into its own
# multiple_ids, which is why the workspace of 2026-09-04 was built that way. The branch needs the
# same, and it needs it without confusing two different things that `preferred_parent` was doing at
# once: what to mine, and which parent wins when a sid sits in two clusters.

test_that("what is mined and which parent resolves duplicates are two different things", {
  # Nothing given: the branch mines its focal clade, as it has since Phase 2
  expect_identical(.bc_mining_taxids(NULL, "3593"), "3593")
  # Given: those are mined, and the preferred parent is left out of it
  six_taxids <- c("107598", "107617", "107583", "3582", "107600", "108056")
  expect_identical(.bc_mining_taxids(six_taxids, "3593"), six_taxids)
})

test_that("a missing workspace is mined with every taxid it was given, not only the first", {
  six_taxids <- c("107598", "107617", "107583", "3582", "107600", "108056")
  setup_args <- NULL
  n_read <- 0L
  obj <- .phylotar_load_or_mine(
    "wd", preferred_parent = "3593", txid = six_taxids, ncbi_dr = "/opt/blast/bin",
    reader   = function(wd) {
      n_read <<- n_read + 1L
      if (n_read == 1L) stop("no workspace") else "PHYLOTA"
    },
    setup_fn = function(...) setup_args <<- list(...),
    run_fn   = function(...) invisible(NULL)
  )
  expect_identical(obj, "PHYLOTA")
  expect_identical(setup_args$txid, six_taxids)
  # And the search terms are the shared ones: the outgroup is not mined with a different criterion
  expect_identical(setup_args$srch_trm, .phylotar_search_terms())
})

test_that("txid defaults to preferred_parent, so nothing that already worked changes", {
  setup_args <- NULL
  n_read <- 0L
  .phylotar_load_or_mine(
    "wd", preferred_parent = "3593",
    reader   = function(wd) {
      n_read <<- n_read + 1L
      if (n_read == 1L) stop("no workspace") else "PHYLOTA"
    },
    setup_fn = function(...) setup_args <<- list(...),
    run_fn   = function(...) invisible(NULL)
  )
  expect_identical(setup_args$txid, "3593")
})

test_that("the branch holds its own copy of the outgroup taxids, and notices if the two drift apart", {
  # Decision of BMM, 2026-09-23: the branch keeps its own copy rather than reading the phylogeny's
  # default, so that it runs without the phylogeny and so that the list is visible where it is used.
  # Nothing of the phylogeny is touched. This test is the drift detector: the day somebody adds a
  # family to one of the two lists, it fails and says which one moved.
  six_taxids <- c("107598", "107617", "107583", "3582", "107600", "108056")
  expect_identical(barcoding_outgroup_taxids(), six_taxids)
  expect_identical(eval(formals(assemble_outgroup_phylotar)$outgroups), barcoding_outgroup_taxids())
  # And the branch can be told to mine them, without that changing what it mines by default
  expect_true("taxids" %in% names(formals(assemble_barcoding_dataset)))
  expect_null(formals(assemble_barcoding_dataset)$taxids)
})

# Phase 10 (decisions L1, L4 and E1 of BMM, 28-09)
.ws_reader_once_missing <- function() {
  n_read <- 0L
  function(wd) { n_read <<- n_read + 1L; if (n_read == 1L) stop("no workspace") else "PHYLOTA" }
}

test_that("the mining keeps records of 100 to 5000 bases (L1)", {
  setup_args <- NULL
  .phylotar_load_or_mine("wd", reader = .ws_reader_once_missing(),
                         setup_fn = function(...) setup_args <<- list(...), run_fn = function(...) NULL)
  expect_identical(setup_args$mnsql, 100L)
  expect_identical(setup_args$mxsql, 5000L)
})

test_that("the outgroup of the phylogeny is mined through the shared function, with the same parameters (E1)", {
  b <- paste(deparse(body(assemble_outgroup_phylotar)), collapse = "\n")
  expect_match(b, ".phylotar_load_or_mine(", fixed = TRUE)
  expect_false(grepl("phylotaR::setup(", b, fixed = TRUE))
})

test_that("a mining sends one notification when it ends, and none when the workspace is only read (L4)", {
  sent <- list()
  nf <- function(subject, body, ...) { sent[[length(sent) + 1L]] <<- list(subject = subject, body = body, ...); invisible(TRUE) }
  .phylotar_load_or_mine("0_phylotaR_raw_Ingroup", notify = TRUE, notify_fn = nf, reader = .ws_reader_once_missing(),
                         setup_fn = function(...) NULL, run_fn = function(...) NULL)
  expect_length(sent, 1L)
  expect_match(sent[[1]]$subject, "finished", fixed = TRUE)
  expect_match(sent[[1]]$subject, "phylotaR mining", fixed = TRUE)
  expect_match(sent[[1]]$subject, "\U0001f335", fixed = TRUE)
  expect_match(sent[[1]]$body, "0_phylotaR_raw_Ingroup", fixed = TRUE)

  sent <- list()
  .phylotar_load_or_mine("wd", notify = TRUE, notify_fn = nf, reader = function(wd) "PHYLOTA",
                         setup_fn = function(...) NULL, run_fn = function(...) NULL)
  expect_length(sent, 0L)
  # notify = FALSE, the default, sends nothing
  .phylotar_load_or_mine("wd", notify_fn = nf, reader = .ws_reader_once_missing(),
                         setup_fn = function(...) NULL, run_fn = function(...) NULL)
  expect_length(sent, 0L)
})

test_that("a mining that fails sends the failure and still stops (L4)", {
  sent <- list()
  nf <- function(subject, body, ...) { sent[[length(sent) + 1L]] <<- list(subject = subject, body = body); invisible(TRUE) }
  expect_error(
    .phylotar_load_or_mine("wd", notify = TRUE, notify_fn = nf, reader = function(wd) stop("no workspace"),
                           setup_fn = function(...) NULL, run_fn = function(...) NULL),
    "Could not load")
  expect_length(sent, 1L)
  expect_match(sent[[1]]$subject, "FAILED", fixed = TRUE)
})

test_that("the three mining functions take the notification arguments (L4)", {
  for (f in list(assemble_ingroup_phylotar, assemble_outgroup_phylotar, assemble_barcoding_dataset)) {
    expect_true(all(c("notify", "notify_to", "notify_credentials") %in% names(formals(f))))
    expect_false(eval(formals(f)$notify))
  }
})

# Phase 10, Y1 (BMM, 29-09): the cut of min_species applies to the locus, the clusters named alike
# pooled, in both branches; a cluster without a locus name is judged alone.
.y1_species <- function(cid, sp) data.frame(cluster_id = cid, species = sp, stringsAsFactors = FALSE)

test_that("clusters of one locus are pooled before the cut of min_species (Y1)", {
  sp <- rbind(.y1_species(1L, paste0("s", 1:30)), .y1_species(2L, paste0("s", 31:60)),
              .y1_species(3L, paste0("t", 1:30)), .y1_species(4L, paste0("t", 21:50)),
              .y1_species(5L, paste0("u", 1:60)), .y1_species(6L, paste0("v", 1:40)),
              .y1_species(7L, paste0("w", 1:40)))
  loc <- data.frame(cluster_id = 1:7, locus = c("psbA-trnH", "psbA-trnH", "rpL16", "rpL16", NA, NA, ""),
                    stringsAsFactors = FALSE)
  out <- .cp_select_clusters_by_locus(sp, loc, min_species = 50)
  # psbA-trnH 60 species in two clusters of 30: kept; rpL16 50 species (overlap): not over 50
  expect_setequal(out$keep, c(1L, 2L, 5L))
  tab <- out$table
  expect_equal(tab$species[tab$locus == "psbA-trnH"], 60L)
  expect_equal(tab$species[tab$locus == "rpL16"], 50L)
  # Unnamed clusters are not pooled with one another
  expect_false(6L %in% out$keep)
  expect_false(7L %in% out$keep)
})

test_that("the locus of every cluster is read from the definition lines of its sequences (Y1)", {
  skip_if_not_installed("phylotaR")
  data("aotus", package = "phylotaR", envir = environment())
  lookup <- data.frame(search_marker = "cytb", Marker_std = "cytb", stringsAsFactors = FALSE)
  gm <- data.frame(search = "cytb", replace = "cytb", stringsAsFactors = FALSE)
  loc <- .cp_cluster_loci(aotus, pattern = cp_build_pattern(c("cytb", "cytochrome b")), genes_map_df = gm,
                          marker_lookup = lookup)
  expect_setequal(names(loc), c("cluster_id", "locus"))
  expect_equal(nrow(loc), length(aotus@cids))
  expect_true(any(loc$locus %in% "cytb"))
})

test_that("both branches select clusters by locus and translate the manual exclusions (Y1, X2)", {
  for (f in list(assemble_ingroup_phylotar, assemble_barcoding_dataset)) {
    b <- paste(deparse(body(f)), collapse = "\n")
    expect_match(b, ".cp_select_clusters_by_locus(", fixed = TRUE)
    expect_match(b, ".cp_exclusion_pairs(", fixed = TRUE)
    expect_false(grepl("ntaxa > min_species", b, fixed = TRUE))
  }
})

# X2 (BMM, 29-09): the curated exclusions are keyed by locus and accession, so that they survive a
# new mining that renumbers the clusters.
test_that("exclusions by (locus, sid) become the (cluster, sid) pairs of the present workspace (X2)", {
  records <- data.frame(cluster_id = c(10L, 11L, 12L, 12L, 20L), sid = c("A.1", "A.1", "B.1", "C.1", "B.1"),
                        stringsAsFactors = FALSE)
  loc <- data.frame(cluster_id = c(10L, 11L, 12L, 20L), locus = c("ITS", "ITS", "matK", "rbcL"),
                    stringsAsFactors = FALSE)
  exc <- data.frame(locus = c("ITS", "matK", "trnL-trnF"), sid = c("A.1", "B.1", "Z.1"),
                    reason = c("not homologous", "chimeric", "gone"), stringsAsFactors = FALSE)
  p <- .cp_exclusion_pairs(exc, records, loc)
  # A.1 leaves both ITS clusters; B.1 leaves matK and stays in rbcL; Z.1 matches nothing
  expect_setequal(paste(p$pairs$cluster_id, p$pairs$sid), c("10 A.1", "11 A.1", "12 B.1"))
  expect_identical(p$pairs$reason[p$pairs$sid == "B.1"], "chimeric")
  expect_equal(p$counts$n_file, 3L)
  expect_equal(p$counts$n_matched, 2L)
  expect_identical(p$unmatched$sid, "Z.1")
})

test_that("the shipped exclusion lists are keyed by locus and accession (X2)", {
  for (f in c("manual_exclusions_ingroup.csv", "manual_exclusions_outgroup.csv")) {
    x <- utils::read.csv(system.file("extdata", f, package = "PhyloCactus"), stringsAsFactors = FALSE)
    expect_true(all(c("locus", "sid", "reason") %in% names(x)), info = f)
    expect_false("cluster_id" %in% names(x), info = f)
    expect_true(all(nzchar(x$locus) & !is.na(x$locus)), info = f)
    expect_false(any(duplicated(paste(x$locus, x$sid))), info = f)
  }
})
