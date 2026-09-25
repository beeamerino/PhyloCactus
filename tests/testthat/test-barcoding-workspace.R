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

  get <- function(id) cmp$estado[cmp$ID == id]
  expect_equal(get("0"), "igual")
  expect_equal(get("11"), "solo_rama")
  expect_equal(get("45"), "solo_rama")
  expect_equal(get("7"), "distinto")
  expect_equal(get("9"), "solo_filogenia")
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
  seis <- c("107598", "107617", "107583", "3582", "107600", "108056")
  expect_identical(.bc_mining_taxids(seis, "3593"), seis)
})

test_that("a missing workspace is mined with every taxid it was given, not only the first", {
  seis <- c("107598", "107617", "107583", "3582", "107600", "108056")
  setup_args <- NULL
  n_read <- 0L
  obj <- .phylotar_load_or_mine(
    "wd", preferred_parent = "3593", txid = seis, ncbi_dr = "/opt/blast/bin",
    reader   = function(wd) {
      n_read <<- n_read + 1L
      if (n_read == 1L) stop("no workspace") else "PHYLOTA"
    },
    setup_fn = function(...) setup_args <<- list(...),
    run_fn   = function(...) invisible(NULL)
  )
  expect_identical(obj, "PHYLOTA")
  expect_identical(setup_args$txid, seis)
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
  seis <- c("107598", "107617", "107583", "3582", "107600", "108056")
  expect_identical(barcoding_outgroup_taxids(), seis)
  expect_identical(eval(formals(assemble_outgroup_phylotar)$outgroups), barcoding_outgroup_taxids())
  # And the branch can be told to mine them, without that changing what it mines by default
  expect_true("taxids" %in% names(formals(assemble_barcoding_dataset)))
  expect_null(formals(assemble_barcoding_dataset)$taxids)
})
