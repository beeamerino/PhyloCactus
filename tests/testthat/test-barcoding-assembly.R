# Tests of the molecular diagnostic branch (11_barcoding/), assembly step.
#
# Written before the functions they test (validation plan, sec. 7; executive document, Phase 2).
# The logic is tested on small tables, without phylotaR S4 objects: assemble_barcoding_dataset()
# only reads the workspace and hands plain data frames to these internal functions.

test_that("a sid in several clusters is kept in the cluster whose parent is the focal ingroup", {
  records <- data.frame(
    cluster_id = c(1L, 7L, 2L, 5L, 9L),
    sid        = c("A1.1", "A1.1", "B2.1", "B2.1", "C3.1"),
    stringsAsFactors = FALSE
  )
  parent_of <- c("1" = "186265", "7" = "3593", "2" = "186270", "5" = "186275", "9" = "3593")

  res <- .bc_resolve_cluster_duplicates(records, parent_of, preferred_parent = "3593")

  # A1.1: cluster 7 has the focal parent, so it wins over the lower identifier.
  # B2.1: no focal parent in either cluster, so the lower identifier wins.
  # C3.1: in one cluster only, untouched.
  expect_setequal(paste(res$kept$cluster_id, res$kept$sid), c("7 A1.1", "2 B2.1", "9 C3.1"))
  expect_setequal(paste(res$removed$cluster_id, res$removed$sid), c("1 A1.1", "5 B2.1"))
  expect_equal(nrow(res$kept) + nrow(res$removed), nrow(records))
})

test_that("manual exclusions act on (cluster, sid) pairs, and a pair that does not exist removes nothing", {
  records <- data.frame(
    cluster_id = c(3L, 3L, 4L),
    sid        = c("X1.1", "X2.1", "X3.1"),
    stringsAsFactors = FALSE
  )
  exclusions <- data.frame(
    cluster_id = c(3L, 3L, 12L),
    sid        = c("X1.1", "X3.1", "X9.1"),
    stringsAsFactors = FALSE
  )

  res <- .bc_apply_manual_exclusions(records, exclusions)

  # (3, X1.1) exists and is removed. (3, X3.1) names an existing sid in the wrong cluster, as
  # (3, HO058696.1) does in the real run, and must not remove X3.1 from cluster 4.
  expect_setequal(res$kept$sid, c("X2.1", "X3.1"))
  expect_equal(res$counts$n_file, 3L)
  expect_equal(res$counts$n_existing, 1L)
  expect_setequal(res$missing_pairs, c("3 X3.1", "12 X9.1"))
})

test_that("no reduction to one sequence per species: every accession of a species is kept", {
  records <- data.frame(
    cluster_id   = c(0L, 0L, 0L, 0L),
    sid          = c("S1.1", "S2.1", "S3.1", "T1.1"),
    locus        = "trnL-trnF",
    species      = c("Opuntia_ficus-indica", "Opuntia_ficus-indica", "Opuntia_ficus-indica", "Opuntia_robusta"),
    genbank_name = c("Opuntia_ficus-indica", "Opuntia_ficus-indica", "Opuntia_ficus-indica", "Opuntia_robusta"),
    stringsAsFactors = FALSE
  )

  registry <- .bc_build_registry(records)

  expect_equal(sum(registry$species == "Opuntia_ficus-indica"), 3L)
  expect_equal(nrow(registry), 4L)
})

test_that("the name rule accepts and rejects exactly what clean_taxonomic_names() does", {
  skip_if_not_installed("Biostrings")

  # GenBank names as cp_clean_species_name() writes them: dots removed, spaces to underscores.
  genbank <- c(
    "Opuntia_ficus-indica",                         # hyphenated epithet, accepted
    "Mammillaria_polyedra",                         # accepted
    "Cylindropuntia_fosbergii",                     # nothospecies: checklist writes the hybrid sign
    "Gymnocalycium_monvillei_subsp_achirasense",    # infraspecific: checklist keeps the dot
    "Lophophora_alberto-vojtechii"                  # not in the checklist
  )
  checklist_names <- c("Opuntia ficus-indica", "Mammillaria polyedra",
                       "Cylindropuntia \u00d7 fosbergii",
                       "Gymnocalycium monvillei subsp. achirasense")

  tmp <- withr::local_tempdir()
  fasta <- file.path(tmp, "locus.fasta")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(stats::setNames(rep(paste(rep("A", 120), collapse = ""), length(genbank)), genbank)),
    fasta
  )
  checklist <- file.path(tmp, "checklist.csv")
  utils::write.csv(data.frame(pureName = checklist_names), checklist, row.names = FALSE)
  out <- clean_taxonomic_names(raw_input_fasta = fasta, checklist_path = checklist,
                               output_clean_dir = file.path(tmp, "clean"))
  kept_phylogeny <- gsub("[[:space:]]+", "_", names(Biostrings::readDNAStringSet(out)))

  matched <- .bc_match_names(genbank, checklist_names)

  expect_setequal(stats::na.omit(matched), kept_phylogeny)
  expect_equal(matched, c("Opuntia_ficus-indica", "Mammillaria_polyedra", NA, NA, NA))
})

test_that("headers are Genus_species|sid, unique, and parse back without losing hyphenated epithets", {
  h <- .bc_header(c("Opuntia_ficus-indica", "Opuntia_ficus-indica", "Mammillaria_polyedra"),
                  c("AB1.1", "AB2.1", "CD3.2"))

  expect_equal(h, c("Opuntia_ficus-indica|AB1.1", "Opuntia_ficus-indica|AB2.1", "Mammillaria_polyedra|CD3.2"))
  expect_false(anyDuplicated(h) > 0)

  parsed <- .bc_parse_header(h)
  expect_equal(parsed$species, c("Opuntia_ficus-indica", "Opuntia_ficus-indica", "Mammillaria_polyedra"))
  expect_equal(parsed$sid, c("AB1.1", "AB2.1", "CD3.2"))
})

test_that("the registry has one row per (locus, sid) and refuses a sid in two loci", {
  base <- data.frame(
    cluster_id = c(0L, 0L, 3L),
    sid        = c("S1.1", "S2.1", "S3.1"),
    locus      = c("matK", "matK", "rbcL"),
    species    = c("Opuntia_robusta", "Opuntia_robusta", "Opuntia_robusta"),
    genbank_name = "Opuntia_robusta",
    stringsAsFactors = FALSE
  )

  registry <- .bc_build_registry(base)
  expect_named(registry, c("sid", "species", "genus", "locus", "cluster_id", "genbank_name"), ignore.order = TRUE)
  expect_equal(registry$genus, rep("Opuntia", 3L))
  expect_equal(as.integer(table(registry$locus)[c("matK", "rbcL")]), c(2L, 1L))

  in_two_loci <- rbind(base, data.frame(cluster_id = 4L, sid = "S1.1", locus = "rbcL",
                                        species = "Opuntia_robusta", genbank_name = "Opuntia_robusta"))
  expect_error(.bc_build_registry(in_two_loci), "more than one locus")
})

test_that("the per-locus summary counts replication within a locus, never across loci", {
  registry <- data.frame(
    sid     = c("a", "b", "c", "d", "e", "f"),
    species = c("Opuntia_robusta", "Opuntia_robusta", "Opuntia_stricta", "Opuntia_stricta", "Cereus_jamacaru", "Cereus_jamacaru"),
    genus   = c("Opuntia", "Opuntia", "Opuntia", "Opuntia", "Cereus", "Cereus"),
    locus   = c("matK", "matK", "matK", "rbcL", "matK", "rbcL"),
    stringsAsFactors = FALSE
  )

  s <- .bc_locus_summary(registry)
  mk <- s[s$locus == "matK", ]
  rb <- s[s$locus == "rbcL", ]

  # Opuntia stricta and Cereus jamacaru have one accession in each locus: no replica in either.
  expect_equal(mk$total_species, 3L)
  expect_equal(mk$species_with_replicate, 1L)
  expect_equal(mk$total_accessions, 4L)
  expect_equal(mk$total_genera, 2L)
  expect_equal(mk$genera_with_2plus_species, 1L)
  expect_equal(rb$species_with_replicate, 0L)
})

test_that("branch functions refuse to write into the phylogeny's directories", {
  tmp <- withr::local_tempdir()
  for (d in c("0_phylotaR_raw_Ingroup", "1_phylotaR_out_Ingroup", "4_Cleaned")) {
    expect_error(.bc_assert_output_dir(file.path(tmp, d, "sub")), "phylogeny")
  }
  expect_silent(.bc_assert_output_dir(file.path(tmp, "11_barcoding")))
})

test_that("writing the branch FASTA files leaves the input directories unchanged", {
  skip_if_not_installed("Biostrings")
  tmp <- withr::local_tempdir()
  input_dir <- file.path(tmp, "1_phylotaR_out_Ingroup")
  dir.create(input_dir)
  writeLines("ID,Marker_std\n0,matK", file.path(input_dir, "TABLE_CLUSTER_MARKER_ASSIGNMENT_INGROUP.csv"))
  md5_before <- tools::md5sum(list.files(input_dir, full.names = TRUE))

  registry <- data.frame(
    sid = c("S1.1", "S2.1"), species = "Opuntia_robusta", genus = "Opuntia",
    locus = "matK", cluster_id = 0L, genbank_name = "Opuntia_robusta", stringsAsFactors = FALSE
  )
  sequences <- c(S1.1 = "ACGTACGTAC", S2.1 = "ACGTACGTAA")
  out_dir <- file.path(tmp, "11_barcoding", "1_assembly")

  .bc_write_locus_fastas(registry, sequences, out_dir)

  written <- Biostrings::readDNAStringSet(file.path(out_dir, "matK.fasta"))
  expect_setequal(names(written), c("Opuntia_robusta|S1.1", "Opuntia_robusta|S2.1"))
  expect_equal(tools::md5sum(list.files(input_dir, full.names = TRUE)), md5_before)
})

# Added on 2026-09-26. When no name is discarded, table() over empty vectors lost both factor columns
# and the table was written with the header "accessions","" (outgroup of CN2, run of 2026-09-23).
test_that("the table of discarded names keeps its three columns when nothing is discarded", {
  empty_kept <- data.frame(genbank_name = character(0), locus = character(0), species = character(0),
                      stringsAsFactors = FALSE)
  d0 <- .bc_discarded_names(empty_kept)
  expect_identical(names(d0), c("genbank_name", "locus", "accessions"))
  expect_equal(nrow(d0), 0L)

  k <- data.frame(genbank_name = c("Opuntia_x", "Opuntia_x", "Cereus_y", "Cereus_z"),
                  locus = c("matK", "matK", "rbcL", "rbcL"),
                  species = c(NA, NA, NA, "Cereus_z"), stringsAsFactors = FALSE)
  d <- .bc_discarded_names(k)
  expect_identical(names(d), c("genbank_name", "locus", "accessions"))
  expect_equal(d$accessions[d$genbank_name == "Opuntia_x" & d$locus == "matK"], 2L)
  expect_equal(d$accessions[d$genbank_name == "Cereus_y"], 1L)
  expect_false("Cereus_z" %in% d$genbank_name)
})
