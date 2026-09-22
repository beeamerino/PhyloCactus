# Tests of locus selection in the molecular diagnostic branch (11_barcoding/). Written before the code.
# Decisions of BMM, 2026-09-22 (Phase 2 proposal, sec. 10): the loci come out of the data, not out of a
# list. Steps 2 and 3 process every assembled locus and the threshold decides; clusters whose gene is not
# recognised in target_genes.txt stay in the funnel as cluster_<id>, flagged; every filter is counted.

test_that("clusters without a recognised gene keep a provisional name instead of being dropped", {
  clusters <- data.frame(ID = c("0", "5", "7"), Marker_std = c("trnL-trnF", NA, ""),
                         Description = c("seed 0 trnL-trnF", "seed 5 hypothetical protein", "seed 7 unknown"),
                         stringsAsFactors = FALSE)

  out <- .bc_name_clusters(clusters)

  expect_equal(out$map$ID, c("0", "5", "7"))
  expect_equal(out$map$Marker_std, c("trnL-trnF", "cluster_5", "cluster_7"))
  expect_equal(out$map$locus_provisional, c(FALSE, TRUE, TRUE))
  expect_equal(out$unnamed$ID, c("5", "7"))
  expect_equal(out$unnamed$locus, c("cluster_5", "cluster_7"))
  expect_equal(out$unnamed$descripcion_semilla, c("seed 5 hypothetical protein", "seed 7 unknown"))
})

test_that("with every cluster named, the unnamed table is empty but keeps its columns", {
  out <- .bc_name_clusters(data.frame(ID = "0", Marker_std = "matK", Description = "d", stringsAsFactors = FALSE))
  expect_equal(nrow(out$unnamed), 0L)
  expect_true(all(c("ID", "locus", "descripcion_semilla") %in% names(out$unnamed)))
  expect_false(out$map$locus_provisional)
})

test_that("provisional loci are recognised by their cluster_<id> name only", {
  tab <- .bc_flag_provisional(data.frame(locus = c("matK", "cluster_12", "cluster_x", "ITS_cluster_3"),
                                         stringsAsFactors = FALSE))
  expect_equal(tab$locus_provisional, c(FALSE, TRUE, FALSE, FALSE))
  expect_equal(.bc_flag_provisional(data.frame(locus = character(0)))$locus_provisional, logical(0))
})

test_that("the cluster funnel of step 1 counts every stage", {
  map <- data.frame(ID = c("0", "1", "5"), Marker_std = c("trnL-trnF", "trnL-trnF", "cluster_5"),
                    locus_provisional = c(FALSE, FALSE, TRUE), stringsAsFactors = FALSE)

  f <- .bc_cluster_funnel(n_workspace = 1371L, n_selected = 72L, n_with_sequences = 32L, marker_map = map)

  expect_equal(f$etapa, c("clusteres_en_workspace", "clusteres_con_mas_de_min_species",
                          "clusteres_con_secuencias_tras_filtros", "clusteres_con_nombre",
                          "clusteres_sin_nombre_conservados", "loci_ensamblados"))
  expect_equal(f$n, c(1371L, 72L, 32L, 2L, 1L, 2L))
})

test_that("the assembly names clusters and writes the funnel through the tested helpers", {
  b <- paste(deparse(body(assemble_barcoding_dataset)), collapse = "\n")
  expect_match(b, ".bc_name_clusters(", fixed = TRUE)
  expect_match(b, ".bc_cluster_funnel(", fixed = TRUE)
})

test_that("curation processes every assembled locus when no list is given", {
  skip_if_not_installed("Biostrings")
  skip_if(!nzchar(Sys.which("mafft")), "MAFFT not installed")

  tmp <- withr::local_tempdir()
  input_dir <- file.path(tmp, "1_assembly")
  dir.create(input_dir)
  s <- paste(rep("ACGT", 60), collapse = "")
  for (l in c("matK", "cluster_7")) {
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(`Opuntia_robusta|A1.1` = s, `Opuntia_robusta|A2.1` = s,
                                                          `Opuntia_stricta|B1.1` = s, `Opuntia_stricta|B2.1` = s)),
                                file.path(input_dir, paste0(l, ".fasta")))
  }

  suppressMessages(curate_barcoding_markers(input_dir, file.path(tmp, "2_curated")))

  aln <- file.path(tmp, "2_curated", "alignments", paste0("ALN_masked_final_", c("matK", "cluster_7"), ".fasta"))
  expect_true(all(file.exists(aln)))
})

test_that("screening screens every curated alignment when no list is given, and flags provisional loci", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("ape")

  tmp <- withr::local_tempdir()
  set.seed(5L)
  base <- paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = "")
  mutate <- function(s, pos) { substr(s, pos, pos) <- ifelse(substr(s, pos, pos) == "A", "C", "A"); s }
  seqs <- c(base, mutate(base, 10), mutate(base, 100), mutate(mutate(base, 100), 200))

  asm_dir <- file.path(tmp, "1_assembly")
  cur_dir <- file.path(tmp, "2_curated")
  dir.create(asm_dir)
  dir.create(file.path(cur_dir, "alignments"), recursive = TRUE)
  reg <- NULL
  for (l in c("matK", "cluster_7")) {
    sids <- paste0(l, "_", 1:4, ".1")
    sp <- rep(c("Opuntia_robusta", "Opuntia_stricta"), each = 2)
    reg <- rbind(reg, data.frame(sid = sids, species = sp, genus = "Opuntia", locus = l, cluster_id = 1L,
                                 nombre_genbank = sp, stringsAsFactors = FALSE))
    Biostrings::writeXStringSet(Biostrings::DNAStringSet(stats::setNames(seqs, paste0(sp, "|", sids))),
                                file.path(cur_dir, "alignments", paste0("ALN_masked_final_", l, ".fasta")))
  }
  utils::write.csv(reg, file.path(asm_dir, "TABLE_barcoding_accession_registry.csv"), row.names = FALSE)

  suppressMessages(screen_barcoding_markers(assembly_dir = asm_dir, curated_dir = cur_dir,
                                            output_dir = file.path(tmp, "3_screening"), min_species_with_replica = 2L))

  for (f in c("TABLE_barcoding_marker_screening.csv", "TABLE_barcoding_identical_sequences.csv")) {
    t <- utils::read.csv(file.path(tmp, "3_screening", f), stringsAsFactors = FALSE)
    expect_setequal(unique(t$locus), c("matK", "cluster_7"))
    expect_true(all(t$locus_provisional[t$locus == "cluster_7"]), info = f)
    expect_false(any(t$locus_provisional[t$locus == "matK"]), info = f)
  }
})

test_that("screening refuses to run when the curated directory holds no alignment", {
  tmp <- withr::local_tempdir()
  dir.create(file.path(tmp, "1_assembly"))
  utils::write.csv(data.frame(sid = "a", species = "Opuntia_robusta", genus = "Opuntia", locus = "matK",
                              cluster_id = 1L, nombre_genbank = "Opuntia_robusta"),
                   file.path(tmp, "1_assembly", "TABLE_barcoding_accession_registry.csv"), row.names = FALSE)
  dir.create(file.path(tmp, "2_curated", "alignments"), recursive = TRUE)
  expect_error(
    screen_barcoding_markers(assembly_dir = file.path(tmp, "1_assembly"), curated_dir = file.path(tmp, "2_curated"),
                             output_dir = file.path(tmp, "3_screening")),
    "curate_barcoding_markers"
  )
})
