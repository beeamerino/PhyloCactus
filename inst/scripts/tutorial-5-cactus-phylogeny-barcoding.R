# -------------------------------------------------------------
# PhyloCactus: Tutorial 5 - Molecular Diagnostic Branch (11_barcoding/)
# -------------------------------------------------------------
# Phases 2 to 6A of PhyloCactus 0.5.0: assembly, curation, screening, reference library,
# validation folds, barcode gap, classification of those folds, negative controls and the
# threshold of remoteness. The controls of step 8 come before any real figure is read.
#
# The branch reads the same phylotaR workspace as the phylogeny (0_phylotaR_raw_Ingroup/)
# and writes only under 11_barcoding/. The first run downloads GenBank metadata for the
# accessions the phylogeny's cache does not hold, so it needs a network connection; later
# runs read the caches in 11_barcoding/cache/.
#
# Every step reads its inputs from the directory of the step before it and writes its own
# directory (1_assembly, 2_curated, 3_screening, 4_library), as in the phylogeny: any step can be
# run on its own in a new session when the directories it reads exist.
# -------------------------------------------------------------
library(PhyloCactus)

tutorial_dir <- "~/Desktop/PhyloCactus_Tutorial"
setwd(tutorial_dir)

# -------------------------------------------------------------
# Switches of the long steps, declared once so the whole script is changed from one place.
# -------------------------------------------------------------
# Email when step 7 ends, whether it finished or failed. FALSE sends nothing and needs no
# configuration. TRUE needs MY_EMAIL in your .Renviron and a blastula credentials file created once
# with blastula::create_smtp_creds_file(). The package reads neither: it hands the path to blastula.
# A notification that cannot be sent is reported and ignored, so it never fails a run.
notify_email <- FALSE

# Folds per locus and scheme for the IdTaxa contrast, sampled with a fixed seed; loci with fewer
# folds run complete. NULL runs all 5928 folds.
#
# Measured on this machine (Apple M2 Pro), seconds per fold, scheme species: trnS-trnG 3.0,
# pepC_like 3.0, phyC 3.1, rbcL 5.5, psbJ-petA 5.5, ITS 6.1, atpB-rbcL 6.1, psbA-trnH 6.3,
# trnL-trnF 11.9, ndhF-rpl32 26.2, matK 28.4, rpL16 47.0. The cost does not follow the size of the
# library: ndhF-rpl32 is one of the smallest and takes eight times what phyC takes.
#
# With 100 the two schemes took 9 h 11 min and wrote 2773 predictions. All 5928 folds project to
# about 34 h at these rates. The nearest neighbour needs no sample: it runs in seconds.
idtaxa_max_folds <- 100

# -------------------------------------------------------------
# Step 1: Assembly. Clusters with more than min_species species; every accession kept; clusters
# whose gene is not recognised are kept as cluster_<id>. Headers Genus_species|sid.
# -------------------------------------------------------------
assemble_barcoding_dataset(
  wd_path = "0_phylotaR_raw_Ingroup",
  output_dir = "11_barcoding",
  target_genes_file = system.file("extdata", "target_genes.txt", package = "PhyloCactus"),
  genes_map_file = system.file("extdata", "genes_map.csv", package = "PhyloCactus"),
  manual_exclusions_file = system.file("extdata", "manual_exclusions_ingroup.csv", package = "PhyloCactus"),
  apply_manual_exclusions = TRUE,
  checklist_path = system.file("extdata", "CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx", package = "PhyloCactus"),
  min_species = 50,
  preferred_parent = "3593",
  force_download = FALSE
)

# -------------------------------------------------------------
# Step 2: Curation of every assembled locus, with the phylogeny's ingroup parameters
# -------------------------------------------------------------
curate_barcoding_markers(
  input_dir = "11_barcoding/1_assembly",
  output_dir = "11_barcoding/2_curated",
  mask_alignment_regions = TRUE,
  min_non_gap_fraction = 0.30,
  max_missing_fraction = 0.30,
  min_masked_alignment_length = 100L,
  preserve_iupac = TRUE,
  fix_strand = TRUE,
  mafft_exec = "mafft",
  mafft_opts = "--auto"
)

# -------------------------------------------------------------
# Step 3: Screening of every curated locus. The replication threshold decides which loci enter;
# saturation and paralogy are flagged, never excluding
# -------------------------------------------------------------
screen_barcoding_markers(
  assembly_dir = "11_barcoding/1_assembly",
  curated_dir = "11_barcoding/2_curated",
  output_dir = "11_barcoding/3_screening",
  min_species_with_replica = 20L,
  saturation_flag_cutoff = 0.3,
  paralog_loci = "pepC_like"   # flagged in every table, never excluded
)

# -------------------------------------------------------------
# Step 4: Final library. Non-homologous sequences out and declared,
# identical sequences collapsed within species, threshold evaluated again on distinct sequences
# -------------------------------------------------------------
finalize_barcoding_library(
  assembly_dir = "11_barcoding/1_assembly",
  curated_dir = "11_barcoding/2_curated",
  screening_dir = "11_barcoding/3_screening",
  output_dir = "11_barcoding/4_library",
  metadata_file = "11_barcoding/cache/CACHE_GENBANK_METADATA_BARCODING.csv",
  min_species_with_replica = 20L,
  paralog_loci = "pepC_like"
)

# -------------------------------------------------------------
# Step 5: Validation folds. Scheme E leaves one sequence out and asks for the species; scheme G
# leaves one whole species out and asks for the genus. Deterministic, and checked for leakage
# -------------------------------------------------------------
build_barcoding_folds(
  library_dir = "11_barcoding/4_library",
  output_dir = "11_barcoding/5_folds"
)

# -------------------------------------------------------------
# Step 6: Barcode gap per locus. Intraspecific distances by pairs and distance to the nearest
# neighbour of another species, in raw and K80; candidate threshold = 95th intraspecific percentile
# -------------------------------------------------------------
analyze_barcode_gap(
  library_dir = "11_barcoding/4_library",
  output_dir = "11_barcoding/6_gap",
  models = c("raw", "K80"),
  min_comparable = 100L,   # pairs resting on fewer positions have no value
  paralog_loci = "pepC_like",
  figures = TRUE
)

# -------------------------------------------------------------
# Step 7: Classification of the folds. Three states from the structure of the tie at the minimum
# distance, never from a threshold: the threshold is swept in Phase 6 over the scores this step
# writes. One row per query and nothing aggregated.
# -------------------------------------------------------------
# The branch classifier of the internal validation: exact leave-one-out over the distance matrix
# of the locus, and seconds of running time.
classify_barcoding_folds(
  library_dir = "11_barcoding/4_library",
  folds_dir = "11_barcoding/5_folds",
  output_dir = "11_barcoding/7_classifier",
  method = "nn",
  model = "raw",
  min_comparable = 100L
)

# The same classifier, with every query measured the way a query of a user is measured: stripped of
# gaps, oriented against the training set of its fold and added to it with MAFFT --add --keeplength,
# one call per query. Writes the tables with the suffix _add and the running time per locus. About
# 8 h on this machine (29 475 s on 2026-09-26); step 9 reads these tables.
classify_barcoding_folds(
  library_dir = "11_barcoding/4_library",
  folds_dir = "11_barcoding/5_folds",
  output_dir = "11_barcoding/7_classifier",
  method = "nn",
  alignment = "add",
  notify = notify_email
)

# The classifier of the real use: IdTaxa needs no alignment, so it is the one that can answer a
# query a user brings. It retrains once per fold, hence the sample declared in the SETUP block.
classify_barcoding_folds(
  library_dir = "11_barcoding/4_library",
  folds_dir = "11_barcoding/5_folds",
  output_dir = "11_barcoding/7_classifier",
  method = "idtaxa",
  max_folds = idtaxa_max_folds,
  seed = 1L,
  notify = notify_email
)

# -------------------------------------------------------------
# Step 8: Negative controls, before any real accuracy is read. CN1 permutes the labels and asks
# whether the folds leak; CN2 asks the library about sequences that are not cacti; CN3 measures the
# resubstitution bias on the same queries as step 7.
# -------------------------------------------------------------
# CN2 needs an outgroup, assembled from the raw phylotaR workspace of the outgroup with the filters
# of step 1. checklist_path = NA applies no checklist: NULL would load the checklist of Cactaceae
# and discard every outgroup name. No manual exclusions, and no minimum of species per cluster.
# The comparison table is written against the outgroup's own cluster-to-locus assignment.
assemble_barcoding_dataset(
  wd_path = "0_phylotaR_raw_Outgroup",
  output_dir = "11_barcoding/8_controls/outgroup",
  target_genes_file = system.file("extdata", "target_genes.txt", package = "PhyloCactus"),
  genes_map_file = system.file("extdata", "genes_map.csv", package = "PhyloCactus"),
  apply_manual_exclusions = FALSE,
  checklist_path = NA,
  min_species = 0,
  preferred_parent = "3593",
  force_download = FALSE,
  phylogeny_map_file = "1_phylotaR_out_Outgroup/TABLE_CLUSTER_MARKER_ASSIGNMENT_OUTGROUP.csv"
)

# CN2 aligns each outgroup query to the library with MAFFT --add --keeplength, one call per query,
# so MAFFT has to be available. The three controls take about 11 minutes.
run_barcoding_controls(
  library_dir = "11_barcoding/4_library",
  folds_dir = "11_barcoding/5_folds",
  output_dir = "11_barcoding/8_controls",
  outgroup_dir = "11_barcoding/8_controls/outgroup/1_assembly",
  method = "nn",
  permutations = 10L,
  seed = 1L
)

# -------------------------------------------------------------
# Step 9: Threshold of remoteness. A query whose distance to its nearest neighbour exceeds the
# threshold is not named (state 3, reason remoteness). The curve reports, for every threshold, what the
# legitimate queries become and how many outgroup queries are rejected; the operating threshold of
# each locus is the quantile 0.99 of scheme G, fixed without looking at the outgroup.
# -------------------------------------------------------------
sweep_barcoding_threshold(
  classifier_dir = "11_barcoding/7_classifier",
  controls_dir = "11_barcoding/8_controls",
  output_dir = "11_barcoding/9_threshold",
  alignment = "add",
  q = 0.99,
  quantile_type = 1L
)

# -------------------------------------------------------------
# Results
# -------------------------------------------------------------
print(utils::read.csv("11_barcoding/3_screening/TABLE_barcoding_marker_screening.csv"))
print(utils::read.csv("11_barcoding/4_library/TABLE_barcoding_library_summary.csv"))
# Accessions and decisions per locus, from assembly to library
print(utils::read.csv("11_barcoding/4_library/TABLE_barcoding_funnel.csv"))
# Species and genera inside and outside the denominator, per locus and scheme
print(utils::read.csv("11_barcoding/5_folds/TABLE_barcoding_folds_summary.csv"))
# Barcode gap per locus and model, with the candidate threshold and its decision
print(utils::read.csv("11_barcoding/6_gap/TABLE_barcoding_gap_summary.csv"))
# Negative controls: leakage verdict per locus and scheme, and the outgroup queries per locus
print(utils::read.csv("11_barcoding/8_controls/TABLE_barcoding_cn1_verdict.csv"))
print(utils::read.csv("11_barcoding/8_controls/TABLE_barcoding_cn2_outgroup.csv"))

# Threshold of remoteness per locus and scheme, with the outgroup rejection where it is measured
print(utils::read.csv("11_barcoding/9_threshold/TABLE_barcoding_threshold_operating.csv"))
