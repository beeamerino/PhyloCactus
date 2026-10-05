# -------------------------------------------------------------
# PhyloCactus: Tutorial 5 - Molecular Diagnostic Branch (11_barcoding/)
# -------------------------------------------------------------
# Phases 2 to 11 of PhyloCactus 0.5.0: assembly (with own data, if any), curation, screening,
# reference library, validation folds, barcode gap, classification of those folds, negative
# controls, the threshold of remoteness and the metrics per class. The controls of step 8 come
# before any real figure is read.
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
# Email when steps 1 and 7 end, whether they finished or failed. FALSE sends nothing and needs no
# configuration. TRUE needs MY_EMAIL in your .Renviron and a blastula credentials file created once
# with blastula::create_smtp_creds_file(). The package reads neither: it hands the path to blastula.
# A notification that cannot be sent is reported and ignored, so it never fails a run.
notify_email <- FALSE

# IdTaxa retrains once per fold, which on the whole library takes tens of core-hours: it is run
# once, on a computing cluster, on the frozen library (generate_barcoding_job_scripts() and
# merge_barcoding_chunks(); CN2 with run_barcoding_controls(method = "idtaxa", controls = "CN2")),
# and its merged tables are the result, since LearnTaxa() differs between machines in floating
# point. This script never trains IdTaxa. NULL leaves IdTaxa out; the path of the folder holding
# the merged tables (TABLE_barcoding_predictions_<scheme>_idtaxa.csv, TABLE_barcoding_timing_idtaxa.csv
# and TABLE_barcoding_cn2_queries_idtaxa.csv) copies them into steps 7 and 8, and steps 9 and 10
# then include IdTaxa.
idtaxa_tables <- NULL

# Parallel processes of this machine for the longest step of the laptop, the nearest neighbour
# through MAFFT --add (step 7). Two cores are left free. Each process writes its log in
# 11_barcoding/7_classifier/chunks/, the progress is in 11_barcoding/7_classifier/PROGRESS_nn_add.txt
# (from a terminal: cat or watch that file), and one email is sent at the end.
n_workers <- max(1L, parallel::detectCores() - 2L)

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
  force_download = FALSE,
  notify = notify_email,
  # Own data (O1): NULL leaves them out; "0_own" adds the records of 0_own/ (SAMPLES.csv and one
  # FASTA per locus, headers >sample). What cannot enter is listed in
  # 11_barcoding/1_assembly/TABLE_barcoding_own_left_out.csv.
  own_dir = NULL
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
# one call per query. Writes the tables with the suffix _add and the running time per locus. The
# folds are split among n_workers processes and merged into the same tables as one process would
# write. Step 9 reads these tables. Measured on 2026-09-30, Apple M2 Pro, 8 workers, 11 417 queries:
# 3 h 12 min, 24.7 hours summed over the processes, matK 61 % of them.
classify_barcoding_folds(
  library_dir = "11_barcoding/4_library",
  folds_dir = "11_barcoding/5_folds",
  output_dir = "11_barcoding/7_classifier",
  method = "nn",
  alignment = "add",
  workers = n_workers,
  notify = notify_email
)

# The classifier of the real use: IdTaxa needs no alignment, so it is the one that can answer a
# query a user brings. Its tables come from the cluster run (see idtaxa_tables in the SETUP block).
if (!is.null(idtaxa_tables)) {
  file.copy(list.files(idtaxa_tables, pattern = "^TABLE_barcoding_(predictions_(species|genus)|timing)_idtaxa\\.csv$",
                       full.names = TRUE),
            "11_barcoding/7_classifier", overwrite = FALSE)
}

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
# CN2 with IdTaxa comes from the cluster run, with the tables of step 7
if (!is.null(idtaxa_tables)) {
  file.copy(file.path(idtaxa_tables, "TABLE_barcoding_cn2_queries_idtaxa.csv"), "11_barcoding/8_controls",
            overwrite = FALSE)
}

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
# IdTaxa: the curve runs over its confidence. The operating threshold is 60, DECIPHER's default
# (decision K7 of 26-09), and step 10 reads that row; the quantile 0.01 of the genus confidence of
# scheme G is reported next to it, not used
if (!is.null(idtaxa_tables)) {
  sweep_barcoding_threshold(
    classifier_dir = "11_barcoding/7_classifier",
    controls_dir = "11_barcoding/8_controls",
    output_dir = "11_barcoding/9_threshold",
    method = "idtaxa",
    quantile_type = 1L
  )
}

# -------------------------------------------------------------
# Step 10: Metrics per class (6D): per locus and scheme, precision and recall per class and the shares
# of correct, honest genus, wrong and unassigned answers, never summed into one error; for the nearest
# neighbour at its operating threshold and, with the cluster tables, for IdTaxa at 60.
# -------------------------------------------------------------
summarise_barcoding_metrics(
  classifier_dir = "11_barcoding/7_classifier",
  threshold_dir = "11_barcoding/9_threshold",
  library_dir = "11_barcoding/4_library",
  output_dir = "11_barcoding/11_metrics",
  methods = if (is.null(idtaxa_tables)) "nn_add" else c("nn_add", "idtaxa")
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
# Metrics per locus and scheme
print(utils::read.csv("11_barcoding/11_metrics/TABLE_barcoding_metrics_summary_nn_add.csv"))
