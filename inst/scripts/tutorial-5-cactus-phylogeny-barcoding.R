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
# scheme G is reported next to it, not used. The sweep also writes
# TABLE_barcoding_species_threshold_idtaxa.csv: for each locus, the quantile 0.95 of the species
# confidence of scheme G (species absent from the library), never under 60. The identification of
# a user query names a species only above it (rule R-c, Phase 11); the validation stays at 60
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
# Step 11: Identification of a sample. Each locus of the library is cut from the query, classified
# with IdTaxa against the library and given one of three states: a species (state 1), a genus
# (state 2) or no answer (state 3). A species is named only above the species threshold of its locus,
# written by step 9 when the IdTaxa tables are present (rule R-c); without them every locus uses 60.
# The example is the plastome of Mammillaria zephyranthoides (GenBank MN517611.1) distributed with
# the package. The table and a self-contained HTML report are written in 11_barcoding/10_identify/.
# -------------------------------------------------------------
identify_barcoding_query(
  system.file("extdata", "barcoding_example_plastome_MN517611.1.gb", package = "PhyloCactus"),
  library_dir = "11_barcoding/4_library",
  metrics_dir = "11_barcoding/11_metrics",
  threshold_dir = "11_barcoding/9_threshold",
  output_dir = "11_barcoding/10_identify",
  run_name = "MN517611.1"
)
print(utils::read.csv("11_barcoding/10_identify/TABLE_barcoding_identify_MN517611.1.csv")[
  , c("locus", "state", "predicted_species", "predicted_genus", "genus_confidence", "species_confidence",
      "species_threshold", "reason")])

# Your own samples take the same call or one of the three below. They need your files, and the reads
# need fastp and GetOrganelle on the PATH (minimap2 and samtools for an assembly), so they are shown
# here and not run.
# Several samples at once, from a sample sheet (sample_id, input_type, file_1 and, for reads, file_2):
# identify_barcoding_samples("my_samples.csv")
# An assembly (FASTA, plain or gzipped), from which the library loci are mapped and cut:
# identify_barcoding_assembly("my_assembly.fasta", run_name = "my_assembly")
# Paired Illumina reads (genome skim or whole-genome sequencing): plastome and nrDNA are assembled,
# then identified locus by locus:
# identify_barcoding_reads("my_run_1.fastq.gz", "my_run_2.fastq.gz", run_name = "my_run")

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

# -------------------------------------------------------------
# Figures of the library and its validation (TB1 to TB7 and TB9), drawn from the files of steps 1
# to 10 into 11_barcoding/figures/, each as PDF and PNG. The schema of the branch (TB8) is drawn once
# and is not produced here. The figures of the folds read the IdTaxa predictions when the cluster
# tables were copied into step 7, and the nearest neighbour through MAFFT --add otherwise; the
# states of the nearest neighbour are those of the tie structure, before the threshold of
# remoteness of step 9.
# -------------------------------------------------------------
library(ggplot2)
fig_dir <- "11_barcoding/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
save_fig <- function(p, name, width, height) {
  for (ext in c("pdf", "png")) {
    ggsave(file.path(fig_dir, paste0(name, ".", ext)), p, width = width, height = height, dpi = 300, limitsize = FALSE)
  }
}
# The sequences of the library: one per accession, identical sequences of a species collapsed
lib_records <- utils::read.csv("11_barcoding/4_library/TABLE_barcoding_collapsed_identical.csv")
lib_loci <- sort(sub("^LIB_(.*)\\.fasta$", "\\1", list.files("11_barcoding/4_library", pattern = "^LIB_.*\\.fasta$")),
                 method = "radix")
lib_records <- lib_records[!as.logical(lib_records$collapsed) & lib_records$locus %in% lib_loci, ]
lib_records$genus <- sub("_.*$", "", lib_records$species)

# Accepted species of the checklist (Korotkova et al. 2021) and the genera with a native species in
# Chile (sheet Facts_Cactaceae)
checklist_file <- system.file("extdata", "CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx", package = "PhyloCactus")
accepted <- as.data.frame(readxl::read_excel(checklist_file, sheet = "Cactaceae"))
accepted <- unique(gsub(" ", "_", accepted$pureName[accepted$RANK %in% "Species"]))
accepted_genus <- table(sub("_.*$", "", accepted))
facts <- as.data.frame(readxl::read_excel(checklist_file, sheet = "Facts_Cactaceae"))
chilean_genera <- unique(sub(" .*$", "", facts$fullName[facts$RANK %in% "Genus" & facts$area %in% "Chile" &
                                                           facts$status %in% "native"]))

# The predictions of the folds, one outcome per query
fold_method <- if (file.exists("11_barcoding/7_classifier/TABLE_barcoding_predictions_species_idtaxa.csv")) "idtaxa" else "nn_add"
folds <- do.call(rbind, lapply(c("species", "genus"), function(scheme) {
  d <- utils::read.csv(sprintf("11_barcoding/7_classifier/TABLE_barcoding_predictions_%s_%s.csv", scheme, fold_method))
  d[, c("locus", "scheme", "sid", "true_species", "true_genus", "state", "predicted_species", "predicted_genus")]
}))
same <- function(a, b) !is.na(a) & !is.na(b) & a == b
folds$outcome <- ifelse(folds$state %in% 1L & same(folds$predicted_species, folds$true_species), "species correct",
                 ifelse(folds$state %in% 1L, "wrong species",
                 ifelse(folds$state %in% 2L & same(folds$predicted_genus, folds$true_genus), "genus correct",
                 ifelse(folds$state %in% 2L, "wrong genus", "abstention"))))
folds$genus_right <- folds$state %in% c(1L, 2L) & same(folds$predicted_genus, folds$true_genus)
folds$scheme <- factor(folds$scheme, levels = c("species", "genus"))
outcome_levels <- c("species correct", "genus correct", "abstention", "wrong species", "wrong genus")
outcome_colours <- c("species correct" = "#1b7837", "genus correct" = "#a6dba0", "abstention" = "#bdbdbd",
                     "wrong species" = "#f4a582", "wrong genus" = "#b2182b")

# TB1. Alignment overview of each locus: one row per sequence, ordered by genus, coloured by base,
# gaps blank, with the occupancy of each column above
for (l in lib_loci) {
  aln <- Biostrings::readDNAStringSet(file.path("11_barcoding/4_library", paste0("LIB_", l, ".fasta")))
  m <- as.matrix(aln)
  ord <- order(sub("_.*$", "", names(aln)), names(aln), method = "radix")
  m <- m[ord, , drop = FALSE]
  cells <- data.frame(row = rep(seq_len(nrow(m)), ncol(m)), column = rep(seq_len(ncol(m)), each = nrow(m)),
                      base = factor(toupper(as.vector(m)), levels = c("A", "C", "G", "T")))
  cells <- cells[!is.na(cells$base), ]
  occupancy <- data.frame(column = seq_len(ncol(m)), occupancy = colMeans(m != "-" & m != "N"))
  p_aln <- ggplot(cells, aes(column, row, fill = base)) + geom_raster() +
    scale_fill_manual(values = c(A = "#4daf4a", C = "#377eb8", G = "#ff7f00", T = "#e41a1c")) +
    scale_y_reverse(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
    labs(x = "Alignment column", y = "Sequences, ordered by genus", fill = "Base") + theme_minimal(base_size = 9)
  p_occ <- ggplot(occupancy, aes(column, occupancy)) + geom_area(fill = "grey60") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = "Occupancy", title = paste0(l, ": ", nrow(m), " sequences, ", ncol(m), " columns")) +
    theme_minimal(base_size = 9)
  save_fig(patchwork::wrap_plots(p_occ, p_aln, ncol = 1, heights = c(1, 5)), paste0("Figure_TB1_alignment_", l), 10, 8)
}

# TB2. Sequences of each species in each locus, species ordered by genus
species_locus <- as.data.frame(table(species = lib_records$species, locus = lib_records$locus))
species_locus <- species_locus[species_locus$Freq > 0, ]
species_locus$species <- factor(species_locus$species, levels = rev(sort(unique(lib_records$species), method = "radix")))
p_tb2 <- ggplot(species_locus, aes(locus, species, fill = Freq)) + geom_tile() +
  scale_fill_viridis_c(trans = "log10", name = "Sequences") +
  labs(x = NULL, y = paste0(length(unique(lib_records$species)), " species, ordered by genus")) +
  theme_minimal(base_size = 9) + theme(axis.text.y = element_blank(), panel.grid = element_blank(),
                                       axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_tb2, "Figure_TB2_species_by_locus", 7, 12)

# TB3. Accessions of each locus from assembly to library
funnel <- utils::read.csv("11_barcoding/4_library/TABLE_barcoding_funnel.csv")
funnel <- funnel[funnel$locus %in% lib_loci, ]
funnel$after_exclusion <- funnel$accessions_curated - funnel$accessions_excluded_nonhomologous
funnel_long <- data.frame(
  locus = rep(funnel$locus, 4),
  stage = factor(rep(c("assembled", "curated", "homologous", "library"), each = nrow(funnel)),
                 levels = c("assembled", "curated", "homologous", "library")),
  accessions = c(funnel$accessions_assembled, funnel$accessions_curated, funnel$after_exclusion,
                 funnel$accessions_after_collapse))
p_tb3 <- ggplot(funnel_long, aes(stage, accessions, group = locus)) + geom_line(colour = "grey40") + geom_point() +
  facet_wrap(~ locus, scales = "free_y") + labs(x = NULL, y = "Accessions") +
  theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_tb3, "Figure_TB3_funnel", 10, 8)

# TB4. Length of the sequences of each locus before and after the cut of records over 2000 bases
# (L2), and sequences per species in the library
cut_records <- utils::read.csv("11_barcoding/1_assembly/TABLE_barcoding_long_records_cut.csv")
lengths <- do.call(rbind, lapply(lib_loci, function(l) {
  s <- Biostrings::readDNAStringSet(file.path("11_barcoding/4_library", paste0("UNMASKED_", l, ".fasta")))
  sid <- sub("^.*\\|", "", names(s))
  after <- Biostrings::width(s)
  before <- after
  hit <- match(paste(l, sid), paste(cut_records$locus, cut_records$sid))
  before[!is.na(hit)] <- cut_records$length[hit[!is.na(hit)]]
  rbind(data.frame(locus = l, stage = "before the cut", length = before),
        data.frame(locus = l, stage = "after the cut", length = after))
}))
lengths$stage <- factor(lengths$stage, levels = c("before the cut", "after the cut"))
p_len <- ggplot(lengths, aes(locus, length, fill = stage)) + geom_boxplot(outlier.size = 0.4) +
  scale_y_log10() + labs(x = NULL, y = "Sequence length (bases)", fill = NULL) +
  theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
replicates <- as.data.frame(table(locus = lib_records$locus, species = lib_records$species))
replicates <- replicates[replicates$Freq > 0, ]
replicates$class <- factor(cut(replicates$Freq, breaks = c(0, 1, 2, 5, Inf), labels = c("1", "2", "3 to 5", "6 or more")),
                           levels = c("6 or more", "3 to 5", "2", "1"))
p_rep <- ggplot(replicates, aes(locus, fill = class)) + geom_bar() +
  scale_fill_viridis_d(direction = -1) +
  labs(x = NULL, y = "Species", fill = "Sequences per species") +
  theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(patchwork::wrap_plots(p_len, p_rep, ncol = 1), "Figure_TB4_lengths_and_replicates", 9, 8)

# TB5. Answers of the folds per genus against its sampling: species of the genus in the library and
# the fraction of its accepted species they cover. For a sample of an accepted species of the genus,
# the expected rate of a correct genus is f x (right genus in scheme E) + (1 - f) x (right genus in
# scheme G), with f the sampled fraction and the right genus counted in states 1 and 2
lib_genus <- tapply(lib_records$species, lib_records$genus, function(x) length(unique(x)))
genus_rates <- do.call(rbind, lapply(split(folds, folds$true_genus), function(d) {
  e <- d[d$scheme == "species", ]
  g <- d[d$scheme == "genus", ]
  data.frame(genus = d$true_genus[1],
             species_correct_E = if (nrow(e)) mean(e$outcome == "species correct") else NA,
             genus_right_E = if (nrow(e)) mean(e$genus_right) else NA,
             genus_right_G = if (nrow(g)) mean(g$genus_right) else NA,
             queries = nrow(d))
}))
genus_rates$species_in_library <- as.integer(lib_genus[genus_rates$genus])
genus_rates$accepted <- as.integer(accepted_genus[genus_rates$genus])
genus_rates$fraction <- pmin(1, genus_rates$species_in_library / genus_rates$accepted)
genus_rates$expected_genus <- with(genus_rates, fraction * genus_right_E + (1 - fraction) * genus_right_G)
p_tb5a <- ggplot(genus_rates[!is.na(genus_rates$species_correct_E), ],
                 aes(species_in_library, species_correct_E, size = queries)) +
  geom_point(alpha = 0.5) + scale_x_log10() + scale_size_area(max_size = 5) +
  labs(x = "Species of the genus in the library", y = "Species correct, scheme E", size = "Queries") +
  theme_minimal(base_size = 9)
p_tb5b <- ggplot(genus_rates[!is.na(genus_rates$expected_genus) & !is.na(genus_rates$fraction), ],
                 aes(fraction, expected_genus, size = queries)) +
  geom_point(alpha = 0.5) + scale_size_area(max_size = 5) +
  labs(x = "Fraction of the accepted species of the genus in the library",
       y = "Expected rate of a correct genus", size = "Queries") +
  theme_minimal(base_size = 9)
save_fig(patchwork::wrap_plots(p_tb5a, p_tb5b, ncol = 2) +
           patchwork::plot_annotation(caption = paste0("Folds classified by ", fold_method, "; one point per genus.")),
         "Figure_TB5_accuracy_against_sampling", 10, 4.5)

# TB6. Answers of the folds against the missing data of the query: the share of the alignment
# columns its sequence covers, in bins
coverage <- do.call(rbind, lapply(lib_loci, function(l) {
  aln <- Biostrings::readDNAStringSet(file.path("11_barcoding/4_library", paste0("LIB_", l, ".fasta")))
  m <- as.matrix(aln)
  data.frame(locus = l, sid = sub("^.*\\|", "", names(aln)), coverage = rowMeans(m != "-" & m != "N"))
}))
folds_cov <- merge(folds, coverage, by = c("locus", "sid"))
folds_cov$coverage_bin <- cut(folds_cov$coverage, breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1), include.lowest = TRUE)
folds_cov$outcome <- factor(folds_cov$outcome, levels = outcome_levels)
p_tb6 <- ggplot(folds_cov, aes(coverage_bin, fill = outcome)) + geom_bar(position = "fill") +
  facet_wrap(~ scheme, labeller = labeller(scheme = c(species = "Scheme E (species present)", genus = "Scheme G (species absent)"))) +
  scale_fill_manual(values = outcome_colours) +
  labs(x = "Share of the alignment columns covered by the query", y = "Share of queries", fill = NULL,
       caption = paste0("Folds classified by ", fold_method, ".")) +
  theme_minimal(base_size = 9)
save_fig(p_tb6, "Figure_TB6_accuracy_against_missing_data", 9, 4.5)

# TB7. Taxonomic coverage against the checklist: species of each locus and of the library over the
# accepted species, and the sampled fraction of each genus
cover_locus <- data.frame(locus = c(lib_loci, "any locus"),
                          species = c(vapply(lib_loci, function(l) length(unique(lib_records$species[lib_records$locus == l])), 1L),
                                      length(unique(lib_records$species))))
cover_locus$fraction <- cover_locus$species / length(accepted)
cover_locus$locus <- factor(cover_locus$locus, levels = cover_locus$locus[order(cover_locus$species)])
p_tb7a <- ggplot(cover_locus, aes(fraction, locus)) + geom_col(fill = "grey40") +
  geom_text(aes(label = species), hjust = -0.2, size = 2.5) + scale_x_continuous(limits = c(0, 1)) +
  labs(x = paste0("Fraction of the ", length(accepted), " accepted species"), y = NULL) + theme_minimal(base_size = 9)
cover_genus <- data.frame(genus = names(accepted_genus), accepted = as.integer(accepted_genus))
cover_genus$in_library <- as.integer(ifelse(is.na(lib_genus[cover_genus$genus]), 0L, lib_genus[cover_genus$genus]))
cover_genus$fraction <- pmin(1, cover_genus$in_library / cover_genus$accepted)
p_tb7b <- ggplot(cover_genus, aes(accepted, fraction)) + geom_point(alpha = 0.5) + scale_x_log10() +
  labs(x = "Accepted species of the genus", y = "Fraction in the library") + theme_minimal(base_size = 9)
save_fig(patchwork::wrap_plots(p_tb7a, p_tb7b, ncol = 2), "Figure_TB7_taxonomic_coverage", 10, 4.5)

# TB9. The genera native to Chile read against two large, well sampled genera
folds$region <- ifelse(folds$true_genus %in% chilean_genera, "Chilean genera",
                       ifelse(folds$true_genus %in% c("Opuntia", "Mammillaria"), folds$true_genus, "other genera"))
folds$region <- factor(folds$region, levels = c("Chilean genera", "Opuntia", "Mammillaria", "other genera"))
folds$outcome <- factor(folds$outcome, levels = outcome_levels)
p_tb9 <- ggplot(folds, aes(region, fill = outcome)) + geom_bar(position = "fill") +
  facet_wrap(~ scheme, labeller = labeller(scheme = c(species = "Scheme E (species present)", genus = "Scheme G (species absent)"))) +
  scale_fill_manual(values = outcome_colours) +
  labs(x = NULL, y = "Share of queries", fill = NULL,
       caption = paste(strwrap(paste0("Folds classified by ", fold_method, "; Chilean genera: ",
                                       paste(sort(chilean_genera), collapse = ", "), "."), 120), collapse = "\n")) +
  theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_tb9, "Figure_TB9_chilean_genera", 9, 4.5)
