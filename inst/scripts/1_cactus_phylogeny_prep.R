# -------------------------------------------------------------
# PhyloCactus: Tutorial 1 - Phylogenetics Pipeline: Data Assembly & Preparation
# -------------------------------------------------------------
# This script covers Stages 1 to 6 of the phylogenetic pipeline:
# Sequence Mining, Alignment, Screening, and Concatenation.
# -------------------------------------------------------------
library(PhyloCactus)

# -------------------------------------------------------------
# Setup: Creating a Clean Workspace
# -------------------------------------------------------------
tutorial_dir <- "~/Desktop/PhyloCactus_Tutorial"
dir.create(tutorial_dir, showWarnings = FALSE)
setwd(tutorial_dir)

dir.create("0_phylotaR_raw_Ingroup", showWarnings = FALSE)
dir.create("0_phylotaR_raw_Outgroup", showWarnings = FALSE)

# -------------------------------------------------------------
# Stage 1: Mine GenBank Clusters and Fetch Metadata
# -------------------------------------------------------------
# Mine ingroup taxonomic database
ingroup_assembly <- assemble_ingroup_phylotar(
  wd_path = "0_phylotaR_raw_Ingroup",
  target_genes_file = system.file("extdata", "target_genes.txt", package = "PhyloCactus"),
  genes_map_file = system.file("extdata", "genes_map.csv", package = "PhyloCactus"),
  manual_exclusions_file = system.file("extdata", "manual_exclusions_ingroup.csv", package = "PhyloCactus"),
  min_species = 50,
  force_download = FALSE
)

# Mine outgroup taxonomic database
outgroup_assembly <- assemble_outgroup_phylotar(
  wd_path = "0_phylotaR_raw_Outgroup",
  target_genes_file = system.file("extdata", "target_genes.txt", package = "PhyloCactus"),
  genes_map_file = system.file("extdata", "genes_map.csv", package = "PhyloCactus"),
  manual_exclusions_file = system.file("extdata", "manual_exclusions_outgroup.csv", package = "PhyloCactus"),
  outgroups = c("107598", "107617", "107583", "3582"),
  force_download = FALSE
)

# -------------------------------------------------------------
# Stage 2: Enforcing Positional Homology and Mitigating Systematic Errors
# -------------------------------------------------------------
ingroup_alignment_manifest <- run_alignment_pipeline(
  input_folder = "1_phylotaR_out_ingroup",
  output_dir = "2_MAFFT_Cactaceae",
  min_non_gap_fraction = 0.30,
  max_missing_fraction = 0.30
)

outgroup_alignment_manifest <- run_alignment_pipeline(
  input_folder = "1_phylotaR_out_outgroup",
  output_dir = "2_MAFFT_Outgroup",
  min_non_gap_fraction = 0.30,
  max_missing_fraction = 0.30
)

# -------------------------------------------------------------
# Stage 3: Evaluating Marker Phylogenetic Utility (Saturation and Anomalies)
# -------------------------------------------------------------
screening_summary <- run_marker_screening(
  fasta_folder = "2_MAFFT_Cactaceae/alignments",
  out_base = "3_Saturation",
  min_cols_to_evaluate = 50,
  min_aln_len_to_retain = 200,
  min_nseq_to_retain = 100,
  saturation_keep_cutoff = 0.5
)

# -------------------------------------------------------------
# Stage 4: Standardizing Nomenclatural Frameworks and Unified Matrices
# -------------------------------------------------------------
integrate_and_clean_markers(
  ingroup_dir = "3_Saturation/filtered_markers",
  outgroup_dir = "2_MAFFT_Outgroup/alignments",
  output_dir = "4_Cleaned",
  accepted_list_file = system.file("extdata", "CactaceaeFullList_accepted.csv", package = "PhyloCactus"),
  metadata_in_file = "1_phylotaR_out_ingroup/TABLE_ACCESSION_OCCUPANCY_INGROUP_CLEAN.csv",
  metadata_out_file = "1_phylotaR_out_outgroup/TABLE_ACCESSION_OCCUPANCY_OUTGROUP_CLEAN.csv"
)

# -------------------------------------------------------------
# Stage 5: Enforcing Absolute Positional Homology Across Aggregated Lineages
# -------------------------------------------------------------
run_joint_realignment(
  input_dir = "4_Cleaned/cleaned_markers",
  output_fasta_dir = "5_MAFFT_Cleaned",
  output_aln_dir = "5_MAFFT_Cleaned/aligned_markers"
)

# -------------------------------------------------------------
# Stage 6: Synthesizing the Unified Multilocus Supermatrix
# -------------------------------------------------------------
run_concatenation_pipeline(
  input_dir = "5_MAFFT_Cleaned/aligned_markers",
  output_dir = "6_Concatenated"
)
