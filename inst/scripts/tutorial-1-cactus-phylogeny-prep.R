# -------------------------------------------------------------
# PhyloCactus: Tutorial 1 - Phylogenetics Pipeline: Data Assembly & Preparation
# -------------------------------------------------------------
# This script covers Stages 1 to 6 of the phylogenetic pipeline:
# Sequence Mining, Alignment, Screening, and Concatenation.
#
# Keep this script in sync with vignettes/tutorial-1-cactus-phylogeny-prep.Rmd:
# every argument below is stated explicitly, including the ones that match the
# package defaults, so that a run can be reproduced from the script alone.
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
  # The exclusion lists are expert curation by the team supporting PhyloCactus: each
  # accession was inspected and removed on grounds not recoverable from GenBank metadata.
  # Stated explicitly so the curation is visible rather than inherited from a default.
  # Set to FALSE to reproduce the uncurated cluster set and measure what curation removes.
  apply_manual_exclusions = TRUE,
  min_species = 50,
  force_download = FALSE
)

# Mine outgroup taxonomic database
outgroup_assembly <- assemble_outgroup_phylotar(
  wd_path = "0_phylotaR_raw_Outgroup",
  target_genes_file = system.file("extdata", "target_genes.txt", package = "PhyloCactus"),
  genes_map_file = system.file("extdata", "genes_map.csv", package = "PhyloCactus"),
  manual_exclusions_file = system.file("extdata", "manual_exclusions_outgroup.csv", package = "PhyloCactus"),
  apply_manual_exclusions = TRUE,
  # Anacampserotaceae: Talinopsis (107598), Grahamia (107617), Anacampseros (107583).
  # Portulacaceae: Portulaca (3582).
  # Talinaceae: Talinum (107600), Talinella (108056).
  #
  # Talinopsis (107598) is Anacampserotaceae, not Talinaceae, despite the name and the adjacent id.
  #
  # Talinaceae was added on 2026-09-04, reversing the rejection of 2026-09-02. The reason for the
  # reversal is not the one considered then. With only Cactaceae, Anacampserotaceae and
  # Portulacaceae in the matrix, and the tree rooted on Portulacaceae, the grouping of Cactaceae
  # with Anacampserotaceae follows from the sampling and the rooting: no fourth lineage is present
  # that would let the data prefer any other resolution. That grouping is the node addressed by the
  # Cactaceae_Anacampserotaceae_stem calibration, so the calibration rests on a relationship this
  # matrix cannot test.
  #
  # The relationship is not settled in the literature either. Ramirez-Barahona et al. (2020) recover
  # Cactaceae with Anacampserotaceae; Zuntini et al. (2024) and Kew Tree of Life release 4.0 recover
  # Anacampserotaceae with Portulacaceae, but with a quartet support of 0.47 and 0.33 on an
  # alternative, so neither resolution is established. Adding Talinaceae supplies the fourth lineage
  # and lets the likelihood decide.
  #
  # Amphipetalum (1835425) is deliberately excluded: GenBank held no nucleotide records for it when
  # checked on 2026-09-04, so requesting it only produces an empty download.
  outgroups = c("107598", "107617", "107583", "3582", "107600", "108056"),
  force_download = FALSE
)

# -------------------------------------------------------------
# Stage 2: Enforcing Positional Homology and Mitigating Systematic Errors
# -------------------------------------------------------------
# Masking policy. DECIPHER::MaskAlignment defines its column set from the sequences
# it is given. The ingroup has dozens to hundreds of accessions per locus, so masking
# behaves reliably. The outgroup has only two to seven highly divergent accessions per
# locus, and masking it on its own produces a column set the ingroup does not share:
# rbcL collapses from 1629 to 67 bp, psbA-trnH from 377 to 1 bp. Those fragments then
# fail the occupancy filter of Stage 5 and disappear, stripping the outgroup from
# precisely the loci needed to root the tree. Masking of the outgroup is therefore
# deferred to Stage 5, where both groups are masked together against a single column set.
#
# Note: with mask_alignment_regions = FALSE the thresholds min_non_gap_fraction,
# max_missing_fraction and min_masked_alignment_length are not applied at this stage
# (they are relative to a post-masking column set that does not exist here). Occupancy
# filtering happens once, in run_joint_realignment() below.
# Ambiguity policy: preserve_iupac = TRUE keeps the IUPAC codes (R, Y, S, W, K, M,
# B, D, H, V) instead of collapsing them to N. RAxML-NG and ModelTest-NG treat them
# as a partial constraint on the state, so an R site means "A or G" while an N means
# "anything": collapsing them discards information for no gain. It is stated
# explicitly here, rather than left to the default, because it is a deliberate
# analytical choice. Set it to FALSE only if an alignment free of ambiguity is
# required for a specific reason.
#
# How much ambiguity there actually is per locus is reported in the manifest columns
# mean_fraction_ambiguous_* and n_sites_ambiguous_*. The raw_input figures are measured
# before any cleaning, so they show what GenBank actually delivered.
#
# Strand policy: fix_strand = TRUE puts every sequence of a marker on the same strand
# before MAFFT sees it. GenBank stores each record on whichever strand the submitter
# deposited, and MAFFT compares only the orientation it is given: a reverse-complemented
# accession is aligned anyway and contributes columns with no positional homology.
# Nothing downstream catches it. The row looks like a very divergent sequence, or fails
# the occupancy filter of Stage 5 and disappears without explanation.
#
# Found 2026-09-01 in the rbcL and matK accessions of Portulaca oleracea and P. pilosa,
# at 0.51 and 0.40 observed divergence from Cactaceae where P. grandiflora, same genus,
# sits at 0.030 and 0.066. P. oleracea was the terminal the acceptance check at the end
# of this script nominated for rooting.
#
# The ingroup was NOT clean, and the first count said it was: measured on 4_Cleaned,
# which is post-filtering, it read 0 of 3207. Measured on the raw input it is 163 --
# 77 in psbA-trnH and 86 in rpL16 -- which the Stage 5 occupancy filter had been
# discarding without saying why, because a reversed sequence aligns to nothing and
# ends with almost no occupancy. Correcting the strand recovered them: psbA_trnH went
# from 296 to 368 ingroup sequences and rpL16 from 600 to 677.
#
# Every sequence is logged in both directions in tables/LOG_STRAND_<marker>.csv. The
# ingroup and outgroup are aligned in separate runs, so neither can see that the two
# pools disagree with each other; homology_check in Stage 4 reports that case.
ingroup_alignment_manifest <- run_alignment_pipeline(
  input_folder = "1_phylotaR_out_Ingroup",
  output_dir = "2_MAFFT_Cactaceae",
  mask_alignment_regions = TRUE,
  min_non_gap_fraction = 0.30,
  max_missing_fraction = 0.30,
  min_masked_alignment_length = 100L,   # floor against masking-degraded alignments
  preserve_iupac = TRUE,                # keep ambiguity codes; see note above
  fix_strand = TRUE                     # one strand per marker before MAFFT; see note below
)

outgroup_alignment_manifest <- run_alignment_pipeline(
  input_folder = "1_phylotaR_out_Outgroup",
  output_dir = "2_MAFFT_Outgroup",
  mask_alignment_regions = FALSE,
  preserve_iupac = TRUE,
  fix_strand = TRUE
)

# Ambiguity actually present per locus, before and after processing. This is the table
# to read before deciding whether any individual marker warrants preserve_iupac = FALSE.
print(ingroup_alignment_manifest[, c("marker",
                                     "n_sites_ambiguous_raw_input",
                                     "mean_fraction_ambiguous_raw_input",
                                     "mean_fraction_ambiguous_final_filtered")])

# Quality check: confirm the masking policy that was actually applied. The outgroup rows
# must read masking_applied = FALSE and status = "OK_NO_MASK". If they read TRUE / "OK",
# the call above did not take effect and everything downstream is the old dataset.
print(outgroup_alignment_manifest[, c("marker", "masking_applied", "status")])
stopifnot(all(!outgroup_alignment_manifest$masking_applied))

# -------------------------------------------------------------
# Stage 3: Evaluating Marker Phylogenetic Utility (Saturation and Anomalies)
# -------------------------------------------------------------
screening_summary <- run_marker_screening(
  fasta_folder = "2_MAFFT_Cactaceae/alignments",
  out_base = "3_Saturation",
  min_cols_to_evaluate = 50,
  min_aln_len_to_retain = 200,
  min_nseq_to_retain = 100,
  # Reporting only, and only for the loci that have an outgroup counterpart. Which outgroup
  # markers have none is Stage 4's business, not this module's.
  # The summary table gains n_outgroup and n_total so a rejected locus can be
  # seen to carry outgroup data. No retention decision here depends on the outgroup.
  outgroup_folder = "2_MAFFT_Outgroup/alignments",
  saturation_keep_cutoff = 0.5
)

# -------------------------------------------------------------
# Stage 4: Standardizing Nomenclatural Frameworks and Decoupling Lineages
# -------------------------------------------------------------
integrate_and_clean_markers(
  ingroup_dir = "3_Saturation/filtered_markers",
  outgroup_dir = "2_MAFFT_Outgroup/alignments",
  output_dir = "4_Cleaned",
  accepted_list_file = system.file("extdata", "CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx", package = "PhyloCactus"),
  metadata_in_file = "1_phylotaR_out_Ingroup/TABLE_ACCESSION_OCCUPANCY_INGROUP_CLEAN.csv",
  metadata_out_file = "1_phylotaR_out_Outgroup/TABLE_ACCESSION_OCCUPANCY_OUTGROUP_CLEAN.csv",
  # The two phylotaR runs cluster independently, so the same region can carry different names.
  # trnL (outgroup, the trnL intron) shares 63% of its 20-mers with trnL-trnF (ingroup), whose
  # amplicon contains that intron. Verified before aliasing: outgroup ndhF shares none with
  # ndhF-rpl32 despite the resemblance, and is deliberately not aliased.
  marker_aliases = c(trnL = "trnL_trnF"),
  # readmit_markers overrides the Module 3 verdict for a locus kept on outgroup grounds rather
  # than ingroup resolution. Left empty deliberately.
  #
  # trnT-psbD was readmitted here on 2026-09-01 and removed the same day. Its counts were the most
  # balanced in the dataset, 49 ingroup against 49 outgroup at 0.998 median coverage, and the
  # sequences are not the same region: 0.039 of outgroup 20-mers occur in the ingroup, against
  # 0.28 to 0.60 for every genuine counterpart, and observed divergence against Cactaceae of 0.440
  # where matK gives 0.077 and rbcL 0.033. Ingroup median 598 bp, outgroup 1347 bp. It was placing
  # 521 non-homologous columns on the terminals that define the root, which is the branch the
  # dating depends on. homology_check below now measures this for every marker.
  readmit_markers = NULL,
  readmit_dir = NULL,
  # Alignment-free homology check, ingroup against outgroup, per marker. Writes
  # 4_Cleaned/tables/TABLE_marker_homology_check.csv and warns below both thresholds.
  homology_check = TRUE
)

# -------------------------------------------------------------
# Stage 5: Enforcing Absolute Positional Homology Across Aggregated Lineages
# -------------------------------------------------------------
# This is the single point at which the outgroup is masked, and the step that can drop
# individual terminals from a locus without dropping the locus itself. Removed terminals
# are listed in 5_MAFFT_Cleaned/LOG_SEQ_FILTER_<marker>.csv with Retained = FALSE.
#
# protect_pattern retains matching terminals regardless of occupancy. protect_markers
# restricts that exemption to named loci ("trnL_trnF"), safeguarding the short Portulaca
# amplicons necessary for root estimation without globally inflating missing data.
run_joint_realignment(
  input_dir = "4_Cleaned/cleaned_markers_joint",
  output_fasta_dir = "5_MAFFT_Cleaned",
  output_aln_dir = "5_MAFFT_Cleaned/aligned_markers",
  min_non_gap_fraction = 0.30,
  max_missing_fraction = 0.30,
  rooting_pattern = "^(Anacampseros|Grahamia|Talinopsis|Portulaca|Talinum|Talinella)_",
  protect_pattern = "^(Anacampseros|Grahamia|Talinopsis|Portulaca|Talinum|Talinella)_",
  protect_markers = "trnL_trnF"
)

# -------------------------------------------------------------
# Stage 6: Synthesizing the Unified Multilocus Supermatrix
# -------------------------------------------------------------
# exclude_markers names the loci left out of THIS supermatrix. The alignments are neither
# modified nor removed, so an excluded locus stays available to every other analysis; the
# exclusions are echoed to the run log and written to logs_and_qc/TABLE_markers_excluded.csv,
# and a name matching no alignment warns instead of passing silently.
#
# phyC is excluded because it carries 167 ingroup terminals and 15 outgroup terminals, all 15 in
# Anacampserotaceae and none in Portulacaceae: in a matrix built to date the divergence between
# those families it supplies columns only one side of the root can occupy. That is a property of
# this matrix, not of the locus, and phyC remains the more discriminating nuclear marker for
# barcoding. Set exclude_markers = NULL to build the full matrix.
run_concatenation_pipeline(
  input_dir = "5_MAFFT_Cleaned/aligned_markers",
  output_dir = "6_Concatenated",
  outgroup_pattern = "^(Anacampseros|Grahamia|Talinopsis|Portulaca|Talinum|Talinella)_",
  exclude_markers = "phyC"
)

# -------------------------------------------------------------
# Acceptance check
# -------------------------------------------------------------
# The number of partitions retaining outgroup is the headline result of this pipeline:
# it is how many of the retained loci can contribute to rooting the tree. Read it from the
# run log rather than assuming it.
cat(readLines("6_Concatenated/logs_and_qc/LOG_concatenation_run.txt", n = 8), sep = "\n")

# Outgroup occupancy per terminal, used in Tutorial 2 to choose the rooting terminal.
final_summary <- utils::read.csv("6_Concatenated/final_tables/TABLE_final_species_alignment_summary.csv")
outgroup_rows <- final_summary[final_summary$species_class == "outgroup", c("species", "retained_markers", "pct_markers")]
print(outgroup_rows[order(-outgroup_rows$retained_markers), ])
