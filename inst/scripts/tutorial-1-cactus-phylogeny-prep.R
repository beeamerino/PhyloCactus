# -------------------------------------------------------------
# PhyloCactus: Tutorial 1 - Phylogenetics Pipeline: Data Assembly & Preparation
# -------------------------------------------------------------
# This script covers Modules 1 to 6 of the phylogenetic pipeline:
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
# Module 1: Mine Orthologous Sequence Clusters and Retrieve Metadata
# -------------------------------------------------------------
# Mine ingroup taxonomic database
ingroup_assembly <- assemble_ingroup_phylotar(
  wd_path = "0_phylotaR_raw_Ingroup",
  target_genes_file = system.file("extdata", "target_genes.txt", package = "PhyloCactus"),
  genes_map_file = system.file("extdata", "genes_map.csv", package = "PhyloCactus"),
  manual_exclusions_file = system.file("extdata", "manual_exclusions_ingroup.csv", package = "PhyloCactus"),
  # The exclusion lists are expert curation by the team supporting PhyloCactus: each
  # accession was inspected and removed on grounds not recoverable from GenBank metadata.
  # Stated explicitly so that the curation is visible.
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
  # Talinaceae roots the tree outside the ACP clade. With only Cactaceae, Anacampserotaceae and
  # Portulacaceae in the matrix and the tree rooted on Portulacaceae, the grouping of Cactaceae with
  # Anacampserotaceae follows from the sampling and the rooting, and the data cannot prefer another
  # resolution. The relationship is not settled in the literature: Ramirez-Barahona et al. (2020)
  # recover Cactaceae with Anacampserotaceae; Zuntini et al. (2024) recover Anacampserotaceae with
  # Portulacaceae, with a quartet support of 0.47 against 0.33 for an alternative. Sampling
  # Talinaceae supplies the fourth lineage, and the likelihood decides among the resolutions.
  #
  # Amphipetalum (1835425) is not requested: GenBank holds no nucleotide records for it, so the
  # request would only produce an empty download.
  outgroups = c("107598", "107617", "107583", "3582", "107600", "108056"),
  force_download = FALSE
)

# -------------------------------------------------------------
# Module 2: Align Sequences and Mask Low Confidence Regions
# -------------------------------------------------------------
# Masking policy. DECIPHER::MaskAlignment defines its column set from the sequences
# it is given. The ingroup has dozens to hundreds of accessions per locus, so masking
# behaves reliably. The outgroup has only two to seven highly divergent accessions per
# locus, and masking it on its own produces a column set the ingroup does not share:
# rbcL collapses from 1629 to 67 bp, psbA-trnH from 377 to 1 bp. Those fragments then
# fail the occupancy filter of Module 5 and disappear, stripping the outgroup from
# precisely the loci needed to root the tree. Masking of the outgroup is therefore
# deferred to Module 5, where both groups are masked together against a single column set.
#
# Note: with mask_alignment_regions = FALSE the thresholds min_non_gap_fraction,
# max_missing_fraction and min_masked_alignment_length are not applied at this stage
# (they are relative to a post-masking column set that does not exist here). Occupancy
# filtering happens once, in run_joint_realignment() below.
# Ambiguity policy: preserve_iupac = TRUE keeps the IUPAC codes (R, Y, S, W, K, M,
# B, D, H, V) instead of collapsing them to N. RAxML-NG and ModelTest-NG treat them
# as a partial constraint on the state, so an R site means "A or G" while an N means
# "anything": collapsing them discards information. It is stated explicitly because
# it is an analytical choice. Set it to FALSE for analyses that require alignments
# without ambiguity codes.
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
# the occupancy filter of Module 5 and disappears without explanation.
#
# In the reference dataset, the rbcL and matK accessions of Portulaca oleracea and P. pilosa are
# deposited reversed, and so are 163 ingroup sequences, 77 of psbA-trnH and 86 of rpL16.
# Uncorrected, they align to nothing, lose occupancy and are removed by the occupancy filter of
# Module 5, so the curated output shows no reversed sequence. Correcting the strand retains them:
# psbA_trnH keeps 368 ingroup sequences in place of 296, and rpL16 677 in place of 600.
#
# Every sequence is logged in both directions in tables/LOG_STRAND_<marker>.csv. The
# ingroup and outgroup are aligned in separate runs, so neither can see that the two
# pools disagree with each other; homology_check in Module 4 reports that case.
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
# Module 3: Evaluate Phylogenetic Signal and Screen Ingroup Markers
# -------------------------------------------------------------
screening_summary <- run_marker_screening(
  fasta_folder = "2_MAFFT_Cactaceae/alignments",
  out_base = "3_Saturation",
  min_cols_to_evaluate = 50,
  min_aln_len_to_retain = 200,
  min_nseq_to_retain = 100,
  # Reporting only, and only for the loci that have an outgroup counterpart. Which outgroup
  # markers have none is Module 4's business, not this module's.
  # The summary table gains n_outgroup and n_total so a rejected locus can be
  # seen to carry outgroup data. No retention decision here depends on the outgroup.
  outgroup_folder = "2_MAFFT_Outgroup/alignments",
  saturation_keep_cutoff = 0.5
)

# -------------------------------------------------------------
# Module 4: Integrate and Decouple Ingroup and Outgroup Markers
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
  # ndhF-rpl32 despite the resemblance, and is not aliased.
  marker_aliases = c(trnL = "trnL_trnF"),
  # readmit_markers overrides the Module 3 verdict for a locus kept on outgroup grounds rather
  # than ingroup resolution. Left empty.
  #
  # trnT-psbD is not readmitted. It has more outgroup sequences than any other marker (52), but
  # the two sets are not the same region: 0.040 of the outgroup 20-mers occur in the ingroup,
  # against 0.26 to 0.59 for trnL_trnF, matK, phyC and rbcL. homology_check below measures
  # this for every marker.
  readmit_markers = NULL,
  readmit_dir = NULL,
  # Alignment-free homology check, ingroup against outgroup, per marker. Writes
  # 4_Cleaned/tables/TABLE_marker_homology_check.csv and warns below both thresholds.
  homology_check = TRUE
)

# -------------------------------------------------------------
# Module 5: Joint Realignment of Curated Markers
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
# Module 6: Construct the Multilocus Supermatrix
# -------------------------------------------------------------
# exclude_markers names the loci left out of THIS supermatrix. The alignments are neither
# modified nor removed, so an excluded locus stays available to every other analysis; the
# exclusions are echoed to the run log and written to logs_and_qc/TABLE_markers_excluded.csv,
# and a name matching no alignment warns instead of passing silently.
#
# phyC is excluded because it carries no Portulacaceae terminal: of its 16 outgroup terminals, 15
# are Anacampserotaceae and one is Talinaceae, against 167 ingroup terminals. In a matrix built to
# date the ACP clade, a locus sampled for Anacampserotaceae and not for Portulacaceae supplies
# characters to only one side of the divergence between the two families. The exclusion is a
# property of this matrix; phyC remains the more discriminating nuclear marker for barcoding.
# Set exclude_markers = NULL to build the full matrix.
#
# outgroup_pattern is a named vector, not a single expression, so the coverage report resolves the
# outgroup by family and writes n_<family>, n_<family>_covered and median_cov_<family>. A combined
# count can look adequate while the family carrying the root is absent from the locus entirely.
run_concatenation_pipeline(
  input_dir = "5_MAFFT_Cleaned/aligned_markers",
  output_dir = "6_Concatenated",
  outgroup_pattern = c(
    Anacampserotaceae = "^(Anacampseros|Grahamia|Talinopsis)_",
    Portulacaceae     = "^Portulaca_",
    Talinaceae        = "^(Talinum|Talinella)_"
  ),
  exclude_markers = "phyC"
)

# -------------------------------------------------------------
# Acceptance check
# -------------------------------------------------------------
# The number of partitions retaining outgroup is the headline result of this pipeline:
# it is how many of the retained loci can contribute to rooting the tree. Read it from the
# run log.
cat(readLines("6_Concatenated/logs_and_qc/LOG_concatenation_run.txt", n = 8), sep = "\n")

# Outgroup occupancy per terminal, used in Tutorial 2 to choose the rooting terminal.
final_summary <- utils::read.csv("6_Concatenated/final_tables/TABLE_final_species_alignment_summary.csv")
outgroup_rows <- final_summary[final_summary$species_class == "outgroup", c("species", "retained_markers", "pct_markers")]
print(outgroup_rows[order(-outgroup_rows$retained_markers), ])
