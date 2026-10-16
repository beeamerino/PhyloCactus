# Package index

## Sequence Assemblies & phylotaR Import

Setup local taxonomic databases, pull orthologous sequence clusters via
`phylotaR`, and retrieve GenBank metadata.

- [`assemble_ingroup_phylotar()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_ingroup_phylotar.md)
  : Assemble Ingroup Sequence Clusters via phylotaR Similarity Mining
- [`assemble_outgroup_phylotar()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_outgroup_phylotar.md)
  : Assemble Outgroup Sequence Clusters via phylotaR
- [`fetch_genbank_metadata()`](https://beeamerino.github.io/PhyloCactus/reference/fetch_genbank_metadata.md)
  : Fetch GenBank Sequence Metadata via NCBI Entrez Utilities

## Taxonomic Cleaning & Reconciliations

Verify nomenclatural spelling, audit duplicates, and reconcile species
names against accepted botanical checklists.

- [`clean_taxonomic_names()`](https://beeamerino.github.io/PhyloCactus/reference/clean_taxonomic_names.md)
  : Reconcile and Validate Taxonomic Nomenclature
- [`integrate_and_clean_markers()`](https://beeamerino.github.io/PhyloCactus/reference/integrate_and_clean_markers.md)
  : Final Marker Integration, Taxonomic Cleaning, and Ingroup-Outgroup
  Partitioning

## Multiple Alignment, DECIPHER Cleaners, & Masking

Infer positional homology using `MAFFT` and refine nucleotide alignments
with `DECIPHER` automated masking.

- [`run_mafft()`](https://beeamerino.github.io/PhyloCactus/reference/run_mafft.md)
  : Infer Positional Homology via MAFFT Alignment
- [`run_alignment_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_alignment_pipeline.md)
  : Execute Complete Alignment and Gap-Masking Pipeline
- [`run_joint_realignment()`](https://beeamerino.github.io/PhyloCactus/reference/run_joint_realignment.md)
  : Perform Joint Realignment Across Integrated Ingroup and Outgroup
  Sequences

## Locus Quality Control & Saturation Screening

Detect sequence length outliers and evaluate substitution saturation
regressions across orthologous markers.

- [`run_marker_screening()`](https://beeamerino.github.io/PhyloCactus/reference/run_marker_screening.md)
  : Screen Locus Alignments for Substitution Saturation and
  Informativeness

## Concatenation and Partitioning Maps

Concatenate locus alignments end-to-end into a multilocus supermatrix
and define partition coordinates.

- [`run_concatenation_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_concatenation_pipeline.md)
  : Concatenate Locus Alignments and Build Partition Coordinate Maps
- [`report_marker_group_coverage()`](https://beeamerino.github.io/PhyloCactus/reference/report_marker_group_coverage.md)
  : Report How Each Marker Covers the Ingroup and the Outgroup

## Phylogenetic Inference & Constrained Search

Evaluate nucleotide substitution models (`ModelTest-NG`) and infer
maximum-likelihood phylogenies (`RAxML-NG`) guided by constraint
scaffolds.

- [`preprocess_partitions()`](https://beeamerino.github.io/PhyloCactus/reference/preprocess_partitions.md)
  : Preprocess Partitions and Validate Alignment Syntax
- [`run_modeltest_ng()`](https://beeamerino.github.io/PhyloCactus/reference/run_modeltest_ng.md)
  : Evaluate Nucleotide Substitution Models via ModelTest-NG
- [`build_constraint_scaffold()`](https://beeamerino.github.io/PhyloCactus/reference/build_constraint_scaffold.md)
  : Synthesize Multifurcating Monophyly Constraint Scaffold
- [`generate_ml_search_script()`](https://beeamerino.github.io/PhyloCactus/reference/generate_ml_search_script.md)
  : Generate a SLURM Batch Script for the Constrained Maximum-Likelihood
  Search
- [`calculate_ml_tree()`](https://beeamerino.github.io/PhyloCactus/reference/calculate_ml_tree.md)
  : Infer Maximum-Likelihood Phylogeny under Constrained Search
- [`calculate_rf_distances()`](https://beeamerino.github.io/PhyloCactus/reference/calculate_rf_distances.md)
  : Compute Robinson-Foulds Distances Across Maximum-Likelihood Trees
- [`generate_bootstrap_script()`](https://beeamerino.github.io/PhyloCactus/reference/generate_bootstrap_script.md)
  : Generate HPC SLURM Batch Script for Parallel Bootstrapping
- [`run_local_bootstraps()`](https://beeamerino.github.io/PhyloCactus/reference/run_local_bootstraps.md)
  : Generate Non-Parametric Bootstrap Trees Locally
- [`collect_bootstraps()`](https://beeamerino.github.io/PhyloCactus/reference/collect_bootstraps.md)
  : Collect and Concatenate Parallel Bootstrap Tree Outputs
- [`check_bs_convergence()`](https://beeamerino.github.io/PhyloCactus/reference/check_bs_convergence.md)
  : Check Bootstrap Convergence Criterion in RAxML-NG
- [`generate_temporal_bootstrap_script()`](https://beeamerino.github.io/PhyloCactus/reference/generate_temporal_bootstrap_script.md)
  : Generate High-Performance Computing (HPC/SLURM) Batch Script for
  Temporal Bootstrapping
- [`calculate_temporal_bootstraps()`](https://beeamerino.github.io/PhyloCactus/reference/calculate_temporal_bootstraps.md)
  : Estimate Temporal Bootstrap Replicates Constrained to Best ML
  Topology
- [`map_branch_supports()`](https://beeamerino.github.io/PhyloCactus/reference/map_branch_supports.md)
  : Map Bootstrap Support Values onto Reference Phylogeny
- [`classify_constrained_nodes()`](https://beeamerino.github.io/PhyloCactus/reference/classify_constrained_nodes.md)
  : Classify Every Internal Node as Constrained or Estimated

## Divergence Time Estimation & treePL Automation

Automate `treePL` penalized-likelihood cross-validation and estimate
divergence time chronograms over temporal bootstrap cohorts.

- [`rescale_tree()`](https://beeamerino.github.io/PhyloCactus/reference/rescale_tree.md)
  : Rescale Branch Lengths of a Phylogenetic Tree
- [`run_treePL_direct()`](https://beeamerino.github.io/PhyloCactus/reference/run_treePL_direct.md)
  : Run treePL Executable Directly
- [`run_treePL_cv()`](https://beeamerino.github.io/PhyloCactus/reference/run_treePL_cv.md)
  : Prime, Cross-Validate and Date a Tree with treePL
- [`automate_treePL()`](https://beeamerino.github.io/PhyloCactus/reference/automate_treePL.md)
  : Automate treePL Divergence Time Estimation Pipeline Across Bootstrap
  Cohorts
- [`send_run_notification()`](https://beeamerino.github.io/PhyloCactus/reference/send_run_notification.md)
  : Send a Run Notification by Email
- [`report_smoothing_sensitivity()`](https://beeamerino.github.io/PhyloCactus/reference/report_smoothing_sensitivity.md)
  : Age of Every Calibrated Node Across a Range of Rate-Smoothing Values
- [`check_calibration_consistency()`](https://beeamerino.github.io/PhyloCactus/reference/check_calibration_consistency.md)
  : Check That a Set of Calibration Bounds Is Internally Satisfiable
- [`audit_calibration_support()`](https://beeamerino.github.io/PhyloCactus/reference/audit_calibration_support.md)
  : Audit the Nodes a Calibration Table Addresses
- [`calibration_tips()`](https://beeamerino.github.io/PhyloCactus/reference/calibration_tips.md)
  : Resolve the Terminals a Calibration Row Addresses
- [`report_bound_adherence()`](https://beeamerino.github.io/PhyloCactus/reference/report_bound_adherence.md)
  : Report Which Calibrated Nodes Came Back Sitting on a Bound

## Topological Validation & Tree Space Projections

Quantify topological congruence across comparative trees using
Robinson-Foulds distances and multidimensional tree-space projections.

- [`validate_phylogenies()`](https://beeamerino.github.io/PhyloCactus/reference/validate_phylogenies.md)
  : Comparative Tree-Space Projection and Topological Validation
  Pipeline

## Publication Integration & Clade Summaries

Map branch support values, collapse poorly supported nodes into soft
polytomies (`integrate_publication_tree`), and export figures.

- [`integrate_publication_tree()`](https://beeamerino.github.io/PhyloCactus/reference/integrate_publication_tree.md)
  : Render Final Publication Figures and Registry

## IUCN Conservation Metadata Integration

Extract threat metrics and habitats via `rredlist` and apply
standardized IUCN color scales to phylogenetic trees.

- [`get_iucn_data()`](https://beeamerino.github.io/PhyloCactus/reference/get_iucn_data.md)
  : Fetch IUCN Red List Metadata via rredlist API
- [`scale_fill_iucn()`](https://beeamerino.github.io/PhyloCactus/reference/scale_fill_iucn.md)
  : Standardized IUCN Category Color Scale for ggplot2
- [`phylocactus_colours()`](https://beeamerino.github.io/PhyloCactus/reference/phylocactus_colours.md)
  : Colours of the Figures of PhyloCactus
- [`phylocactus_panel_tags()`](https://beeamerino.github.io/PhyloCactus/reference/phylocactus_panel_tags.md)
  : Panel Tags of the Composite Figures of PhyloCactus
- [`theme_phylocactus()`](https://beeamerino.github.io/PhyloCactus/reference/theme_phylocactus.md)
  : Theme of the Figures of PhyloCactus

## Barcoding Library

Build the reference library of the barcoding branch from the phylotaR
mining of the phylogeny and, optionally, own sequences (Tutorial 5,
steps 1 to 4).

- [`assemble_barcoding_dataset()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_barcoding_dataset.md)
  : Assemble the Reference Library of the Molecular Diagnostic Branch
- [`curate_barcoding_markers()`](https://beeamerino.github.io/PhyloCactus/reference/curate_barcoding_markers.md)
  : Curate the Loci of the Molecular Diagnostic Branch
- [`screen_barcoding_markers()`](https://beeamerino.github.io/PhyloCactus/reference/screen_barcoding_markers.md)
  : Screen the Curated Loci of the Molecular Diagnostic Branch
- [`finalize_barcoding_library()`](https://beeamerino.github.io/PhyloCactus/reference/finalize_barcoding_library.md)
  : Build the Final Reference Library of the Molecular Diagnostic Branch
- [`extract_barcoding_plastome_loci()`](https://beeamerino.github.io/PhyloCactus/reference/extract_barcoding_plastome_loci.md)
  : Extra library records cropped from GenBank plastomes (step 1)
- [`barcoding_outgroup_taxids()`](https://beeamerino.github.io/PhyloCactus/reference/barcoding_outgroup_taxids.md)
  : Taxonomic IDs of the Outgroup of the Molecular Diagnostic Branch

## Barcoding Validation

Validation folds, barcode gap, classification of the folds by nearest
neighbour and IdTaxa (on one machine or a computing cluster), negative
controls, thresholds and metrics (Tutorial 5, steps 5 to 10).

- [`build_barcoding_folds()`](https://beeamerino.github.io/PhyloCactus/reference/build_barcoding_folds.md)
  : Build the Validation Folds of the Molecular Diagnostic Branch
- [`analyze_barcode_gap()`](https://beeamerino.github.io/PhyloCactus/reference/analyze_barcode_gap.md)
  : Barcode Gap of the Molecular Diagnostic Branch
- [`classify_barcoding_folds()`](https://beeamerino.github.io/PhyloCactus/reference/classify_barcoding_folds.md)
  : Classify the Validation Folds of the Molecular Diagnostic Branch
- [`generate_barcoding_job_scripts()`](https://beeamerino.github.io/PhyloCactus/reference/generate_barcoding_job_scripts.md)
  : Generate the SLURM Job Scripts of a Chunked Run of Step 7
- [`merge_barcoding_chunks()`](https://beeamerino.github.io/PhyloCactus/reference/merge_barcoding_chunks.md)
  : Merge the Chunks of a Chunked Run of Step 7
- [`run_barcoding_controls()`](https://beeamerino.github.io/PhyloCactus/reference/run_barcoding_controls.md)
  : Negative Controls of the Molecular Diagnostic Branch
- [`sweep_barcoding_threshold()`](https://beeamerino.github.io/PhyloCactus/reference/sweep_barcoding_threshold.md)
  : Sweep the Threshold of Remoteness of the Molecular Diagnostic Branch
- [`summarise_barcoding_metrics()`](https://beeamerino.github.io/PhyloCactus/reference/summarise_barcoding_metrics.md)
  : Summarise the Per-Class Metrics of the Molecular Diagnostic Branch
- [`summarise_barcoding_genus_discrimination()`](https://beeamerino.github.io/PhyloCactus/reference/summarise_barcoding_genus_discrimination.md)
  : Per-genus discrimination of each locus in the validation folds

## Molecular Identification

Identify a sample against the library, locus by locus, from sequences, a
sample sheet, an assembly or Illumina reads, write its report, and draw
its answer (Tutorials 5 and 6).

- [`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md)
  : Identify user sequences against the barcoding library (step 11 of
  Tutorial 5)
- [`identify_barcoding_samples()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_samples.md)
  : Identify the samples of a sample sheet
- [`identify_barcoding_assembly()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_assembly.md)
  : Identify a genome assembly against the library
- [`identify_barcoding_reads()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_reads.md)
  : Identify paired Illumina reads through route A
- [`report_barcoding_identification()`](https://beeamerino.github.io/PhyloCactus/reference/report_barcoding_identification.md)
  : Build the identification report of one run
- [`plot_barcoding_identification()`](https://beeamerino.github.io/PhyloCactus/reference/plot_barcoding_identification.md)
  : Draw the answer of an identification
- [`plot_barcoding_identification_summary()`](https://beeamerino.github.io/PhyloCactus/reference/plot_barcoding_identification_summary.md)
  : Draw the answers of several identifications side by side

## Internal Tools and Visualization Helpers

Internal functions for phylogenetic visualization and data preparation.

- [`resolve_run_paths()`](https://beeamerino.github.io/PhyloCactus/reference/resolve_run_paths.md)
  [`print(`*`<cactus_run_paths>`*`)`](https://beeamerino.github.io/PhyloCactus/reference/resolve_run_paths.md)
  : Resolve Every File Path of a Phylogenetic Run from Its Naming
  Convention
- [`resolve_rooting_outgroup()`](https://beeamerino.github.io/PhyloCactus/reference/resolve_rooting_outgroup.md)
  : Resolve the Set of Rooting Terminals from a Tip Label Vector
- [`root_on_clade()`](https://beeamerino.github.io/PhyloCactus/reference/root_on_clade.md)
  : Root a Tree on a Clade of Terminals, Tolerating a Basal Polytomy
- [`add_clade_labels_by_level()`](https://beeamerino.github.io/PhyloCactus/reference/add_clade_labels_by_level.md)
  : Add Clade Labels to ggtree Plot by Taxonomic Level
- [`apply_clade_label_layers()`](https://beeamerino.github.io/PhyloCactus/reference/apply_clade_label_layers.md)
  : Apply Multiple Taxonomic Clade Label Layers onto ggtree Plot
- [`assert_constraint_columns()`](https://beeamerino.github.io/PhyloCactus/reference/assert_constraint_columns.md)
  : Assert Mandatory Constraint Table Columns
- [`augment_treedata_with_registry()`](https://beeamerino.github.io/PhyloCactus/reference/augment_treedata_with_registry.md)
  : Augment treedata Object with Registry Metadata
- [`build_annotation_registry()`](https://beeamerino.github.io/PhyloCactus/reference/build_annotation_registry.md)
  : Build Taxonomic Annotation Registry across Tree Levels
- [`compute_group_nodes()`](https://beeamerino.github.io/PhyloCactus/reference/compute_group_nodes.md)
  : Compute Most Recent Common Ancestor (MRCA) Nodes for Taxon Groups
- [`extract_species_binomial()`](https://beeamerino.github.io/PhyloCactus/reference/extract_species_binomial.md)
  : Extract Standardized Species Binomial from Sequence Header
- [`get_annotation_nodes()`](https://beeamerino.github.io/PhyloCactus/reference/get_annotation_nodes.md)
  : Extract Monophyletic Internal Annotation Nodes
- [`prepare_tip_annotation()`](https://beeamerino.github.io/PhyloCactus/reference/prepare_tip_annotation.md)
  : Prepare Tip Annotations Matched to Tree Leaves
- [`standardize_taxon()`](https://beeamerino.github.io/PhyloCactus/reference/standardize_taxon.md)
  : Standardize Taxon Name Strings
- [`get_tree_max_depth()`](https://beeamerino.github.io/PhyloCactus/reference/get_tree_max_depth.md)
  : Compute Maximum Depth of a Dated Chronogram Tree
- [`format_age_labels()`](https://beeamerino.github.io/PhyloCactus/reference/format_age_labels.md)
  : Format Age Labels for Chronogram Axis Ticks
- [`add_chronogram_axis()`](https://beeamerino.github.io/PhyloCactus/reference/add_chronogram_axis.md)
  : Add Geological Time Scale Axis to a Chronogram ggtree Plot
