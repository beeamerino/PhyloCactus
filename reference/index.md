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

## DNA Barcoding (section under construction)

Molecular diagnostic section. No implementation is shipped while its
analytical strategy is being defined; see the section manifest in
Tutorial 5.

- [`barcoding`](https://beeamerino.github.io/PhyloCactus/reference/barcoding.md)
  : DNA Barcoding Evaluation (Section Under Construction)
- [`evaluate_dna_barcoding()`](https://beeamerino.github.io/PhyloCactus/reference/evaluate_dna_barcoding.md)
  : Evaluate DNA Barcoding Resolution and Taxonomic Classification Power

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
- [`infer_gene_trees()`](https://beeamerino.github.io/PhyloCactus/reference/infer_gene_trees.md)
  : Infer Single-Locus Maximum-Likelihood Gene Trees and Assess Species
  Monophyly

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
