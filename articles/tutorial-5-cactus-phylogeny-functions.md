# Tutorial 5: Function Reference & Dictionary

## Overview and Methodological Dictionary

Following the step-by-step analytical modules presented in Tutorials 1
through 4, this reference vignette establishes a complete dictionary of
the core R functions implemented in `PhyloCactus`.

### Narrative Framework of the Reference Suite

- **Where are we?** Having assembled supermatrices, inferred
  maximum-likelihood phylogenies, estimated divergence times with
  `treePL`, and validated topologies against external backbones, we
  require an organized reference dictionary detailing function
  signatures and biological motivations.
- **Why are we here?** Complex computational pipelines require explicit
  documentation linking software parameters directly to underlying
  evolutionary mechanics, such as mitigating long-branch attraction
  (LBA), enforcing positional homology, and accounting for incomplete
  lineage sorting (ILS).
- **What will we achieve?** A complete, reproducible function reference
  mapping sequence mining, alignment quality control, substitution model
  evaluation, maximum-likelihood search, temporal calibration, and
  visualization utilities.
- **What comes next?** Users can query individual function manual pages
  (`?function_name`) or incorporate specific functions into customized
  plant phylogenomic pipelines.

## Function Reference Dictionary

### Sequence Assembly & Mining

Functions in this category retrieve orthologous sequence clusters from
public repositories using sequence similarity clustering (`phylotaR`),
avoiding nomenclature ambiguities across public databases.

| Function | Biological & Methodological Application | Canonical Citation |
|:---|:---|:---|
| [`assemble_ingroup_phylotar()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_ingroup_phylotar.md) | Mines orthologous DNA sequence clusters for a focal ingroup clade (e.g., Cactaceae; Guerrero et al. 2019) from GenBank using `phylotaR`. Applies strict taxonomic filters, evaluates cluster occupancy, and resolves nomenclature inconsistencies. | Bennett et al. (2018) |
| [`assemble_outgroup_phylotar()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_outgroup_phylotar.md) | Retrieves orthologous sequences for specified outgroup lineages to establish accurate root positions for downstream likelihood search and chronogram estimation. | Bennett et al. (2018) |
| [`fetch_genbank_metadata()`](https://beeamerino.github.io/PhyloCactus/reference/fetch_genbank_metadata.md) | Interrogates NCBI Entrez utilities to extract sequence lengths, organism taxonomy, publication titles, and accession IDs for downloaded sequence clusters. | Camacho et al. (2009) |

### Alignment & Quality Control

Functions in this category establish hypotheses of positional homology
(`MAFFT`) and apply automated masking of non-homologous or ambiguous
regions (`DECIPHER`).

| Function | Biological & Methodological Application | Canonical Citation |
|:---|:---|:---|
| [`run_alignment_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_alignment_pipeline.md) | Generates multiple sequence alignments using `MAFFT` and masks poorly aligned, ambiguous, or hypervariable regions using `DECIPHER`, eliminating systematic noise while retaining phylogenetically informative sites. | Katoh & Standley (2013); Wright (2024) |
| [`run_marker_screening()`](https://beeamerino.github.io/PhyloCactus/reference/run_marker_screening.md) | Evaluates phylogenetic informativeness, sequence coverage, and substitution saturation across individual locus alignments, excluding loci exhibiting high saturation slope erosion. | Xia et al. (2003) |
| [`clean_taxonomic_names()`](https://beeamerino.github.io/PhyloCactus/reference/clean_taxonomic_names.md) | Standardizes tip labels against authoritative botanical checklists (e.g., Caryophyllales.org checklist), resolving synonymies and orthographic variants. | Korotkova et al. (2021) |
| [`run_joint_realignment()`](https://beeamerino.github.io/PhyloCactus/reference/run_joint_realignment.md) | Re-estimates positional homology alignments (`MAFFT`) after aggregating ingroup and outgroup sequence markers to refine alignment hypotheses. | Katoh & Standley (2013) |

### Concatenation & Supermatrix Assembly

| Function | Biological & Methodological Application | Canonical Citation |
|:---|:---|:---|
| [`run_concatenation_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_concatenation_pipeline.md) | Concatenates individual orthologous locus alignments into a unified multilocus supermatrix and generates partition coordinate maps compatible with `RAxML-NG`, `IQ-TREE`, and `MrBayes`. | Kozlov et al. (2019) |

### Likelihood Inference & Model Evaluation

Functions in this category select optimal nucleotide substitution models
per partition and infer maximum-likelihood phylogenies under topological
constraints.

| Function | Biological & Methodological Application | Canonical Citation |
|:---|:---|:---|
| [`preprocess_partitions()`](https://beeamerino.github.io/PhyloCactus/reference/preprocess_partitions.md) | Validates structural integrity and syntax of the supermatrix and its corresponding partition definitions prior to likelihood evaluation using `RAxML-NG`. | Kozlov et al. (2019) |
| [`run_modeltest_ng()`](https://beeamerino.github.io/PhyloCactus/reference/run_modeltest_ng.md) | Identifies the optimal substitution model per partition (`ModelTest-NG`) using AICc criteria, controlling for mutational rate heterogeneity and mitigating long-branch attraction (LBA). | Darriba et al. (2020) |
| [`build_constraint_scaffold()`](https://beeamerino.github.io/PhyloCactus/reference/build_constraint_scaffold.md) | Synthesizes a Newick multifurcating constraint scaffold enforcing monophyletic relationships of established higher taxonomic ranks. | Guerrero et al. (2019) |
| [`calculate_ml_tree()`](https://beeamerino.github.io/PhyloCactus/reference/calculate_ml_tree.md) | Infers the maximum-likelihood evolutionary hypothesis explaining the supermatrix under selected partition models and constraint scaffolds (`RAxML-NG`). | Kozlov et al. (2019) |
| [`map_branch_supports()`](https://beeamerino.github.io/PhyloCactus/reference/map_branch_supports.md) | Maps Transfer Bootstrap Expectation (TBE) support values onto the best maximum-likelihood tree, providing robust clade support for large phylogenies. | Lemoine et al. (2018) |

### Divergence Time Estimation & treePL Automation

| Function | Biological & Methodological Application | Canonical Citation |
|:---|:---|:---|
| [`rescale_tree()`](https://beeamerino.github.io/PhyloCactus/reference/rescale_tree.md) | Multiplies tree edge lengths by a scaling factor prior to penalized likelihood dating (`treePL`) to prevent numerical underflow. | Sanderson (2002) |
| [`run_treePL()`](https://beeamerino.github.io/PhyloCactus/reference/run_treePL.md) | Interfacing with `treePL` via a shell wrapper script to estimate ultrametric chronograms under penalized likelihood. | Smith & O’Meara (2012) |
| [`automate_treePL()`](https://beeamerino.github.io/PhyloCactus/reference/automate_treePL.md) | Automates cross-validation parameter optimization, rate smoothing selection, and chronogram estimation across temporal bootstrap replicates. | Maurin (2020) |

### Comparative Validation & IUCN Utilities

| Function | Biological & Methodological Application | Canonical Citation |
|:---|:---|:---|
| [`validate_phylogenies()`](https://beeamerino.github.io/PhyloCactus/reference/validate_phylogenies.md) | Quantifies topological congruence between the focal supermatrix tree and external published backbones using Robinson-Foulds distances and multidimensional scaling (MDS). | Robinson & Foulds (1981); Smith (2020) |
| [`integrate_publication_tree()`](https://beeamerino.github.io/PhyloCactus/reference/integrate_publication_tree.md) | Maps statistical support onto chronograms, collapsing weakly supported nodes (TBE \< 0.70) into soft polytomies to represent analytical uncertainty. | Lemoine et al. (2018) |
| [`get_iucn_data()`](https://beeamerino.github.io/PhyloCactus/reference/get_iucn_data.md) | Interrogates the IUCN Red List API via `rredlist` to retrieve extinction threat categories, population trends, and geographic endemism metrics. | IUCN (2024) |
| [`scale_fill_iucn()`](https://beeamerino.github.io/PhyloCactus/reference/scale_fill_iucn.md) | Provides a standardized `ggplot2` manual fill scale mapping official IUCN threat categories onto comparative figures. | IUCN (2024) |

### Reproducible Execution Example

``` r

library(PhyloCactus)

# Query help documentation for any pipeline function
?run_alignment_pipeline
?calculate_ml_tree
?automate_treePL
```
