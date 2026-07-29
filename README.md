# PhyloCactus <img src="man/figures/logo.png" align="right" height="300"/>

<!-- badges: start -->

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)[![R-CMD-check](https://github.com/beeamerino/PhyloCactus/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/beeamerino/PhyloCactus/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

**Reproducible R workflows for multilocus data assembly, evolutionary modeling, maximum-likelihood inference, and divergence time estimation in plant radiations**

## Overview & Biological Motivation

Constructing reliable multilocus molecular sequence matrices directly from public repositories such as GenBank presents a fundamental computational challenge in systematic biology. Sequence records deposited across decades frequently exhibit inconsistent locus annotations, duplicated accessions, unvouchered identifications, orthographic variants, and massive nomenclatural synonymies. Transforming these uncurated records into high-quality phylogenetic matrices requires intensive manual curation, taxonomic reconciliation against authoritative botanical checklists, and rigorous alignment quality control.

This computational problem becomes exponentially more complex when targeting plant lineages with intricate evolutionary histories. Clades characterized by recent explosive adaptive radiations, low plastid sequence divergence, incomplete lineage sorting (ILS), and ancient reticulate evolution amplify the risk of misidentifying orthology, accumulating systematic alignment noise, and producing biased topological reconstructions.

The family **Cactaceae** serves as the prime empirical exemplar of these combined computational and biological hurdles. Comprising one of the largest succulent plant radiations in the Neotropics (Arakaki *et al*. 2011; Hernández-Hernández *et al*. 2014; Guerrero *et al*. 2019), cactus phylogenetics requires extensive data curation to resolve persistent gene tree discordance and handle heterogeneous molecular datasets.

`PhyloCactus` was developed to transform this complex data preparation and analytical process into a fully automated, transparent, and reproducible R workflow. Rather than relying on manual ad hoc scripts, `PhyloCactus` orchestrates established bioinformatics software into a standardized 4-stage pipeline. The package automates orthology-based sequence retrieval (`phylotaR`), alignment of positional homology (`MAFFT`), objective quality control masking (`DECIPHER`), mutational saturation screening, substitution model evaluation (`ModelTest-NG`), constrained maximum-likelihood topology inference (`RAxML-NG`), Transfer Bootstrap Expectation support mapping (`TBE`), penalized likelihood chronogram estimation (`treePL`), IUCN Red List metadata enrichment (`rredlist`), and multispecies coalescent tree-space validation (`ASTRAL-III`).

## Four-Stage Analytical Architecture

The `PhyloCactus` pipeline is organized into **four analytical stages** comprising **thirteen interoperable modules**:

![PhyloCactus Workflow Architecture](man/figures/Fig1_PhyloCactus_Workflow.png)

1. **Stage 1: Data Assembly and Preparation (Modules 1 to 6)**  
   Retrieves orthologous sequence clusters using similarity-based clustering via `phylotaR`, bypassing GenBank annotation errors. Establishes positional homology with `MAFFT`, applies objective quality control masking with `DECIPHER`, screens locus alignments for mutational saturation slope erosion (Xia test), reconciles species nomenclature against the Caryophyllales.org checklist (`CactaceaeFullList_accepted.csv`; Korotkova *et al*. 2021), executes joint alignment of ingroup (**Cactaceae**) and outgroup (*Portulaca*, *Anacampseros*, *Talinopsis*, *Grahamia*) markers, and concatenates locus alignments into partitioned supermatrices with explicit coordinate boundaries.

2. **Stage 2: Phylogenetic Inference and Dating (Modules 7 to 10)**  
   Evaluates partition-specific nucleotide substitution models using `ModelTest-NG` under corrected Akaike Information Criteria (AICc) to control mutational rate heterogeneity and mitigate Long-Branch Attraction (LBA). Infers maximum-likelihood topologies using `RAxML-NG` under topological constraint scaffolds (`cactus_constraints.csv`) enforcing higher-level clade monophyly (**Cactoideae**, **Opuntioideae**, **Leuenbergeria**, **Pereskia**). Maps clade support via Transfer Bootstrap Expectation (`TBE`) to handle missing data, generates topologically constrained temporal bootstrap alignments, and estimates ultrametric chronograms using penalized likelihood in `treePL` guided by secondary calibration boundaries (`calibrations_bounds.csv`).

3. **Stage 3: Visualization and Metadata Integration (Modules 11 and 12)**  
   Enriches the phylogenetic framework with species-level conservation attributes retrieved automatically from the IUCN Red List database via `rredlist`. Collapses weakly supported internal nodes (TBE < 0.70) into soft politomies to prevent over-interpreting unresolved rapid radiation nodes, and renders publication-ready phylogenetic figures, chronograms, and conservation cadastres using `ggplot2`.

4. **Stage 4: Validation and Sub-tree Comparisons (Module 13)**  
   Quantifies topological congruence between the focal supermatrix tree and external published backbones (e.g., Amaral *et al*. 2022; Thompson *et al*. 2024; de Vos *et al*. 2025) using Robinson-Foulds distances and multidimensional scaling (MDS) tree-space projections (`validate_phylogenies`). Infers summary species trees using `ASTRAL-III` under multispecies coalescent theory to evaluate biological gene tree discordance driven by incomplete lineage sorting.

## Installation

Install the development version of `PhyloCactus` directly from GitHub:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

remotes::install_github("beeamerino/PhyloCactus")
```

Verify package installation:

``` r
library(PhyloCactus)

packageVersion("PhyloCactus")
```

## R Dependencies

`PhyloCactus` integrates tools from CRAN, Bioconductor, and GitHub:

| Repository | Main packages |
|:---|:---|
| **remotes** | `phylotaR` |
| **Bioconductor** | `DECIPHER`, `Biostrings` |
| **CRAN** | `ape`, `ggplot2`, `dplyr`, `tidyr`, `readr`, `stringr`, `purrr`, `rredlist`, `forcats`, `scales`, `RColorBrewer` |

To install the Bioconductor dependencies:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DECIPHER", "Biostrings"))
```

## External Software Binaries

Execution of `PhyloCactus` requires the following external command-line binaries, which must be installed and available in your system `$PATH`:

| Software | Purpose | Citation |
|:---|:---|:---|
| [`BLAST+`](https://www.ncbi.nlm.nih.gov/books/NBK279690/) | Local sequence similarity searching & cluster identification | Camacho *et al*. (2009) *BMC Bioinformatics* |
| [`MAFFT`](https://mafft.cbrc.jp/alignment/software/) | Multiple sequence alignment & positional homology | Katoh & Standley (2013) *Mol. Biol. Evol.* |
| [`ModelTest-NG`](https://github.com/ddarriba/modeltest) | Partitioned substitution model selection under AICc | Darriba *et al*. (2020) *Mol. Biol. Evol.* |
| [`RAxML-NG`](https://github.com/amkozlov/raxml-ng) | Constrained maximum-likelihood tree search & TBE support | Kozlov *et al*. (2019) *Bioinformatics* |
| [`treePL`](https://github.com/blackrim/treePL) | Penalized likelihood divergence time estimation | Smith & O’Meara (2012) *Bioinformatics* |

Ensure these executables are accessible by adding their installation path to your `~/.Renviron` file (`usethis::edit_r_environ()`).

## Documentation & Tutorials

Learn `PhyloCactus` through the sequential tutorial suite:

| Tutorial | Focus & Scope |
|:---|:---|
| [**Get Started**](https://beeamerino.github.io/PhyloCactus/articles/PhyloCactus.html) | Package overview, design principles, and dependency configuration |
| [**Tutorial 1**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-1-cactus-phylogeny-prep.html) | Stage 1: Data Assembly, Alignment, Saturation Screening, & Supermatrix Concatenation |
| [**Tutorial 2**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-2-cactus-phylogeny-inference.html) | Stage 2: Substitution Modeling, Constrained ML Search, TBE Support, & treePL Dating |
| [**Tutorial 3**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-3-cactus-phylogeny-visualization.html) | Stage 3: IUCN Red List Enrichment, Soft Polytomy Collapsing, & Figure Rendering |
| [**Tutorial 4**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-4-cactus-phylogeny-validation.html) | Stage 4: ASTRAL-III Coalescence, Robinson-Foulds Distances, & MDS Tree-Space Validation |
| [**Function Reference**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-5-cactus-phylogeny-functions.html) | Comprehensive dictionary of package functions and analytical signatures |

## Citation

If you use `PhyloCactus` in your research, please cite:

``` r
citation("PhyloCactus")
```

- Arakaki *et al*. 2011. Contemporaneous and recent radiations of the world’s major succulent plant lineages. *PNAS*, 108(20), 8379–8384.
- Guerrero *et al*. 2019. Phylogenetic Relationships and Evolutionary Trends in the Cactus Family. *Journal of Heredity*, 110(1), 4–21.
- Hernández-Hernández *et al*. 2014. Beyond aridification: Multiple explanations for the elevated diversification of cacti in the New World Succulent Biome. *New Phytologist*, 202(4), 1382–1397.
- Korotkova *et al*. 2021. Cactaceae at Caryophyllales.org - A dynamic online species-level taxonomic backbone for the family. *Willdenowia*, 51(2), 251–270. <https://doi.org/10.3372/wi.51.51208>

## License

`PhyloCactus` is licensed under the **GNU General Public License v3.0 (GPL-3)**.


