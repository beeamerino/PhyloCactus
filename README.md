# PhyloCactus <img src="man/figures/logo.png" align="right" height="300"/>

<!-- badges: start -->

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)[![R-CMD-check](https://github.com/beeamerino/PhyloCactus/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/beeamerino/PhyloCactus/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

**Reproducible phylogenetic workflows for data assembly, inference, and divergence time estimation**

Reconstructing evolutionary relationships within plant clades characterized by rapid adaptive radiations—such as the family **Cactaceae**—presents major methodological challenges. Low plastid sequence divergence, persistent gene tree discordance driven by incomplete lineage sorting (ILS) and reticulate evolution, and pervasive taxonomic synonymies in public sequence repositories hinder computational reproducibility.

`PhyloCactus` provides a standardized, reproducible R workflow designed to address these biological challenges. By integrating sequence retrieval, taxonomic reconciliation, alignment quality control, substitution model evaluation, maximum-likelihood inference, divergence time estimation, and comparative tree validation into a unified pipeline, `PhyloCactus` transforms heterogeneous genomic sequence data into statistically evaluated phylogenetic hypotheses. While optimized for Cactaceae, the modular architecture applies to any plant lineage requiring rigorous phylogenetic data preparation and macroevolutionary analysis.

## Main Features

The package provides tools for:

- Automating the retrieval of orthologous DNA sequence clusters from GenBank using similarity clustering via `phylotaR` (Bennett *et al.* 2018) and `BLAST+` (Camacho *et al.* 2009).
- Inferring positional homology alignments using `MAFFT` (Katoh & Standley 2013).
- Mitigating systematic alignment errors through automated masking of non-homologous segments with `DECIPHER` (Wright 2024).
- Standardizing nomenclatural frameworks to assemble curated, unified multilocus datasets against dynamic taxonomic checklists.
- Synthesizing concatenated supermatrices and defining structural partition bounds for partitioned phylogenetic inference.
- Statistically controlling for mutational heterogeneity and long-branch attraction by evaluating nucleotide substitution models via `ModelTest-NG` (Darriba *et al.* 2020).
- Inferring maximum-likelihood evolutionary hypotheses using `RAxML-NG` (Kozlov *et al.* 2019).
- Mapping Transfer Bootstrap Expectation (`TBE`; Lemoine *et al.* 2018) onto maximum-likelihood topologies to evaluate clade support under missing data.
- Estimating ultrametric chronograms to accommodate evolutionary rate heterogeneity via penalized likelihood with `treePL` (Sanderson 2002; Smith & O’Meara 2012).
- Integrating complex taxonomic and conservation metadata, including statuses from the **IUCN Red List**.
- Producing publication-quality phylogenetic figures and chronological cadastres.
- Quantifying topological congruence and validating alternative phylogenetic hypotheses via multidimensional tree-space projections.

## Workflow

The complete workflow is organized into four analytical stages comprising thirteen interoperable modules.

![PhyloCactus Workflow Architecture](man/figures/Fig1_PhyloCactus_Workflow.png)

> Each stage can be executed independently or combined into a fully reproducible phylogenetic workflow.

## Installation

Install the development version directly from GitHub:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

remotes::install_github("beeamerino/PhyloCactus")
```

## R Dependencies

`PhyloCactus` relies on packages available from both CRAN and Bioconductor:

| Repository | Main packages |
|:---|:---|
| **remotes** | `phylotaR` |
| **Bioconductor** | `DECIPHER`, `Biostrings` |
| **CRAN** | `ape`, `ggplot2`, `dplyr`, `tidyr`, `readr`, `stringr`, `purrr`, `rredlist`, `forcats`, `scales`, `RColorBrewer` |

## External Software

Several analyses performed by `PhyloCactus` rely on external command-line software that must be installed separately and available from your system `PATH`:

| Software | Purpose | Citation |
|:---|:---|:---|
| [`BLAST+`](https://www.ncbi.nlm.nih.gov/books/NBK279690/) | Local sequence alignment and database searching | Camacho, C. *et al*. (2009). BLAST+: Architecture and applications. *BMC Bioinformatics*, 10. <https://doi.org/10.1186/1471-2105-10-421> |
| [`MAFFT`](https://mafft.cbrc.jp/alignment/software/) | Multiple sequence alignment | Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: Improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772–780. <https://doi.org/10.1093/molbev/mst010> |
| [`ModelTest-NG`](https://github.com/ddarriba/modeltest) | Model selection | Darriba, D. *et al*. (2020). ModelTest-NG: a new and scalable tool for the selection of DNA and protein evolutionary models. *Molecular Biology and Evolution*, 37(1), 291-294. <https://doi.org/10.1093/molbev/msz189> |
| [`RAxML-NG`](https://github.com/amkozlov/raxml-ng) | Maximum-likelihood phylogenetic inference | Kozlov, A. M. *et al*. (2019). RAxML-NG: A fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics*, 35(21), 4453–4455. <https://doi.org/10.1093/bioinformatics/btz305> |
| [`treePL`](https://github.com/blackrim/treePL) | Divergence time estimation | Smith, S. A., & O’Meara, B. C. (2012). TreePL: Divergence time estimation using penalized likelihood for large phylogenies. *Bioinformatics*, 28(20), 2689–2690. <https://doi.org/10.1093/bioinformatics/bts492> |

## Documentation

The recommended way to learn `PhyloCactus` is by following the tutorials in sequence:

| Tutorial | Description |
|:---|:---|
| [**Get Started**](https://beeamerino.github.io/PhyloCactus/articles/PhyloCactus.html) | Package overview and installation |
| [**Tutorial 1**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-1-cactus-phylogeny-prep.html) | Data Assembly and Preparation |
| [**Tutorial 2**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-2-cactus-phylogeny-inference.html) | Phylogenetic Inference and Divergence Time Estimation |
| [**Tutorial 3**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-3-cactus-phylogeny-visualization.html) | Visualization and Metadata Integration |
| [**Tutorial 4**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-4-cactus-phylogeny-validation.html) | Phylogenetic Validation and Comparative Analyses |
| [**Function Reference**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-5-cactus-phylogeny-functions.html) | Complete reference of package functions |

## Citation

If you use `PhyloCactus` in your research, please cite the package:

``` r
citation("PhyloCactus")
```

If a software paper is available, please cite both the package and the associated publication.

## License

`PhyloCactus` is released under the **GPL-3 License**.

> Bug reports, feature requests, and source code are available through the project's GitHub repository.
