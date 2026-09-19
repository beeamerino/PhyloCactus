# PhyloCactus

**Reproducible R workflows for multilocus data assembly, evolutionary
modeling, maximum-likelihood inference, and divergence time estimation
in plant radiations**

## Overview

Multilocus sequence matrices assembled from public repositories such as
GenBank inherit the heterogeneity of decades of deposition: inconsistent
locus annotations, duplicated accessions, records without vouchers,
orthographic variants and nomenclatural synonyms. Converting these
records into a phylogenetic matrix requires curation, taxonomic
reconciliation against an accepted checklist, and quality control of
each alignment.

The difficulty increases in plant lineages that diversified recently and
rapidly. Low plastid divergence, incomplete lineage sorting (ILS) and
reticulate evolution increase the risk of erroneous orthology assessment
and of alignment error, and both can bias the inferred topology.

**Cactaceae** combine these properties. The family is one of the largest
succulent radiations of the Neotropics (Arakaki *et al*. 2011;
Hernández-Hernández *et al*. 2014; Guerrero *et al*. 2019), and its
GenBank record is uneven in taxonomic and locus coverage.

`PhyloCactus` is an R package that implements this assembly and analysis
as a reproducible workflow. It links established software in four
stages: orthology-based sequence retrieval (`phylotaR`), multiple
sequence alignment (`MAFFT`), alignment masking (`DECIPHER`), a
saturation screen, substitution model selection (`ModelTest-NG`),
constrained maximum-likelihood inference (`RAxML-NG`), Felsenstein
bootstrap proportions (FBP), penalized likelihood dating (`treePL`),
retrieval of IUCN Red List assessments (`rredlist`), and comparison with
published phylogenies, one of which can be summarized under the
multispecies coalescent with `ASTRAL-III`.

## Workflow

The workflow comprises thirteen modules grouped in four stages.

![PhyloCactus Workflow
Architecture](reference/figures/Fig1_PhyloCactus_Workflow.png)

PhyloCactus Workflow Architecture

1.  **Stage 1: Data Assembly and Preparation (Modules 1 to 6)**  
    Retrieves orthologous sequence clusters from GenBank by
    similarity-based clustering with `phylotaR`, which reduces the
    dependence on locus annotations. Aligns each locus with `MAFFT` and
    masks unreliable alignment regions with `DECIPHER`. Screens each
    locus with a regression of uncorrected p-distance against a
    Gamma-corrected K80 distance (Kimura 1980; Jin & Nei 1990), a
    criterion implemented in this package, and rejects loci whose slope
    falls below the retention threshold. Reconciles species names
    against the Caryophyllales.org checklist
    (`CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx`; Korotkova *et
    al*. 2021), realigns ingroup (Cactaceae) and outgroup (*Portulaca*,
    *Anacampseros*, *Talinopsis*, *Grahamia*, *Talinum*, *Talinella*)
    sequences jointly, and concatenates the loci into a partitioned
    supermatrix.

2.  **Stage 2: Phylogenetic Inference and Dating (Modules 7 to 10)**  
    Selects a substitution model for each partition with `ModelTest-NG`
    under the corrected Akaike information criterion (AICc). Infers the
    maximum-likelihood tree with `RAxML-NG` under a topological
    constraint (`cactus_constraints.csv`) that fixes 34 of the 1021
    internal bipartitions of the reference tree, all of them among its
    deepest. The constraint imposes the monophyly of 21 named clades and
    the branching order among them within Cactoideae and within
    Opuntioideae; the remaining 987 bipartitions are estimated.
    Bootstrap values on imposed bipartitions reproduce the constraint
    and are reported as constrained
    ([`classify_constrained_nodes()`](https://beeamerino.github.io/PhyloCactus/reference/classify_constrained_nodes.md)).
    Support is summarized as FBP, with the transfer bootstrap
    expectation (TBE) as a secondary measure. Divergence times are
    estimated by penalized likelihood in `treePL` with secondary
    calibrations (`calibrations_bounds.csv`), and their confidence
    intervals from temporal bootstrap replicates.

3.  **Stage 3: Visualization and Metadata Integration (Modules 11 and
    12)**  
    Adds species-level IUCN Red List assessments retrieved with
    `rredlist`, collapses internal nodes with FBP \< 0.70
    (`collapse_cutoff`) into polytomies, and produces the tree,
    chronogram and conservation figures with `ggplot2` and `ggtree`.

4.  **Stage 4: Validation (Module 13)**  
    Compares the focal tree with published phylogenies (Amaral *et al*.
    2022; Thompson *et al*. 2024; Zuntini *et al*. 2024; de Vos *et
    al*. 2025) by Robinson-Foulds and information-theoretic distances
    and by multidimensional scaling (MDS) of tree space
    ([`validate_phylogenies()`](https://beeamerino.github.io/PhyloCactus/reference/validate_phylogenies.md)).
    The published gene trees of de Vos *et al*. (2025) can optionally be
    summarized into a species tree with `ASTRAL-III`, which then enters
    the comparison as a further reference; by default the distributed
    summary tree is used. `PhyloCactus` does not infer gene trees from
    its own supermatrix and does not quantify gene tree discordance.

## Installation

Install the development version from GitHub:

``` r

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

remotes::install_github("beeamerino/PhyloCactus")
```

Check the installed version:

``` r

library(PhyloCactus)

packageVersion("PhyloCactus")
```

## R Dependencies

| Repository | Main packages |
|:---|:---|
| **GitHub** | `phylotaR` |
| **Bioconductor** | `DECIPHER`, `Biostrings`, `ggtree` |
| **CRAN** | `ape`, `phangorn`, `TreeDist`, `ggplot2`, `dplyr`, `tidyr`, `readr`, `stringr`, `purrr`, `rredlist` |

The Bioconductor dependencies are installed with:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DECIPHER", "Biostrings", "ggtree"))
```

## External Software

The following programs must be installed and available on the system
`PATH`, or declared in `~/.Renviron` (`usethis::edit_r_environ()`):

| Software | Purpose | Citation |
|:---|:---|:---|
| [`BLAST+`](https://www.ncbi.nlm.nih.gov/books/NBK279690/) | Sequence similarity search for cluster identification | Camacho *et al*. (2009) *BMC Bioinformatics* |
| [`MAFFT`](https://mafft.cbrc.jp/alignment/software/) | Multiple sequence alignment | Katoh & Standley (2013) *Mol. Biol. Evol.* |
| [`ModelTest-NG`](https://github.com/ddarriba/modeltest) | Substitution model selection for partitioned data | Darriba *et al*. (2020) *Mol. Biol. Evol.* |
| [`RAxML-NG`](https://github.com/amkozlov/raxml-ng) | Constrained maximum-likelihood inference and bootstrap support | Kozlov *et al*. (2019) *Bioinformatics* |
| [`treePL`](https://github.com/blackrim/treePL) | Penalized likelihood divergence time estimation | Smith & O’Meara (2012) *Bioinformatics* |

## Documentation

| Tutorial | Content |
|:---|:---|
| [**Get Started**](https://beeamerino.github.io/PhyloCactus/articles/PhyloCactus.html) | Package overview, design and dependency configuration |
| [**Tutorial 1**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-1-cactus-phylogeny-prep.html) | Stage 1: data assembly, alignment, saturation screen and supermatrix concatenation |
| [**Tutorial 2**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-2-cactus-phylogeny-inference.html) | Stage 2: substitution models, constrained maximum-likelihood search, bootstrap support and `treePL` dating |
| [**Tutorial 3**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-3-cactus-phylogeny-visualization.html) | Stage 3: IUCN Red List data, collapse of weakly supported nodes and figures |
| [**Tutorial 4**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-4-cactus-phylogeny-validation.html) | Stage 4: tree distances, MDS of tree space and the optional `ASTRAL-III` reference summary |
| [**Tutorial 5**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-5-cactus-phylogeny-barcoding.html) | Section manifest: scope, decisions taken, open questions and reopening criteria for the molecular diagnostic section |
| [**Function Reference**](https://beeamerino.github.io/PhyloCactus/articles/tutorial-6-cactus-phylogeny-functions.html) | Package functions and the methods they implement |

## Citation

To cite `PhyloCactus`, use the entry returned by:

``` r

citation("PhyloCactus")
```

## References

- Amaral *et al*. 2022. Spatial patterns of evolutionary diversity in
  Cactaceae show low ecological representation within protected areas.
  *Biological Conservation*, 273, 109677.
  <https://doi.org/10.1016/j.biocon.2022.109677>
- Arakaki *et al*. 2011. Contemporaneous and recent radiations of the
  world’s major succulent plant lineages. *PNAS*, 108(20), 8379-8384.
  <https://doi.org/10.1073/pnas.1100628108>
- de Vos *et al*. 2025. Phylogenomics and classification of Cactaceae
  based on hundreds of nuclear genes. *Plant Systematics and Evolution*,
  311(5), 28. <https://doi.org/10.1007/s00606-025-01948-z>
- Guerrero *et al*. 2019. Phylogenetic relationships and evolutionary
  trends in the cactus family. *Journal of Heredity*, 110(1), 4-21.
  <https://doi.org/10.1093/jhered/esy064>
- Hernández-Hernández *et al*. 2014. Beyond aridification: multiple
  explanations for the elevated diversification of cacti in the New
  World Succulent Biome. *New Phytologist*, 202(4), 1382-1397.
  <https://doi.org/10.1111/nph.12752>
- Jin, L., & Nei, M. 1990. Limitations of the evolutionary parsimony
  method of phylogenetic analysis. *Molecular Biology and Evolution*,
  7(1), 82-102. <https://doi.org/10.1093/oxfordjournals.molbev.a040588>
- Kimura, M. 1980. A simple method for estimating evolutionary rates of
  base substitutions through comparative studies of nucleotide
  sequences. *Journal of Molecular Evolution*, 16(2), 111-120.
  <https://doi.org/10.1007/BF01731581>
- Korotkova *et al*. 2021. Cactaceae at Caryophyllales.org, a dynamic
  online species-level taxonomic backbone for the family. *Willdenowia*,
  51(2), 251-270. <https://doi.org/10.3372/wi.51.51208>
- Thompson *et al*. 2024. Identifying the multiple drivers of cactus
  diversification. *Nature Communications*, 15(1), 7114.
  <https://doi.org/10.1038/s41467-024-51666-2>
- Zuntini *et al*. 2024. Phylogenomics and the rise of the angiosperms.
  *Nature*, 629, 843-850. <https://doi.org/10.1038/s41586-024-07324-0>

## License

`PhyloCactus` is distributed under the GNU General Public License v3.0
(GPL-3).
