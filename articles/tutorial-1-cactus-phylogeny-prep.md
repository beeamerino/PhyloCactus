# Tutorial 1: Phylogenetic Pipeline: Data Assembly and Preparation

## Abstract

Constructing reliable multilocus molecular sequence matrices directly
from public repositories such as GenBank presents a fundamental
computational challenge in systematic biology. Sequence records
deposited across decades frequently exhibit inconsistent locus
annotations, duplicated accessions, unvouchered identifications,
orthographic variants, and massive nomenclatural synonymies.
Transforming these uncurated records into high-quality phylogenetic
matrices requires intensive manual curation, taxonomic reconciliation
against authoritative botanical checklists, and rigorous alignment
quality control.

This computational problem becomes exponentially more complex when
targeting plant lineages with intricate evolutionary histories. Clades
characterized by recent explosive adaptive radiations, low plastid
sequence divergence, incomplete lineage sorting (ILS), and ancient
reticulate evolution amplify the risk of misidentifying orthology,
accumulating systematic alignment noise, and producing biased
topological reconstructions.

The family **Cactaceae** serves as the prime empirical exemplar of these
combined computational and biological hurdles (Guerrero *et al*. 2019).
Comprising one of the largest succulent plant radiations in the
Neotropics, cactus phylogenetics requires extensive data curation to
resolve persistent gene tree discordance and handle heterogeneous
molecular datasets.

`PhyloCactus` is an R package designed to automate the assembly and
curation of multilocus phylogenetic datasets. The workflow builds upon
the orthology-based sequence mining strategy implemented in `phylotaR`,
extending it with comprehensive taxonomic reconciliation against the
Caryophyllales.org checklist
(`CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx`; Korotkova *et al*.
2021), automated quality control procedures using `DECIPHER`,
substitution saturation filtering, multiple sequence alignment using
`MAFFT` to infer positional homology, and concatenation of curated loci
into partitioned supermatrices.

The complete `PhyloCactus` workflow is organized into thirteen
interoperable modules distributed across four analytical stages. This
vignette introduces the first stage of the pipeline (Modules 1 to 6),
guiding users through orthology-based sequence retrieval, taxonomic
standardization, sequence quality assessment, alignment, and matrix
assembly. The resulting curated multilocus datasets constitute the
starting point for all subsequent phylogenetic, biogeographic, and
macroevolutionary analyses performed within the `PhyloCactus` framework.

## Complete Pipeline Execution Workflow

Below is a demonstration of how the `PhyloCactus` package functions
align end-to-end to assemble the initial multilocus sequences.

### Setup: Creating a Clean Workspace

Before beginning the pipeline, it is highly recommended to create a
dedicated folder for your tutorial outputs.

``` r

# Create a dedicated tutorial folder on your Desktop (or any preferred location)
tutorial_dir <- "~/Desktop/PhyloCactus_Tutorial"
dir.create(tutorial_dir, showWarnings = FALSE)

# Set the tutorial folder as your working directory
setwd(tutorial_dir)

# You can also copy the entire tutorial R script to this folder for easy execution
file.copy(
  system.file("scripts", "tutorial-1-cactus-phylogeny-prep.R", package = "PhyloCactus"),
  file.path(tutorial_dir, "tutorial-1-cactus-phylogeny-prep.R")
)
```

### Module 1: Mine Orthologous Sequence Clusters and Retrieve Metadata

The first stage of the workflow assembles the molecular dataset by
retrieving orthologous sequence clusters identified by `phylotaR`.
Rather than relying on heterogeneous gene annotations deposited in
GenBank, `phylotaR` identifies homologous sequences using sequence
similarity, providing a reproducible starting point for multilocus
phylogenetic analyses.

For the default workflow, the ingroup is defined as the family
**Cactaceae** (NCBI Taxonomy ID: **3593**), allowing all descendant taxa
within this lineage to be automatically retrieved from the `phylotaR`
database. Outgroup sampling spans six genera across three families of
the Cactineae suborder within Caryophyllales: **Anacampserotaceae**,
with *Anacampseros* (NCBI Taxonomy ID: **107583**), *Talinopsis*
(**107598**), and *Grahamia* (**107617**); **Portulacaceae**, with
*Portulaca* (**3582**); and **Talinaceae**, with *Talinum* (**107600**)
and *Talinella* (**108056**).

Sampling **Talinaceae** provides the essential outgroup lineage required
to root the tree and evaluate the topological placement of the ACP clade
(Anacampserotaceae, Cactaceae, Portulacaceae). In a matrix containing
only Cactaceae, Anacampserotaceae, and Portulacaceae rooted on
Portulacaceae, the grouping of Cactaceae with Anacampserotaceae is
constrained by the rooting rather than tested by empirical signal.
Including Talinaceae establishes the root outside the core trio, placing
the ACP crown as an internal node and enabling maximum-likelihood
inference to test alternative topological resolutions among the three
families (Ramírez-Barahona *et al.*, 2020; Zuntini *et al.*, 2024).
*Amphipetalum* (**1835425**) is excluded because no nucleotide
accessions were available in GenBank.

Note that *Talinopsis* (**107598**) belongs to Anacampserotaceae and not
to Talinaceae, despite its name and adjacent taxonomic identifier.
Although these taxonomic identifiers constitute the default
configuration distributed with the package, users may specify
alternative ingroup or outgroup NCBI Taxonomy IDs to adapt the workflow
to other evolutionary systems.

To ensure reproducible taxonomic and molecular data curation,
`PhyloCactus` relies on a collection of editable reference files
distributed with the package (`inst/extdata`) and accessed internally
using the [`system.file()`](https://rdrr.io/r/base/system.file.html)
function. Taxonomic standardization is based on a curated taxonomic
backbone (`CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx`) derived
from the **Caryophyllales.org** project (Korotkova *et al.*, 2021). This
checklist contains the currently accepted species of **Cactaceae**,
**Portulacaceae**, and **Anacampserotaceae** and serves as the
authoritative reference for reconciling species names throughout the
pipeline. By relying on an external taxonomic backbone rather than
GenBank nomenclature alone, `PhyloCactus` minimizes inconsistencies
arising from outdated names, spelling variants, and unresolved
synonymies.

Locus selection is controlled by `target_genes.txt`, which defines the
set of molecular markers targeted during sequence retrieval. The default
file contains a curated list of loci widely used in **Cactaceae**
phylogenetics, although users may freely modify or replace this list to
accommodate alternative taxonomic groups or marker sets. Gene
nomenclature is standardized using `genes_map.csv`, a curated dictionary
of synonymous gene names commonly encountered in GenBank. During data
processing, heterogeneous locus annotations are harmonized into
standardized marker names, ensuring that equivalent loci are
consistently recognized despite differences in nomenclature among
independent sequencing projects. Together, these editable configuration
files separate biological customization from analytical code, allowing
the same analytical workflow to be readily adapted to different
evolutionary systems without modifying the package source code.

The exclusion lists deserve a note of their own, because they encode a
decision rather than a threshold. `manual_exclusions_ingroup.csv` and
`manual_exclusions_outgroup.csv` remove 73 unique accessions across 87
records, spanning 16 ingroup clusters and 24 outgroup clusters. These
are the product of manual curation by the expert team supporting
`PhyloCactus`: each accession was inspected and removed on taxonomic or
sequence-quality grounds that are not recoverable from GenBank metadata
alone. The `apply_manual_exclusions` argument makes that decision
visible and reversible. It defaults to `TRUE`, and is stated explicitly
in the call below rather than inherited, so that a reader sees the
curation exists. Setting it to `FALSE` reproduces the uncurated cluster
set, which is the way to quantify what the curation actually removes.
Either way the applied list is written to
`TABLE_MANUAL_EXCLUSIONS_*.csv` in the output directory, so every run
documents its own curation state.

The
[`assemble_ingroup_phylotar()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_ingroup_phylotar.md)
function retrieves all orthologous sequence clusters associated with the
focal taxonomic ingroup from the `phylotaR` database. Cluster metadata
are parsed, taxonomic names are reconciled against the curated taxonomic
backbone, gene names are standardized using `genes_map.csv`, and only
loci matching the predefined marker list in `target_genes.txt` are
retained. The resulting orthologous sequence clusters are exported as
individual FASTA files within the `1_phylotaR_out_ingroup` directory,
together with an occupancy table summarizing taxonomic representation,
cluster characteristics, and associated metadata.

The
[`assemble_outgroup_phylotar()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_outgroup_phylotar.md)
function performs the same procedure for the selected outgroup taxa. To
maximize comparability across the final multilocus dataset, only
orthologous loci corresponding to the standardized marker set recovered
for the ingroup are retained. The resulting FASTA files and occupancy
metadata are written to the `1_phylotaR_out_outgroup` directory,
producing a harmonized collection of orthologous loci that serves as the
input for downstream analyses.

``` r

library(PhyloCactus)

# Define and create working directories for `phylotaR`
dir.create("0_phylotaR_raw_Ingroup", showWarnings = FALSE)
dir.create("0_phylotaR_raw_Outgroup", showWarnings = FALSE)

# Mine ingroup taxonomic database
ingroup_assembly <- assemble_ingroup_phylotar(
  wd_path = "0_phylotaR_raw_Ingroup",
  target_genes_file = system.file("extdata", "target_genes.txt", package = "PhyloCactus"),
  genes_map_file = system.file("extdata", "genes_map.csv", package = "PhyloCactus"),
  manual_exclusions_file = system.file("extdata", "manual_exclusions_ingroup.csv", package = "PhyloCactus"),
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
  # Talinopsis (107598) is Anacampserotaceae, not Talinaceae, despite the name and adjacent ID.
  outgroups = c("107598", "107617", "107583", "3582", "107600", "108056"),
  force_download = FALSE
)

# Load the accepted species checklist
checklist_path <- system.file("extdata", "CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx", package = "PhyloCactus")
if (file.exists(checklist_path)) {
  sheets <- readxl::excel_sheets(checklist_path)
  checklist_sheets <- sheets[!grepl("^facts", sheets, ignore.case = TRUE)]
  list_df <- lapply(checklist_sheets, function(sh) {
    readxl::read_excel(checklist_path, sheet = sh)
  })
  cactaceae_checklist <- dplyr::bind_rows(list_df)
  message(paste("Loaded", nrow(cactaceae_checklist), "taxa from the accepted checklist. \U0001f335"))
}
```

### Module 2: Align Sequences and Mask Low Confidence Regions

Multiple sequence alignment represents the fundamental hypothesis of
positional homology upon which all subsequent phylogenetic analyses
depend. Errors introduced during sequence alignment can propagate
throughout the analytical pipeline, biasing branch length estimation,
reducing nodal support, and potentially leading to incorrect
phylogenetic inference. Consequently, objective alignment and quality
assessment are essential before constructing concatenated multilocus
datasets.

The
[`run_alignment_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_alignment_pipeline.md)
function processes the orthologous sequence clusters generated in Module
1 to produce high quality multiple sequence alignments (MSAs). Primary
alignments are inferred using **MAFFT** (Katoh & Standley, 2013), a fast
and accurate multiple sequence alignment algorithm widely adopted in
molecular phylogenetics. Following alignment, `PhyloCactus` applies the
**DECIPHER** framework (Wright, 2024) to identify and mask poorly
aligned regions, ambiguous nucleotide positions, long insertions or
deletions, and other alignment segments with low confidence that are
unlikely to represent reliable positional homology. This automated
masking procedure minimizes the influence of alignment uncertainty while
preserving informative phylogenetic signal.

Finally, whenever masking is enabled the pipeline filters sequences and
nucleotide sites according to user defined occupancy thresholds,
removing loci or taxa with excessive missing data and retaining only
well supported alignment columns suitable for downstream analyses. The
`min_masked_alignment_length` argument sets an absolute floor on the
number of columns that must survive masking: because
`min_non_gap_fraction` is evaluated relative to the post-masking width,
a locus collapsed to a handful of columns would otherwise pass the
filter and be exported as a near-empty alignment.

Two of the arguments above encode analytical decisions rather than mere
thresholds, and both are stated explicitly in the call instead of being
left to the default.

`fix_strand = TRUE` puts every sequence of a marker on the same strand
before `MAFFT` sees it. GenBank stores each record on whichever strand
the submitter deposited, and `MAFFT` compares only the orientation it is
given: a reverse-complemented accession is aligned anyway, and the
resulting row carries no positional homology to the rest of the block.
Nothing downstream detects this. The row simply looks like a very
divergent sequence, or fails the occupancy filter of Module 5 and
disappears without explanation.

The 2026-09-01 audit found it in the `rbcL` and `matK` accessions of
*Portulaca oleracea* and *P. pilosa*, both deposited reversed. Their
aligned `rbcL` sat at 0.51 observed divergence from Cactaceae, and their
`matK` at 0.40, where *P. grandiflora* of the same genus sits at 0.030
and 0.066. *P. oleracea* was at that moment the terminal with the most
retained markers in the outgroup, and therefore the one the acceptance
check at the end of this tutorial nominated for rooting the tree. The
ingroup was not clean either, and the first measurement said it was.
Counted on the curated output of Module 4 it read 0 reversed of 3207;
counted on the raw input it is 163, of which 77 in `psbA-trnH` and 86 in
`rpL16`. The difference is the point: a reversed sequence aligns to
nothing, ends with almost no occupancy, and is discarded by the
occupancy filter of Module 5, so by the time the data reach a curated
directory the evidence of the defect has already been removed along with
the sequences. Correcting the strand recovered them, taking `psbA_trnH`
from 296 to 368 ingroup sequences and `rpL16` from 600 to 677.

Orientation is decided against each marker’s own majority rather than an
external reference, because the pipeline mines whatever GenBank holds
for a locus and no reference is guaranteed to exist. Every sequence is
logged with its k-mer match in both directions in
`tables/LOG_STRAND_<marker>.csv`, and a sequence matching the marker in
neither direction is left untouched and reported separately: that is a
homology problem rather than an orientation one, and flipping it would
only hide it. `mafft --adjustdirection` solves the same problem and is
not used here because it renames the sequences it flips with an `_R_`
prefix, which would then have to be undone in the sequence filter log,
the name crosswalk and the concatenation.

One case this cannot see: the ingroup and the outgroup are aligned in
separate runs, so each is oriented within its own pool and neither can
tell that the two pools disagree with each other. That comparison
belongs to Module 4, where the two sets are joined, and `homology_check`
reports it there as `reversed_relative_to_ingroup`.

`preserve_iupac = TRUE` retains the IUPAC ambiguity codes (`R`, `Y`,
`S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`) instead of collapsing them to
`N`. The tools downstream all accept them, and they are not equivalent
to missing data: `RAxML-NG` and `ModelTest-NG` treat an ambiguity code
as a partial constraint on the state, so an `R` site restricts the
possibilities to A or G whereas an `N` restricts nothing. Collapsing the
codes therefore discards a real constraint for no analytical gain, and
it also erases genuine heterozygous signal in multicopy nuclear markers
such as `ITS`. The alternative remains available for anyone who needs an
alignment free of ambiguity, but the pipeline does not impose it.

Ambiguity is measured rather than assumed away. The manifest reports
`mean_fraction_ambiguous_*` for each processing stage and
`n_sites_ambiguous_*` for the raw input and the final alignment,
alongside the existing gap and missing-data columns. The `raw_input`
figures are computed before any cleaning is applied, so they record what
the source records actually contained irrespective of the policy in
force. That is what makes a per-locus judgement possible: a marker whose
ambiguity is concentrated can be examined on its own evidence rather
than subjected to a global rule.

`min_masked_alignment_length = 100L` sets the floor discussed above. Its
purpose is to catch alignments that masking has degraded to the point of
being uninformative, not to arbitrate between loci of different lengths:
no marker in this dataset is shorter than 100 columns, so a masked
alignment below that value is degenerate rather than merely short.
Raising the floor can only reject markers, never admit them, and any
marker it rejects is reported with a `decision_reason` naming the
threshold, so its effect is always visible in the screening table.

Note that `min_non_gap_fraction`, `max_missing_fraction` and
`min_masked_alignment_length` are all evaluated against the post-masking
column set, so when `mask_alignment_regions = FALSE` none of them is
applied at this stage. That is intentional: with masking deferred,
occupancy filtering happens exactly once, on the joint ingroup plus
outgroup alignment of Module 5. The chosen policy and the value of every
threshold are recorded in `2_MAFFT_*/logs/LOG_alignment_run_info.txt`
and stamped into each `TABLE_marker_alignment_summary_<marker>.csv`, so
a given output directory always documents the call that produced it.
Those same stamped values drive the caching behaviour: a marker is
reused from a previous run only if all five parameters match, so
changing `mask_alignment_regions` on an already populated output
directory forces reprocessing rather than silently returning the earlier
alignments.

Masking is applied to the ingroup, where each locus is represented by
dozens to hundreds of sequences and
[`DECIPHER::MaskAlignment`](https://rdrr.io/pkg/DECIPHER/man/MaskAlignment.html)
has enough signal to behave reliably. It is **deliberately disabled for
the outgroup** (`mask_alignment_regions = FALSE`). With only two to
seven highly divergent accessions per locus, masking the outgroup on its
own defines a column set that the ingroup does not share and can erode
an outgroup alignment to a few base pairs; those fragments then fail the
occupancy filter during the joint realignment of Module 5 and vanish
from the supermatrix, stripping the outgroup precisely from the loci
needed to root the tree. Deferring masking to Module 5, where ingroup
and outgroup are masked together against a single column set, preserves
outgroup coverage without weakening quality control.

``` r

# Align ingroup sequences and mask alignments using the automated batch pipeline
ingroup_alignment_manifest <- run_alignment_pipeline(
  input_folder = "1_phylotaR_out_Ingroup",
  output_dir = "2_MAFFT_Cactaceae",
  mask_alignment_regions = TRUE,
  min_non_gap_fraction = 0.30,
  max_missing_fraction = 0.30,
  min_masked_alignment_length = 100L,
  preserve_iupac = TRUE,
  fix_strand = TRUE
)

# Align outgroup sequences, deferring masking to the joint realignment of Module 5.
# The occupancy thresholds are omitted deliberately: they are not applied when masking
# is disabled, so passing them here would only suggest a filtering that does not happen.
outgroup_alignment_manifest <- run_alignment_pipeline(
  input_folder = "1_phylotaR_out_Outgroup",
  output_dir = "2_MAFFT_Outgroup",
  mask_alignment_regions = FALSE,
  preserve_iupac = TRUE,
  fix_strand = TRUE
)

# Confirm the masking policy that was actually applied. The outgroup rows must read
# masking_applied = FALSE and status = "OK_NO_MASK". If they read TRUE and "OK", the
# argument did not take effect and every downstream table is the previous dataset.
print(outgroup_alignment_manifest[, c("marker", "masking_applied", "status")])
stopifnot(all(!outgroup_alignment_manifest$masking_applied))
```

### Module 3: Evaluate Phylogenetic Signal and Screen Ingroup Markers

Although orthologous loci can be accurately aligned, not all markers
contribute equally to phylogenetic inference. Differences in
evolutionary rate, substitution saturation, sequence completeness, and
structural anomalies may reduce the phylogenetic information contained
within individual loci or introduce systematic bias into concatenated
analyses. Consequently, each marker should be evaluated before inclusion
in the final multilocus dataset.

The
[`run_marker_screening()`](https://beeamerino.github.io/PhyloCactus/reference/run_marker_screening.md)
function performs an automated quality assessment of the ingroup
alignments generated in Module 2. Passing `outgroup_folder` records
`n_outgroup` and `n_total` in the summary table for diagnostic
visibility, though retention decisions remain strictly based on ingroup
thresholds to prevent distant outgroup divergence from distorting
saturation regressions. For example, `trnT-psbD` possessed 50 ingroup
sequences and was filtered out by `min_nseq_to_retain = 100`, despite
carrying extensive outgroup representation. First, substitution
saturation is evaluated using regression-based statistics to identify
loci in which multiple substitutions may have eroded the underlying
phylogenetic signal. The function then detects sequence length outliers
using the interquartile range (IQR), identifying sequences that may
represent incomplete assemblies, sequencing artifacts, or annotation
errors.

Based on user defined thresholds, loci and sequences that fail these
quality criteria are excluded from subsequent analyses. The function
produces a curated collection of phylogenetically informative markers
together with a comprehensive diagnostic table summarizing saturation
statistics, sequence length distributions, outlier detection, and
filtering decisions for each locus. These curated markers constitute the
final set of ingroup loci that will be combined with the corresponding
outgroup sequences for downstream concatenation and phylogenetic
inference.

``` r

# Sub-select markers with stable saturation regression slopes (slope > 0.5)
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
```

### Module 4: Integrate and Decouple Ingroup and Outgroup Markers

The previous modules independently retrieve, align, and evaluate ingroup
and outgroup sequences. Before downstream concatenation and evolutionary
modeling, these datasets must be integrated into a standardized
taxonomic framework to reconcile species binomials against the
**Caryophyllales.org** taxonomic backbone
(`CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx`; Korotkova *et al*.
2021). Taxonomic synonyms, obsolete names, infraspecific designations,
and orthographic variants are resolved to their accepted species names
across **Cactaceae** and **Anacampserotaceae**, while root outgroups
(**Portulacaceae** and **Talinaceae**) are preserved for phylogenetic
rooting.

To prevent outgroup divergence from inflating molecular informativeness
metrics (such as parsimony informative sites and alignment lengths) of
the focal radiation, the
[`integrate_and_clean_markers()`](https://beeamerino.github.io/PhyloCactus/reference/integrate_and_clean_markers.md)
function organizes all curated outputs into structured subdirectories:

1.  `4_Cleaned/cleaned_markers_ingroup`: Contains curated sequences
    exclusively for accepted **Cactaceae** taxa, preserved for
    single-locus analyses. These are the inputs the molecular diagnostic
    section will consume once its strategy is defined; see the section
    manifest in Tutorial 5.
2.  `4_Cleaned/cleaned_markers_outgroup`: Contains corresponding
    outgroup sequences (*Portulaca*, *Anacampseros*, *Talinopsis*,
    *Grahamia*, *Talinum*, *Talinella*).
3.  `4_Cleaned/cleaned_markers_joint`: Contains the integrated ingroup
    and outgroup sequences prepared for joint realignment in Module 5.
4.  `4_Cleaned/tables`: Contains the complete taxonomic synthesis and
    curation table suite:
    - `TABLE_taxonomic_backbone_molecular_matrix.csv`: Master taxonomic
      matrix anchored to all accepted species in the checklist
      (**Cactaceae** and **Anacampserotaceae**) plus outgroups. Features
      dynamic GenBank accession columns (`sid_<marker>`),
      Caryophyllales.org UUIDs, and standardized placeholder columns for
      phylogenetic inclusion (`retained_in_phylogeny`) and conservation
      status (`iucn_redlist_category`).
    - `TABLE_dataset_species_summary.csv`: Global dataset-level
      synthesis reporting total species in the checklist, raw GenBank
      species recovered, accepted vs. rejected species counts, outgroups
      retained, total unique taxa in the joint dataset, and overall
      checklist recovery percentage, plus the number of standardized
      loci mined per taxonomic group (**Cactaceae**,
      **Anacampserotaceae**, **Portulacaceae**, **Talinaceae**) at
      Stage 1. Every row that depends on what survives the joint
      realignment, both the “Final loci retained” counts and the
      outgroup and joint **species** counts, is written here as a
      provisional (pre-realignment) value and overwritten in place by
      [`run_concatenation_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_concatenation_pipeline.md)
      (Module 6) once the true post-realignment supermatrix exists. This
      single file therefore holds the definitive counts after Module 6
      runs, and any row still labelled “Provisional” in its `details`
      column indicates that Module 6 has not been executed against it.
      Do not quote figures from this table before Module 6 has run.
    - `TABLE_species_curation_status.csv`: Taxon-level master registry
      detailing acceptance status, source classification
      (`accepted_ingroup`, `rejected_ingroup`, `accepted_outgroup`), and
      retained multi-marker occupancy for every species retrieved.
    - `TABLE_marker_metrics_ingroup.csv`,
      `TABLE_marker_metrics_outgroup.csv`,
      `TABLE_marker_metrics_joint.csv`: Decoupled molecular alignment
      metrics (length, gaps, GC content, variable sites, parsimony
      informative sites) calculated independently for the accepted
      ingroup, outgroups, and joint datasets.
    - `TABLE_marker_summary.csv`: Per-marker comparative summary
      calculating the parsimony informative site inflation ratio
      (`pis_inflation_ratio = pis_joint / pis_ingroup`) and taxon
      retention percentages.
    - `TABLE_species_marker_sid_matrix.csv`,
      `TABLE_sequence_registry_with_acceptance.csv`,
      `TABLE_duplicate_resolution_species_marker.csv`: Full
      sequence-level traceability and deduplication audit registries.
5.  `4_Cleaned/logs`: Contains the comprehensive integration and
    curation audit log (`LOG_clean_integration_summary.txt`).

``` r

# Integrate datasets, resolve botanical synonymies, and export decoupled directories
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
  # than ingroup resolution. Left empty deliberately; see the note below.
  readmit_markers = NULL,
  readmit_dir = NULL,
  # Alignment-free homology check, ingroup against outgroup, per marker.
  homology_check = TRUE
)
```

`trnT-psbD` was readmitted here on 2026-09-01 and removed the same day,
and the reason is worth stating because the counts argued the other way.
It had the most balanced coverage in the dataset, 49 ingroup terminals
against 49 outgroup at 0.998 median occupancy, where the next best locus
offered six outgroup terminals. The sequences are not the same region:
0.039 of the outgroup 20-mers occur anywhere in the ingroup sequences,
against 0.28 to 0.60 for every genuine counterpart, and observed
divergence against Cactaceae is 0.440 where `matK` gives 0.077 and
`rbcL` 0.033. The ingroup median length is 598 bp and the outgroup 1347.
It was contributing 521 columns of non-homologous characters to the
terminals that define the root, which is the branch every calibration is
estimated across.

`homology_check` exists because counts cannot see this and neither can a
coverage table. It measures, for each marker present on both sides, the
fraction of outgroup k-mers occurring in the ingroup sequences of the
same name, writes `4_Cleaned/tables/TABLE_marker_homology_check.csv`,
and warns for the markers below both thresholds. The measure is
alignment-free, so a low value cannot be fixed by aligning better: it
separates sequences that are hard to align from sequences that are not
the same region. It reports and does not block, because a low share on a
fast-evolving locus is a reason to look rather than a verdict.

### Module 5: Joint Realignment of Curated Markers

Although ingroup and outgroup markers were aligned independently during
earlier stages, integrating both datasets introduces homologous sequence
variation across divergent lineages. Independent alignment can introduce
insertion/deletion (indel) boundary shifts at the junction between
divergent lineages, potentially biasing downstream branch length
estimation if concatenated directly.

The
[`run_joint_realignment()`](https://beeamerino.github.io/PhyloCactus/reference/run_joint_realignment.md)
function executes multiple sequence alignment using `MAFFT` on the joint
marker FASTA files, followed by alignment quality masking with
`DECIPHER`
([`DECIPHER::MaskAlignment`](https://rdrr.io/pkg/DECIPHER/man/MaskAlignment.html)).
This establishes positional homology across the complete taxonomic
sampling before supermatrix assembly, and it is the single point at
which masking is applied to the outgroup.

This module also applies a per-sequence occupancy filter
(`min_non_gap_fraction`, `max_missing_fraction`), and it is the step
that can remove individual taxa from a locus without removing the locus
itself. Sequences that arrive as short fragments fail the filter and
disappear from the supermatrix; the affected terminals are recorded in
`5_MAFFT_Cleaned/LOG_SEQ_FILTER_<marker>.csv` with `Retained = FALSE`.

The threshold is a fraction of the alignment width, and the alignment
width is set by the longest sequences, which are the ingroup ones. An
outgroup accession covering a shorter amplicon of the same region
therefore fails on length alone. That is how the five *Portulaca*
sequences of `trnL_trnF`, 297 bp against a threshold of about 353, left
the matrix: the branch subtending the outgroup fell from roughly 600
expected substitutions to 0.03, and the two deepest calibrated nodes
stopped being estimated at all.

Three arguments address this. `rooting_pattern` exempts nothing and
names, in a warning, any rooting terminal the filter removes, so the
loss is visible where it happens rather than four modules downstream.
`protect_pattern` retains matching terminals regardless of occupancy.
`protect_markers` restricts that exemption to named loci, which matters
because a matrix-wide exemption retains every short fragment of every
protected terminal in every locus, raising the gap fraction of the
supermatrix and destabilising the affected terminals during inference.
Scoping the exemption to the loci that actually lose rooting terminals
buys the separation without that cost.

``` r

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
```

### Module 6: Construct the Multilocus Supermatrix

After individual locus alignments have been curated, realigned, and
evaluated, they are concatenated into a unified multilocus supermatrix
for partitioned maximum-likelihood phylogenetic inference.

Passing `outgroup_pattern` runs \[report_marker_group_coverage()\]
first, which counts how many ingroup and how many outgroup terminals
carry real sequence in each locus and warns when either side is empty. A
locus sampled almost entirely on one side of the root contributes
columns the other side cannot share, and the branch lengths spanning
that bipartition are then estimated from the loci that remain. The
report does not block: whether such a locus belongs in a matrix depends
on the analysis the matrix is for, and a nuclear locus with no outgroup
coverage is unusable for dating and valuable for species discrimination.

The
[`run_concatenation_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_concatenation_pipeline.md)
function combines the aligned markers, matches taxa across loci, inserts
missing data characters (`-`) where accessions are absent, and
calculates the exact nucleotide coordinates defining every gene
partition. The module generates concatenated alignments in FASTA, PHYLIP
(`.phy`), and NEXUS (`.nex`) formats, exports partition schemes
formatted for `RAxML-NG` and `IQ-TREE`, and writes
`TABLE_partition_informativeness_comparison.csv` comparing ingroup
versus joint parsimony informative sites across all partitions. It also
updates the “Final loci retained” rows of
`4_Cleaned/tables/TABLE_dataset_species_summary.csv` (written
provisionally by Module 4) in place, replacing them with the true
post-realignment counts now that Module 5’s joint realignment may have
dropped individual taxa from a locus.

`exclude_markers` names the loci to leave out of this particular
supermatrix. The alignments are neither modified nor removed, so a locus
excluded here remains available to every other analysis; only this
matrix is built without it. The excluded names are echoed to the run log
and written to `logs_and_qc/TABLE_markers_excluded.csv`, and a name
matching no alignment raises a warning rather than passing silently, so
a typo cannot leave a locus in the matrix while the script reads as
though it had been dropped.

Set `exclude_markers = NULL` to build the full matrix. `phyC` is
excluded below for the reason the coverage report gives above: it
carries 167 ingroup terminals and 15 outgroup terminals, all 15 in
Anacampserotaceae and none in Portulacaceae, so in a matrix built to
date the divergence between those families it supplies columns that only
one side of the root can occupy. That is a property of this matrix, not
of the locus, and `phyC` remains the more discriminating of the two
nuclear markers for barcoding.

``` r

run_concatenation_pipeline(
  input_dir = "5_MAFFT_Cleaned/aligned_markers",
  output_dir = "6_Concatenated",
  outgroup_pattern = "^(Anacampseros|Grahamia|Talinopsis|Portulaca|Talinum|Talinella)_",
  exclude_markers = "phyC"   # NULL to retain all loci
)
```

## Conclusion

At this stage, the molecular dataset has been fully assembled, curated,
and prepared for phylogenetic analysis. Orthologous loci have been
retrieved, taxonomic names standardized, multiple sequence alignments
refined, phylogenetically informative markers selected, and the final
multilocus supermatrix constructed together with its corresponding
partition scheme. These outputs constitute the complete analytical
dataset required for evolutionary inference.

The next stage of the `PhyloCactus` workflow focuses on reconstructing
evolutionary relationships and estimating the temporal framework of
diversification. In the following tutorial, we will infer a maximum
likelihood phylogeny using `RAxML-NG`, evaluate branch support, and
estimate divergence times under a penalized likelihood framework using
`treePL`. The resulting time calibrated phylogeny will serve as the
foundation for historical biogeographic and macroevolutionary analyses.

[Continue to Tutorial 2: Inference and Divergence Time
Estimation](https://beeamerino.github.io/PhyloCactus/articles/tutorial-2-cactus-phylogeny-inference.html)

## References

- Guerrero *et al*. 2019. Phylogenetic Relationships and Evolutionary
  Trends in the Cactus Family. *Journal of Heredity*, 110(1), 4–21.
  <https://doi.org/10.1093/jhered/esy064>
- Katoh, K., & Standley, D. M. 2013. MAFFT multiple sequence alignment
  software version 7: Improvements in performance and usability.
  *Molecular Biology and Evolution*, 30(4), 772–780.
  <https://doi.org/10.1093/molbev/mst010>
- Korotkova *et al*. 2021. Cactaceae at Caryophyllales.org - A dynamic
  online species-level taxonomic backbone for the family. *Willdenowia*,
  51(2), 251–270. <https://doi.org/10.3372/wi.51.51208>
- Wright, E. 2024. Fast and Flexible Search for Homologous Biological
  Sequences with DECIPHER v3. *The R Journal*, 16(2), 191-200.
  <https://doi.org/10.18129/B9.bioc.DECIPHER>
