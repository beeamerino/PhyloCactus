# Tutorial 1: Data Assembly and Preparation

## Abstract

Multilocus matrices assembled from GenBank inherit the heterogeneity of
decades of deposition: inconsistent locus annotations, duplicated
accessions, records without vouchers, orthographic variants and
nomenclatural synonyms. In plant lineages that diversified recently and
rapidly, such as **Cactaceae** (Guerrero *et al*. 2019), low plastid
divergence, incomplete lineage sorting (ILS) and reticulate evolution
add to these problems the risk of erroneous orthology assessment and of
alignment error.

This tutorial covers the first stage of the `PhyloCactus` workflow
(Modules 1 to 6): orthology-based sequence retrieval with `phylotaR`,
taxonomic reconciliation against the Caryophyllales.org checklist
(`CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx`; Korotkova *et al*.
2021), alignment with `MAFFT`, masking with `DECIPHER`, a saturation
screen, and concatenation of the curated loci into a partitioned
supermatrix. The supermatrix and its partition scheme are the input of
Tutorial 2.

## Complete Pipeline Execution Workflow

The code below runs Modules 1 to 6 in order.

### Setup: Creating a Clean Workspace

Create a dedicated folder for the outputs of the tutorial.

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

The first module retrieves orthologous sequence clusters identified by
`phylotaR`, which groups sequences by similarity and does not depend on
the gene annotations of GenBank records.

For the default workflow, the ingroup is defined as the family
**Cactaceae** (NCBI Taxonomy ID: **3593**), allowing all descendant taxa
within this lineage to be automatically retrieved from the `phylotaR`
database. Outgroup sampling spans six genera across three families of
the Cactineae suborder within Caryophyllales: **Anacampserotaceae**,
with *Anacampseros* (NCBI Taxonomy ID: **107583**), *Talinopsis*
(**107598**), and *Grahamia* (**107617**); **Portulacaceae**, with
*Portulaca* (**3582**); and **Talinaceae**, with *Talinum* (**107600**)
and *Talinella* (**108056**).

**Talinaceae** are sampled to root the tree outside the ACP clade
(Anacampserotaceae, Cactaceae, Portulacaceae). In a matrix containing
only Cactaceae, Anacampserotaceae and Portulacaceae and rooted on
Portulacaceae, the grouping of Cactaceae with Anacampserotaceae is
imposed by the rooting and not tested by the data. With Talinaceae
sampled, the root lies outside the three families, the crown of the ACP
clade is an internal node, and maximum-likelihood inference can compare
alternative resolutions among the three families (Ramírez-Barahona *et
al.*, 2020; Zuntini *et al.*, 2024). *Amphipetalum* (**1835425**) is not
sampled because GenBank holds no nucleotide records for it.

*Talinopsis* (**107598**) belongs to Anacampserotaceae, not to
Talinaceae, despite its name and adjacent taxonomic identifier. These
identifiers are the default configuration; other ingroup or outgroup
NCBI Taxonomy IDs can be supplied to apply the workflow to other groups.

Taxonomic and molecular curation relies on reference files distributed
with the package in `inst/extdata` and read with
[`system.file()`](https://rdrr.io/r/base/system.file.html). Species
names are reconciled against a taxonomic backbone
(`CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx`) derived from the
Caryophyllales.org project (Korotkova *et al.*, 2021), which lists the
accepted species of Cactaceae, Portulacaceae and Anacampserotaceae.
Reconciling against this checklist, and not relying on the names in
GenBank records, reduces inconsistencies from outdated names, spelling
variants and unresolved synonymy.

Locus selection is controlled by `target_genes.txt`, which lists the
markers targeted during retrieval; the default list contains loci
commonly used in Cactaceae phylogenetics and can be edited or replaced.
Gene names are standardized with `genes_map.csv`, a dictionary of
synonymous gene names found in GenBank, so that the same locus is
recognized under the different names used by independent studies. Both
files are separate from the package code, so the workflow can be adapted
to other groups or marker sets without modifying the package.

The exclusion lists record a curatorial decision.
`manual_exclusions_ingroup.csv` and `manual_exclusions_outgroup.csv`
remove 99 unique accessions across 113 records, in 16 ingroup clusters
and 24 outgroup clusters. Each accession was inspected and removed by
the team supporting `PhyloCactus` on taxonomic or sequence-quality
grounds that cannot be recovered from GenBank metadata. The
`apply_manual_exclusions` argument defaults to `TRUE` and is stated
explicitly in the call below so that the curation is visible; setting it
to `FALSE` reproduces the uncurated cluster set and allows the effect of
the curation to be measured. The applied list is written to
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
function applies the same procedure to the outgroup taxa and retains
only the loci of the standardized marker set recovered for the ingroup.
The FASTA files and occupancy metadata are written to
`1_phylotaR_out_outgroup` and are the input of the next module.

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

A multiple sequence alignment is the hypothesis of positional homology
on which all later analyses depend. Alignment errors propagate to branch
length estimation and nodal support and can lead to incorrect
topologies, so each alignment is checked before the loci are
concatenated.

The
[`run_alignment_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_alignment_pipeline.md)
function aligns the clusters produced in Module 1 with `MAFFT` (Katoh &
Standley, 2013) and masks poorly aligned regions, ambiguous positions,
long insertions or deletions and other segments of doubtful positional
homology with `DECIPHER` (Wright, 2024).

When masking is enabled, the pipeline then filters sequences and sites
by occupancy, removing sequences with excessive missing data and poorly
occupied columns. The `min_masked_alignment_length` argument sets an
absolute floor on the number of columns that must remain after masking:
because `min_non_gap_fraction` is evaluated relative to the post-masking
width, a locus reduced to a handful of columns would otherwise pass the
filter and be exported as a near-empty alignment.

Two arguments of the call below set analytical choices, and both are
stated explicitly.

`fix_strand = TRUE` puts every sequence of a marker on the same strand
before `MAFFT` sees it. GenBank stores each record on whichever strand
the submitter deposited, and `MAFFT` compares only the orientation it is
given: a reverse-complemented accession is aligned anyway, and the
resulting row carries no positional homology to the rest of the block.
Nothing downstream detects this. The row simply looks like a very
divergent sequence, or fails the occupancy filter of Module 5 and
disappears without explanation.

In the reference dataset, the `rbcL` and `matK` accessions of *Portulaca
oleracea* and *P. pilosa* are deposited reversed, and so are 163 ingroup
sequences, 77 of `psbA-trnH` and 86 of `rpL16`. Uncorrected, these
sequences align to nothing, lose occupancy and are removed by the
occupancy filter of Module 5, so no reversed sequence remains visible in
the curated output. Correcting the strand retains them: `psbA_trnH`
keeps 368 ingroup sequences in place of 296, and `rpL16` 677 in place of
600.

Orientation is decided against the majority orientation of each marker,
because no external reference is guaranteed to exist for a mined locus.
Every sequence is logged with its k-mer match in both directions in
`tables/LOG_STRAND_<marker>.csv`. A sequence matching the marker in
neither direction is left unchanged and reported separately: that is a
homology problem, which reversing the sequence would conceal.
`mafft --adjustdirection` addresses the same problem but renames the
sequences it reverses with an `_R_` prefix, which would have to be
undone in the sequence filter log, the name crosswalk and the
concatenation.

The ingroup and the outgroup are aligned in separate runs, so each is
oriented within its own pool and a disagreement between the two pools is
not detected here. That comparison is made in Module 4, where the two
sets are joined, and `homology_check` reports it as
`reversed_relative_to_ingroup`.

`preserve_iupac = TRUE` retains the IUPAC ambiguity codes (`R`, `Y`,
`S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`) instead of collapsing them to
`N`. `RAxML-NG` and `ModelTest-NG` treat an ambiguity code as a partial
constraint on the state: an `R` site restricts the state to A or G,
whereas an `N` does not restrict it. Collapsing the codes discards that
information and erases heterozygous signal in multicopy nuclear markers
such as `ITS`. Setting `preserve_iupac = FALSE` collapses them for
analyses that require alignments without ambiguity codes.

The amount of ambiguity is reported. The manifest gives
`mean_fraction_ambiguous_*` for each processing stage and
`n_sites_ambiguous_*` for the raw input and the final alignment,
alongside the gap and missing-data columns. The `raw_input` values are
computed before any cleaning, so they record what the source records
contain under either setting, and a locus with concentrated ambiguity
can be examined individually.

`min_masked_alignment_length = 100L` sets the floor described above. It
detects alignments that masking has reduced to an uninformative remnant;
it does not select between loci of different lengths, because no marker
in this dataset is shorter than 100 columns. Raising the floor can only
reject markers, and a rejected marker is reported with a
`decision_reason` naming the threshold.

`min_non_gap_fraction`, `max_missing_fraction` and
`min_masked_alignment_length` are evaluated on the post-masking column
set, so none of them is applied at this stage when
`mask_alignment_regions = FALSE`. With masking deferred, occupancy
filtering takes place once, on the joint ingroup and outgroup alignment
of Module 5. The masking policy and the value of every threshold are
recorded in `2_MAFFT_*/logs/LOG_alignment_run_info.txt` and in each
`TABLE_marker_alignment_summary_<marker>.csv`. The same values control
caching: a marker is reused from a previous run only if all five
parameters match, so changing `mask_alignment_regions` in a populated
output directory forces reprocessing.

Masking is applied to the ingroup, where each locus is represented by
dozens to hundreds of sequences and
[`DECIPHER::MaskAlignment`](https://rdrr.io/pkg/DECIPHER/man/MaskAlignment.html)
has enough signal to behave reliably. It is disabled for the outgroup
(`mask_alignment_regions = FALSE`). With two to seven divergent
accessions per locus, masking the outgroup separately defines a column
set that the ingroup does not share and can reduce an outgroup alignment
to a few base pairs; those fragments then fail the occupancy filter in
the joint realignment of Module 5 and are lost from the supermatrix,
together with the characters needed to root the tree. Masking is
therefore deferred to Module 5, where ingroup and outgroup are masked
together against a single column set.

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
# The occupancy thresholds are omitted: they are not applied when masking
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

Loci and sequences that fail the thresholds are excluded from the
following modules. The function writes the retained markers and a
diagnostic table with the saturation statistics, sequence length
distributions, outliers and filtering decision for each locus. The
retained markers are combined with the corresponding outgroup sequences
in Module 4.

``` r

# Sub-select markers with stable saturation regression slopes (slope > 0.5)
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
```

### Module 4: Integrate and Decouple Ingroup and Outgroup Markers

The previous modules retrieve, align and evaluate ingroup and outgroup
sequences separately. Before concatenation, the two datasets are
reconciled against the Caryophyllales.org taxonomic backbone
(`CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx`; Korotkova *et al*.
2021): synonyms, obsolete names, infraspecific designations and
orthographic variants are resolved to accepted species names in
Cactaceae and Anacampserotaceae, and the Portulacaceae and Talinaceae
terminals are retained as outgroups.

So that outgroup divergence does not inflate the informativeness metrics
of the ingroup (parsimony informative sites, alignment length),
[`integrate_and_clean_markers()`](https://beeamerino.github.io/PhyloCactus/reference/integrate_and_clean_markers.md)
writes its outputs to separate subdirectories:

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
      ingroup, outgroup and joint **species** counts, is written here as
      a provisional (pre-realignment) value and overwritten in place by
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
5.  `4_Cleaned/logs`: Contains the integration and curation log
    (`LOG_clean_integration_summary.txt`).

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
  # ndhF-rpl32 despite the resemblance, and is not aliased.
  marker_aliases = c(trnL = "trnL_trnF"),
  # readmit_markers overrides the Module 3 verdict for a locus kept on outgroup grounds rather
  # than ingroup resolution. Left empty; see the note below.
  readmit_markers = NULL,
  readmit_dir = NULL,
  # Alignment-free homology check, ingroup against outgroup, per marker.
  homology_check = TRUE
)
```

`trnT-psbD` shows why a locus is not readmitted on coverage alone. It
has more outgroup sequences than any other marker (52 sequences of 52
*Portulaca* species, against 50 ingroup sequences), but the two sets are
not the same region: 0.040 of the outgroup 20-mers occur in the ingroup
sequences, against 0.26 to 0.59 for `trnL_trnF`, `matK`, `phyC` and
`rbcL`. The median length is 605 bp in the ingroup and 1347 bp in the
outgroup.

`homology_check` detects this case, which neither counts nor a coverage
table can show. For each marker present on both sides, it measures the
fraction of outgroup k-mers that occur in the ingroup sequences of the
same name, writes `4_Cleaned/tables/TABLE_marker_homology_check.csv`,
and warns for the markers below both thresholds. The measure is
alignment-free, so a low value cannot be corrected by a better
alignment: it separates sequences that are difficult to align from
sequences of different regions. It does not block, because a low value
on a fast-evolving locus calls for inspection and does not by itself
establish non-homology.

### Module 5: Joint Realignment of Curated Markers

Ingroup and outgroup markers were aligned separately in Module 2, and
concatenating separate alignments can shift indel boundaries at the
junction between divergent lineages and bias branch lengths.

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

The threshold is a fraction of the alignment width, which is set by the
longest sequences, those of the ingroup. An outgroup accession that
covers a shorter amplicon of the same region therefore fails on length
alone. In the reference dataset this affects the *Portulaca* sequences
of `trnL_trnF`; losing them removes most of the characters on the branch
that subtends the outgroup, and with them the information that dates the
deepest calibrated nodes.

Three arguments address this. `rooting_pattern` exempts nothing but
names, in a warning, any rooting terminal removed by the filter, so the
loss is reported where it occurs. `protect_pattern` retains matching
terminals regardless of occupancy. `protect_markers` restricts that
exemption to the named loci: an exemption applied to the whole matrix
retains every short fragment of every protected terminal in every locus,
which raises the gap fraction of the supermatrix and destabilizes those
terminals during inference. Restricting it to the loci that lose rooting
terminals keeps them without that cost.

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

The curated and realigned loci are concatenated into a multilocus
supermatrix for partitioned maximum-likelihood inference.

Passing `outgroup_pattern` runs
[`report_marker_group_coverage()`](https://beeamerino.github.io/PhyloCactus/reference/report_marker_group_coverage.md)
on the loci that enter the matrix. It counts the ingroup and outgroup
terminals carrying sequence in each locus and warns when either side is
empty. A locus sampled almost entirely on one side of the root
contributes characters that the other side cannot share, and the branch
lengths across that bipartition are then estimated from the remaining
loci. The report does not block, because whether such a locus belongs in
a matrix depends on the purpose of the matrix: a nuclear locus without
outgroup coverage is of no use for dating and can be informative for
species discrimination.

Passing a named vector of expressions, one per family, resolves the
outgroup by family and adds `n_<family>`, `n_<family>_covered` and
`median_cov_<family>` to
`logs_and_qc/SUPP_TABLE_marker_group_coverage.csv`. The three outgroup
families are not sampled alike, and a combined count can look adequate
while one family, in particular the one that carries the root, is absent
from the locus.

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

`exclude_markers` names the loci left out of this supermatrix. The
alignments are not modified or removed, so an excluded locus remains
available to other analyses. The excluded names are written to the run
log and to `logs_and_qc/TABLE_markers_excluded.csv`, and a name matching
no alignment raises a warning, so a misspelled name cannot leave a locus
in the matrix unnoticed.

Set `exclude_markers = NULL` to build the full matrix. `phyC` is
excluded below because it carries no Portulacaceae terminal: of its 16
outgroup terminals, 15 are Anacampserotaceae and one is Talinaceae,
against 167 ingroup terminals (counted in
`5_MAFFT_Cleaned/aligned_markers/phyC.fasta`; the coverage report lists
only the loci that enter the supermatrix). In a matrix built to date the
ACP clade, a locus sampled for Anacampserotaceae and not for
Portulacaceae supplies characters to only one side of the divergence
between the two families. The exclusion is a property of this matrix;
`phyC` remains the more discriminating of the two nuclear markers for
barcoding.

``` r

run_concatenation_pipeline(
  input_dir = "5_MAFFT_Cleaned/aligned_markers",
  output_dir = "6_Concatenated",
  outgroup_pattern = c(
    Anacampserotaceae = "^(Anacampseros|Grahamia|Talinopsis)_",
    Portulacaceae     = "^Portulaca_",
    Talinaceae        = "^(Talinum|Talinella)_"
  ),
  exclude_markers = "phyC"   # NULL to retain all loci
)
```

## Conclusion

Stage 1 ends with the curated loci, the partitioned supermatrix and its
partition scheme, and the tables documenting each curation step.

Tutorial 2 infers the maximum-likelihood phylogeny with `RAxML-NG`,
estimates branch support, and dates the tree by penalized likelihood
with `treePL`.

[Continue to Tutorial 2: Phylogenetic Inference and Divergence Time
Estimation](https://beeamerino.github.io/PhyloCactus/articles/tutorial-2-cactus-phylogeny-inference.html)

## References

- Guerrero *et al*. 2019. Phylogenetic relationships and evolutionary
  trends in the cactus family. *Journal of Heredity*, 110(1), 4–21.
  <https://doi.org/10.1093/jhered/esy064>
- Katoh, K., & Standley, D. M. 2013. MAFFT multiple sequence alignment
  software version 7: improvements in performance and usability.
  *Molecular Biology and Evolution*, 30(4), 772–780.
  <https://doi.org/10.1093/molbev/mst010>
- Korotkova *et al*. 2021. Cactaceae at Caryophyllales.org, a dynamic
  online species-level taxonomic backbone for the family. *Willdenowia*,
  51(2), 251–270. <https://doi.org/10.3372/wi.51.51208>
- Wright, E. 2024. Fast and flexible search for homologous biological
  sequences with DECIPHER v3. *The R Journal*, 16(2), 191-200.
  <https://doi.org/10.18129/B9.bioc.DECIPHER>
