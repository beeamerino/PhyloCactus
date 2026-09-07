# Final Marker Integration, Taxonomic Cleaning, and Ingroup-Outgroup Partitioning

Integrates independently curated ingroup and outgroup sequence datasets,
standardizes species binomials against the authoritative taxonomic
checklist (Korotkova et al. 2021), and exports decoupled FASTA sequence
directories (`cleaned_markers_ingroup/`, `cleaned_markers_outgroup/`,
`cleaned_markers_joint/`). Computes isolated molecular informativeness
metrics (tips, length, variable sites, parsimony informative sites, GC
content, and missingness) for the focal ingroup radiation to avoid
outgroup-driven inflation artifacts.

## Usage

``` r
integrate_and_clean_markers(
  ingroup_dir,
  outgroup_dir,
  output_dir,
  accepted_list_file,
  metadata_in_file,
  metadata_out_file,
  marker_aliases = NULL,
  readmit_markers = NULL,
  readmit_dir = NULL,
  drop_report = TRUE,
  homology_check = TRUE,
  homology_k20_min = 0.05,
  homology_k10_min = 0.15
)
```

## Arguments

- ingroup_dir:

  Character. Path to directory containing filtered ingroup FASTA
  alignments.

- outgroup_dir:

  Character. Path to directory containing filtered outgroup FASTA
  alignments.

- output_dir:

  Character. Root destination directory for cleaned FASTA outputs and
  diagnostic tables.

- accepted_list_file:

  Character. Path to accepted botanical checklist CSV or Excel file
  (e.g., `CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx`).

- metadata_in_file:

  Character. Path to ingroup accession occupancy metadata CSV file.

- metadata_out_file:

  Character. Path to outgroup accession occupancy metadata CSV file.

- marker_aliases:

  Named character vector or `NULL`. Renames outgroup markers onto their
  ingroup counterpart before the two sets are joined, as
  `c(trnL = "trnL_trnF")`. Names and values are normalised the same way
  the filenames are, so either spelling works.

  The two `phylotaR` runs cluster independently, so the same region can
  be named differently in each: the outgroup carries two `trnL`
  accessions of the trnL intron, which is the first half of the ingroup
  `trnL-trnF` amplicon, and they share 63% of their 20-mers with it.
  Aliasing is a homology claim and has to be evidenced, not assumed from
  the name: outgroup `ndhF` shares no 20-mer with ingroup `ndhF-rpl32`
  despite the obvious resemblance, because the two datasets cover
  different parts of the gene. Defaults to `NULL`.

- readmit_markers:

  Character vector or `NULL`. Markers that
  [`run_marker_screening()`](https://beeamerino.github.io/PhyloCactus/reference/run_marker_screening.md)
  rejected and that are brought back into the joint dataset because they
  connect the outgroup to the ingroup, read from `readmit_dir`. This is
  a different question from the one Module 3 answers, and it is stated
  as such: `trnT-psbD` resolves the ingroup poorly, with 50 sequences
  over 1008 Cactaceae terminals, and carries the best outgroup coverage
  of the dataset, 49 *Portulaca* accessions of 1347 bp against the five
  or six terminals of overlap the retained loci provide. Retaining it
  for ingroup resolution would be wrong; retaining it for rooting is the
  point. Defaults to `NULL`.

- readmit_dir:

  Character or `NULL`. Directory holding the pre-screening ingroup
  alignments, normally `2_MAFFT_Cactaceae/alignments`. Required when
  `readmit_markers` is supplied. Defaults to `NULL`.

- drop_report:

  Logical. Write `tables/TABLE_outgroup_markers_dropped.csv` listing the
  outgroup markers that have no ingroup counterpart, with their sequence
  counts, and warn about them. The joint marker set is the ingroup's, so
  those markers leave the analysis at this point; a run takes hours and
  a console warning alone scrolls away. Defaults to `TRUE`.

- homology_check:

  Logical. For every marker present on both sides, measure the fraction
  of outgroup k-mers that occur in the ingroup sequences of the same
  name, write `tables/TABLE_marker_homology_check.csv`, and warn about
  the markers that share almost none. The measure is alignment-free, so
  it distinguishes sequences that are hard to align from sequences that
  are not the same region. Reported, never blocking. Defaults to `TRUE`.

- homology_k20_min, homology_k10_min:

  Numeric. A marker is flagged `no_detectable_homology` when it falls
  below **both**, at `k = 20` and `k = 10` respectively. Genuine
  counterparts in this dataset return 0.28 to 0.60 at `k = 20`; the two
  false pairs found on 2026-09-01 returned 0.000 (`pepC`, two PEPC
  paralogues) and 0.039 (`trnT-psbD`, ingroup median 598 bp against
  outgroup 1347 bp). Default to `0.05` and `0.15`.

## Value

A data frame containing the comprehensive marker summary with decoupled
ingroup and joint metrics.

## References

Korotkova, N., Aquino, D., Arias, S., Eggli, U., Franck, A.,
Gómez-Hinostrosa, C., Guerrero, P. C., Hernández, H. M., Kohlbecker, A.,
Köhler, M., Luna, R., Machado, M., Merclinger, M., Nyffeler, R.,
Salvador-Montiel, S., Sánchez, D., Schlumpberger, B. O., & Berendsohn,
W. G. (2021). Cactaceae at Caryophyllales.org - a dynamic online
species-level taxonomic backbone for the family. *Willdenowia*, 51(2),
251–270. [doi:10.3372/wi.51.51208](https://doi.org/10.3372/wi.51.51208)
