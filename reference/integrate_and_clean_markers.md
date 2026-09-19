# Final Marker Integration, Taxonomic Cleaning, and Ingroup-Outgroup Partitioning

Integrates independently curated ingroup and outgroup sequence datasets,
standardizes species binomials against the accepted checklist (Korotkova
et al. 2021), and exports decoupled FASTA sequence directories
(`cleaned_markers_ingroup/`, `cleaned_markers_outgroup/`,
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
  `c(trnL = "trnL_trnF")`. Names and values are normalized the same way
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
  a different question from the one Module 3 answers. A marker brought
  back for this reason must be homologous on both sides; check
  `TABLE_marker_homology_check.csv` (`homology_check = TRUE`) before
  readmitting it. Defaults to `NULL`.

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
  below **both**, at `k = 20` and `k = 10` respectively. In the
  reference dataset `trnL_trnF`, `matK`, `phyC` and `rbcL` return 0.26
  to 0.59 at `k = 20`; `trnT-psbD` returns 0.040 (ingroup median 605 bp
  against outgroup 1347 bp). Default to `0.05` and `0.15`.

## Value

A data frame with the marker summary, with separate ingroup and joint
metrics.

## Alias collisions

`marker_aliases` can point two source files at the same `marker_key`.
When both carry a record for the same species, that species appears
twice under one FASTA header in the exported marker, and the downstream
matrix assembly resolves the duplication by file order, silently and in
favor of whichever record sorts first, which is not necessarily the more
informative one. The duplicates are therefore resolved here by an
explicit rule: the record with the greatest number of non-gap sites is
retained, ties are broken on the source file name so that the outcome
does not depend on the order in which the directory was listed, every
competing record is written to
`tables/TABLE_alias_collisions_resolved.csv` with its source file, its
length and whether it was retained, and a warning names the affected
headers.

## References

Korotkova, N., Aquino, D., Arias, S., Eggli, U., Franck, A.,
Gómez-Hinostrosa, C., Guerrero, P. C., Hernández, H. M., Kohlbecker, A.,
Köhler, M., Luther, K., Majure, L. C., Müller, A., Metzing, D.,
Nyffeler, R., Sánchez, D., Schlumpberger, B. O., & Berendsohn, W. G.
(2021). Cactaceae at Caryophyllales.org - a dynamic online species-level
taxonomic backbone for the family. *Willdenowia*, 51(2), 251–270.
[doi:10.3372/wi.51.51208](https://doi.org/10.3372/wi.51.51208)
