# Assemble Ingroup Sequence Clusters via phylotaR Similarity Mining

Retrieves orthologous sequence clusters for a focal taxonomic ingroup
(e.g., family **Cactaceae**, NCBI Taxonomy ID: 3593) directly from
GenBank using similarity clustering via `phylotaR` (Bennett *et al.*,
2018). Clusters are formed by sequence similarity, not by GenBank locus
annotations, so sequences are not lost to gene-name synonyms or
mislabeled records.

## Usage

``` r
assemble_ingroup_phylotar(
  wd_path,
  target_genes_file = NULL,
  genes_map_file = NULL,
  manual_exclusions_file = NULL,
  apply_manual_exclusions = TRUE,
  min_species = 50,
  preferred_parent = "3593",
  ncbi_dr = NULL,
  force_download = FALSE,
  out_dir = "1_phylotaR_out_Ingroup",
  mnsql = 100L,
  mxsql = 5000L,
  notify = FALSE,
  notify_to = NULL,
  notify_credentials = NULL
)
```

## Arguments

- wd_path:

  Character. Path to the `phylotaR` workspace directory storing local
  parameters and database caches.

- target_genes_file:

  Character. Path to the target locus list text file. If `NULL`,
  defaults to package `inst/extdata/target_genes.txt`.

- genes_map_file:

  Character. Path to the gene synonymy mapping CSV file. If `NULL`,
  defaults to package `inst/extdata/genes_map.csv`.

- manual_exclusions_file:

  Character. Path to the accession exclusion CSV file. If `NULL`,
  defaults to package `inst/extdata/manual_exclusions_ingroup.csv`.

- apply_manual_exclusions:

  Logical. Apply the curated accession exclusion list? **Defaults to
  `TRUE`.** The lists shipped with the package are the product of manual
  curation by the expert team supporting `PhyloCactus`: each excluded
  accession was inspected and removed on taxonomic or sequence-quality
  grounds that are not recoverable from GenBank metadata alone. The
  lists are keyed by locus and accession (`locus`, `sid`, `reason`,
  `source`), so that they survive a new mining that renumbers the
  clusters: an accession is removed from every cluster of the locus
  named. Across both lists they hold 102 accessions, 60 in 5 ingroup
  loci and 42 in 13 outgroup loci, and they live in
  `inst/extdata/manual_exclusions_ingroup.csv` and
  `inst/extdata/manual_exclusions_outgroup.csv`. Setting this to `FALSE`
  reproduces the uncurated cluster set, which is the way to quantify
  what the curation actually removes; the applied list is always written
  to `TABLE_MANUAL_EXCLUSIONS_*.csv` in the output directory, so any run
  documents its own curation state.

- min_species:

  Integer. Species-richness threshold for cluster retention. Clusters
  are retained when they contain **strictly more than** `min_species`
  distinct species, so `min_species = 50` keeps clusters with 51 species
  or more. Defaults to `50`.

- preferred_parent:

  Character. NCBI Taxonomy ID of the focal ingroup parent node. Defaults
  to `"3593"` (Cactaceae).

- ncbi_dr:

  Character. Path to local `BLAST+` binaries directory. If `NULL`,
  attempts system environment auto-detection.

- force_download:

  Logical. Force fresh database retrieval instead of using local cache?
  Defaults to `FALSE`.

- out_dir:

  Character. Output directory path to save cluster summaries, FASTA
  sequence matrices, and log reports.

- mnsql, mxsql:

  Integer. Shortest and longest sequence, in bases, that a new mining
  keeps
  ([`phylotaR::setup()`](https://docs.ropensci.org/phylotaR/reference/setup.html));
  100 and 5000 by default. Used only when GenBank is mined.

- notify:

  Logical. Send an email through
  [`send_run_notification()`](https://beeamerino.github.io/PhyloCactus/reference/send_run_notification.md)
  when a mining ends or fails; reading an existing workspace sends
  nothing. Defaults to `FALSE`.

- notify_to, notify_credentials:

  Passed to
  [`send_run_notification()`](https://beeamerino.github.io/PhyloCactus/reference/send_run_notification.md)
  as `to` and `credentials`.

## Value

A data frame summarizing sequence occupancy, taxon representation, and
cluster characteristics across retained loci.

## Details

Records of more than 2000 bases (for example `matK` records that carry
the flanking `trnK` intron) are cut to the extent that the records of
2000 bases or fewer of the same marker cover, before one record per
species is kept, so that the records of a partition span the same
region. A long record that does not overlap that extent is left out, and
a marker with no record of 2000 bases or fewer gets no FASTA file. Every
long record, with what was done to it, is written to
`TABLE_LONG_RECORDS_CUT_INGROUP.csv` in `out_dir`, and the log gives the
counts.

## References

Bennett, D. J., Hettling, H., Silvestro, D., Zizka, A., Bacon, C. D.,
Faurby, S., ... & Antonelli, A. (2018). phylotaR: An automated pipeline
for retrieving orthologous DNA sequences from GenBank in R. *Life*,
8(2), 20. [doi:10.3390/life8020020](https://doi.org/10.3390/life8020020)

## Examples

``` r
if (FALSE) { # \dontrun{
assemble_ingroup_phylotar(
  wd_path = "0_phylotaR_raw_Ingroup",
  min_species = 50,
  preferred_parent = "3593"
)
} # }
```
