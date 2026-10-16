# Assemble the Reference Library of the Molecular Diagnostic Branch

Builds the multi-accession reference library of the molecular diagnostic
section from the `phylotaR` workspace shared with the phylogeny. The
workspace is read, or created with the same `phylotaR` call as
[`assemble_ingroup_phylotar()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_ingroup_phylotar.md)
when it does not exist, so either branch can run first. Clusters are
selected with the phylogeny's criterion (more than `min_species`
species). Sequences present in several clusters are resolved with the
phylogeny's rule and the curated accession exclusions are applied, but
none of the reductions to one sequence per species is: every accession
of every species is kept. Each cluster is assigned to a locus with the
same code as the phylogeny, from GenBank metadata cached under
`output_dir/cache/`. Names are reconciled against the checklist with the
rule of
[`clean_taxonomic_names()`](https://beeamerino.github.io/PhyloCactus/reference/clean_taxonomic_names.md);
names that do not match are dropped and reported. Synonyms are not
resolved.

## Usage

``` r
assemble_barcoding_dataset(
  wd_path,
  output_dir = "11_barcoding",
  target_genes_file = NULL,
  genes_map_file = NULL,
  manual_exclusions_file = NULL,
  apply_manual_exclusions = TRUE,
  checklist_path = NULL,
  min_species = 50,
  preferred_parent = "3593",
  taxids = NULL,
  ncbi_dr = NULL,
  force_download = FALSE,
  mnsql = 100L,
  mxsql = 5000L,
  phylogeny_map_file = NULL,
  extra_records_dir = NULL,
  own_dir = NULL,
  notify = FALSE,
  notify_to = NULL,
  notify_credentials = NULL
)
```

## Arguments

- wd_path:

  Character. The `phylotaR` workspace, shared with the phylogeny (e.g.
  `0_phylotaR_raw_Ingroup`).

- output_dir:

  Character. Root of the branch. Refused if it lies inside a directory
  of the phylogeny. Defaults to `"11_barcoding"`.

- target_genes_file, genes_map_file:

  Character. Marker definitions; default to the files distributed with
  the package, as in
  [`assemble_ingroup_phylotar()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_ingroup_phylotar.md).

- manual_exclusions_file:

  Character. Curated accession exclusions; defaults to
  `inst/extdata/manual_exclusions_ingroup.csv`.

- apply_manual_exclusions:

  Logical. Apply the curated exclusions? Defaults to `TRUE`.

- checklist_path:

  Character. Accepted checklist (Korotkova *et al*. 2021); `NULL`, the
  default, uses the checklist distributed with the package. `NA` applies
  no checklist and keeps every GenBank name, cleaned to `Genus_species`:
  that is what the outgroup of the outgroup control needs, since a
  checklist of Cactaceae would discard every outgroup name.

- min_species:

  Integer. A cluster is retained when it holds more than this number of
  species. Defaults to `50`, as in the phylogeny.

- preferred_parent:

  Character. NCBI Taxonomy ID of the focal ingroup, used to decide which
  cluster keeps a sid that sits in two. Defaults to `"3593"`.

- taxids:

  Character vector or `NULL`. What to mine from GenBank when the
  workspace of `wd_path` does not exist yet, which is what lets the
  branch run without the phylogeny having run first. `NULL` mines
  `preferred_parent`, the behaviour of earlier versions. Several IDs are
  mined together, as
  [`barcoding_outgroup_taxids()`](https://beeamerino.github.io/PhyloCactus/reference/barcoding_outgroup_taxids.md)
  returns them for the outgroup control.

- ncbi_dr:

  Character or `NULL`. Directory of the BLAST+ binaries, used only when
  the workspace has to be created.

- force_download:

  Logical. Mine GenBank again even if the workspace exists. Defaults to
  `FALSE`.

- mnsql, mxsql:

  Integer. Shortest and longest sequence, in bases, that a new mining
  keeps
  ([`phylotaR::setup()`](https://docs.ropensci.org/phylotaR/reference/setup.html));
  100 and 5000 by default. Used only when GenBank is mined.

- phylogeny_map_file:

  Character or `NULL`. Cluster-to-locus table of the phylogeny, compared
  with the branch's own when it exists. Defaults to
  `1_phylotaR_out_Ingroup/TABLE_CLUSTER_MARKER_ASSIGNMENT_INGROUP.csv`
  next to `wd_path`.

- extra_records_dir:

  Character or `NULL`. A folder written by
  [`extract_barcoding_plastome_loci()`](https://beeamerino.github.io/PhyloCactus/reference/extract_barcoding_plastome_loci.md);
  its records are appended to their loci and the registry gains a
  `source` column (`phylotaR` or `genbank_plastome`). `NULL`, the
  default, leaves the step as it was.

- own_dir:

  Character or `NULL`. A folder of own data (`0_own/`): `SAMPLES.csv`
  (`sample`, `species`, `voucher`, `data_type`, `published`, `notes`)
  and one FASTA per library locus, named as the locus, with headers
  `>sample`. The records join their loci before the cut of long records,
  with the sid `own:<sample>:<locus>` and `source = "own"`; a record
  whose sample is not in the sheet, whose species is not in the
  checklist or whose locus is not in the registry is left out and listed
  in `TABLE_barcoding_own_left_out.csv`. `NULL`, the default, leaves the
  step as it was.

- notify:

  Logical. Send an email through
  [`send_run_notification()`](https://beeamerino.github.io/PhyloCactus/reference/send_run_notification.md)
  when step 1 ends, finished or failed, and another when a mining ends
  or fails (reading an existing workspace mines nothing). A notification
  that cannot be sent is reported and never fails the step. Defaults to
  `FALSE`.

- notify_to, notify_credentials:

  Passed to
  [`send_run_notification()`](https://beeamerino.github.io/PhyloCactus/reference/send_run_notification.md)
  as `to` and `credentials`.

## Value

Invisibly, a list with `registry`, `summary`, `marker_map`, `comparison`
and `counts`. Files written under `output_dir/1_assembly/`: one FASTA
per locus, `TABLE_barcoding_accession_registry.csv`,
`TABLE_barcoding_cluster_marker_assignment.csv`,
`TABLE_barcoding_marker_assignment_comparison.csv` (when the phylogeny's
table exists), `TABLE_barcoding_discarded_names.csv` and
`TABLE_barcoding_replication_summary_assembly.csv`.

## References

Korotkova, N., Aquino, D., Arias, S., Eggli, U., Franck, A.,
Gómez-Hinostrosa, C., Guerrero, P. C., Hernández, H. M., Kohlbecker, A.,
Köhler, M., Luther, K., Majure, L. C., Müller, A., Metzing, D.,
Nyffeler, R., Sánchez, D., Schlumpberger, B. O., & Berendsohn, W. G.
(2021). Cactaceae at Caryophyllales.org - a dynamic online species-level
taxonomic backbone for the family. *Willdenowia*, 51(2), 251-270.
[doi:10.3372/wi.51.51208](https://doi.org/10.3372/wi.51.51208)

## Examples

``` r
if (FALSE) { # \dontrun{
assemble_barcoding_dataset(
  wd_path = "0_phylotaR_raw_Ingroup",
  output_dir = "11_barcoding"
)
} # }
```
