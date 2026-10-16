# Screen the Curated Loci of the Molecular Diagnostic Branch

For each curated locus, recomputes the five context columns (species,
species with a replicate, accessions, genera, genera with two or more
species) on the accessions that survived curation, applies the branch's
replication threshold to them, flags saturation (which never excludes a
locus), and counts identical aligned sequences per species. Inputs are
read from the directories of the earlier steps, so the screening runs on
its own whenever those directories exist.

## Usage

``` r
screen_barcoding_markers(
  assembly_dir = file.path("11_barcoding", "1_assembly"),
  curated_dir = file.path("11_barcoding", "2_curated"),
  loci = NULL,
  output_dir = file.path("11_barcoding", "3_screening"),
  min_species_with_replica = 20L,
  saturation_flag_cutoff = 0.3,
  paralog_loci = "pepC_like"
)
```

## Arguments

- assembly_dir:

  Character. Output directory of
  [`assemble_barcoding_dataset()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_barcoding_dataset.md);
  the accession registry `TABLE_barcoding_accession_registry.csv` is
  read from it.

- curated_dir:

  Character. Output directory of
  [`curate_barcoding_markers()`](https://beeamerino.github.io/PhyloCactus/reference/curate_barcoding_markers.md).

- loci:

  Character vector or `NULL`. Loci to screen. `NULL` (the default)
  screens every curated alignment in `curated_dir`; the threshold then
  decides which loci enter (decision of 2026-09-22).

- output_dir:

  Character. Destination of the screening tables.

- min_species_with_replica:

  Integer. A locus enters the branch when at least this number of
  species keep two or more accessions after curation. Defaults to `20L`
  (decision of 2026-09-21).

- saturation_flag_cutoff:

  Numeric. Slope threshold of the saturation proxy, as in
  [`run_marker_screening()`](https://beeamerino.github.io/PhyloCactus/reference/run_marker_screening.md).
  Defaults to `0.3`.

- paralog_loci:

  Character vector. Loci under paralogy surveillance, flagged in the
  column `possible_paralog` of every table; the flag never excludes a
  locus. Defaults to `"pepC_like"`.

## Value

Invisibly, a list with `screening` and `identical`. Writes
`TABLE_barcoding_marker_screening.csv` and
`TABLE_barcoding_identical_sequences.csv`.
