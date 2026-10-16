# Build the Final Reference Library of the Molecular Diagnostic Branch

Takes the curated alignments written by
[`curate_barcoding_markers()`](https://beeamerino.github.io/PhyloCactus/reference/curate_barcoding_markers.md)
and writes the reference library used by the validation. For each locus
it (1) removes the sequences that the alignment pipeline logged as
matching the locus in neither direction (`no_match_either_direction` in
`LOG_STRAND_<locus>.csv`), which are paralogs or other non-homologous
copies, not intraspecific replication; (2) collapses accessions with the
same aligned sequence within each species to the alphabetically smallest
sid; (3) evaluates the replication threshold again, counting a species
as replicated when it keeps two or more distinct sequences. Identical
sequences of different species are kept. The columns of the curated
alignment are not modified, so removing sequences can leave columns made
only of gaps.

## Usage

``` r
finalize_barcoding_library(
  assembly_dir = file.path("11_barcoding", "1_assembly"),
  curated_dir = file.path("11_barcoding", "2_curated"),
  screening_dir = file.path("11_barcoding", "3_screening"),
  output_dir = file.path("11_barcoding", "4_library"),
  metadata_file = file.path("11_barcoding", "cache",
    "CACHE_GENBANK_METADATA_BARCODING.csv"),
  min_species_with_replica = 20L,
  paralog_loci = "pepC_like"
)
```

## Arguments

- assembly_dir:

  Character. Output directory of
  [`assemble_barcoding_dataset()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_barcoding_dataset.md);
  the accession registry is read from it.

- curated_dir:

  Character. Output directory of
  [`curate_barcoding_markers()`](https://beeamerino.github.io/PhyloCactus/reference/curate_barcoding_markers.md);
  the curated alignments and the `LOG_STRAND_<locus>.csv` files are read
  from it.

- screening_dir:

  Character. Output directory of
  [`screen_barcoding_markers()`](https://beeamerino.github.io/PhyloCactus/reference/screen_barcoding_markers.md);
  the loci that entered the screening are read from
  `TABLE_barcoding_marker_screening.csv`.

- output_dir:

  Character. Destination; refused inside a directory of the phylogeny.

- metadata_file:

  Character. GenBank metadata cache written by
  [`assemble_barcoding_dataset()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_barcoding_dataset.md),
  read for the description of each excluded sequence.

- min_species_with_replica:

  Integer. A locus enters the library when at least this number of
  species keep two or more distinct sequences. Defaults to `20L`
  (decision of 2026-09-21).

- paralog_loci:

  Character vector. Loci under paralogy surveillance, flagged in the
  column `possible_paralog` of every table; the flag never excludes a
  locus. Defaults to `"pepC_like"`.

## Value

Invisibly, a list with `summary`, `excluded` and `collapsed`. Writes
`LIB_<locus>.fasta` for every locus that enters,
`UNMASKED_<locus>.fasta` beside it (the sequence of step 1, before
masking, of every accession the locus represents, collapsed ones
included; used to recognise a query that is itself a library accession),
`TABLE_barcoding_library_summary.csv`,
`TABLE_barcoding_excluded_nonhomologous.csv`,
`TABLE_barcoding_collapsed_identical.csv` and
`TABLE_barcoding_funnel.csv` (per locus: accessions assembled, curated,
excluded, collapsed and kept, species with replica, and whether the
locus entered the screening and the library).

## Details

Every input is read from the directory of the step that wrote it, so the
library runs on its own whenever those directories exist. The loci are
those marked `enters` in the screening table.

## Examples

``` r
if (FALSE) { # \dontrun{
lib <- finalize_barcoding_library(
  assembly_dir = "11_barcoding/1_assembly",
  curated_dir = "11_barcoding/2_curated",
  screening_dir = "11_barcoding/3_screening"
)
} # }
```
