# Extra library records cropped from GenBank plastomes (step 1)

Downloads complete plastomes of Cactaceae from GenBank, crops from each
the loci of the library and writes them as extra records that
[`assemble_barcoding_dataset()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_barcoding_dataset.md)
can take with `extra_records_dir`. Nothing is classified: a record takes
the name of its plastome.

## Usage

``` r
extract_barcoding_plastome_loci(
  query = paste0("Cactaceae[Organism] AND 100000:170000[SLEN] AND ",
    "(chloroplast OR plastid) AND \"complete genome\" NOT refseq[filter]"),
  library_dir = file.path("11_barcoding", "4_library"),
  output_dir = file.path("0_genomic_raw", "plastomes"),
  checklist_path = system.file("extdata",
    "CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx", package = "PhyloCactus"),
  min_region = 100L,
  keep_genera = character(0),
  fetch = .bc_fetch_genbank
)
```

## Arguments

- query:

  Character. GenBank query of the plastomes.

- library_dir:

  Character. Directory of step 4, with the `LIB_<locus>.fasta` files.

- output_dir:

  Character. Where the downloads and the extra records are written.

- checklist_path:

  Character. The checklist (a path, or a vector of accepted names).

- min_region:

  Integer. Shortest region written. Defaults to 100.

- keep_genera:

  Character. Genera outside the checklist whose plastomes are kept under
  the NCBI binomial (for the outgroup: Portulacaceae and Talinaceae). A
  name without an epithet (`sp.`, `cf.`, `aff.`) is still left out.
  Defaults to none.

- fetch:

  Function `(query, ids)` returning GenBank flat text; defaults to a
  fetcher through `rentrez`.

## Value

Invisibly, a list with `records`, `skipped` and `unmatched`.

## Details

- **Names** are matched to the checklist with the rule of step 1 and cut
  to the binomial, the rank of the library; unmatched names are left out
  and listed in `TABLE_unmatched_names.csv`.

- **Extraction**, per locus: the crop of
  [`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md)
  (the window with the most 20-mer hits of the locus) and the cut to the
  core of the locus (columns covered by half or more of its library
  sequences), so a record has the extent of the library; regions shorter
  than `min_region` are not written. Where the crop finds more than one
  window (an inverted repeat), the best one is kept and `other_windows`
  counts the rest.

- **Duplicates:** RefSeq is left out of the query, and a plastome whose
  sequence (on either strand) repeats one already downloaded adds
  nothing; it is listed in `TABLE_duplicate_plastomes.csv`.

- **Identifier:** `sid = <accession>__<locus>`, so a sid stays unique
  across loci.

- **Same specimen:** a plastome whose `specimen_voucher` matches that of
  a library accession of the same species in the same locus (same
  collector surname and number, or the same normalised string) does not
  add that locus; the pair is listed in
  `TABLE_skipped_same_specimen.csv`.

Writes to `output_dir`: `downloads/<accession>.gb`, `MANIFEST.csv`
(accession, organism, length, md5, query, date), `extra_records.fasta`
(headers `Genus_species|sid`), `TABLE_genomic_extra_records.csv` (sid,
accession, species, organism, locus, region_start, region_end, strand,
region_length, other_windows, md5 of the plastome file, library_md5 of
the `LIB_<locus>.fasta` that set the extent),
`TABLE_skipped_same_specimen.csv` and `TABLE_unmatched_names.csv`.

## See also

[`assemble_barcoding_dataset()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_barcoding_dataset.md),
[`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md).

## Examples

``` r
if (FALSE) { # \dontrun{
extract_barcoding_plastome_loci()
} # }
```
