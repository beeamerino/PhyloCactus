# Curate the Loci of the Molecular Diagnostic Branch

Aligns and curates the per-locus FASTA files written by
[`assemble_barcoding_dataset()`](https://beeamerino.github.io/PhyloCactus/reference/assemble_barcoding_dataset.md)
with the phylogeny's alignment pipeline,
[`run_alignment_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_alignment_pipeline.md),
called unchanged and with the parameters the phylogeny uses for the
ingroup. Only the loci named in `loci` are processed. The pipeline works
sequence by sequence and keeps the `Genus_species|sid` headers, so every
accession remains identified after curation.

## Usage

``` r
curate_barcoding_markers(
  input_dir = file.path("11_barcoding", "1_assembly"),
  output_dir = file.path("11_barcoding", "2_curated"),
  loci = NULL,
  mask_alignment_regions = TRUE,
  min_non_gap_fraction = 0.3,
  max_missing_fraction = 0.3,
  min_masked_alignment_length = 100L,
  preserve_iupac = TRUE,
  fix_strand = TRUE,
  mafft_exec = "mafft",
  mafft_opts = "--auto"
)
```

## Arguments

- input_dir:

  Character. Directory with one FASTA per locus
  (`11_barcoding/1_assembly`).

- output_dir:

  Character. Destination; refused inside a directory of the phylogeny.

- loci:

  Character vector or `NULL`. Loci to curate, as FASTA file names
  without extension. `NULL` (the default) curates every FASTA in
  `input_dir`: the loci come out of the data, not out of a list.

- mask_alignment_regions:

  Logical. Apply automated alignment masking via `DECIPHER`? Defaults to
  `TRUE`.

- min_non_gap_fraction:

  Numeric. Minimum allowable proportion of non-gap characters required
  to retain a site column. Defaults to `0.30`. Applied after masking and
  relative to the post-masking alignment width.

- max_missing_fraction:

  Numeric. Maximum allowable proportion of missing or ambiguous
  characters (`N`) allowed per sequence. Defaults to `0.30`.

- min_masked_alignment_length:

  Integer. Absolute minimum number of alignment columns that must
  survive masking for the locus to be retained. Loci falling below this
  floor are reported as `ZERO_RETAINED` and are not exported. **Defaults
  to `100L`.** The floor exists to catch alignments that masking has
  degraded to the point of being uninformative, not to arbitrate between
  loci of different lengths: no marker in the reference Cactaceae
  dataset is shorter than 100 columns, so any masked alignment below
  that value is degenerate, not merely short. Raising the floor can only
  reject markers, never admit them; a marker rejected by it is reported
  with `decision_reason` naming the threshold, so the effect is always
  visible in the screening table.

  **Scope.** This floor is evaluated only in the masking branch, that is
  when `mask_alignment_regions = TRUE`. With masking deferred (`FALSE`,
  the configuration used for the outgroup) it is not applied, because it
  is defined against a post-masking column set that does not exist in
  that branch. Outgroup terminals are filtered instead by per-sequence
  occupancy in
  [`run_joint_realignment()`](https://beeamerino.github.io/PhyloCactus/reference/run_joint_realignment.md)
  (Module 5). Do not assume this parameter protects both branches.

- preserve_iupac:

  Logical. Retain IUPAC ambiguity codes (`R`, `Y`, `S`, `W`, `K`, `M`,
  `B`, `D`, `H`, `V`) instead of collapsing them to `N` before
  alignment. **Defaults to `TRUE`.** `RAxML-NG` and `ModelTest-NG`
  incorporate ambiguity into the likelihood as a partial constraint, so
  an `R` site restricts the state to A or G whereas an `N` restricts
  nothing: collapsing the codes discards real information for no
  analytical gain. Set to `FALSE` only when an alignment free of
  ambiguity is explicitly required. The amount of ambiguity present at
  every stage is reported in the `mean_fraction_ambiguous_*` and
  `n_sites_ambiguous_*` columns of the marker summary, so the decision
  can be revisited per locus with data. See `@details`.

- fix_strand:

  Logical. Put every sequence of a marker on the same strand before
  alignment, deciding orientation against the marker's own majority.
  GenBank stores each record on whichever strand the submitter
  deposited, and `MAFFT` compares only the orientation it is given: a
  reverse-complemented accession is aligned anyway and contributes
  columns with no positional homology. Each sequence is logged with its
  match in both directions in `tables/LOG_STRAND_<marker>.csv`. Defaults
  to `TRUE`.

- mafft_exec:

  Character. System command or full path to the executable `MAFFT`
  binary. Defaults to `"mafft"`.

- mafft_opts:

  Character. Command-line parameters passed directly to `MAFFT`.
  Defaults to `"--auto"`.

## Value

Invisibly, the summary table returned by
[`run_alignment_pipeline()`](https://beeamerino.github.io/PhyloCactus/reference/run_alignment_pipeline.md).

## Examples

``` r
if (FALSE) { # \dontrun{
curate_barcoding_markers(
  input_dir = "11_barcoding/1_assembly",
  output_dir = "11_barcoding/2_curated"
)
} # }
```
