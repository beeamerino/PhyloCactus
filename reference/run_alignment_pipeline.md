# Execute Complete Alignment and Gap-Masking Pipeline

Orchestrates multiple sequence alignment (MSA) and automated
quality-control masking across orthologous sequence clusters. Primary
alignment hypotheses are inferred using `MAFFT` (Katoh & Standley,
2013). Subsequently, ambiguous sites, poorly aligned terminal fragments,
and non-homologous insertions are masked using the `DECIPHER` framework
(Wright, 2024), eliminating systematic noise while retaining
phylogenetically informative nucleotide positions for downstream
supermatrix assembly.

## Usage

``` r
run_alignment_pipeline(
  input_folder,
  output_dir,
  fasta_pattern = "\\.fasta$",
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

- input_folder:

  Character. Directory path containing raw unaligned orthologous FASTA
  files.

- output_dir:

  Character. Root destination directory for output subfolders
  (`alignments/`, `tables/`, `logs/`).

- fasta_pattern:

  Character. Regular expression pattern matching target FASTA files.
  Defaults to `"\\.fasta$"`.

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
  floor are reported as `ZERO_RETAINED` rather than exported as
  near-empty alignments. **Defaults to `100L`.** The floor exists to
  catch alignments that masking has degraded to the point of being
  uninformative, not to arbitrate between loci of different lengths: no
  marker in the reference Cactaceae dataset is shorter than 100 columns,
  so any masked alignment falling below that value is degenerate rather
  than merely short. Raising the floor can only reject markers, never
  admit them; a marker rejected by it is reported with `decision_reason`
  naming the threshold, so the effect is always visible in the screening
  table.

  **Scope.** This floor is evaluated only in the masking branch, that is
  when `mask_alignment_regions = TRUE`. With masking deferred (`FALSE`,
  the configuration used for the outgroup) it is deliberately not
  applied, because it is defined against a post-masking column set that
  does not exist in that branch. Outgroup terminals are filtered instead
  by per-sequence occupancy in
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
  columns with no positional homology. Found on 2026-09-01 in the `rbcL`
  and `matK` accessions of `Portulaca oleracea` and `P. pilosa`, which
  sat at 0.51 and 0.40 observed divergence from Cactaceae where
  `P. grandiflora` sits at 0.030 and 0.066. Each sequence is logged with
  its match in both directions in `tables/LOG_STRAND_<marker>.csv`.
  Defaults to `TRUE`.

- mafft_exec:

  Character. System command or full path to the executable `MAFFT`
  binary. Defaults to `"mafft"`.

- mafft_opts:

  Character. Command-line parameters passed directly to `MAFFT`.
  Defaults to `"--auto"`.

## Value

A data frame containing site length, missingness, and sequence retention
statistics across processed loci.

## Details

**Masking policy.** With `mask_alignment_regions = TRUE` (the default)
the `MAFFT` alignment is passed through
[`DECIPHER::RemoveGaps()`](https://rdrr.io/pkg/DECIPHER/man/RemoveGaps.html)
and
[`DECIPHER::MaskAlignment()`](https://rdrr.io/pkg/DECIPHER/man/MaskAlignment.html),
and the resulting post-masking column set is what
`min_non_gap_fraction`, `max_missing_fraction` and
`min_masked_alignment_length` are evaluated against. With
`mask_alignment_regions = FALSE` the `MAFFT` alignment is exported
unmodified and **those three thresholds are not applied at this stage**:
they are relative to a masked column set that does not exist in that
branch. Occupancy filtering then happens once, on the joint ingroup plus
outgroup alignment, in
[`run_joint_realignment()`](https://beeamerino.github.io/PhyloCactus/reference/run_joint_realignment.md).
This is the recommended setting for a sparsely sampled outgroup: masking
two to seven highly divergent accessions on their own defines a column
set the ingroup does not share and can erode an outgroup alignment to a
few base pairs, which then fails the occupancy filter of the joint
realignment and removes the outgroup precisely from the loci needed to
root the tree. Both the chosen policy and the parameter values are
written to `logs/LOG_alignment_run_info.txt` and stamped into every
per-marker summary table.

**Caching.** A marker is reused from a previous run only when its
exported alignment, summary table and filter log all exist *and* the
five parameters stamped into the cached summary table match the current
call. Any parameter change invalidates the cache and the marker is
reprocessed, so changing `mask_alignment_regions` on an already
populated output directory takes effect instead of silently returning
the previous alignments.

**Ambiguity policy.** IUPAC ambiguity codes are retained by default
(`preserve_iupac = TRUE`). The downstream tools all accept them, and
they carry information that `N` does not: an `R` site constrains the
state to A or G, an `N` constrains nothing. There is no analytical
reason to discard that constraint, so the pipeline does not.

Ambiguity is nevertheless measured rather than assumed away. The marker
summary reports `mean_fraction_ambiguous_*` for each of the five
processing stages, and `n_sites_ambiguous_*` for the raw input and the
final alignment. The `raw_input` figures are computed before
`clean_ambiguous()` is applied, so they record what the source records
actually contained regardless of the policy in force. That is what makes
a per-locus decision possible: a marker whose ambiguity is concentrated
rather than diffuse can be examined on its own evidence instead of being
subjected to a global rule.

## References

Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment
software version 7: Improvements in performance and usability.
*Molecular Biology and Evolution*, 30(4), 772–780.
[doi:10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)

Wright, E. S. (2024). Fast and Flexible Search for Homologous Biological
Sequences with DECIPHER v3. *The R Journal*, 16(2), 191-200.
[doi:10.18129/B9.bioc.DECIPHER](https://doi.org/10.18129/B9.bioc.DECIPHER)

## Examples

``` r
if (FALSE) { # \dontrun{
run_alignment_pipeline(
  input_folder = "1_phylotaR_out_ingroup",
  output_dir = "2_MAFFT_Cactaceae",
  min_non_gap_fraction = 0.30,
  max_missing_fraction = 0.30
)
} # }
```
