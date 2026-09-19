# Reconcile and Validate Taxonomic Nomenclature

Reconciles sequence names with an accepted checklist (e.g., the
Caryophyllales.org checklist; Korotkova et al., 2021). Keeps the
sequences whose name matches a name in the checklist, ignoring
differences in spaces, hyphens and underscores, renames them to the
spelling of the checklist, and keeps one sequence per name. Sequences
whose name is not in the checklist are dropped.

## Usage

``` r
clean_taxonomic_names(
  raw_input_fasta,
  checklist_path,
  output_clean_dir,
  force_process = FALSE
)
```

## Arguments

- raw_input_fasta:

  Character. Path to input FASTA file containing raw GenBank sequence
  accessions.

- checklist_path:

  Character. Path to accepted taxonomic checklist CSV or Excel file.

- output_clean_dir:

  Character. Directory path where standardized FASTA sequence output
  will be saved.

- force_process:

  Logical. Force reprocessing if output cached file exists? Defaults to
  `FALSE`.

## Value

Character string path to the generated FASTA file with standardized
species binomials.

## Examples

``` r
if (FALSE) { # \dontrun{
clean_taxonomic_names(
  raw_input_fasta = "marker_rbcl.fasta",
  checklist_path = "CactaceaeFullList_2026_07_01_Beatriz_Merino.xlsx",
  output_clean_dir = "cleaned_names"
)
} # }
```
