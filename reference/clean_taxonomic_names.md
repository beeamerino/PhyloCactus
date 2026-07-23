# Reconcile and Validate Taxonomic Nomenclature

Reconciles sequence tip labels against authoritative botanical
checklists (e.g., Caryophyllales.org checklist). Resolves taxonomic
synonymies, infraspecific variants, and orthographic errors,
guaranteeing nomenclatural stability across public GenBank sequence
downloads.

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
  \`FALSE\`.

## Value

Character string path to the generated FASTA file with standardized
species binomials.

## Examples

``` r
if (FALSE) { # \dontrun{
clean_taxonomic_names(
  raw_input_fasta = "marker_rbcl.fasta",
  checklist_path = "CactaceaeFullList_accepted.csv",
  output_clean_dir = "cleaned_names"
)
} # }
```
