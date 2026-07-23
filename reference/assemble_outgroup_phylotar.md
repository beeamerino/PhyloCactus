# Assemble Outgroup Sequence Clusters via phylotaR

Retrieves orthologous sequence clusters for specified outgroup lineages
(e.g., \*Portulaca\*, \*Anacampseros\*, \*Talinopsis\*, \*Grahamia\*)
matching the locus target constraints defined for the focal ingroup.
Outer reference sampling provides phylogenetically informative root
positions necessary for maximum-likelihood tree search and divergence
time estimation.

## Usage

``` r
assemble_outgroup_phylotar(
  wd_path,
  target_genes_file = NULL,
  genes_map_file = NULL,
  manual_exclusions_file = NULL,
  outgroups = c("107598", "107617", "107583", "3582"),
  force_download = FALSE,
  out_dir = "1_phylotaR_out_Outgroup"
)
```

## Arguments

- wd_path:

  Character. Path to the \`phylotaR\` workspace directory.

- target_genes_file:

  Character. Path to the target locus list text file. If \`NULL\`,
  defaults to package \`inst/extdata/target_genes.txt\`.

- genes_map_file:

  Character. Path to the gene synonymy mapping CSV file. If \`NULL\`,
  defaults to package \`inst/extdata/genes_map.csv\`.

- manual_exclusions_file:

  Character. Path to outgroup exclusions CSV file. If \`NULL\`, defaults
  to package \`inst/extdata/manual_exclusions_outgroup.csv\`.

- outgroups:

  Character vector of NCBI Taxonomy IDs for outgroup lineages. Defaults
  to \`c("107598", "107617", "107583", "3582")\`.

- force_download:

  Logical. Force fresh database retrieval instead of using local cache?
  Defaults to \`FALSE\`.

- out_dir:

  Character. Output directory path to save outgroup cluster tables and
  FASTA sequence files.

## Value

A list containing the processed outgroup cluster objects and retained
cluster IDs.

## References

Bennett, D. J., Hettling, H., Silvestro, D., Zizka, A., Bacon, C. D.,
Faurby, S., ... & Antonelli, A. (2018). phylotaR: An automated pipeline
for retrieving orthologous DNA sequences from GenBank in R. \*Life\*,
8(2), 20. [doi:10.3390/life8020020](https://doi.org/10.3390/life8020020)

## Examples

``` r
if (FALSE) { # \dontrun{
assemble_outgroup_phylotar(
  wd_path = "0_phylotaR_raw_Outgroup",
  outgroups = c("3582", "107583")
)
} # }
```
