# Draw the answers of several identifications side by side

One row per identification and one column per locus found in any of
them, each cell coloured by the state of that locus (1 species, 2 genus,
3 no answer); a cross marks a locus named at genus rank only by rule
(reason `genus_only_locus`, argument `genus_only_loci` of
[`identify_barcoding_query()`](https://beeamerino.github.io/PhyloCactus/reference/identify_barcoding_query.md)).
The answer of each identification (one species, one genus, a conflict,
or no answer) is written next to its name. When a locus was found in
more than one region of a query (contigs of reads), the cell shows its
most informative row.

## Usage

``` r
plot_barcoding_identification_summary(answers)
```

## Arguments

- answers:

  Named list. Each element is the answer table of one identification
  (`REPORT_<run_name>_answer.csv` or its path, as in
  [`plot_barcoding_identification()`](https://beeamerino.github.io/PhyloCactus/reference/plot_barcoding_identification.md));
  the names label the rows.

## Value

A `ggplot` object; its data hold one row per identification and locus
found, with the columns `example`, `locus`, `state`, `genus_only` and
`answer`. Nothing is written.

## See also

[`plot_barcoding_identification()`](https://beeamerino.github.io/PhyloCactus/reference/plot_barcoding_identification.md).

## Examples

``` r
if (FALSE) { # \dontrun{
plot_barcoding_identification_summary(list(
  "MN517611.1" = "11_barcoding/10_identify/REPORT_MN517611.1_answer.csv",
  "LT745702.1" = "11_barcoding/10_identify/REPORT_LT745702.1_answer.csv"))
} # }
```
