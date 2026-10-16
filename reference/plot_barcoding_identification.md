# Draw the answer of an identification

Draws, from the tables the identification report writes, the answer per
locus of an identification: for each locus found, the IdTaxa confidence
at genus rank (striped bar) and at species rank (solid bar), coloured by
state (1 species, 2 genus, 3 no answer), with the threshold as a dashed
line. When the sampling table is given, the sampling of the named genera
in the library (sequences per species and locus) is drawn below. The
figure is the one of section 2 and section 6 of the report
([`report_barcoding_identification()`](https://beeamerino.github.io/PhyloCactus/reference/report_barcoding_identification.md)),
drawn without the library.

## Usage

``` r
plot_barcoding_identification(
  answer,
  sampling = NULL,
  threshold = 60,
  max_species = 40L
)
```

## Arguments

- answer:

  Data frame or path. `REPORT_<run_name>_answer.csv`, or a table with
  its columns `query`, `locus`, `compartment`, `state`,
  `predicted_species`, `predicted_genus`, `genus_confidence` and
  `species_confidence` (`reason`, if present, leaves out the loci not
  found).

- sampling:

  Data frame, path or `NULL`. `REPORT_<run_name>_sampling.csv`, or a
  table with its columns `genus`, `species`, `locus` and `sequences`;
  restricted to the genera named in `answer`.

- threshold:

  Numeric. The IdTaxa threshold drawn as a dashed line. Defaults to 60.

- max_species:

  Integer. Above this number of species in the sampling, the species not
  named are gathered into one row per genus, whose cells count the
  species with sequences at each locus. Defaults to 40.

## Value

A `ggplot` object, or a `patchwork` object with the sampling below when
`sampling` is given. Nothing is written.

## See also

[`plot_barcoding_identification_summary()`](https://beeamerino.github.io/PhyloCactus/reference/plot_barcoding_identification_summary.md),
[`report_barcoding_identification()`](https://beeamerino.github.io/PhyloCactus/reference/report_barcoding_identification.md).

## Examples

``` r
if (FALSE) { # \dontrun{
plot_barcoding_identification("11_barcoding/10_identify/REPORT_MN517611.1_answer.csv",
                              sampling = "11_barcoding/10_identify/REPORT_MN517611.1_sampling.csv")
} # }
```
