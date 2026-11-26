# Fetch GenBank Sequence Metadata via NCBI Entrez Utilities

Retrieves the organism name and the definition line of each GenBank
sequence identifier (SID) with
[`ape::read.GenBank()`](https://rdrr.io/pkg/ape/man/read.GenBank.html),
in batches with retries, and caches the result.

## Usage

``` r
fetch_genbank_metadata(
  sids,
  cache_file,
  batch_size = 200,
  sleep_time = 0.5,
  max_retries = 5,
  force_download = FALSE,
  max_failed_fraction = 0.05
)
```

## Arguments

- sids:

  Character vector of GenBank Sequence Identifiers (SIDs) to query.

- cache_file:

  Character. File path to store and load cached metadata tables.

- batch_size:

  Integer. Number of sequence IDs requested per HTTP batch query.
  Defaults to `200`.

- sleep_time:

  Numeric. Pause duration in seconds between consecutive batch requests
  to respect NCBI rate limits. Defaults to `0.5`.

- max_retries:

  Integer. Maximum retry attempts permitted per batch before failing.
  Defaults to `5`.

- force_download:

  Deprecated and ignored: records already in `cache_file` are always
  reused. Defaults to `FALSE`; `TRUE` gives a warning. Delete
  `cache_file` to download every record again.

- max_failed_fraction:

  Numeric between 0 and 1. A batch whose `max_retries` attempts all fail
  is logged and its accessions counted; when the failed accessions
  exceed this fraction of those requested, the function stops instead of
  returning a table in which they would lack a species and description.
  Defaults to `0.05`.

## Value

A data frame with columns `sid`, `Species_gb` and `Description_gb`, with
the accessions of failed batches in its attribute `failed_sids`.

## Examples

``` r
if (FALSE) { # \dontrun{
fetch_genbank_metadata(
  sids = c("AY123456", "AY123457"),
  cache_file = "cache_metadata.csv"
)
} # }
```
