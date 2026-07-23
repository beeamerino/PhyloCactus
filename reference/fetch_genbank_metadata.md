# Fetch GenBank Sequence Metadata via NCBI Entrez Utilities

Queries NCBI Entrez Utilities to retrieve sequence lengths, organism
taxonomy, publication titles, and accession IDs for a collection of
GenBank sequence identifiers (SIDs). Metadata retrieval enriches raw
sequence clusters with verifiable audit data.

## Usage

``` r
fetch_genbank_metadata(
  sids,
  cache_file,
  batch_size = 200,
  sleep_time = 0.5,
  max_retries = 5,
  force_download = FALSE
)
```

## Arguments

- sids:

  Character vector of GenBank Sequence Identifiers (SIDs) to query.

- cache_file:

  Character. File path to store and load cached metadata tables.

- batch_size:

  Integer. Number of sequence IDs requested per HTTP batch query.
  Defaults to \`200\`.

- sleep_time:

  Numeric. Pause duration in seconds between consecutive batch requests
  to respect NCBI rate limits. Defaults to \`0.5\`.

- max_retries:

  Integer. Maximum retry attempts permitted per batch before failing.
  Defaults to \`5\`.

- force_download:

  Logical. Force fresh Entrez queries instead of loading local cache?
  Defaults to \`FALSE\`.

## Value

A data frame containing fetched GenBank metadata fields for each
requested sequence ID.

## Examples

``` r
if (FALSE) { # \dontrun{
fetch_genbank_metadata(
  sids = c("AY123456", "AY123457"),
  cache_file = "cache_metadata.csv"
)
} # }
```
