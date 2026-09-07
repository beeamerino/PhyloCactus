# Run treePL Executable Directly

Executes the `treePL` binary directly to estimate divergence times from
a prepared configuration file.

## Usage

``` r
run_treePL_direct(cfg_file, label, cwd = NULL)
```

## Arguments

- cfg_file:

  Character. Path to `treePL` configuration file.

- label:

  Character. Descriptive run label identifier.

- cwd:

  Character. Optional working directory context. Defaults to `NULL`.

## Value

Invisible NULL upon system command execution.

## Details

The tree path is read from the `treefile` line of `cfg_file` and checked
for rootedness before the binary is invoked, because `treePL` fails
opaquely on an unrooted tree. When the configuration declares no
`treefile`, or the declared path cannot be resolved, the check is
skipped with a warning rather than blocking the run.

## References

Smith, S. A., & O’Meara, B. C. (2012). treePL: divergence time
estimation using penalized likelihood for large phylogenies.
*Bioinformatics*, 28(20), 2689-2690.
[doi:10.1093/bioinformatics/bts492](https://doi.org/10.1093/bioinformatics/bts492)
