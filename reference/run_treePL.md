# Estimate Divergence Times via treePL Wrapper

Interfacing with \`treePL\` (Sanderson, 2002; Smith & O'Meara, 2012) via
a shell wrapper script to estimate ultrametric chronograms under
penalized likelihood. Integrates molecular branch lengths with temporal
calibration bounds.

## Usage

``` r
run_treePL(cfg, treefile, label, wrapper_sh)
```

## Arguments

- cfg:

  Character. Path to \`treePL\` configuration text file specifying
  calibrations and parameters.

- treefile:

  Character. Path to input input phylogenetic tree file (Newick format).

- label:

  Character. Descriptive run label identifier.

- wrapper_sh:

  Character. Path to the \`treePL\` shell wrapper script
  (\`treepl_wrapper_v1.sh\`).

## Value

Invisible NULL upon system command execution.

## References

Sanderson, M. J. (2002). Estimating absolute rates of molecular
evolution and divergence times: a penalized likelihood approach.
\*Molecular Biology and Evolution\*, 19(1), 101-109.
[doi:10.1093/oxfordjournals.molbev.a003974](https://doi.org/10.1093/oxfordjournals.molbev.a003974)

Smith, S. A., & O’Meara, B. C. (2012). treePL: divergence time
estimation using penalized likelihood for large phylogenies.
\*Bioinformatics\*, 28(20), 2689-2690.
[doi:10.1093/bioinformatics/bts492](https://doi.org/10.1093/bioinformatics/bts492)
