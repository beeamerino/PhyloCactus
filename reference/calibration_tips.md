# Resolve the Terminals a Calibration Row Addresses

`treePL` locates a calibrated node by the MRCA of the terminals declared
for it, so the taxon set, not the row label, is what determines where a
bound lands. `value` accepts a semicolon-separated list of names, which
the two stem calibrations require: a stem node is shared by the families
it subtends and cannot be addressed by a single family name.

## Usage

``` r
calibration_tips(constraints, column, value, tip_labels)
```

## Arguments

- constraints:

  Data frame of the taxonomic table (`cactus_constraints.csv`), carrying
  a `Specie_name` column and one column per rank used by the calibration
  table.

- column:

  Character. Name of the rank column the calibration row keys on.

- value:

  Character. One name, or several separated by `;`.

- tip_labels:

  Character vector of terminals present in the tree.

## Value

Character vector of terminals, intersected with `tip_labels`.

## Examples

``` r
constraints <- data.frame(
  Specie_name = c("Opuntia_stricta", "Pereskia_aculeata", "Anacampseros_retusa"),
  Family = c("Cactaceae", "Cactaceae", "Anacampserotaceae")
)
calibration_tips(constraints, "Family", "Cactaceae;Anacampserotaceae",
                 constraints$Specie_name)
#> [1] "Opuntia_stricta"     "Pereskia_aculeata"   "Anacampseros_retusa"
```
