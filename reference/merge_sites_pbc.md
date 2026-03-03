# Merge Close Atoms (Pymatgen Style)

Merges atomic sites that are within a specified distance tolerance,
accounting for Periodic Boundary Conditions (PBC). This mimics
pymatgen's `Structure.merge_sites` logic: it performs hierarchical
clustering and then averages the coordinates of the merged sites, taking
care of PBC wrapping.

## Usage

``` r
merge_sites_pbc(atomic_coordinates, unit_cell_metrics, tol = 1e-04)
```

## Arguments

- atomic_coordinates:

  A data.table with columns `Label`, `x_a`, `y_b`, `z_c`.

- unit_cell_metrics:

  A data.table containing cell lengths and angles.

- tol:

  Numeric. Distance tolerance in Angstroms. Sites closer than this will
  be merged. Default is 1e-4.

## Value

A deduplicated data.table with coordinates averaged.

## Examples

``` r
ac <- data.table::data.table(Label = c("A", "B"),
                             x_a = c(0, 0.00001),
                             y_b = c(0, 0),
                             z_c = c(0, 0))
uc <- data.table::data.table(`_cell_length_a` = 10, `_cell_length_b` = 10,
                             `_cell_length_c` = 10, `_cell_angle_alpha` = 90,
                             `_cell_angle_beta` = 90, `_cell_angle_gamma` = 90)
merge_sites_pbc(ac, uc, tol = 1e-4)
#>     Label   x_a   y_b   z_c
#>    <char> <num> <num> <num>
#> 1:      A     0     0     0
```
