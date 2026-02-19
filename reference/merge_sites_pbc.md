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
