# Apply Symmetry Operations to Generate a Full Unit Cell

Generates all symmetry-equivalent atomic positions within the unit cell.

## Usage

``` r
apply_symmetry_operations(
  atomic_coordinates,
  symmetry_operations,
  unit_cell_metrics,
  tolerance = 1e-04
)
```

## Arguments

- atomic_coordinates:

  A `data.table` of asymmetric atoms.

- symmetry_operations:

  A `data.table` of operations.

- unit_cell_metrics:

  A `data.table` of cell parameters (needed for distance merging).

- tolerance:

  Numeric. The distance tolerance in Angstroms for merging close atoms.
  Default is 1e-4.

## Value

A `data.table` containing all unique atomic positions in the unit cell.

## See also

Other coordinate processors:
[`calculate_expansion_factors()`](https://prabhulab.github.io/ml-crystals/reference/calculate_expansion_factors.md),
[`expand_transformed_coords()`](https://prabhulab.github.io/ml-crystals/reference/expand_transformed_coords.md)
