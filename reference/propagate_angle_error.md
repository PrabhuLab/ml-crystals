# Propagate Angle Error

Calculates the standard uncertainty for each bond angle.

## Usage

``` r
propagate_angle_error(
  bond_angles,
  atomic_coordinates,
  expanded_coords,
  unit_cell_metrics
)
```

## Arguments

- bond_angles:

  Data.table of calculated bond angles.

- atomic_coordinates:

  Data.table with fractional coordinates and errors.

- expanded_coords:

  Data.table of supercell atom coordinates.

- unit_cell_metrics:

  Data.table with unit cell parameters and errors.

## Value

The input `bond_angles` data.table with a new 'AngleError' column.

## See also

Other error propagators:
[`propagate_distance_error()`](https://prabhulab.github.io/ml-crystals/reference/propagate_distance_error.md)
