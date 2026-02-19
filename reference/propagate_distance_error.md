# Propagate Distance Error

Calculates the standard uncertainty for each interatomic distance.

## Usage

``` r
propagate_distance_error(bonded_pairs, atomic_coordinates, unit_cell_metrics)
```

## Arguments

- bonded_pairs:

  Data.table of bonded atoms with their distances.

- atomic_coordinates:

  Data.table with fractional coordinates and errors.

- unit_cell_metrics:

  Data.table with unit cell parameters and errors.

## Value

The input `bonded_pairs` data.table with a new 'DistanceError' column.

## See also

Other error propagators:
[`propagate_angle_error()`](https://prabhulab.github.io/ml-crystals/reference/propagate_angle_error.md)
