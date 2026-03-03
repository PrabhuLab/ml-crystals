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

## Examples

``` r
bp <- data.table::data.table(Atom1 = "Si1", Atom2 = "O1_1", Distance = 1.6,
                             DeltaX = 0.1, DeltaY = 0.1, DeltaZ = 0)
ac <- data.table::data.table(Label = c("Si1", "O1"),
                             x_error = c(0.01, 0.01),
                             y_error = c(0.01, 0.01),
                             z_error = c(0.01, 0.01))
uc <- data.table::data.table(`_cell_length_a` = 10, `_cell_length_a_error` = 0.1,
                             `_cell_length_b` = 10, `_cell_length_b_error` = 0.1,
                             `_cell_length_c` = 10, `_cell_length_c_error` = 0.1,
                             `_cell_angle_alpha` = 90, `_cell_angle_alpha_error` = 0,
                             `_cell_angle_beta` = 90, `_cell_angle_beta_error` = 0,
                             `_cell_angle_gamma` = 90, `_cell_angle_gamma_error` = 0)
propagate_distance_error(bp, ac, uc)
#> Key: <Atom1, Atom2, Distance>
#>     Atom1  Atom2 Distance DeltaX DeltaY DeltaZ DistanceError
#>    <char> <char>    <num>  <num>  <num>  <num>         <num>
#> 1:    Si1   O1_1      1.6    0.1    0.1      0     0.1253121
```
