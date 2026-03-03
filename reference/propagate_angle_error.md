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

## Examples

``` r
# 1. Create dummy bond angles
ba <- data.table::data.table(
  CentralAtom = "Si1", Neighbor1 = "O1", Neighbor2 = "O2", Angle = 109.5
)

# 2. Create dummy atomic coordinates with errors
ac <- data.table::data.table(
  Label = c("Si1", "O1", "O2"),
  x_a = c(0, 0.1, 0), y_b = c(0, 0, 0.1), z_c = c(0, 0, 0),
  x_error = c(0.01, 0.01, 0.01),
  y_error = c(0.01, 0.01, 0.01),
  z_error = c(0.01, 0.01, 0.01)
)

# 3. Create dummy expanded coordinates
ec <- data.table::data.table(
  Label = c("O1", "O2"),
  x_a = c(0.1, 0), y_b = c(0, 0.1), z_c = c(0, 0)
)

# 4. Create dummy unit cell metrics
uc <- data.table::data.table(
  `_cell_length_a` = 10, `_cell_length_a_error` = 0.1,
  `_cell_length_b` = 10, `_cell_length_b_error` = 0.1,
  `_cell_length_c` = 10, `_cell_length_c_error` = 0.1,
  `_cell_angle_alpha` = 90, `_cell_angle_alpha_error` = 0,
  `_cell_angle_beta` = 90, `_cell_angle_beta_error` = 0,
  `_cell_angle_gamma` = 90, `_cell_angle_gamma_error` = 0
)

# 5. Run the error propagation
propagate_angle_error(ba, ac, ec, uc)
#> Key: <CentralAtom, Neighbor1, Neighbor2>
#>    CentralAtom Neighbor1 Neighbor2 Angle AngleError
#>         <char>    <char>    <char> <num>      <num>
#> 1:         Si1        O1        O2 109.5   11.45916
```
