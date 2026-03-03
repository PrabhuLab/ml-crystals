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

## Examples

``` r
ac <- data.table::data.table(Label = "A", x_a = 0, y_b = 0, z_c = 0)
ops <- data.table::data.table(x = c("x", "-x"), y = c("y", "-y"), z = c("z", "-z"))
uc <- data.table::data.table(`_cell_length_a` = 10, `_cell_length_b` = 10,
                             `_cell_length_c` = 10, `_cell_angle_alpha` = 90,
                             `_cell_angle_beta` = 90, `_cell_angle_gamma` = 90)
apply_symmetry_operations(ac, ops, uc)
#>     Label SymOp   x_a   y_b   z_c
#>    <char> <int> <num> <num> <num>
#> 1:    A_1     1     0     0     0
```
