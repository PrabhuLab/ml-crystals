# Calculate Supercell Expansion Factors

Computes the number of unit cell repetitions required in each direction
(a, b, c) to encompass a sphere of a given radius.

## Usage

``` r
calculate_expansion_factors(unit_cell_metrics, radius)
```

## Arguments

- unit_cell_metrics:

  A `data.table` with cell parameters.

- radius:

  Numeric. The cutoff radius (e.g., 13.0 for Voronoi).

## Value

A named numeric vector `c(a=..., b=..., c=...)`.

## See also

Other coordinate processors:
[`apply_symmetry_operations()`](https://prabhulab.github.io/ml-crystals/reference/apply_symmetry_operations.md),
[`expand_transformed_coords()`](https://prabhulab.github.io/ml-crystals/reference/expand_transformed_coords.md)

## Examples

``` r
uc <- data.table::data.table(`_cell_length_a` = 10, `_cell_length_b` = 10,
                             `_cell_length_c` = 10, `_cell_angle_alpha` = 90,
                             `_cell_angle_beta` = 90, `_cell_angle_gamma` = 90)
calculate_expansion_factors(uc, radius = 13.0)
#> a b c 
#> 2 2 2 
```
