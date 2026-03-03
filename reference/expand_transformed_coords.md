# Expand Coordinates into a Supercell

Takes a set of atomic coordinates within a single unit cell and
replicates them into a supercell grid defined by expansion factors.

## Usage

``` r
expand_transformed_coords(transformed_coords, expansion_factors = c(1, 1, 1))
```

## Arguments

- transformed_coords:

  A `data.table` of atom positions within one unit cell.

- expansion_factors:

  A numeric vector `c(a, b, c)` specifying the number of images in each
  direction (+/-). Default `c(1,1,1)` creates a 3x3x3 grid.

## Value

A `data.table` containing all atomic positions in the expanded
supercell.

## See also

Other coordinate processors:
[`apply_symmetry_operations()`](https://prabhulab.github.io/ml-crystals/reference/apply_symmetry_operations.md),
[`calculate_expansion_factors()`](https://prabhulab.github.io/ml-crystals/reference/calculate_expansion_factors.md)

## Examples

``` r
tc <- data.table::data.table(Label = "A", x_a = 0, y_b = 0, z_c = 0)
expand_transformed_coords(tc, expansion_factors = c(1, 1, 1))
#>          Label   x_a   y_b   z_c    Tx    Ty    Tz
#>         <char> <num> <num> <num> <int> <int> <int>
#>  1: A_-1_-1_-1    -1    -1    -1    -1    -1    -1
#>  2:  A_0_-1_-1     0    -1    -1     0    -1    -1
#>  3:  A_1_-1_-1     1    -1    -1     1    -1    -1
#>  4:  A_-1_0_-1    -1     0    -1    -1     0    -1
#>  5:   A_0_0_-1     0     0    -1     0     0    -1
#>  6:   A_1_0_-1     1     0    -1     1     0    -1
#>  7:  A_-1_1_-1    -1     1    -1    -1     1    -1
#>  8:   A_0_1_-1     0     1    -1     0     1    -1
#>  9:   A_1_1_-1     1     1    -1     1     1    -1
#> 10:  A_-1_-1_0    -1    -1     0    -1    -1     0
#> 11:   A_0_-1_0     0    -1     0     0    -1     0
#> 12:   A_1_-1_0     1    -1     0     1    -1     0
#> 13:   A_-1_0_0    -1     0     0    -1     0     0
#> 14:    A_0_0_0     0     0     0     0     0     0
#> 15:    A_1_0_0     1     0     0     1     0     0
#> 16:   A_-1_1_0    -1     1     0    -1     1     0
#> 17:    A_0_1_0     0     1     0     0     1     0
#> 18:    A_1_1_0     1     1     0     1     1     0
#> 19:  A_-1_-1_1    -1    -1     1    -1    -1     1
#> 20:   A_0_-1_1     0    -1     1     0    -1     1
#> 21:   A_1_-1_1     1    -1     1     1    -1     1
#> 22:   A_-1_0_1    -1     0     1    -1     0     1
#> 23:    A_0_0_1     0     0     1     0     0     1
#> 24:    A_1_0_1     1     0     1     1     0     1
#> 25:   A_-1_1_1    -1     1     1    -1     1     1
#> 26:    A_0_1_1     0     1     1     0     1     1
#> 27:    A_1_1_1     1     1     1     1     1     1
#>          Label   x_a   y_b   z_c    Tx    Ty    Tz
#>         <char> <num> <num> <num> <int> <int> <int>
```
