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
