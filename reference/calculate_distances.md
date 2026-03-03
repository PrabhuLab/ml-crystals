# Calculate Interatomic Distances

Computes the distances between a central set of atoms and an expanded
set, using the metric tensor for accuracy.

## Usage

``` r
calculate_distances(
  atomic_coordinates,
  expanded_coords,
  unit_cell_metrics,
  tolerance = 1e-06
)
```

## Arguments

- atomic_coordinates:

  A `data.table` of the primary (asymmetric) atom set.

- expanded_coords:

  A `data.table` of atoms in the expanded supercell.

- unit_cell_metrics:

  A `data.table` with cell parameters.

- tolerance:

  A numeric cutoff (default 1e-6). Distances smaller than this value are
  considered floating-point noise (overlapping atoms) and are filtered
  out.

## Value

A `data.table` of all non-zero distances.

## See also

Other property calculators:
[`calculate_angles()`](https://prabhulab.github.io/ml-crystals/reference/calculate_angles.md),
[`calculate_neighbor_counts()`](https://prabhulab.github.io/ml-crystals/reference/calculate_neighbor_counts.md),
[`calculate_weighted_neighbor_counts()`](https://prabhulab.github.io/ml-crystals/reference/calculate_weighted_neighbor_counts.md),
[`filter_atoms_by_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_atoms_by_symbol.md),
[`filter_by_elements()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_elements.md),
[`filter_by_wyckoff_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_wyckoff_symbol.md)

## Examples

``` r
ac <- data.table::data.table(Label = "Si1", x_a = 0, y_b = 0, z_c = 0)
ec <- data.table::data.table(Label = "O1_1", x_a = 0.5, y_b = 0.5, z_c = 0.5)
uc <- data.table::data.table(`_cell_length_a` = 10, `_cell_length_b` = 10,
                             `_cell_length_c` = 10, `_cell_angle_alpha` = 90,
                             `_cell_angle_beta` = 90, `_cell_angle_gamma` = 90)
calculate_distances(ac, ec, uc)
#>     Atom1  Atom2 Distance DeltaX DeltaY DeltaZ
#>    <char> <char>    <num>  <num>  <num>  <num>
#> 1:    Si1   O1_1 8.660254   -0.5   -0.5   -0.5
```
