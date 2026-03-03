# Calculate Bond Angles

Calculates all bond angles centered on each atom, formed by pairs of its
bonded neighbors.

## Usage

``` r
calculate_angles(
  bonded_pairs,
  atomic_coordinates,
  expanded_coords,
  unit_cell_metrics
)
```

## Arguments

- bonded_pairs:

  Data.table of bonded atoms.

- atomic_coordinates:

  Data.table of asymmetric atom coordinates.

- expanded_coords:

  Data.table of supercell atom coordinates.

- unit_cell_metrics:

  Data.table with unit cell parameters.

## Value

A `data.table` of all unique bond angles.

## See also

Other property calculators:
[`calculate_distances()`](https://prabhulab.github.io/ml-crystals/reference/calculate_distances.md),
[`calculate_neighbor_counts()`](https://prabhulab.github.io/ml-crystals/reference/calculate_neighbor_counts.md),
[`calculate_weighted_neighbor_counts()`](https://prabhulab.github.io/ml-crystals/reference/calculate_weighted_neighbor_counts.md),
[`filter_atoms_by_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_atoms_by_symbol.md),
[`filter_by_elements()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_elements.md),
[`filter_by_wyckoff_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_wyckoff_symbol.md)

## Examples

``` r
bp <- data.table::data.table(Atom1 = c("Si1", "Si1"), Atom2 = c("O1", "O2"))
ac <- data.table::data.table(Label = "Si1", x_a = 0, y_b = 0, z_c = 0)
ec <- data.table::data.table(Label = c("O1", "O2"),
                             x_a = c(1, 0), y_b = c(0, 1), z_c = c(0, 0))
uc <- data.table::data.table(`_cell_length_a` = 10, `_cell_length_b` = 10,
                             `_cell_length_c` = 10, `_cell_angle_alpha` = 90,
                             `_cell_angle_beta` = 90, `_cell_angle_gamma` = 90)
calculate_angles(bp, ac, ec, uc)
#> Key: <CentralAtom>
#>    CentralAtom Neighbor1 Neighbor2 Angle
#>         <char>    <char>    <char> <num>
#> 1:         Si1        O1        O2    90
```
