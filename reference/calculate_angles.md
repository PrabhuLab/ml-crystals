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
[`filter_atoms_by_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_atoms_by_symbol.md),
[`filter_by_elements()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_elements.md),
[`filter_by_wyckoff_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_wyckoff_symbol.md)
