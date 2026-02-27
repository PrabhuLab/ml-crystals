# Calculate Coordination Numbers

Counts the number of nearest neighbors for each central atom based on a
table of bonded pairs.

## Usage

``` r
calculate_neighbor_counts(bonded_pairs_table)
```

## Arguments

- bonded_pairs_table:

  A `data.table` of bonded pairs.

## Value

A `data.table` with columns 'Atom' and 'CoordinationNumber'.

## See also

Other property calculators:
[`calculate_angles()`](https://prabhulab.github.io/ml-crystals/reference/calculate_angles.md),
[`calculate_distances()`](https://prabhulab.github.io/ml-crystals/reference/calculate_distances.md),
[`calculate_weighted_neighbor_counts()`](https://prabhulab.github.io/ml-crystals/reference/calculate_weighted_neighbor_counts.md),
[`filter_atoms_by_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_atoms_by_symbol.md),
[`filter_by_elements()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_elements.md),
[`filter_by_wyckoff_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_wyckoff_symbol.md)
