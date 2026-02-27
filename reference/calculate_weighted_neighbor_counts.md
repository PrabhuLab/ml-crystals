# Calculate Weighted Coordination Numbers

Counts the weighted number of nearest neighbors for each central atom
based on a table of bonded pairs, accounting for the fractional
occupancy of the neighbor sites. This is particularly useful for
analyzing disordered crystal structures.

If all occupancies are 1.0, this efficiently returns the standard
coordination number. If the fractional occupancies on any single
crystallographic site sum to greater than 1.0, this indicates an invalid
or physically impossible structure; the function will throw a warning
and return `NULL` so the file can be discarded by the analysis pipeline.

## Usage

``` r
calculate_weighted_neighbor_counts(bonded_pairs_table, atomic_coordinates)
```

## Arguments

- bonded_pairs_table:

  A `data.table` of bonded pairs (e.g., from `minimum_distance`).

- atomic_coordinates:

  A `data.table` of asymmetric atoms from `extract_atomic_coordinates`.

## Value

A `data.table` with columns 'Atom', 'CoordinationNumber', and
'WeightedCoordinationNumber'. Returns `NULL` if the occupancies sum to
\> 1 on any shared site.

## See also

Other property calculators:
[`calculate_angles()`](https://prabhulab.github.io/ml-crystals/reference/calculate_angles.md),
[`calculate_distances()`](https://prabhulab.github.io/ml-crystals/reference/calculate_distances.md),
[`calculate_neighbor_counts()`](https://prabhulab.github.io/ml-crystals/reference/calculate_neighbor_counts.md),
[`filter_atoms_by_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_atoms_by_symbol.md),
[`filter_by_elements()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_elements.md),
[`filter_by_wyckoff_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_wyckoff_symbol.md)
