# Calculate Weighted Average Network Bond Distance

Computes a single, representative bond length for a specified atomic
network. This function precisely implements the validated logic that
accounts for site multiplicity and occupancy of the central atoms.

## Usage

``` r
calculate_weighted_average_network_distance(
  distances,
  atomic_coordinates,
  wyckoff_symbols
)
```

## Arguments

- distances:

  A `data.table` of interatomic distances filtered to include **only
  bonded pairs** (e.g., from `minimum_distance`).

- atomic_coordinates:

  A `data.table` of asymmetric atoms from `extract_atomic_coordinates`.

- wyckoff_symbols:

  A character vector of Wyckoff symbols defining the atomic network
  (e.g., `c("6c", "16i", "24k")`). Must be the full symbol.

## Value

A single numeric value representing the weighted average bond distance.

## See also

Other post-processing:
[`filter_ghost_distances()`](https://prabhulab.github.io/ml-crystals/reference/filter_ghost_distances.md),
[`set_radii_data()`](https://prabhulab.github.io/ml-crystals/reference/set_radii_data.md)
