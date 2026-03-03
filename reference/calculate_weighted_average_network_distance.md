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

## Examples

``` r
dists <- data.table::data.table(Atom1 = "Si1", Atom2 = "O1", Distance = 1.6)
ac <- data.table::data.table(Label = c("Si1", "O1"),
                             WyckoffMultiplicity = c(4, 4),
                             WyckoffSymbol = c("c", "c"),
                             Occupancy = c(1, 1))
calculate_weighted_average_network_distance(dists, ac, wyckoff_symbols = "4c")
#> [1] 1.6
```
