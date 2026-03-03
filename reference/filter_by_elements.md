# Filter Distances by Element Symbols

Removes any distance pair where at least one of the atoms corresponds to
a specified list of element symbols.

## Usage

``` r
filter_by_elements(distances, atomic_coordinates, elements_to_exclude)
```

## Arguments

- distances:

  A `data.table` of interatomic distances.

- atomic_coordinates:

  A `data.table` of atomic coordinates.

- elements_to_exclude:

  A character vector of element symbols to exclude.

## Value

A `data.table` of distances with the specified elements removed.

## See also

Other property calculators:
[`calculate_angles()`](https://prabhulab.github.io/ml-crystals/reference/calculate_angles.md),
[`calculate_distances()`](https://prabhulab.github.io/ml-crystals/reference/calculate_distances.md),
[`calculate_neighbor_counts()`](https://prabhulab.github.io/ml-crystals/reference/calculate_neighbor_counts.md),
[`calculate_weighted_neighbor_counts()`](https://prabhulab.github.io/ml-crystals/reference/calculate_weighted_neighbor_counts.md),
[`filter_atoms_by_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_atoms_by_symbol.md),
[`filter_by_wyckoff_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_wyckoff_symbol.md)

## Examples

``` r
dists <- data.table::data.table(Atom1 = "Si1_1_0_0", Atom2 = "O1_1_0_0", Distance = 1.6)
ac <- data.table::data.table(Label = c("Si1", "O1"))
filter_by_elements(dists, ac, elements_to_exclude = "O")
#> Empty data.table (0 rows and 3 cols): Atom1,Atom2,Distance
```
