# Filter Ghost Distances Using Atomic Radii

Cleans a distance table by removing physically implausible distances. It
uses a table of atomic radii to establish a plausible bond length range
for each atom pair. Any calculated distance falling outside this range
(defined by a `margin`) is considered a "ghost" distance and is removed.
This is particularly useful for cleaning data from disordered crystal
structures.

## Usage

``` r
filter_ghost_distances(
  distances,
  atomic_coordinates,
  margin = 0.1,
  radii_type = "covalent"
)
```

## Arguments

- distances:

  A `data.table` of interatomic distances, typically from
  `calculate_distances`. Must contain `Atom1`, `Atom2`, and `Distance`.

- atomic_coordinates:

  A `data.table` of asymmetric atoms from `extract_atomic_coordinates`.
  Used to link atom labels to element types.

- margin:

  A numeric value (default 0.1) specifying the tolerance. A distance `d`
  between atoms with radii `r1` and `r2` is kept if
  `(r1+r2)*(1-margin) <= d <= (r1+r2)*(1+margin)`.

- radii_type:

  A character string specifying the type of radius to use for the
  calculation (e.g., `"covalent"`, `"ionic"`). This value must
  correspond to an entry in the `Type` column of the active radii table.
  Defaults to `"covalent"`. The radii table can be customized for the
  session using
  [`set_radii_data()`](https://prabhulab.github.io/ml-crystals/reference/set_radii_data.md).

## Value

A list containing two `data.table`s:

- kept:

  The distances considered physically plausible.

- removed:

  The "ghost" distances that were filtered out, with columns explaining
  the reason for removal.

## See also

Other post-processing:
[`calculate_weighted_average_network_distance()`](https://prabhulab.github.io/ml-crystals/reference/calculate_weighted_average_network_distance.md),
[`set_radii_data()`](https://prabhulab.github.io/ml-crystals/reference/set_radii_data.md)

## Examples

``` r
# Create minimal dummy data for demonstration
distances <- data.table::data.table(
  Atom1 = c("Si1_1_0_0", "O1_1_0_0"),
  Atom2 = c("O1_1_0_0", "Si1_1_0_0"),
  Distance = c(1.6, 0.5) # 0.5 is implausibly short
)

atomic_coords <- data.table::data.table(
  Label = c("Si1", "O1")
)

# Run the filter
result <- filter_ghost_distances(distances, atomic_coords, margin = 0.1)

print(result$kept)
#> Empty data.table (0 rows and 3 cols): Atom1,Atom2,Distance
print(result$removed)
#>        Atom1     Atom2 Distance expected_dist lower_bound upper_bound
#>       <char>    <char>    <num>         <num>       <num>       <num>
#> 1: Si1_1_0_0  O1_1_0_0      1.6          1.83       1.647       2.013
#> 2:  O1_1_0_0 Si1_1_0_0      0.5          1.83       1.647       2.013
#>                   Reason
#>                   <char>
#> 1: Distance is TOO SHORT
#> 2: Distance is TOO SHORT
```
