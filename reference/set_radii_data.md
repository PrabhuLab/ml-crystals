# Set or Reset a Custom Atomic Radii Table

Allows the user to provide their own table of atomic radii for the
current R session. This custom table will be used by functions like
`filter_ghost_distances`.

## Usage

``` r
set_radii_data(radii_data = NULL)
```

## Arguments

- radii_data:

  A `data.table` or `data.frame` containing the custom radii data, or
  `NULL` to reset to the package default.

## Value

Invisibly returns `NULL`. A message is printed to the console confirming
the action.

## Details

The provided `data.table` or `data.frame` must contain at least three
columns:

- `Symbol` (character): The chemical symbol of the element (e.g., "Si").

- `Radius` (numeric): The atomic radius in Angstroms (Å).

- `Type` (character): A descriptor for the radius type (e.g.,
  "covalent", "ionic").

This allows for storing multiple types of radii in the same table, which
can be selected using the `radii_type` argument in relevant functions.

To revert to using the package's default radii table, call the function
with `NULL` or without any arguments.

## See also

Other post-processing:
[`calculate_weighted_average_network_distance()`](https://prabhulab.github.io/ml-crystals/reference/calculate_weighted_average_network_distance.md),
[`filter_ghost_distances()`](https://prabhulab.github.io/ml-crystals/reference/filter_ghost_distances.md)

## Examples

``` r
# 1. Create a custom radii table with both covalent and ionic radii
my_radii <- data.table::data.table(
  Symbol = c("Si", "O", "O"),
  Radius = c(1.11, 0.66, 1.40),
  Type   = c("covalent", "covalent", "ionic")
)

# 2. Set this table for the current session
set_radii_data(my_radii)
#> Custom radii table has been set for the current R session.

# 3. Reset to the package's default covalent radii table
set_radii_data(NULL)
#> Session radii data reset to package default.
```
