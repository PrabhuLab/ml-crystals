# Filter Data by Wyckoff Symbol

Filters a data table (e.g., bonds or angles) to include only entries
where a specified atom occupies one of the given Wyckoff sites.

## Usage

``` r
filter_by_wyckoff_symbol(
  data_table,
  atomic_coordinates,
  atom_col,
  wyckoff_symbols
)
```

## Arguments

- data_table:

  A `data.table` object, such as one produced by `minimum_distance` or
  `calculate_angles`.

- atomic_coordinates:

  A `data.table` from `extract_atomic_coordinates` containing
  `WyckoffMultiplicity` and `WyckoffSymbol` columns.

- atom_col:

  A character string specifying the column in `data_table` that contains
  the atom labels to filter by (e.g., "Atom1", "CentralAtom").

- wyckoff_symbols:

  A character vector of the full Wyckoff symbols to keep (e.g.,
  `c("4c")` or `c("6c", "16i", "24k")`).

## Value

A `data.table` filtered to include only rows where the specified atom
occupies one of the desired Wyckoff sites.

## Details

This function is designed to facilitate analysis based on
crystallographic site symmetry. It is particularly useful for analyzing
complex structures like clathrates where different atoms perform
distinct structural roles based on their site symmetry.

The function works by:

1.  Creating a full Wyckoff label (e.g., "4c", "24k") by combining the
    `WyckoffMultiplicity` and `WyckoffSymbol` columns from the
    `atomic_coordinates` table.

2.  Identifying the parent atom for each entry in the `data_table`.

3.  Merging this Wyckoff information into the results table.

4.  Filtering to keep only rows where the atom's full Wyckoff label
    matches one of the symbols provided in the `wyckoff_symbols` vector.

## See also

Other property calculators:
[`calculate_angles()`](https://prabhulab.github.io/ml-crystals/reference/calculate_angles.md),
[`calculate_distances()`](https://prabhulab.github.io/ml-crystals/reference/calculate_distances.md),
[`calculate_neighbor_counts()`](https://prabhulab.github.io/ml-crystals/reference/calculate_neighbor_counts.md),
[`filter_atoms_by_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_atoms_by_symbol.md),
[`filter_by_elements()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_elements.md)

## Examples

``` r
cif_file <- system.file("extdata", "ICSD422.cif", package = "crystract")
if (file.exists(cif_file)) {
  # 1. Perform a standard analysis to get bond and coordinate tables
  cif_content <- read_cif_files(cif_file)[[1]]
  atoms <- extract_atomic_coordinates(cif_content)
  metrics <- extract_unit_cell_metrics(cif_content)
  sym_ops <- extract_symmetry_operations(cif_content)
  full_cell <- apply_symmetry_operations(atoms, sym_ops)
  super_cell <- expand_transformed_coords(full_cell)
  dists <- calculate_distances(atoms, super_cell, metrics)
  bonds <- minimum_distance(dists)

  # 2. In ICSD422, all atoms are on the "4c" site.
  print("Original atomic coordinates showing Wyckoff sites:")
  print(atoms[, .(Label, WyckoffSymbol, WyckoffMultiplicity)])

  filtered_bonds <- filter_by_wyckoff_symbol(
    data_table = bonds,
    atomic_coordinates = atoms,
    atom_col = "Atom1",
    wyckoff_symbols = "4c"
  )

  cat("\nNumber of bonds in original table:", nrow(bonds), "\n")
  cat("Number of bonds after filtering for '4c' site:", nrow(filtered_bonds), "\n")
}
```
