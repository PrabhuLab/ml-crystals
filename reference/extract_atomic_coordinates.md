# Extract Atomic Coordinates

Parses atomic site info, including labels, fractional coordinates,
occupancies, oxidation states, and thermal parameters. It includes logic
to extract oxidation states from both site tags and type loops.

## Usage

``` r
extract_atomic_coordinates(cif_content, chemical_formula = NA)
```

## Arguments

- cif_content:

  A `data.table` containing the lines of a CIF file.

- chemical_formula:

  The chemical formula string for validation.

## Value

A `data.table` with atomic coordinate data, or `NULL` if not found.

## See also

Other extractors:
[`extract_chemical_formula()`](https://prabhulab.github.io/ml-crystals/reference/extract_chemical_formula.md),
[`extract_database_code()`](https://prabhulab.github.io/ml-crystals/reference/extract_database_code.md),
[`extract_space_group_name()`](https://prabhulab.github.io/ml-crystals/reference/extract_space_group_name.md),
[`extract_space_group_number()`](https://prabhulab.github.io/ml-crystals/reference/extract_space_group_number.md),
[`extract_structure_type()`](https://prabhulab.github.io/ml-crystals/reference/extract_structure_type.md),
[`extract_symmetry_operations()`](https://prabhulab.github.io/ml-crystals/reference/extract_symmetry_operations.md),
[`extract_unit_cell_metrics()`](https://prabhulab.github.io/ml-crystals/reference/extract_unit_cell_metrics.md)
