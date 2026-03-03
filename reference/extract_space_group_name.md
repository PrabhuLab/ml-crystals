# Extract Space Group Name from CIF Content

Extracts the Hermann-Mauguin space group name (e.g., from
`_space_group_name_H-M_alt`).

## Usage

``` r
extract_space_group_name(cif_content)
```

## Arguments

- cif_content:

  A `data.table` containing the lines of a CIF file.

## Value

A character string of the space group name, or `NA` if not found.

## See also

Other extractors:
[`extract_atomic_coordinates()`](https://prabhulab.github.io/ml-crystals/reference/extract_atomic_coordinates.md),
[`extract_chemical_formula()`](https://prabhulab.github.io/ml-crystals/reference/extract_chemical_formula.md),
[`extract_database_code()`](https://prabhulab.github.io/ml-crystals/reference/extract_database_code.md),
[`extract_space_group_number()`](https://prabhulab.github.io/ml-crystals/reference/extract_space_group_number.md),
[`extract_structure_type()`](https://prabhulab.github.io/ml-crystals/reference/extract_structure_type.md),
[`extract_symmetry_operations()`](https://prabhulab.github.io/ml-crystals/reference/extract_symmetry_operations.md),
[`extract_unit_cell_metrics()`](https://prabhulab.github.io/ml-crystals/reference/extract_unit_cell_metrics.md)

## Examples

``` r
cif_path <- system.file("extdata", "1590946.cif", package = "crystract")
if (file.exists(cif_path)) {
  cif_content <- read_cif_files(cif_path)[[1]]
  extract_space_group_name(cif_content)
}
#> [1] "P n m a"
```
