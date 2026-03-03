# Extract Database Code from CIF Content

Extracts the database code identifier (e.g., from `_database_code_`)
from the CIF.

## Usage

``` r
extract_database_code(cif_content)
```

## Arguments

- cif_content:

  A `data.table` containing the lines of a CIF file.

## Value

A character string of the database code, or `NA` if not found.

## See also

Other extractors:
[`extract_atomic_coordinates()`](https://prabhulab.github.io/ml-crystals/reference/extract_atomic_coordinates.md),
[`extract_chemical_formula()`](https://prabhulab.github.io/ml-crystals/reference/extract_chemical_formula.md),
[`extract_space_group_name()`](https://prabhulab.github.io/ml-crystals/reference/extract_space_group_name.md),
[`extract_space_group_number()`](https://prabhulab.github.io/ml-crystals/reference/extract_space_group_number.md),
[`extract_structure_type()`](https://prabhulab.github.io/ml-crystals/reference/extract_structure_type.md),
[`extract_symmetry_operations()`](https://prabhulab.github.io/ml-crystals/reference/extract_symmetry_operations.md),
[`extract_unit_cell_metrics()`](https://prabhulab.github.io/ml-crystals/reference/extract_unit_cell_metrics.md)

## Examples

``` r
cif_path <- system.file("extdata", "1590946.cif", package = "crystract")
if (file.exists(cif_path)) {
  cif_content <- read_cif_files(cif_path)[[1]]
  extract_database_code(cif_content)
}
#> [1] "depnum_ccdc_archive CCDC 1590946"
```
