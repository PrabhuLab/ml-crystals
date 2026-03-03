# Extract Symmetry Operations

Parses the symmetry operation definitions from the CIF content.

## Usage

``` r
extract_symmetry_operations(cif_content)
```

## Arguments

- cif_content:

  A `data.table` containing the lines of a CIF file.

## Value

A `data.table` with symmetry operations. Returns `NULL` if not found.

## See also

Other extractors:
[`extract_atomic_coordinates()`](https://prabhulab.github.io/ml-crystals/reference/extract_atomic_coordinates.md),
[`extract_chemical_formula()`](https://prabhulab.github.io/ml-crystals/reference/extract_chemical_formula.md),
[`extract_database_code()`](https://prabhulab.github.io/ml-crystals/reference/extract_database_code.md),
[`extract_space_group_name()`](https://prabhulab.github.io/ml-crystals/reference/extract_space_group_name.md),
[`extract_space_group_number()`](https://prabhulab.github.io/ml-crystals/reference/extract_space_group_number.md),
[`extract_structure_type()`](https://prabhulab.github.io/ml-crystals/reference/extract_structure_type.md),
[`extract_unit_cell_metrics()`](https://prabhulab.github.io/ml-crystals/reference/extract_unit_cell_metrics.md)

## Examples

``` r
cif_path <- system.file("extdata", "1590946.cif", package = "crystract")
if (file.exists(cif_path)) {
  cif_content <- read_cif_files(cif_path)[[1]]
  extract_symmetry_operations(cif_content)
}
#>         x      y      z
#>    <char> <char> <char>
#> 1:  x+1/2      y -z+1/2
#> 2:      x -y+1/2      z
#> 3: -x+1/2  y+1/2  z+1/2
#> 4:     -x     -y     -z
#> 5: -x+1/2     -y  z+1/2
#> 6:     -x  y+1/2     -z
#> 7:  x+1/2 -y+1/2 -z+1/2
#> 8:      x      y      z
```
