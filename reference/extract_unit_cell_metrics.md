# Extract Unit Cell Metrics

Parses the unit cell parameters (lengths a, b, c and angles alpha, beta,
gamma) and their associated standard uncertainties from CIF content.

## Usage

``` r
extract_unit_cell_metrics(cif_content)
```

## Arguments

- cif_content:

  A `data.table` containing the lines of a CIF file.

## Value

A one-row `data.table` with columns for each cell parameter and its
error.

## See also

Other extractors:
[`extract_atomic_coordinates()`](https://prabhulab.github.io/ml-crystals/reference/extract_atomic_coordinates.md),
[`extract_chemical_formula()`](https://prabhulab.github.io/ml-crystals/reference/extract_chemical_formula.md),
[`extract_database_code()`](https://prabhulab.github.io/ml-crystals/reference/extract_database_code.md),
[`extract_space_group_name()`](https://prabhulab.github.io/ml-crystals/reference/extract_space_group_name.md),
[`extract_space_group_number()`](https://prabhulab.github.io/ml-crystals/reference/extract_space_group_number.md),
[`extract_structure_type()`](https://prabhulab.github.io/ml-crystals/reference/extract_structure_type.md),
[`extract_symmetry_operations()`](https://prabhulab.github.io/ml-crystals/reference/extract_symmetry_operations.md)

## Examples

``` r
cif_path <- system.file("extdata", "1590946.cif", package = "crystract")
if (file.exists(cif_path)) {
  cif_content <- read_cif_files(cif_path)[[1]]
  extract_unit_cell_metrics(cif_content)
}
#>    _cell_length_a _cell_length_b _cell_length_c _cell_angle_alpha
#>             <num>          <num>          <num>             <num>
#> 1:           8.11           5.15           9.54                90
#>    _cell_angle_beta _cell_angle_gamma _cell_length_a_error _cell_length_b_error
#>               <num>             <num>                <num>                <num>
#> 1:               90                90                   NA                   NA
#>    _cell_length_c_error _cell_angle_alpha_error _cell_angle_beta_error
#>                   <num>                   <num>                  <num>
#> 1:                   NA                      NA                     NA
#>    _cell_angle_gamma_error
#>                      <num>
#> 1:                      NA
```
