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

## Examples

``` r
cif_path <- system.file("extdata", "1590946.cif", package = "crystract")
if (file.exists(cif_path)) {
  cif_content <- read_cif_files(cif_path)[[1]]
  extract_atomic_coordinates(cif_content)
}
#>     Label TypeSymbol WyckoffSymbol WyckoffMultiplicity OxidationState
#>    <char>     <char>        <char>               <num>          <num>
#> 1:    Sr1       Sr0+          <NA>                  NA             NA
#> 2:    Sr2       Sr0+          <NA>                  NA             NA
#> 3:    Si1       Si0+          <NA>                  NA             NA
#>    ThermalParam Occupancy OccupancyError    x_a   y_b    z_c x_error y_error
#>           <num>     <num>          <num>  <num> <num>  <num>   <num>   <num>
#> 1:           NA         1             NA 0.6529  0.25 0.0769       0       0
#> 2:           NA         1             NA 0.5192  0.25 0.6748       0       0
#> 3:           NA         1             NA 0.2539  0.25 0.1028       0       0
#>    z_error
#>      <num>
#> 1:       0
#> 2:       0
#> 3:       0
```
