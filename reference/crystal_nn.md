# Identify Atomic Bonds using CrystalNN

Port of Pymatgen's `CrystalNN` weighting logic. Includes Porosity,
Electronegativity, and Distance penalties. Uses specific Decision Tree
logic for radii (Shannon -\> Covalent -\> Atomic) per pymatgen's
fallback logic.

## Usage

``` r
crystal_nn(
  distances,
  atomic_coordinates,
  expanded_coords,
  unit_cell_metrics,
  cutoff_length = 7,
  x_diff_weight = 3,
  porosity_adjustment = TRUE,
  distance_cutoffs = c(0.5, 1)
)
```

## Arguments

- distances:

  Ignored.

- atomic_coordinates:

  Primary atom set.

- expanded_coords:

  Expanded supercell.

- unit_cell_metrics:

  Cell parameters.

- cutoff_length:

  Voronoi cutoff. Default 7.0 matches Pymatgen CrystalNN default.

- x_diff_weight:

  Electronegativity weight (default 3.0).

- porosity_adjustment:

  Logical (default TRUE).

- distance_cutoffs:

  Vector c(0.5, 1.0).

## Value

A `data.table` of bonded pairs.

## See also

Other bonding algorithms:
[`brunner_nn_reciprocal()`](https://prabhulab.github.io/ml-crystals/reference/brunner_nn_reciprocal.md),
[`econ_nn()`](https://prabhulab.github.io/ml-crystals/reference/econ_nn.md),
[`minimum_distance()`](https://prabhulab.github.io/ml-crystals/reference/minimum_distance.md),
[`voronoi_nn()`](https://prabhulab.github.io/ml-crystals/reference/voronoi_nn.md)

## Examples

``` r
# \donttest{
cif_path <- system.file("extdata", "1590946.cif", package = "crystract")
if (file.exists(cif_path)) {
  res <- analyze_single_cif(cif_path, bonding_algorithms = "crystal_nn")
  print(head(res$bonds_crystal_nn[[1]]))
}
#> Key: <Atom1, Atom2, Distance>
#>     Atom1         Atom2 Distance    DeltaX        DeltaY    DeltaZ Weight
#>    <char>        <char>    <num>     <num>         <num>     <num>  <num>
#> 1:  Si1_2   Sr1_1_0_0_0 3.163544 -0.819110  2.220446e-16  3.055662  1.000
#> 2:  Si1_2   Sr1_2_0_0_0 3.245310  3.235890  0.000000e+00 -0.247086  0.914
#> 3:  Si1_2 Sr1_4_0_-1_-1 3.184477  0.755852 -2.575000e+00 -1.714338  0.970
#> 4:  Si1_2  Sr1_4_0_0_-1 3.184477  0.755852  2.575000e+00 -1.714338  0.970
#> 5:  Si1_2  Sr2_1_0_0_-1 3.261366 -1.903417  0.000000e+00 -2.648304  0.908
#> 6:  Si1_2 Sr2_3_-1_-1_0 3.465249 -2.214841 -2.575000e+00  0.686880  0.729
#>    DistanceError
#>            <num>
#> 1:             0
#> 2:             0
#> 3:             0
#> 4:             0
#> 5:             0
#> 6:             0
# }
```
