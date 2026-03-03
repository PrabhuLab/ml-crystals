# Identify Atomic Bonds using Voronoi Tessellation

Performs 3D Voronoi analysis on the supercell.

## Usage

``` r
voronoi_nn(
  atomic_coordinates,
  expanded_coords,
  unit_cell_metrics,
  cutoff = 13,
  tol = 0
)
```

## Arguments

- atomic_coordinates:

  A `data.table` of the primary (asymmetric) atom set.

- expanded_coords:

  A `data.table` of atoms in the expanded supercell.

- unit_cell_metrics:

  A `data.table` with cell parameters.

- cutoff:

  Distance cutoff (default 13.0).

- tol:

  Tolerance for solid angle weights (default 0).

## Value

A `data.table` of bonded pairs.

## See also

Other bonding algorithms:
[`brunner_nn_reciprocal()`](https://prabhulab.github.io/ml-crystals/reference/brunner_nn_reciprocal.md),
[`crystal_nn()`](https://prabhulab.github.io/ml-crystals/reference/crystal_nn.md),
[`econ_nn()`](https://prabhulab.github.io/ml-crystals/reference/econ_nn.md),
[`minimum_distance()`](https://prabhulab.github.io/ml-crystals/reference/minimum_distance.md)

## Examples

``` r
# \donttest{
cif_path <- system.file("extdata", "1590946.cif", package = "crystract")
if (file.exists(cif_path)) {
  res <- analyze_single_cif(cif_path, bonding_algorithms = "voronoi")
  print(head(res$bonds_voronoi[[1]]))
}
#> Key: <Atom1, Atom2, Distance>
#>     Atom1         Atom2 Distance  SolidAngle       Area    DeltaX        DeltaY
#>    <char>        <char>    <num>       <num>      <num>     <num>         <num>
#> 1:  Si1_2  Si1_1_-1_0_0 4.932659 0.004328448 0.02851541 -4.055000  2.220446e-16
#> 2:  Si1_2   Si1_1_0_0_0 4.932659 0.004328448 0.02851541  4.055000  2.220446e-16
#> 3:  Si1_2   Sr1_1_0_0_0 3.163544 1.687714706 6.98543858 -0.819110  2.220446e-16
#> 4:  Si1_2  Sr1_2_-1_0_0 4.880369 0.051103880 0.31256672 -4.874110  0.000000e+00
#> 5:  Si1_2   Sr1_2_0_0_0 3.245310 1.546340541 6.41452789  3.235890  0.000000e+00
#> 6:  Si1_2 Sr1_4_0_-1_-1 3.184477 1.609248649 6.54596788  0.755852 -2.575000e+00
#>       DeltaZ    MaxSA     Weight DistanceError
#>        <num>    <num>      <num>         <num>
#> 1:  2.808576 1.687715 0.00256468             0
#> 2:  2.808576 1.687715 0.00256468             0
#> 3:  3.055662 1.687715 1.00000000             0
#> 4: -0.247086 1.687715 0.03027993             0
#> 5: -0.247086 1.687715 0.91623337             0
#> 6: -1.714338 1.687715 0.95350751             0
# }
```
