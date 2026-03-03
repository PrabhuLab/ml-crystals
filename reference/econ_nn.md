# Identify Atomic Bonds using Hoppe's EconNN Method

Uses Effective Coordination Numbers (ECoN) and Mean Fictive Ionic Radii
(MEFIR).

## Usage

``` r
econ_nn(distances, atomic_coordinates, tol = 0.2, use_fictive_radius = FALSE)
```

## Arguments

- distances:

  A `data.table` of interatomic distances.

- atomic_coordinates:

  A `data.table` of atomic coordinates used to map species to radii.

- tol:

  A bond strength cutoff (default 0.2).

- use_fictive_radius:

  Logical. If `TRUE`, calculates Hoppe's fictive ionic radius.

## Value

A `data.table` of bonded pairs.

## References

Hoppe, R. (1979). "Effective Coordination Numbers (ECoN) and Mean
Fictive Ionic Radii (MEFIR)." *Zeitschrift Für Kristallographie*,
150(1–4), 23–52.
[doi:10.1524/zkri.1979.150.14.23](https://doi.org/10.1524/zkri.1979.150.14.23)

Shannon, R. D. (1976). "Revised Effective Ionic Radii and Systematic
Studies of Interatomic Distances in Halides and Chalcogenides." *Acta
Crystallographica Section A*, 32(5), 751–767.
[doi:10.1107/S0567739476001551](https://doi.org/10.1107/S0567739476001551)

## See also

Other bonding algorithms:
[`brunner_nn_reciprocal()`](https://prabhulab.github.io/ml-crystals/reference/brunner_nn_reciprocal.md),
[`crystal_nn()`](https://prabhulab.github.io/ml-crystals/reference/crystal_nn.md),
[`minimum_distance()`](https://prabhulab.github.io/ml-crystals/reference/minimum_distance.md),
[`voronoi_nn()`](https://prabhulab.github.io/ml-crystals/reference/voronoi_nn.md)

## Examples

``` r
dists <- data.table::data.table(Atom1 = c("Si1", "Si1"), Atom2 = c("O1", "O2"),
                                Distance = c(1.6, 2.0),
                                DeltaX = c(1, 0), DeltaY = c(0, 1), DeltaZ = c(0, 0))
ac <- data.table::data.table(Label = c("Si1", "O1", "O2"),
                             OxidationState = c(4, -2, -2))
econ_nn(dists, ac)
#>     Atom1  Atom2 Distance DeltaX DeltaY DeltaZ   Weight
#>    <char> <char>    <num>  <num>  <num>  <num>    <num>
#> 1:    Si1     O1      1.6      1      0      0 1.110372
```
