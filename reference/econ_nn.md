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

## See also

Other bonding algorithms:
[`brunner_nn_reciprocal()`](https://prabhulab.github.io/ml-crystals/reference/brunner_nn_reciprocal.md),
[`crystal_nn()`](https://prabhulab.github.io/ml-crystals/reference/crystal_nn.md),
[`minimum_distance()`](https://prabhulab.github.io/ml-crystals/reference/minimum_distance.md),
[`voronoi_nn()`](https://prabhulab.github.io/ml-crystals/reference/voronoi_nn.md)
