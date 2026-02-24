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
