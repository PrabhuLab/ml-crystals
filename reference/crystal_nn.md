# Identify Atomic Bonds using CrystalNN

Rebuild of Pymatgen's `CrystalNN` algorithm from the ground up.

## Usage

``` r
crystal_nn(
  distances,
  atomic_coordinates,
  expanded_coords,
  unit_cell_metrics,
  cutoff_length = 7,
  x_diff_weight = 3,
  porous_adjustment = TRUE,
  distance_cutoffs = c(0.5, 1),
  cation_anion = FALSE,
  weighted_cn = FALSE
)
```

## Arguments

- distances:

  Ignored (CrystalNN generates its own Voronoi basis).

- atomic_coordinates:

  Primary atom set.

- expanded_coords:

  Expanded supercell.

- unit_cell_metrics:

  Cell parameters.

- cutoff_length:

  Numeric. Cutoff in Angstroms for initial neighbor search.

- x_diff_weight:

  Numeric. Electronegativity difference weight.

- porous_adjustment:

  Logical. If TRUE, adjusts Voronoi weights.

- distance_cutoffs:

  Numeric vector. Penalizes neighbor distances greater than sum of
  radii.

- cation_anion:

  Logical. Restrictions targets to opposite charge.

- weighted_cn:

  Logical. Return fractional probabilities vs strict max probability.

## Value

A `data.table` of bonded pairs.

## See also

Other bonding algorithms:
[`brunner_nn_reciprocal()`](https://prabhulab.github.io/ml-crystals/reference/brunner_nn_reciprocal.md),
[`econ_nn()`](https://prabhulab.github.io/ml-crystals/reference/econ_nn.md),
[`minimum_distance()`](https://prabhulab.github.io/ml-crystals/reference/minimum_distance.md),
[`voronoi_nn()`](https://prabhulab.github.io/ml-crystals/reference/voronoi_nn.md)
