# Identify Atomic Bonds using CrystalNN

Poer of Pymatgen's `CrystalNN` algorithm from the ground up. It uses a
Voronoi-based algorithm and solid angle weights to determine the
probability of various coordination environments. It modifies
probability using smooth distance cutoffs and Pauling electronegativity
differences. The output is either the most probable coordination
environment (`weighted_cn = FALSE`) or a weighted list of coordination
environments (`weighted_cn = TRUE`).

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

  Numeric. Cutoff in Angstroms for initial neighbor search (default
  7.0).

- x_diff_weight:

  Numeric. Electronegativity difference weight (default 3.0).

- porous_adjustment:

  Logical. If TRUE, adjusts Voronoi weights to better describe layered /
  porous structures (default TRUE).

- distance_cutoffs:

  Numeric vector. Penalizes neighbor distances greater than sum of
  covalent radii plus these cutoffs. Set to NULL to disable.

- cation_anion:

  Logical. If TRUE, restricts bonding targets to sites with opposite or
  zero charge (requires oxidation states).

- weighted_cn:

  Logical. If FALSE (default), returns neighbors for the most probable
  coordination environment with weight 1.0. If TRUE, returns all
  potential neighbors with fractional probabilities.

## Value

A `data.table` of bonded pairs.

## See also

Other bonding algorithms:
[`brunner_nn_reciprocal()`](https://prabhulab.github.io/ml-crystals/reference/brunner_nn_reciprocal.md),
[`econ_nn()`](https://prabhulab.github.io/ml-crystals/reference/econ_nn.md),
[`minimum_distance()`](https://prabhulab.github.io/ml-crystals/reference/minimum_distance.md),
[`voronoi_nn()`](https://prabhulab.github.io/ml-crystals/reference/voronoi_nn.md)
