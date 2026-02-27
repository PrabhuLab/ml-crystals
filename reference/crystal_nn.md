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
