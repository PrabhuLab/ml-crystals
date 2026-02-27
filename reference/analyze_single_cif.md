# Analyze the Content of a Single CIF File

The core worker function that orchestrates the analysis pipeline for a
single crystal structure's data. It dynamically adjusts supercell size
and coordinate sets based on the requested bonding algorithms. For
geometric algorithms (Voronoi/CrystalNN), it automatically merges atoms
occupying the same site to ensure mathematical validity.

## Usage

``` r
analyze_single_cif(
  cif_content,
  file_name = NULL,
  perform_extraction = TRUE,
  perform_calcs_and_transforms = TRUE,
  bonding_algorithms = c("minimum_distance"),
  calculate_bond_angles = TRUE,
  perform_error_propagation = TRUE,
  tolerance = 1e-04,
  minimum_distance_delta = 0.1
)
```

## Arguments

- cif_content:

  Either a `data.table` containing the lines of a CIF file, OR a
  character string specifying the file path to a CIF file.

- file_name:

  The name of the original CIF file, used for labeling output.

- perform_extraction:

  Logical. If `TRUE`, extracts all metadata.

- perform_calcs_and_transforms:

  Logical. If `TRUE`, generates unit cell and supercell.

- bonding_algorithms:

  Character vector. Options: `"minimum_distance"`, `"brunner"`,
  `"econ"`, `"crystal_nn"`, `"voronoi"`.

- calculate_bond_angles:

  Logical.

- perform_error_propagation:

  Logical.

- tolerance:

  Numeric. Cutoff for floating-point noise and atom merging (default
  1e-4).

- minimum_distance_delta:

  Numeric. Tolerance for min-dist algorithm.

## Value

A one-row `data.table` with results.
