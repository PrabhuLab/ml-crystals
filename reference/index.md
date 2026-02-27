# Package index

## Core Analysis

The main entry points for processing CIF files.

- [`analyze_cif_files()`](https://prabhulab.github.io/ml-crystals/reference/analyze_cif_files.md)
  : Analyze a Batch of CIF Files
- [`analyze_single_cif()`](https://prabhulab.github.io/ml-crystals/reference/analyze_single_cif.md)
  : Analyze the Content of a Single CIF File
- [`export_analysis_to_csv()`](https://prabhulab.github.io/ml-crystals/reference/export_analysis_to_csv.md)
  : Export Analysis Results to a Directory of CSVs

## Calculators

Functions for geometry and crystallography calculations.

- [`calculate_distances()`](https://prabhulab.github.io/ml-crystals/reference/calculate_distances.md)
  : Calculate Interatomic Distances
- [`calculate_angles()`](https://prabhulab.github.io/ml-crystals/reference/calculate_angles.md)
  : Calculate Bond Angles
- [`calculate_neighbor_counts()`](https://prabhulab.github.io/ml-crystals/reference/calculate_neighbor_counts.md)
  : Calculate Coordination Numbers
- [`calculate_weighted_neighbor_counts()`](https://prabhulab.github.io/ml-crystals/reference/calculate_weighted_neighbor_counts.md)
  : Calculate Weighted Coordination Numbers
- [`calculate_weighted_average_network_distance()`](https://prabhulab.github.io/ml-crystals/reference/calculate_weighted_average_network_distance.md)
  : Calculate Weighted Average Network Bond Distance

## Error Propagation

Functions to propagate uncertainties through calculations.

- [`propagate_distance_error()`](https://prabhulab.github.io/ml-crystals/reference/propagate_distance_error.md)
  : Propagate Distance Error
- [`propagate_angle_error()`](https://prabhulab.github.io/ml-crystals/reference/propagate_angle_error.md)
  : Propagate Angle Error

## Bonding Algorithms

Methods for determining atomic connectivity.

- [`minimum_distance()`](https://prabhulab.github.io/ml-crystals/reference/minimum_distance.md)
  : Identify Atomic Bonds using the Minimum Distance Method
- [`brunner_nn_reciprocal()`](https://prabhulab.github.io/ml-crystals/reference/brunner_nn_reciprocal.md)
  : Identify Atomic Bonds using Brunner's Method (Reciprocal)
- [`econ_nn()`](https://prabhulab.github.io/ml-crystals/reference/econ_nn.md)
  : Identify Atomic Bonds using Hoppe's EconNN Method
- [`voronoi_nn()`](https://prabhulab.github.io/ml-crystals/reference/voronoi_nn.md)
  : Identify Atomic Bonds using Voronoi Tessellation
- [`crystal_nn()`](https://prabhulab.github.io/ml-crystals/reference/crystal_nn.md)
  : Identify Atomic Bonds using CrystalNN

## Filtering & Post-Processing

Tools to clean and refine data.

- [`filter_atoms_by_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_atoms_by_symbol.md)
  : Filter Data by Atom Symbol Interactively
- [`filter_by_wyckoff_symbol()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_wyckoff_symbol.md)
  : Filter Data by Wyckoff Symbol
- [`filter_by_elements()`](https://prabhulab.github.io/ml-crystals/reference/filter_by_elements.md)
  : Filter Distances by Element Symbols
- [`filter_ghost_distances()`](https://prabhulab.github.io/ml-crystals/reference/filter_ghost_distances.md)
  : Filter Ghost Distances Using Atomic Radii
- [`aggregate_batch_results()`](https://prabhulab.github.io/ml-crystals/reference/aggregate_batch_results.md)
  : Aggregate Batched Analysis Results

## Coordinate Transformations

Tools for symmetry application and supercell expansion.

- [`apply_symmetry_operations()`](https://prabhulab.github.io/ml-crystals/reference/apply_symmetry_operations.md)
  : Apply Symmetry Operations to Generate a Full Unit Cell
- [`expand_transformed_coords()`](https://prabhulab.github.io/ml-crystals/reference/expand_transformed_coords.md)
  : Expand Coordinates into a Supercell
- [`calculate_expansion_factors()`](https://prabhulab.github.io/ml-crystals/reference/calculate_expansion_factors.md)
  : Calculate Supercell Expansion Factors
- [`merge_sites_pbc()`](https://prabhulab.github.io/ml-crystals/reference/merge_sites_pbc.md)
  : Merge Close Atoms (Pymatgen Style)

## Data Extraction (Low-Level)

Parsers for specific CIF data blocks.

- [`read_cif_files()`](https://prabhulab.github.io/ml-crystals/reference/read_cif_files.md)
  : Read CIF Files into Memory
- [`extract_atomic_coordinates()`](https://prabhulab.github.io/ml-crystals/reference/extract_atomic_coordinates.md)
  : Extract Atomic Coordinates
- [`extract_unit_cell_metrics()`](https://prabhulab.github.io/ml-crystals/reference/extract_unit_cell_metrics.md)
  : Extract Unit Cell Metrics
- [`extract_symmetry_operations()`](https://prabhulab.github.io/ml-crystals/reference/extract_symmetry_operations.md)
  : Extract Symmetry Operations
- [`extract_database_code()`](https://prabhulab.github.io/ml-crystals/reference/extract_database_code.md)
  : Extract Database Code from CIF Content
- [`extract_chemical_formula()`](https://prabhulab.github.io/ml-crystals/reference/extract_chemical_formula.md)
  : Extract Chemical Formula from CIF Content
- [`extract_structure_type()`](https://prabhulab.github.io/ml-crystals/reference/extract_structure_type.md)
  : Extract Structure Type from CIF Content
- [`extract_space_group_name()`](https://prabhulab.github.io/ml-crystals/reference/extract_space_group_name.md)
  : Extract Space Group Name from CIF Content
- [`extract_space_group_number()`](https://prabhulab.github.io/ml-crystals/reference/extract_space_group_number.md)
  : Extract Space Group Number from CIF Content

## Data

- [`covalent_radii`](https://prabhulab.github.io/ml-crystals/reference/covalent_radii.md)
  : Atomic Radii Data for Bond-Length Estimation
- [`set_radii_data()`](https://prabhulab.github.io/ml-crystals/reference/set_radii_data.md)
  : Set or Reset a Custom Atomic Radii Table
