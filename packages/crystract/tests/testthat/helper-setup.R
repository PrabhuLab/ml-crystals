# This helper script prepares data for running tests locally.
# NOTE: Due to licensing restrictions, the CCDC/ICSD CIF files are not included in this package.
# To run these tests, you must:
#   1. Download the file for CCDC entry 1590946 (ICSD 422).
#   2. Save it as "1590946.cif" inside the 'tests/testthat/' directory.

cif_path_1590946 <- "1590946.cif"

# Only run the analysis if the test file has been provided by the user.
if (file.exists(cif_path_1590946)) {
  # Load the single CIF file content
  cif_content_1590946 <- read_cif_files(cif_path_1590946)[[1]]

  # Run the full analysis pipeline with multiple algorithms to get the final results object
  results_1590946 <- analyze_cif_files(
    cif_path_1590946,
    bonding_algorithms = c(
      "minimum_distance",
      "brunner",
      "econ",
      "crystal_nn",
      "voronoi"
    )
  )

  # For more granular testing, create the intermediate building blocks
  atoms_1590946 <- extract_atomic_coordinates(cif_content_1590946)
  metrics_1590946 <- extract_unit_cell_metrics(cif_content_1590946)
  sym_ops_1590946 <- extract_symmetry_operations(cif_content_1590946)

  # Coordinate processing results
  # Note: apply_symmetry_operations now requires unit_cell_metrics for PBC distance calculation
  full_cell_1590946 <- apply_symmetry_operations(atoms_1590946, sym_ops_1590946, metrics_1590946)
  super_cell_1590946 <- expand_transformed_coords(full_cell_1590946)
  distances_1590946 <- calculate_distances(atoms_1590946, super_cell_1590946, metrics_1590946)

} else {
  # Skip all tests that require the data if the file is not present.
  # A message will be printed during testing to inform the developer.
  testthat::skip("Skipping tests: '1590946.cif' not found. See 'helper-setup.R' for instructions.")
}
