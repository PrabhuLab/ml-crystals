# This helper script prepares data for running tests locally.
# NOTE: Due to licensing restrictions, the ICSD CIF files are not included in this package.
# To run these tests, you must:
#   1. Consult the 'required_ICSD_files.csv' file in the root of the 'ml-crystals' repository.
#   2. Download the file for ICSD entry #422 from a licensed ICSD provider.
#   3. Save it as "ICSD_422.cif" inside the 'tests/testthat/' directory.

cif_path_422 <- "ICSD_422.cif"

# Only run the analysis if the test file has been provided by the user.
if (file.exists(cif_path_422)) {
  # Load the single CIF file content
  cif_content_422 <- read_cif_files(cif_path_422)[[1]]

  # Run the full analysis pipeline to get the final results object
  results_422 <- analyze_cif_files(cif_path_422)

  # For more granular testing, create the intermediate building blocks
  atoms_422 <- extract_atomic_coordinates(cif_content_422)
  metrics_422 <- extract_unit_cell_metrics(cif_content_422)
  sym_ops_422 <- extract_symmetry_operations(cif_content_422)

  # And the coordinate processing results
  full_cell_422 <- apply_symmetry_operations(atoms_422, sym_ops_422)
  super_cell_422 <- expand_transformed_coords(full_cell_422)
} else {
  # Skip all tests that require the data if the file is not present.
  # A message will be printed during testing to inform the developer.
  testthat::skip("Skipping tests: 'ICSD_422.cif' not found. See 'helper-setup.R' for instructions.")
}
