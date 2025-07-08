# This helper script runs the full analysis on the ICSD422.cif file once.
# The resulting objects are then available to all test files, making tests
# faster and less repetitive.

# Load the single CIF file content
cif_path_422 <- system.file("extdata", "ICSD422.cif", package = "crysmal")
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
