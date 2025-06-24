# Load the example CIF file content
cif_file_path <- test_path("testdata", "ICSD422.cif")
cif_content_ICSD422 <- data.table::fread(cif_file_path, sep = "\n", header = FALSE, strip.white = FALSE, quote = "", data.table = TRUE)

# Expected data for ICSD422.cif (from your provided examples)
# These will be used as ground truth in the tests.

# --- From extract_cif_data.R ---
expected_database_code_ICSD422 <- "ICSD 422"
expected_chemical_formula_ICSD422 <- "Si1 Sr2"
expected_structure_type_ICSD422 <- "TiNiSi#MgSrSi"
expected_sg_name_ICSD422 <- "P n m a"
expected_sg_number_ICSD422 <- "62"
expected_unit_cell_metrics_ICSD422 <- data.table::data.table(
  `_cell_length_a` = 8.11, `_cell_length_b` = 5.15, `_cell_length_c` = 9.54,
  `_cell_angle_alpha` = 90, `_cell_angle_beta` = 90, `_cell_angle_gamma` = 90
)

# --- From process_coordinates.R ---
expected_atomic_coordinates_ICSD422 <- data.table::data.table(
  Label = c("Sr1", "Sr2", "Si1"),
  x_a = c(0.6529, 0.5192, 0.2539),
  y_b = c(0.25, 0.25, 0.25),
  z_c = c(0.0769, 0.6748, 0.1028),
  x_error = c(0.0006, 0.0006, 0.0016),
  y_error = c(NA_real_, NA_real_, NA_real_),
  z_error = c(0.0005, 0.0005, 0.0014)
)

expected_symmetry_operations_ICSD422 <- data.table::data.table(
  x = c("x+1/2", "x", "-x+1/2", "-x", "-x+1/2", "-x", "x+1/2", "x"),
  y = c("y", "-y+1/2", "y+1/2", "-y", "-y", "y+1/2", "-y+1/2", "y"),
  z = c("-z+1/2", "z", "z+1/2", "-z", "z+1/2", "-z", "-z+1/2", "z")
)
# Order matters for comparison, so ensure it matches function output if not reordered
expected_symmetry_operations_ICSD422 <- expected_symmetry_operations_ICSD422[order(x, y, z)]


expected_transformed_coords_ICSD422 <- data.table::data.table(
  Label = c("Sr1_sym1", "Sr1_sym2", "Sr1_sym3", "Sr1_sym4", "Sr1_sym5", "Sr1_sym6", "Sr1_sym7", "Sr1_sym8",
            "Sr2_sym1", "Sr2_sym2", "Sr2_sym3", "Sr2_sym4", "Sr2_sym5", "Sr2_sym6", "Sr2_sym7", "Sr2_sym8",
            "Si1_sym1", "Si1_sym2", "Si1_sym3", "Si1_sym4", "Si1_sym5", "Si1_sym6", "Si1_sym7", "Si1_sym8"),
  x_a = c(0.1529, 0.6529, 0.8471, 0.3471, 0.8471, 0.3471, 0.1529, 0.6529,
          0.0192, 0.5192, 0.9808, 0.4808, 0.9808, 0.4808, 0.0192, 0.5192,
          0.7539, 0.2539, 0.2461, 0.7461, 0.2461, 0.7461, 0.7539, 0.2539),
  y_b = c(0.25, 0.75, 0.75, 0.25, 0.75, 0.25, 0.25, 0.75,
          0.25, 0.75, 0.75, 0.25, 0.75, 0.25, 0.25, 0.75,
          0.25, 0.75, 0.75, 0.25, 0.75, 0.25, 0.25, 0.75),
  z_c = c(0.4231, 0.0769, 0.5769, 0.9231, 0.5769, 0.9231, 0.4231, 0.0769,
          0.8252, 0.6748, 0.1748, 0.3252, 0.1748, 0.3252, 0.8252, 0.6748,
          0.3972, 0.1028, 0.6028, 0.8972, 0.6028, 0.8972, 0.3972, 0.1028)
)
expected_bonded_pairs_min_dist_ICSD422_dims <- c(14, 5) # 14 rows, 5 columns
expected_bonded_pairs_min_dist_ICSD422_names <- c("Atom1", "Atom2", "Distance", "dcut", "dmin")

expected_neighbor_counts_ICSD422 <- data.table::data.table(
  Atom = c("Si1", "Sr1", "Sr2"), # Order might vary, so test content
  NeighborCount = c(7L, 4L, 3L)
)

expected_bond_angles_ICSD422_dims <- c(30, 4) # 30 rows, 4 columns
expected_bond_angles_ICSD422_names <- c("CentralAtom", "Neighbor1", "Neighbor2", "Angle")
