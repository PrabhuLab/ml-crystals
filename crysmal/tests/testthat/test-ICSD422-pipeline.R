context("Full pipeline processing for ICSD422.cif")

# --- Setup: Run the full analysis once ---
cif_path <- system.file("extdata", "ICSD422.cif", package = "crysmal")
results <- analyze_cif_files(cif_path)

test_that("Pipeline returns a single valid result row", {
  expect_equal(nrow(results), 1)
  expect_true(is.data.frame(results))
})

test_that("Metadata is extracted correctly", {
  expect_equal(results$database_code, "ICSD 422")
  expect_equal(results$chemical_formula, "Si1 Sr2")
  expect_equal(results$structure_type, "TiNiSi#MgSrSi")
  expect_equal(results$space_group_name, "P n m a")
  expect_equal(results$space_group_number, "62")
})

test_that("Unit cell metrics are correct", {
  metrics <- results$unit_cell_metrics[[1]]
  expect_equal(metrics$`_cell_length_a`, 8.11)
  expect_equal(metrics$`_cell_length_b`, 5.15)
  expect_equal(metrics$`_cell_length_c`, 9.54)
  expect_equal(metrics$`_cell_angle_alpha`, 90)
  expect_true(is.na(metrics$`_cell_length_a_error`))
})

test_that("Asymmetric coordinates and symmetry ops are correct", {
  coords <- results$atomic_coordinates[[1]]
  sym_ops <- results$symmetry_operations[[1]]

  expect_equal(nrow(coords), 3)
  expect_equal(coords[Label == "Sr1", x_a], 0.6529)
  expect_equal(coords[Label == "Si1", x_error], 0.0016)
  expect_true(is.na(coords[Label == "Sr1", y_error]))

  expect_equal(nrow(sym_ops), 8)
})

test_that("Symmetry application and expansion are correct", {
  t_coords <- results$transformed_coords[[1]]
  e_coords <- results$expanded_coords[[1]]

  # Ground truth from notebook is 12 atoms in full cell
  expect_equal(nrow(t_coords), 12)
  # Ground truth from notebook is 324 atoms in supercell
  expect_equal(nrow(e_coords), 324)
})

test_that("Bonding analysis results match ground truth", {
  bonds <- results$bonded_pairs[[1]]

  # Ground truth from notebook is 14 bonds using minimum_distance
  expect_equal(nrow(bonds), 14)

  # Check a specific bond: Si1 to Sr1_1_0_0_0
  si1_sr1_bond <- bonds[Atom1 == "Si1" & Atom2 == "Sr1_1_0_0_0"]
  expect_equal(si1_sr1_bond$Distance, 3.163544, tolerance = 1e-6)
  expect_equal(si1_sr1_bond$DistanceError, 0.014160750, tolerance = 1e-6)
})

test_that("Neighbor counts are correct", {
  counts <- results$neighbor_counts[[1]]

  expect_equal(nrow(counts), 3)
  expect_equal(counts[Atom == "Sr1", NeighborCount], 4)
  expect_equal(counts[Atom == "Sr2", NeighborCount], 3)
  expect_equal(counts[Atom == "Si1", NeighborCount], 7)
})

test_that("Bond angle calculations match ground truth", {
  angles <- results$bond_angles[[1]]

  # Ground truth from notebook is 30 angles
  expect_equal(nrow(angles), 30)

  # Check the first angle from the notebook's sorted list
  first_angle_row <- angles[CentralAtom == "Si1" &
                              ((Neighbor1 == "Sr1_4_0_-1_-1" & Neighbor2 == "Sr1_4_0_0_-1") |
                                 (Neighbor2 == "Sr1_4_0_-1_-1" & Neighbor1 == "Sr1_4_0_0_-1"))]

  expect_equal(first_angle_row$Angle, 107.92071, tolerance = 1e-5)
  expect_equal(first_angle_row$AngleError, 0.3991814, tolerance = 1e-5)

  # Check the last angle from the notebook's sorted list
  last_angle_row <- angles[CentralAtom == "Sr2" &
                             ((Neighbor1 == "Si1_1_0_0_0" & Neighbor2 == "Si1_3_0_0_0") |
                                (Neighbor2 == "Si1_1_0_0_0" & Neighbor1 == "Si1_3_0_0_0"))]

  expect_equal(last_angle_row$Angle, 102.24369, tolerance = 1e-5)
  expect_equal(last_angle_row$AngleError, 0.2819465, tolerance = 1e-5)
})
