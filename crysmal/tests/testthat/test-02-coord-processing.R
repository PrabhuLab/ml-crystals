context("Coordinate Processing Functions")

test_that("apply_symmetry_operations generates correct full cell for ICSD422", {
  # The full_cell_422 object is created in the helper file

  # Ground truth from notebook is 12 atoms in the full unit cell
  expect_equal(nrow(full_cell_422), 12)

  # Check a specific transformed coordinate from the notebook output
  # Sr1_1: (0.1529, 0.25, 0.4231)
  sr1_1 <- full_cell_422[Label == "Sr1_1"]
  expect_equal(sr1_1$x_a, 0.1529, tolerance = 1e-6)
  expect_equal(sr1_1$y_b, 0.2500, tolerance = 1e-6)
  expect_equal(sr1_1$z_c, 0.4231, tolerance = 1e-6)
})

test_that("expand_transformed_coords generates correct supercell for ICSD422", {
  # The super_cell_422 object is created in the helper file

  # Ground truth from notebook is 324 atoms in the 3x3x3 supercell
  expect_equal(nrow(super_cell_422), 324)
})
