test_that("Metadata extraction works for ICSD422", {
  expect_equal(extract_database_code(cif_content_422), "ICSD 422")
  expect_equal(extract_chemical_formula(cif_content_422), "Si1 Sr2")
  expect_equal(extract_space_group_name(cif_content_422), "P n m a")
})

test_that("Unit cell metrics extraction works for ICSD422", {
  expect_s3_class(metrics_422, "data.table")
  expect_equal(nrow(metrics_422), 1)
  expect_equal(metrics_422$`_cell_length_a`, 8.11)
  expect_true(is.na(metrics_422$`_cell_length_a_error`))
})

test_that("Atomic coordinate extraction works for ICSD422", {
  expect_s3_class(atoms_422, "data.table")
  expect_equal(nrow(atoms_422), 3)
  si1 <- atoms_422[Label == "Si1"]
  expect_equal(si1$x_a, 0.2539)
  expect_equal(si1$x_error, 0.0016)
})

test_that("Symmetry operations extraction works for ICSD422", {
  expect_s3_class(sym_ops_422, "data.table")
  expect_equal(nrow(sym_ops_422), 8)
  # The first operation in this specific file is 'x+1/2, y, -z+1/2'
  expect_equal(sym_ops_422$x[1], "x+1/2")
  expect_equal(sym_ops_422$y[1], "y")
  expect_equal(sym_ops_422$z[1], "-z+1/2")
  # The identity operation is last in this file
  expect_equal(sym_ops_422$x[8], "x")
})
