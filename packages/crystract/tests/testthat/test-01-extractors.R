test_that("Metadata extraction works for 1590946", {
  # Check that the extracted database code string contains 1590946
  db_code <- extract_database_code(cif_content_1590946)
  expect_true(grepl("1590946", db_code))

  expect_equal(extract_chemical_formula(cif_content_1590946), "Si1 Sr2")
  expect_equal(extract_space_group_name(cif_content_1590946), "P n m a")
})

test_that("Unit cell metrics extraction works for 1590946", {
  expect_s3_class(metrics_1590946, "data.table")
  expect_equal(nrow(metrics_1590946), 1)
  expect_equal(metrics_1590946$`_cell_length_a`, 8.11)
  expect_equal(metrics_1590946$`_cell_length_b`, 5.15)
  expect_equal(metrics_1590946$`_cell_length_c`, 9.54)
})

test_that("Atomic coordinate extraction works for 1590946", {
  expect_s3_class(atoms_1590946, "data.table")
  expect_equal(nrow(atoms_1590946), 3)

  sr1 <- atoms_1590946[Label == "Sr1"]
  expect_equal(sr1$x_a, 0.6529)
  expect_equal(sr1$y_b, 0.25)

  # Ensure parsing didn't fail on missing error parameters (should default to 0)
  expect_equal(sr1$x_error, 0)
})

test_that("Symmetry operations extraction works for 1590946", {
  expect_s3_class(sym_ops_1590946, "data.table")
  expect_equal(nrow(sym_ops_1590946), 8)
  # The first operation in this specific file is 'x+1/2, y, -z+1/2'
  expect_equal(sym_ops_1590946$x[1], "x+1/2")
  expect_equal(sym_ops_1590946$y[1], "y")
  expect_equal(sym_ops_1590946$z[1], "-z+1/2")
  # The identity operation is last in this file
  expect_equal(sym_ops_1590946$x[8], "x")
})
