context("High-Level Wrapper Functions")

test_that("analyze_cif_files produces a valid, self-consistent result", {
  # The 'results_422' object is created in the helper file

  # Check top-level structure
  expect_s3_class(results_422, "data.table")
  expect_equal(nrow(results_422), 1)

  # Check that list-columns contain data.frames and can be extracted
  expect_type(results_422$bonded_pairs, "list")

  # This syntax correctly extracts the data.table from the list-column
  bonded_pairs_table <- results_422$bonded_pairs[[1]]
  expect_s3_class(bonded_pairs_table, "data.frame")

  # Verify that error propagation ran by checking for the error column
  expect_true("DistanceError" %in% names(bonded_pairs_table))
  expect_true("AngleError" %in% names(results_422$bond_angles[[1]]))
})

test_that("analyze_cif_files handles multiple files", {
  # Create a temporary directory and copy the file twice
  temp_dir <- tempfile()
  dir.create(temp_dir)
  file.copy(cif_path_422, file.path(temp_dir, "file1.cif"))
  file.copy(cif_path_422, file.path(temp_dir, "file2.cif"))

  multi_results <- analyze_cif_files(Sys.glob(file.path(temp_dir, "*.cif")))

  # Expect two rows in the final result
  expect_equal(nrow(multi_results), 2)

  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("analyze_cif_files handles errors gracefully", {
  # Create a bad CIF file
  bad_cif_path <- tempfile(fileext = ".cif")
  writeLines("this is not a valid cif file", bad_cif_path)

  # Expect a warning, but no error, and an empty data.table
  expect_warning(
    results <- analyze_cif_files(bad_cif_path),
    "Could not process CIF due to missing essential data"
  )
  expect_s3_class(results, "data.table")
  expect_equal(nrow(results), 0)

  # Clean up
  unlink(bad_cif_path)
})
