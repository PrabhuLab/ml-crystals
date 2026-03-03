test_that(
  "analyze_cif_files produces a valid, self-consistent result with multiple algorithms",
  {
    expect_s3_class(results_1590946, "data.table")
    expect_equal(nrow(results_1590946), 1)

    # Check that requested algorithm columns were generated
    expect_true("bonds_minimum_distance" %in% names(results_1590946))
    expect_true("bonds_crystal_nn" %in% names(results_1590946))
    expect_true("bonds_voronoi" %in% names(results_1590946))
    expect_true("bonds_econ" %in% names(results_1590946))

    # Validate structure of a specific result table
    bonded_pairs_table <- results_1590946$bonds_minimum_distance[[1]]
    expect_s3_class(bonded_pairs_table, "data.frame")
    expect_true("DistanceError" %in% names(bonded_pairs_table))
  }
)

test_that("analyze_cif_files handles multiple files and batching", {
  temp_dir <- tempfile("test-multi-")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  # Copy file twice to simulate batch
  file.copy(cif_path_1590946, file.path(temp_dir, "file1.cif"))
  file.copy(cif_path_1590946, file.path(temp_dir, "file2.cif"))

  multi_results <- analyze_cif_files(list.files(temp_dir, full.names = TRUE))
  expect_equal(nrow(multi_results), 2)
})

test_that("analyze_cif_files handles errors gracefully", {
  bad_cif_path <- tempfile(fileext = ".cif")
  on.exit(unlink(bad_cif_path), add = TRUE)
  writeLines("this is not a valid cif file", bad_cif_path)

  # Expect a warning that processing failed.
  # The function should return an EMPTY data.table if no files succeed.
  expect_warning(results <- analyze_cif_files(bad_cif_path))
  expect_s3_class(results, "data.table")
  expect_equal(nrow(results), 0)
})
