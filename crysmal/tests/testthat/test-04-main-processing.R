context("Main CIF Processing Workflow")

test_that("read_cif_files works", {
  # Use the ICSD422.cif file saved in testdata
  cif_path <- test_path("testdata", "ICSD422.cif")

  # Single file
  contents_list <- read_cif_files(cif_path)
  expect_type(contents_list, "list")
  expect_length(contents_list, 1)
  expect_named(contents_list, basename(cif_path))
  expect_s3_class(contents_list[[1]], "data.table")
  expect_gt(nrow(contents_list[[1]]), 0) # Check it's not empty

  # Multiple files (by duplicating the path)
  contents_list_multi <- read_cif_files(c(cif_path, cif_path))
  expect_length(contents_list_multi, 2)
  # Names might be made unique by R if basenames are identical, or might be overwritten.
  # If paths are identical, names will be identical. list() allows this.
  expect_equal(names(contents_list_multi), rep(basename(cif_path), 2))

  # Non-existent file
  expect_warning(contents_bad_path <- read_cif_files("non_existent_file.cif"), "Failed to read or process CIF file")
  expect_length(contents_bad_path, 0)

  # Empty input
  expect_warning(contents_empty <- read_cif_files(character(0)), "No CIF file paths provided.")
  expect_length(contents_empty, 0)
})

test_that("process_single_cif_data works for ICSD422", {
  # Use default bonding_method = "min_dist"
  result <- process_single_cif_data(cif_content_ICSD422)

  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), 1)

  expect_equal(result$database_code, expected_database_code_ICSD422)
  expect_equal(result$chemical_formula, expected_chemical_formula_ICSD422)

  # Check list columns
  expect_s3_class(result$unit_cell_metrics[[1]], "data.table")
  expect_equal(result$unit_cell_metrics[[1]], expected_unit_cell_metrics_ICSD422, tolerance = 1e-6)

  expect_s3_class(result$atomic_coordinates_input[[1]], "data.table")
  # Order for comparison
  aci_res <- data.table::copy(result$atomic_coordinates_input[[1]]); data.table::setorder(aci_res, Label)
  aci_exp <- data.table::copy(expected_atomic_coordinates_ICSD422); data.table::setorder(aci_exp, Label)
  expect_equal(aci_res, aci_exp, tolerance = 1e-6)

  expect_s3_class(result$bonded_pairs_identified[[1]], "data.table")
  expect_equal(dim(result$bonded_pairs_identified[[1]]), expected_bonded_pairs_min_dist_ICSD422_dims)

  expect_s3_class(result$neighbor_counts_calculated[[1]], "data.table")
  nc_res <- data.table::copy(result$neighbor_counts_calculated[[1]]); data.table::setorder(nc_res, Atom)
  nc_exp <- data.table::copy(expected_neighbor_counts_ICSD422); data.table::setorder(nc_exp, Atom)
  expect_equal(nc_res, nc_exp)

  expect_s3_class(result$bond_angles_calculated[[1]], "data.table")
  expect_equal(dim(result$bond_angles_calculated[[1]]), expected_bond_angles_ICSD422_dims)

  expect_true(is.na(result$error_message))

  # Test with a different bonding method (Brunner)
  result_brunner <- process_single_cif_data(cif_content_ICSD422, bonding_method = "brunner")
  expect_s3_class(result_brunner$bonded_pairs_identified[[1]], "data.table") # or NULL

  # Test missing essential data: unit cell metrics
  cif_no_cell <- data.table::copy(cif_content_ICSD422)
  cif_no_cell <- cif_no_cell[! (V1 %like% "_cell_length_" | V1 %like% "_cell_angle_")]
  expect_warning(result_no_cell <- process_single_cif_data(cif_no_cell), "Essential unit cell metrics missing")
  expect_false(is.na(result_no_cell$error_message))
  expect_true(grepl("Missing unit cell metrics", result_no_cell$error_message))
  expect_null(result_no_cell$bonded_pairs_identified[[1]]) # Should not proceed to bonding

  # Test missing essential data: atomic coordinates
  cif_no_coords <- data.table::copy(cif_content_ICSD422)
  cif_no_coords <- cif_no_coords[! (V1 %like% "_atom_site_fract_")]
  expect_warning(result_no_coords <- process_single_cif_data(cif_no_coords), "Atomic coordinates missing")
  expect_false(is.na(result_no_coords$error_message))
  expect_true(grepl("Missing atomic coordinates", result_no_coords$error_message))
  expect_null(result_no_coords$bonded_pairs_identified[[1]])
})

test_that("analyze_cif_files works", {
  # Create a temporary directory and place a copy of ICSD422.cif in it
  temp_dir <- withr::local_tempdir()
  temp_cif_path <- file.path(temp_dir, "ICSD422_copy.cif")
  file.copy(test_path("testdata", "ICSD422.cif"), temp_cif_path)

  # Test with folder input
  analysis_result_folder <- analyze_cif_files(temp_dir, pattern = "*.cif")
  expect_s3_class(analysis_result_folder, "data.table")
  expect_equal(nrow(analysis_result_folder), 1)
  expect_equal(analysis_result_folder$source_file, basename(temp_cif_path))
  expect_equal(analysis_result_folder$database_code, expected_database_code_ICSD422)

  # Test with direct file path input
  analysis_result_path <- analyze_cif_files(temp_cif_path)
  expect_s3_class(analysis_result_path, "data.table")
  expect_equal(nrow(analysis_result_path), 1)
  expect_equal(analysis_result_path$source_file, basename(temp_cif_path))

  # Test with non-existent folder
  expect_error(analyze_cif_files("non_existent_folder_path"), "cif_input must be a valid folder path or a vector of valid file paths.")

  # Test with folder having no matching files
  temp_dir_empty <- withr::local_tempdir()
  expect_warning(res_empty_folder <- analyze_cif_files(temp_dir_empty, pattern = "*.cif"), "No files matching pattern")
  expect_equal(nrow(res_empty_folder), 0)
})
