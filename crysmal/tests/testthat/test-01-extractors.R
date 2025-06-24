context("CIF Data Extraction")

# Test internal extract_value function
test_that("extract_value works correctly", {
  dt <- data.table::data.table(V1 = c("key value", "key 'quoted value'", "other line"))
  expect_equal(extract_value(dt, "key"), "value")
  expect_equal(extract_value(dt, "key", remove_pattern = FALSE), "key value")
  expect_equal(extract_value(dt, "key '"), "quoted value") # Pattern "key '" removes that part
  expect_equal(extract_value(dt, "nonexistent"), NA_character_)
  expect_equal(extract_value(data.table::data.table(X1 = "foo"), "foo"), NA_character_) # Wrong col name
  expect_equal(extract_value(NULL, "foo"), NA_character_)
  dt_spaces <- data.table::data.table(V1 = c("_tag_  extra_spaces_val  "))
  expect_equal(extract_value(dt_spaces, "_tag_"), "extra_spaces_val")
})

test_that("Database code extraction works for ICSD422", {
  expect_equal(extract_database_code(cif_content_ICSD422), expected_database_code_ICSD422)
})

test_that("Chemical formula extraction works for ICSD422", {
  expect_equal(extract_chemical_formula(cif_content_ICSD422), expected_chemical_formula_ICSD422)
})

test_that("Structure type extraction works for ICSD422", {
  expect_equal(extract_structure_type(cif_content_ICSD422), expected_structure_type_ICSD422)
})

test_that("Space group name extraction works for ICSD422", {
  expect_equal(extract_space_group_name(cif_content_ICSD422), expected_sg_name_ICSD422)
})

test_that("Space group number extraction works for ICSD422", {
  expect_equal(extract_space_group_number(cif_content_ICSD422), expected_sg_number_ICSD422)
})

test_that("Unit cell metrics extraction works for ICSD422", {
  metrics <- extract_unit_cell_metrics(cif_content_ICSD422)
  expect_s3_class(metrics, "data.table")
  expect_equal(nrow(metrics), 1)
  # Convert to numeric for comparison if helper is character
  expected_num <- data.table::copy(expected_unit_cell_metrics_ICSD422)
  for(j in seq_along(expected_num)) data.table::set(expected_num, j=j, value=as.numeric(expected_num[[j]]))

  expect_equal(metrics, expected_num, tolerance = 1e-6)
})

test_that("Unit cell metrics handles missing values and parentheses", {
  cif_missing <- data.table::data.table(V1 = c("_cell_length_a 10.0(1)", "_cell_length_b 5"))
  metrics <- extract_unit_cell_metrics(cif_missing)
  expect_equal(metrics$`_cell_length_a`, 10.0)
  expect_equal(metrics$`_cell_length_b`, 5.0)
  expect_true(is.na(metrics$`_cell_length_c`))
  expect_equal(nrow(metrics), 1) # Should still return one row of data

  cif_empty <- data.table::data.table(V1 = character(0))
  metrics_empty <- extract_unit_cell_metrics(cif_empty)
  expect_true(is.data.table(metrics_empty))
  expect_equal(nrow(metrics_empty), 1) # Now returns 1 row with NAs
  expect_true(all(is.na(metrics_empty)))

  cif_no_cell_info <- data.table::data.table(V1 = c("some other data"))
  metrics_no_cell <- extract_unit_cell_metrics(cif_no_cell_info)
  expect_true(is.data.table(metrics_no_cell))
  expect_equal(nrow(metrics_no_cell), 1)
  expect_true(all(is.na(metrics_no_cell)))
})
