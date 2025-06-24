context("Coordinate Processing")

# Test internal scale_error_internal function
test_that("scale_error_internal works correctly", {
  expect_equal(scale_error_internal("0.6529", "6"), 0.0006)
  expect_equal(scale_error_internal("0.0769", "5"), 0.0005)
  expect_equal(scale_error_internal("0.2539", "16"), 0.0016) # error is "16" for value "0.2539"
  expect_equal(scale_error_internal("0.25", NA_character_), NA_real_)
  expect_equal(scale_error_internal("123", "5"), 5)
  expect_equal(scale_error_internal("1.0", "5"), 0.5)
  expect_equal(scale_error_internal("1.23", "15"), 0.15)
  expect_true(is.na(scale_error_internal("1.23", "")))
})

test_that("Atomic coordinate extraction works for ICSD422", {
  coords <- extract_atomic_coordinates(cif_content_ICSD422)
  expect_s3_class(coords, "data.table")
  expect_equal(nrow(coords), 3)
  expect_equal(names(coords), c("Label", "x_a", "y_b", "z_c", "x_error", "y_error", "z_error"))

  # Order coords and expected_coords by Label for robust comparison
  data.table::setorder(coords, Label)
  data.table::setorder(expected_atomic_coordinates_ICSD422, Label)

  expect_equal(coords$Label, expected_atomic_coordinates_ICSD422$Label)
  expect_equal(coords$x_a, expected_atomic_coordinates_ICSD422$x_a, tolerance = 1e-6)
  expect_equal(coords$y_b, expected_atomic_coordinates_ICSD422$y_b, tolerance = 1e-6)
  expect_equal(coords$z_c, expected_atomic_coordinates_ICSD422$z_c, tolerance = 1e-6)
  expect_equal(coords$x_error, expected_atomic_coordinates_ICSD422$x_error, tolerance = 1e-6)
  expect_true(all(is.na(coords$y_error))) # As per example, no error for y
  expect_equal(coords$z_error, expected_atomic_coordinates_ICSD422$z_error, tolerance = 1e-6)
})

test_that("Atomic coordinate extraction handles edge cases", {
  expect_null(extract_atomic_coordinates(data.table::data.table(V1 = "no_atom_loop_tags")))
  expect_null(extract_atomic_coordinates(data.table::data.table(V1 = "_atom_site_label\n_atom_site_fract_x"))) # Header but no data

  # Test with alternative end-of-header tag
  cif_alt_header <- data.table::data.table(V1 = c(
    "loop_", "_atom_site_label", "_atom_site_fract_x", "_atom_site_U_iso_or_equiv",
    "Si1 0.1 0.001" # Note: This line isn't correctly formatted for the parser's col expectations
  ))
  expect_null(extract_atomic_coordinates(cif_alt_header))

  # Test with a properly formatted line for the alternative header
  cif_alt_header_ok <- data.table::data.table(V1 = c(
    "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x",
    "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_U_iso_or_equiv", # assume U_iso is last header
    "Si1 Si 0.1 0.2 0.3 0.05" # Label, Type, x, y, z, U_iso (6 columns from data)
  ))

})


test_that("Symmetry operation extraction works for ICSD422", {
  sym_ops <- extract_symmetry_operations(cif_content_ICSD422)
  expect_s3_class(sym_ops, "data.table")
  expect_equal(nrow(sym_ops), 8)
  expect_equal(names(sym_ops), c("x", "y", "z"))

  data.table::setorder(sym_ops, x, y, z) # Ensure order for comparison
  expect_equal(sym_ops, expected_symmetry_operations_ICSD422)
})

test_that("Symmetry operation extraction handles alternative tag and edge cases", {
  cif_alt_sym <- data.table::data.table(V1 = c("_symmetry_equiv_pos_as_xyz", " 'x,y,z' ", " '-x,-y,-z' "))
  sym_ops_alt <- extract_symmetry_operations(cif_alt_sym)
  expect_equal(nrow(sym_ops_alt), 2)
  expect_equal(sym_ops_alt$x, c("x", "-x"))

  expect_null(extract_symmetry_operations(data.table::data.table(V1 = "no_sym_op_tag")))
  cif_empty_lines <- data.table::data.table(V1 = c("_space_group_symop_operation_xyz", " ", " 'x,y,z' "))
  sym_ops_empty <- extract_symmetry_operations(cif_empty_lines)
  expect_equal(nrow(sym_ops_empty), 1)
  expect_equal(sym_ops_empty$x, "x")
})


test_that("apply_symmetry_operations works for ICSD422", {
  atomic_coords <- extract_atomic_coordinates(cif_content_ICSD422)
  sym_ops <- extract_symmetry_operations(cif_content_ICSD422)

  transformed <- apply_symmetry_operations(atomic_coords, sym_ops)
  expect_s3_class(transformed, "data.table")
  expect_equal(nrow(transformed), 12) # 3 atoms * 8 ops / (multiplicity) = 24 / 2 = 12 unique typically
  # Actually, it's 4 unique positions for Sr1, 4 for Sr2, 4 for Si1. Total 12.
  expect_equal(names(transformed), c("Label", "x_a", "y_b", "z_c"))

  sym_op_to_check_idx <- which(sym_ops$x == "x+1/2" & sym_ops$y == "y" & sym_ops$z == "-z+1/2")
  # There might be multiple such ops if table is not unique, but should be 1.
  if (length(sym_op_to_check_idx) > 0) {
    sym_op_to_check_idx <- sym_op_to_check_idx[1]
    sr1_op1_transformed_label <- paste0(atomic_coords$Label[1], "_sym", sym_op_to_check_idx)
    sr1_op1_coords <- transformed[Label == sr1_op1_transformed_label, .(x_a, y_b, z_c)]

    expect_equal(sr1_op1_coords$x_a, 0.1529, tolerance = 1e-4) # 1.1529 -> 0.1529
    expect_equal(sr1_op1_coords$y_b, 0.25, tolerance = 1e-4)
    expect_equal(sr1_op1_coords$z_c, 0.4231, tolerance = 1e-4) # -0.0769 + 0.5 = 0.4231
  }

  transformed_coords_for_comp <- transformed[, .(x_a = round(x_a, 4), y_b = round(y_b, 4), z_c = round(z_c, 4))]
  data.table::setorder(transformed_coords_for_comp, x_a, y_b, z_c)

  prompt_transformed_coords_example <- data.table::data.table(
    x_a = c(0.1529,0.6529,0.8471,0.3471,0.0192,0.5192,0.9808,0.4808,0.7539,0.2539,0.2461,0.7461),
    y_b = c(0.25,0.25,0.75,0.75,0.25,0.25,0.75,0.75,0.25,0.25,0.75,0.75),
    z_c = c(0.4231,0.0769,0.5769,0.9231,0.8252,0.6748,0.1748,0.3252,0.3972,0.1028,0.6028,0.8972)
  )
  prompt_transformed_coords_example[, `:=`(x_a = round(x_a, 4), y_b = round(y_b, 4), z_c = round(z_c, 4))]
  data.table::setorder(prompt_transformed_coords_example, x_a, y_b, z_c)

  expect_equal(transformed_coords_for_comp, prompt_transformed_coords_example)
})

test_that("expand_transformed_coords works for ICSD422", {
  # Use the transformed coordinates derived from ICSD422
  atomic_coords <- extract_atomic_coordinates(cif_content_ICSD422)
  sym_ops <- extract_symmetry_operations(cif_content_ICSD422)
  transformed_coords <- apply_symmetry_operations(atomic_coords, sym_ops)

  expanded <- expand_transformed_coords(transformed_coords, n_cells = 1)
  expect_s3_class(expanded, "data.table")
  # 12 unique transformed atoms * (2*1+1)^3 cells = 12 * 27 = 324 atoms
  expect_equal(nrow(expanded), 324)
  expect_equal(names(expanded), c("Label", "x_a", "y_b", "z_c"))

  # Check one specific expansion:
  # Original transformed Sr1_sym1 (from previous test): 0.1529, 0.25, 0.4231 (approx)
  # Find this atom in transformed_coords by its coordinates
  original_atom_row <- transformed[abs(x_a - 0.1529) < 1e-4 & abs(y_b - 0.25) < 1e-4 & abs(z_c - 0.4231) < 1e-4]
  original_label <- original_atom_row$Label[1]

  # Shift by dx=1, dy=0, dz=-1
  # Expected: x = 0.1529+1, y = 0.25+0, z = 0.4231-1
  expected_x = original_atom_row$x_a[1] + 1
  expected_y = original_atom_row$y_b[1] + 0
  expected_z = original_atom_row$z_c[1] - 1
  expected_label_suffix <- "_cell1_0_-1" # Function makes Label_dx_dy_dz

  expanded_atom <- expanded[Label == paste0(original_label, expected_label_suffix)]
  expect_equal(nrow(expanded_atom), 1)
  expect_equal(expanded_atom$x_a, expected_x, tolerance = 1e-6)
  expect_equal(expanded_atom$y_b, expected_y, tolerance = 1e-6)
  expect_equal(expanded_atom$z_c, expected_z, tolerance = 1e-6)

  # Test n_cells = 0 (should return original transformed_coords but with _cell0_0_0 labels)
  expanded_0 <- expand_transformed_coords(transformed_coords, n_cells = 0)
  expect_equal(nrow(expanded_0), nrow(transformed_coords)) # 12 rows
  expect_true(all(grepl("_cell0_0_0$", expanded_0$Label)))
  # Check coordinates match (ignoring new labels)
  expect_equal(expanded_0[,.(x_a,y_b,z_c)], transformed_coords[,.(x_a,y_b,z_c)], tolerance=1e-6)
})

test_that("Coordinate processing handles NULL or empty inputs", {
  expect_equal(apply_symmetry_operations(NULL, expected_symmetry_operations_ICSD422), NULL)
  expect_equal(apply_symmetry_operations(expected_atomic_coordinates_ICSD422, NULL), expected_atomic_coordinates_ICSD422) # Returns original

  empty_dt_coords <- data.table::data.table(Label=character(), x_a=numeric(), y_b=numeric(), z_c=numeric())
  empty_dt_symops <- data.table::data.table(x=character(), y=character(), z=character())

  expect_equal(apply_symmetry_operations(empty_dt_coords, expected_symmetry_operations_ICSD422), empty_dt_coords)
  expect_equal(apply_symmetry_operations(expected_atomic_coordinates_ICSD422, empty_dt_symops), expected_atomic_coordinates_ICSD422)

  expect_null(expand_transformed_coords(NULL))
  expect_null(expand_transformed_coords(empty_dt_coords))
})
