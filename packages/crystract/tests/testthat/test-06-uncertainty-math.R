test_that("cell transformation derivatives match central differences", {
  cell <- data.table::data.table(
    `_cell_length_a` = 7, `_cell_length_b` = 8, `_cell_length_c` = 9,
    `_cell_angle_alpha` = 79, `_cell_angle_beta` = 83, `_cell_angle_gamma` = 74
  )
  geometry <- uncertainty_cell(cell)
  for (name in names(geometry$dM)) {
    column <- if (name %in% c("a", "b", "c"))
      paste0("_cell_length_", name) else paste0("_cell_angle_", name)
    step <- 1e-5
    shifted <- data.table::copy(cell)
    increment <- if (name %in% c("a", "b", "c")) step else step * 180 / pi
    data.table::set(shifted, 1L, column, cell[[column]] + increment)
    plus <- uncertainty_cell(shifted)$M
    data.table::set(shifted, 1L, column, cell[[column]] - increment)
    minus <- uncertainty_cell(shifted)$M
    expect_equal(geometry$dM[[name]], (plus - minus) / (2 * step), tolerance = 1e-7)
  }
})

test_that("distance propagation combines symmetry-related inputs before squaring", {
  cell <- data.table::data.table(
    `_cell_length_a` = 10, `_cell_length_b` = 10, `_cell_length_c` = 10,
    `_cell_angle_alpha` = 90, `_cell_angle_beta` = 90, `_cell_angle_gamma` = 90,
    `_cell_length_a_error` = 0, `_cell_length_b_error` = 0,
    `_cell_length_c_error` = 0, `_cell_angle_alpha_error` = 0,
    `_cell_angle_beta_error` = 0, `_cell_angle_gamma_error` = 0
  )
  atoms <- data.table::data.table(Label = "A", x_a = 0.1, y_b = 0, z_c = 0,
                                  x_error = 0.01, y_error = 0, z_error = 0)
  ops <- data.table::data.table(x = c("x", "-x"), y = c("y", "-y"), z = c("z", "-z"))
  generated <- apply_symmetry_operations(atoms, ops, cell)
  images <- expand_transformed_coords(generated)
  bonds <- calculate_distances(generated[Label == "A_1"], images, cell)
  same_image <- bonds[Atom2 == "A_1_1_0_0"]
  inversion <- bonds[Atom2 == "A_2_0_0_0"]
  expect_equal(propagate_distance_error(same_image, generated[Label == "A_1"],
                                        cell, images, atoms)$DistanceError, 0, tolerance = 1e-10)
  expect_equal(propagate_distance_error(inversion, generated[Label == "A_1"],
                                        cell, images, atoms)$DistanceError, 0.2, tolerance = 1e-9)
})

test_that("angle propagation uses independent inputs and marks singular angles", {
  cell <- data.table::data.table(
    `_cell_length_a` = 10, `_cell_length_b` = 10, `_cell_length_c` = 10,
    `_cell_angle_alpha` = 90, `_cell_angle_beta` = 90, `_cell_angle_gamma` = 90,
    `_cell_length_a_error` = 0, `_cell_length_b_error` = 0,
    `_cell_length_c_error` = 0, `_cell_angle_alpha_error` = 0,
    `_cell_angle_beta_error` = 0, `_cell_angle_gamma_error` = 0
  )
  atoms <- data.table::data.table(
    Label = c("A", "B", "C"), x_a = c(0, 0.1, 0),
    y_b = c(0, 0, 0.1), z_c = 0,
    x_error = c(0, 0, 0.01), y_error = c(0, 0.01, 0), z_error = 0
  )
  angles <- data.table::data.table(CentralAtom = "A", Neighbor1 = "B",
                                   Neighbor2 = "C", Angle = 90)
  expect_equal(propagate_angle_error(angles, atoms, atoms, cell)$AngleError,
               sqrt(2) * 0.1 * 180 / pi, tolerance = 1e-9)
  atoms[Label == "C", `:=`(x_a = -0.1, y_b = 0)]
  expect_warning(expect_true(is.na(propagate_angle_error(angles, atoms, atoms, cell)$AngleError)),
                 "singular geometry")
})

test_that("full covariance changes the distance variance in the expected direction", {
  cell <- data.table::data.table(
    `_cell_length_a` = 10, `_cell_length_b` = 10, `_cell_length_c` = 10,
    `_cell_angle_alpha` = 90, `_cell_angle_beta` = 90, `_cell_angle_gamma` = 90,
    `_cell_length_a_error` = 0, `_cell_length_b_error` = 0,
    `_cell_length_c_error` = 0, `_cell_angle_alpha_error` = 0,
    `_cell_angle_beta_error` = 0, `_cell_angle_gamma_error` = 0
  )
  atoms <- data.table::data.table(Label = c("A", "B"),
                                  x_a = c(0, 0.1), y_b = 0, z_c = 0,
                                  x_error = 0.01, y_error = 0, z_error = 0)
  bond <- data.table::data.table(Atom1 = "A", Atom2 = "B", Distance = 1)
  keys <- c("cell_a", "A:x", "B:x")
  Sigma <- matrix(0, 3, 3, dimnames = list(keys, keys))
  Sigma[c("A:x", "B:x"), c("A:x", "B:x")] <- 0.01^2
  expect_equal(propagate_distance_error(bond, atoms, cell, covariance = Sigma)$DistanceError,
               0, tolerance = 1e-10)
})

test_that("missing relevant uncertainty is reported as unavailable", {
  cell <- data.table::data.table(
    `_cell_length_a` = 10, `_cell_length_b` = 10, `_cell_length_c` = 10,
    `_cell_angle_alpha` = 90, `_cell_angle_beta` = 90, `_cell_angle_gamma` = 90
  )
  atoms <- data.table::data.table(Label = c("A", "B"), x_a = c(0, 0.1),
                                  y_b = 0, z_c = 0, x_error = NA_real_,
                                  y_error = NA_real_, z_error = NA_real_)
  bond <- data.table::data.table(Atom1 = "A", Atom2 = "B", Distance = 1)
  expect_warning(result <- propagate_distance_error(bond, atoms, cell),
                 "unavailable")
  expect_true(is.na(result$DistanceError))
})

test_that("triclinic propagated uncertainties agree with Monte Carlo", {
  cell <- data.table::data.table(
    `_cell_length_a` = 7, `_cell_length_b` = 8, `_cell_length_c` = 9,
    `_cell_angle_alpha` = 79, `_cell_angle_beta` = 83, `_cell_angle_gamma` = 74,
    `_cell_length_a_error` = 0.002, `_cell_length_b_error` = 0.002,
    `_cell_length_c_error` = 0.002, `_cell_angle_alpha_error` = 0.01,
    `_cell_angle_beta_error` = 0.01, `_cell_angle_gamma_error` = 0.01
  )
  atoms <- data.table::data.table(
    Label = c("A", "B"), x_a = c(0.1, 0.2), y_b = c(0.05, 0.12),
    z_c = c(0.07, 0.3), x_error = c(0.001, 0.001),
    y_error = c(0.001, 0.001), z_error = c(0.001, 0.001)
  )
  ops <- data.table::data.table(x = c("x", "-x"), y = c("y", "-y"), z = c("z", "-z"))
  images <- expand_transformed_coords(apply_symmetry_operations(atoms, ops, cell), c(0, 0, 0))
  bond <- calculate_distances(atoms[Label == "A"], images[Label == "B_2_0_0_0"], cell)
  angle <- data.table::data.table(CentralAtom = "A", Neighbor1 = "B_1_0_0_0",
                                  Neighbor2 = "B_2_0_0_0", Angle = 0)
  analytical_d <- propagate_distance_error(bond, atoms, cell, images)$DistanceError
  analytical_angle <- propagate_angle_error(angle, atoms, images, cell)$AngleError

  mean_input <- c(7, 8, 9, c(79, 83, 74) * pi / 180,
                  0.1, 0.05, 0.07, 0.2, 0.12, 0.3)
  sd_input <- c(rep(0.002, 3), rep(0.01 * pi / 180, 3), rep(0.001, 6))
  set.seed(10)
  draws <- matrix(stats::rnorm(4000 * 12, mean_input, sd_input), ncol = 12, byrow = TRUE)
  simulated <- t(apply(draws, 1, function(p) {
    ca <- cos(p[4]); cb <- cos(p[5]); cg <- cos(p[6]); sg <- sin(p[6])
    D <- 1 - ca^2 - cb^2 - cg^2 + 2 * ca * cb * cg
    M <- rbind(c(p[1], p[2] * cg, p[3] * cb),
               c(0, p[2] * sg, p[3] * (ca - cb * cg) / sg),
               c(0, 0, p[3] * sqrt(D) / sg))
    center <- as.numeric(M %*% p[7:9])
    first <- as.numeric(M %*% p[10:12]) - center
    second <- as.numeric(M %*% (1 - p[10:12])) - center
    c(distance = sqrt(sum(second^2)),
      angle = acos(sum(first * second) / sqrt(sum(first^2) * sum(second^2))) * 180 / pi)
  }))
  expect_lt(abs(analytical_d / stats::sd(simulated[, "distance"]) - 1), 0.06)
  expect_lt(abs(analytical_angle / stats::sd(simulated[, "angle"]) - 1), 0.06)
})
