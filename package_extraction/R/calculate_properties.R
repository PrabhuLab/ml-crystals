#' @title Calculate Interatomic Distances
#' @description Computes the distances between a central set of atoms and an
#'   expanded set, using the metric tensor for accuracy.
#' @param atomic_coordinates A `data.table` of the primary (asymmetric) atom set.
#' @param expanded_coords A `data.table` of atoms in the expanded supercell.
#' @param unit_cell_metrics A `data.table` with cell parameters.
#' @return A `data.table` of all non-zero distances.
#' @family property calculators
#' @export
#' @examples
#' # Example setup
#' cif_file <- system.file("extdata", "ICSD422.cif", package = "crysmal")
#' if (file.exists(cif_file)) {
#'   cif_content <- read_cif_files(cif_file)[[1]]
#'   atoms <- extract_atomic_coordinates(cif_content)
#'   metrics <- extract_unit_cell_metrics(cif_content)
#'   sym_ops <- extract_symmetry_operations(cif_content)
#'   full_cell <- apply_symmetry_operations(atoms, sym_ops)
#'   super_cell <- expand_transformed_coords(full_cell)
#'
#'   dists <- calculate_distances(atoms, super_cell, metrics)
#'   print(head(dists))
#' }
calculate_distances <- function(atomic_coordinates,
                                expanded_coords,
                                unit_cell_metrics) {
  a <- unit_cell_metrics$`_cell_length_a`
  b <- unit_cell_metrics$`_cell_length_b`
  c <- unit_cell_metrics$`_cell_length_c`
  alpha <- unit_cell_metrics$`_cell_angle_alpha` * pi / 180
  beta <- unit_cell_metrics$`_cell_angle_beta` * pi / 180
  gamma <- unit_cell_metrics$`_cell_angle_gamma` * pi / 180

  coords_atomic <- as.matrix(atomic_coordinates[, .(x_a, y_b, z_c)])
  coords_expanded <- as.matrix(expanded_coords[, .(x_a, y_b, z_c)])
  labels_atomic <- atomic_coordinates$Label
  labels_expanded <- expanded_coords$Label

  delta_x <- outer(coords_atomic[, 1], coords_expanded[, 1], "-")
  delta_y <- outer(coords_atomic[, 2], coords_expanded[, 2], "-")
  delta_z <- outer(coords_atomic[, 3], coords_expanded[, 3], "-")

  cos_alpha <- cos(alpha)
  cos_beta <- cos(beta)
  cos_gamma <- cos(gamma)

  r2 <- (a^2 * delta_x^2) + (b^2 * delta_y^2) + (c^2 * delta_z^2) +
    (2 * b * c * cos_alpha * delta_y * delta_z) +
    (2 * c * a * cos_beta * delta_z * delta_x) +
    (2 * a * b * cos_gamma * delta_x * delta_y)

  r <- sqrt(r2)

  atom_pairs <- expand.grid(Atom1 = labels_atomic,
                            Atom2 = labels_expanded,
                            KEEP.OUT.ATTRS = TRUE)
  distances <- data.table(
    Atom1 = atom_pairs$Atom1,
    Atom2 = atom_pairs$Atom2,
    Distance = as.vector(r),
    DeltaX = as.vector(delta_x),
    DeltaY = as.vector(delta_y),
    DeltaZ = as.vector(delta_z)
  )

  distances <- distances[Distance > 1e-6]

  return(distances)
}

#' @title Calculate Neighbor Counts
#' @description Counts the number of nearest neighbors for each central atom based on
#'   a table of bonded pairs.
#' @param bonded_pairs_table A `data.table` of bonded pairs.
#' @return A `data.table` with columns 'Atom' and 'NeighborCount'.
#' @family property calculators
#' @export
calculate_neighbor_counts <- function(bonded_pairs_table) {
  if (is.null(bonded_pairs_table) || nrow(bonded_pairs_table) == 0) {
    return(data.table(Atom = character(), NeighborCount = integer()))
  }

  neighbor_counts <- bonded_pairs_table[, .(NeighborCount = .N), by = .(Atom1)]

  setnames(neighbor_counts, "Atom1", "Atom")

  return(neighbor_counts)
}

#' @title Calculate Bond Angles
#' @description Calculates all bond angles centered on each atom, formed by pairs
#'   of its bonded neighbors.
#' @param bonded_pairs Data.table of bonded atoms.
#' @param atomic_coordinates Data.table of asymmetric atom coordinates.
#' @param expanded_coords Data.table of supercell atom coordinates.
#' @param unit_cell_metrics Data.table with unit cell parameters.
#' @return A `data.table` of all unique bond angles.
#' @family property calculators
#' @export
calculate_angles <- function(bonded_pairs,
                             atomic_coordinates,
                             expanded_coords,
                             unit_cell_metrics) {
  a <- unit_cell_metrics$`_cell_length_a`
  b <- unit_cell_metrics$`_cell_length_b`
  c <- unit_cell_metrics$`_cell_length_c`
  alpha_rad <- unit_cell_metrics$`_cell_angle_alpha` * pi / 180
  beta_rad <- unit_cell_metrics$`_cell_angle_beta` * pi / 180
  gamma_rad <- unit_cell_metrics$`_cell_angle_gamma` * pi / 180

  cos_alpha <- cos(alpha_rad)
  cos_beta <- cos(beta_rad)
  cos_gamma <- cos(gamma_rad)

  all_coords <- unique(rbind(atomic_coordinates[, .(Label, x_a, y_b, z_c)], expanded_coords[, .(Label, x_a, y_b, z_c)], fill = TRUE),
                       by = "Label")
  setkey(all_coords, Label)

  calculate_angle_metric <- function(atom1_label, atom2_label, atom3_label) {
    coord1 <- all_coords[atom1_label, .(x_a, y_b, z_c)]
    coord2 <- all_coords[atom2_label, .(x_a, y_b, z_c)]
    coord3 <- all_coords[atom3_label, .(x_a, y_b, z_c)]

    if (anyNA(coord1) ||
        anyNA(coord2) || anyNA(coord3))
      return(NA_real_)

    v1_frac <- as.numeric(coord2 - coord1)
    v2_frac <- as.numeric(coord3 - coord1)
    xf1 <- v1_frac[1]
    yf1 <- v1_frac[2]
    zf1 <- v1_frac[3]
    xf2 <- v2_frac[1]
    yf2 <- v2_frac[2]
    zf2 <- v2_frac[3]

    dot_product <- (
      xf1 * xf2 * a^2 + yf1 * yf2 * b^2 + zf1 * zf2 * c^2 + (xf1 * yf2 + yf1 *
                                                               xf2) * a * b * cos_gamma + (xf1 * zf2 + zf1 * xf2) * a * c * cos_beta + (yf1 *
                                                                                                                                          zf2 + zf1 * yf2) * b * c * cos_alpha
    )

    mag_sq1 <- (
      xf1^2 * a^2 + yf1^2 * b^2 + zf1^2 * c^2 + 2 * xf1 * yf1 * a * b * cos_gamma + 2 *
        xf1 * zf1 * a * c * cos_beta + 2 * yf1 * zf1 * b * c * cos_alpha
    )
    mag_sq2 <- (
      xf2^2 * a^2 + yf2^2 * b^2 + zf2^2 * c^2 + 2 * xf2 * yf2 * a * b * cos_gamma + 2 *
        xf2 * zf2 * a * c * cos_beta + 2 * yf2 * zf2 * b * c * cos_alpha
    )

    if (mag_sq1 <= 1e-10 || mag_sq2 <= 1e-10)
      return(NA_real_)

    mag1 <- sqrt(mag_sq1)
    mag2 <- sqrt(mag_sq2)
    cos_theta <- min(max(dot_product / (mag1 * mag2), -1.0), 1.0)
    return(acos(cos_theta) * 180 / pi)
  }

  angle_list <- list()
  unique_central_atoms <- unique(bonded_pairs$Atom1)

  for (central_atom in unique_central_atoms) {
    bonded_neighbors <- bonded_pairs[Atom1 == central_atom, Atom2]
    if (length(bonded_neighbors) >= 2) {
      neighbor_combinations <- combn(bonded_neighbors, 2, simplify = FALSE)
      for (pair in neighbor_combinations) {
        angle <- calculate_angle_metric(central_atom, pair[1], pair[2])
        if (!is.na(angle)) {
          angle_list[[length(angle_list) + 1]] <- data.table(
            CentralAtom = central_atom,
            Neighbor1 = pair[1],
            Neighbor2 = pair[2],
            Angle = angle
          )
        }
      }
    }
  }

  if (length(angle_list) > 0) {
    angle_results <- rbindlist(angle_list)
    return(angle_results[order(CentralAtom, Neighbor1, Neighbor2)])
  } else {
    return(
      data.table(
        CentralAtom = character(),
        Neighbor1 = character(),
        Neighbor2 = character(),
        Angle = numeric()
      )
    )
  }
}

#' @title Propagate Distance Error
#' @description Calculates the standard uncertainty for each interatomic distance.
#' @param bonded_pairs Data.table of bonded atoms with their distances.
#' @param atomic_coordinates Data.table with fractional coordinates and errors.
#' @param unit_cell_metrics Data.table with unit cell parameters and errors.
#' @return The input `bonded_pairs` data.table with a new 'DistanceError' column.
#' @family error propagators
#' @export
propagate_distance_error <- function(bonded_pairs,
                                     atomic_coordinates,
                                     unit_cell_metrics) {
  if (is.null(bonded_pairs) ||
      nrow(bonded_pairs) == 0) {
    return(bonded_pairs[, DistanceError := NA_real_])
  }

  a <- unit_cell_metrics$`_cell_length_a`
  s_a <- unit_cell_metrics$`_cell_length_a_error`
  b <- unit_cell_metrics$`_cell_length_b`
  s_b <- unit_cell_metrics$`_cell_length_b_error`
  c <- unit_cell_metrics$`_cell_length_c`
  s_c <- unit_cell_metrics$`_cell_length_c_error`
  alpha_rad <- unit_cell_metrics$`_cell_angle_alpha` * pi / 180
  s_alpha_rad <- unit_cell_metrics$`_cell_angle_alpha_error` * pi / 180
  beta_rad <- unit_cell_metrics$`_cell_angle_beta` * pi / 180
  s_beta_rad <- unit_cell_metrics$`_cell_angle_beta_error` * pi / 180
  gamma_rad <- unit_cell_metrics$`_cell_angle_gamma` * pi / 180
  s_gamma_rad <- unit_cell_metrics$`_cell_angle_gamma_error` * pi / 180

  cos_a <- cos(alpha_rad)
  sin_a <- sin(alpha_rad)
  cos_b <- cos(beta_rad)
  sin_b <- sin(beta_rad)
  cos_g <- cos(gamma_rad)
  sin_g <- sin(gamma_rad)

  s_a <- ifelse(is.na(s_a), 0, s_a)
  s_b <- ifelse(is.na(s_b), 0, s_b)
  s_c <- ifelse(is.na(s_c), 0, s_c)
  s_alpha_rad <- ifelse(is.na(s_alpha_rad), 0, s_alpha_rad)
  s_beta_rad <- ifelse(is.na(s_beta_rad), 0, s_beta_rad)
  s_gamma_rad <- ifelse(is.na(s_gamma_rad), 0, s_gamma_rad)

  work_dt <- copy(bonded_pairs)
  setnames(work_dt, c("DeltaX", "DeltaY", "DeltaZ"), c("dx", "dy", "dz"))
  atom_errors <- atomic_coordinates[, .(Label,
                                        s_xf = x_error,
                                        s_yf = y_error,
                                        s_zf = z_error)]

  work_dt <- merge(
    work_dt,
    atom_errors,
    by.x = "Atom1",
    by.y = "Label",
    all.x = TRUE
  )
  setnames(work_dt,
           c("s_xf", "s_yf", "s_zf"),
           c("s_xf1", "s_yf1", "s_zf1"))
  work_dt[, Original_Atom2 := sub("_.*", "", Atom2)]
  work_dt <- merge(
    work_dt,
    atom_errors,
    by.x = "Original_Atom2",
    by.y = "Label",
    all.x = TRUE
  )
  setnames(work_dt,
           c("s_xf", "s_yf", "s_zf"),
           c("s_xf2", "s_yf2", "s_zf2"))
  work_dt[, Original_Atom2 := NULL]

  coord_err_cols <- c("s_xf1", "s_yf1", "s_zf1", "s_xf2", "s_yf2", "s_zf2")
  for (col in coord_err_cols)
    set(work_dt, which(is.na(work_dt[[col]])), col, 0)

  work_dt[, term_common := 1 / (2 * Distance)]
  work_dt[, `:=`(
    pd_a = (term_common * (
      2 * a * dx^2 + 2 * b * dx * dy * cos_g + 2 * c * dx * dz * cos_b
    )),
    pd_b = (term_common * (
      2 * b * dy^2 + 2 * a * dx * dy * cos_g + 2 * c * dy * dz * cos_a
    )),
    pd_c = (term_common * (
      2 * c * dz^2 + 2 * a * dx * dz * cos_b + 2 * b * dy * dz * cos_a
    )),
    pd_alpha = (term_common * (-2 * b * c * dy * dz * sin_a)),
    pd_beta = (term_common * (-2 * a * c * dx * dz * sin_b)),
    pd_gamma = (term_common * (-2 * a * b * dx * dy * sin_g))
  )]
  work_dt[, `:=`(
    pd_xf1 = (term_common * (
      2 * a^2 * dx + 2 * a * b * dy * cos_g + 2 * a * c * dz * cos_b
    )),
    pd_yf1 = (term_common * (
      2 * b^2 * dy + 2 * a * b * dx * cos_g + 2 * b * c * dz * cos_a
    )),
    pd_zf1 = (term_common * (
      2 * c^2 * dz + 2 * a * c * dx * cos_b + 2 * b * c * dy * cos_a
    ))
  )][, `:=`(
    pd_xf2 = -pd_xf1,
    pd_yf2 = -pd_yf1,
    pd_zf2 = -pd_zf1
  )]
  work_dt[, variance := (pd_a * s_a)^2 + (pd_b * s_b)^2 + (pd_c * s_c)^2 + (pd_alpha *
                                                                              s_alpha_rad)^2 + (pd_beta * s_beta_rad)^2 + (pd_gamma * s_gamma_rad)^2 + (pd_xf1 *
                                                                                                                                                          s_xf1)^2 + (pd_yf1 * s_yf1)^2 + (pd_zf1 * s_zf1)^2 + (pd_xf2 * s_xf2)^2 + (pd_yf2 *
                                                                                                                                                                                                                                       s_yf2)^2 + (pd_zf2 * s_zf2)^2]
  work_dt[, DistanceError := sqrt(variance)]
  return(merge(
    bonded_pairs,
    work_dt[, .(Atom1, Atom2, Distance, DistanceError)],
    by = c("Atom1", "Atom2", "Distance"),
    all.x = TRUE
  ))
}

#' @title Propagate Angle Error
#' @description Calculates the standard uncertainty for each bond angle.
#' @param bond_angles Data.table of calculated bond angles.
#' @param atomic_coordinates Data.table with fractional coordinates and errors.
#' @param expanded_coords Data.table of supercell atom coordinates.
#' @param unit_cell_metrics Data.table with unit cell parameters and errors.
#' @return The input `bond_angles` data.table with a new 'AngleError' column.
#' @family error propagators
#' @export
propagate_angle_error <- function(bond_angles,
                                  atomic_coordinates,
                                  expanded_coords,
                                  unit_cell_metrics) {
  if (is.null(bond_angles) ||
      nrow(bond_angles) == 0) {
    return(bond_angles[, AngleError := NA_real_])
  }
  a <- unit_cell_metrics$`_cell_length_a`
  b <- unit_cell_metrics$`_cell_length_b`
  c <- unit_cell_metrics$`_cell_length_c`
  s_a <- ifelse(
    is.na(unit_cell_metrics$`_cell_length_a_error`),
    0,
    unit_cell_metrics$`_cell_length_a_error`
  )
  s_b <- ifelse(
    is.na(unit_cell_metrics$`_cell_length_b_error`),
    0,
    unit_cell_metrics$`_cell_length_b_error`
  )
  s_c <- ifelse(
    is.na(unit_cell_metrics$`_cell_length_c_error`),
    0,
    unit_cell_metrics$`_cell_length_c_error`
  )
  alpha_rad <- unit_cell_metrics$`_cell_angle_alpha` * pi / 180
  beta_rad <- unit_cell_metrics$`_cell_angle_beta` * pi / 180
  gamma_rad <- unit_cell_metrics$`_cell_angle_gamma` * pi / 180
  s_alpha_rad <- ifelse(
    is.na(unit_cell_metrics$`_cell_angle_alpha_error`),
    0,
    unit_cell_metrics$`_cell_angle_alpha_error` * pi / 180
  )
  s_beta_rad <- ifelse(
    is.na(unit_cell_metrics$`_cell_angle_beta_error`),
    0,
    unit_cell_metrics$`_cell_angle_beta_error` * pi / 180
  )
  s_gamma_rad <- ifelse(
    is.na(unit_cell_metrics$`_cell_angle_gamma_error`),
    0,
    unit_cell_metrics$`_cell_angle_gamma_error` * pi / 180
  )
  cos_a <- cos(alpha_rad)
  sin_a <- sin(alpha_rad)
  cos_b <- cos(beta_rad)
  sin_b <- sin(beta_rad)
  cos_g <- cos(gamma_rad)
  sin_g <- sin(gamma_rad)
  v_sq <- 1 - cos_a^2 - cos_b^2 - cos_g^2 + 2 * cos_a * cos_b * cos_g
  v <- sqrt(v_sq)
  cart_errors <- copy(atomic_coordinates)
  for (col in c("x_error", "y_error", "z_error"))
    set(cart_errors, which(is.na(cart_errors[[col]])), col, 0)
  cart_errors[, `:=`(
    p_xc_a = x_a,
    p_xc_b = y_b * cos_g,
    p_xc_c = z_c * cos_b,
    p_xc_alpha = 0,
    p_xc_beta = -c * z_c * sin_b,
    p_xc_gamma = -b * y_b * sin_g,
    p_xc_xf = a,
    p_xc_yf = b * cos_g,
    p_xc_zf = c * cos_b
  )]
  cart_errors[, `:=`(
    p_yc_a = 0,
    p_yc_b = y_b * sin_g,
    p_yc_c = z_c * (cos_a - cos_b * cos_g) / sin_g,
    p_yc_alpha = -c * z_c * sin_a / sin_g,
    p_yc_beta = c * z_c * sin_b * cos_g / sin_g,
    p_yc_gamma = b * y_b * cos_g + c * z_c * (cos_b - cos_a * cos_g) / sin_g^2,
    p_yc_xf = 0,
    p_yc_yf = b * sin_g,
    p_yc_zf = c * (cos_a - cos_b * cos_g) / sin_g
  )]
  cart_errors[, `:=`(
    p_zc_a = 0,
    p_zc_b = 0,
    p_zc_c = z_c * v / sin_g,
    p_zc_alpha = c * z_c * (cos_b * cos_g - cos_a) / (sin_g * v),
    p_zc_beta = c * z_c * (cos_a * cos_g - cos_b) / (sin_g * v),
    p_zc_gamma = -c * z_c * (v_sq * cos_g + sin_g^2 * (cos_a * cos_b - cos_g)) /
      (sin_g^2 * v),
    p_zc_xf = 0,
    p_zc_yf = 0,
    p_zc_zf = c * v / sin_g
  )]
  cart_errors[, s_xc_sq := (p_xc_a * s_a)^2 + (p_xc_b * s_b)^2 + (p_xc_c *
                                                                    s_c)^2 + (p_xc_alpha * s_alpha_rad)^2 + (p_xc_beta * s_beta_rad)^2 + (p_xc_gamma *
                                                                                                                                            s_gamma_rad)^2 + (p_xc_xf * x_error)^2 + (p_xc_yf * y_error)^2 + (p_xc_zf *
                                                                                                                                                                                                                z_error)^2]
  cart_errors[, s_yc_sq := (p_yc_a * s_a)^2 + (p_yc_b * s_b)^2 + (p_yc_c *
                                                                    s_c)^2 + (p_yc_alpha * s_alpha_rad)^2 + (p_yc_beta * s_beta_rad)^2 + (p_yc_gamma *
                                                                                                                                            s_gamma_rad)^2 + (p_yc_xf * x_error)^2 + (p_yc_yf * y_error)^2 + (p_yc_zf *
                                                                                                                                                                                                                z_error)^2]
  cart_errors[, s_zc_sq := (p_zc_a * s_a)^2 + (p_zc_b * s_b)^2 + (p_zc_c *
                                                                    s_c)^2 + (p_zc_alpha * s_alpha_rad)^2 + (p_zc_beta * s_beta_rad)^2 + (p_zc_gamma *
                                                                                                                                            s_gamma_rad)^2 + (p_zc_xf * x_error)^2 + (p_zc_yf * y_error)^2 + (p_zc_zf *
                                                                                                                                                                                                                z_error)^2]
  all_frac_coords <- unique(rbind(atomic_coordinates[, .(Label, x_a, y_b, z_c)], expanded_coords), by =
                              "Label")
  all_cart_coords <- all_frac_coords[, `:=`(
    xc = a * x_a + b * y_b * cos_g + c * z_c * cos_b,
    yc = b * y_b * sin_g + c * z_c * (cos_a - cos_b * cos_g) / sin_g,
    zc = c * z_c * v / sin_g
  )][, .(Label, xc, yc, zc)]
  setkey(all_cart_coords, Label)
  error_subset <- cart_errors[, .(Label, s_xc_sq, s_yc_sq, s_zc_sq)]
  setkey(error_subset, Label)
  work_dt <- copy(bond_angles)
  work_dt[all_cart_coords, on = .(CentralAtom = Label), `:=`(xc1 = i.xc, yc1 =
                                                               i.yc, zc1 = i.zc)]
  work_dt[all_cart_coords, on = .(Neighbor1 = Label), `:=`(xc2 = i.xc, yc2 =
                                                             i.yc, zc2 = i.zc)]
  work_dt[all_cart_coords, on = .(Neighbor2 = Label), `:=`(xc3 = i.xc, yc3 =
                                                             i.yc, zc3 = i.zc)]
  work_dt[, Parent_C := sub("_.*", "", CentralAtom)][, Parent_N1 := sub("_.*", "", Neighbor1)][, Parent_N2 :=
                                                                                                 sub("_.*", "", Neighbor2)]
  work_dt[error_subset, on = .(Parent_C = Label), `:=`(s_xc1_sq = i.s_xc_sq,
                                                       s_yc1_sq = i.s_yc_sq,
                                                       s_zc1_sq = i.s_zc_sq)]
  work_dt[error_subset, on = .(Parent_N1 = Label), `:=`(s_xc2_sq = i.s_xc_sq,
                                                        s_yc2_sq = i.s_yc_sq,
                                                        s_zc2_sq = i.s_zc_sq)]
  work_dt[error_subset, on = .(Parent_N2 = Label), `:=`(s_xc3_sq = i.s_xc_sq,
                                                        s_yc3_sq = i.s_yc_sq,
                                                        s_zc3_sq = i.s_zc_sq)]
  work_dt[, `:=`(
    a_vx = xc2 - xc1,
    a_vy = yc2 - yc1,
    a_vz = zc2 - zc1,
    b_vx = xc3 - xc1,
    b_vy = yc3 - yc1,
    b_vz = zc3 - zc1
  )]
  work_dt[, `:=`(
    mag_a = sqrt(a_vx^2 + a_vy^2 + a_vz^2),
    mag_b = sqrt(b_vx^2 + b_vy^2 + b_vz^2)
  )]
  work_dt <- work_dt[mag_a > 1e-9 & mag_b > 1e-9]
  work_dt[, C_val := pmin(pmax((a_vx * b_vx + a_vy * b_vy + a_vz * b_vz) /
                                 (mag_a * mag_b), -1.0), 1.0)]
  work_dt[, `:=`(
    p_C_xc1 = -(((
      b_vx / mag_b - C_val * a_vx / mag_a
    ) / mag_a) + ((
      a_vx / mag_a - C_val * b_vx / mag_b
    ) / mag_b)),
    p_C_yc1 = -(((
      b_vy / mag_b - C_val * a_vy / mag_a
    ) / mag_a) + ((
      a_vy / mag_a - C_val * b_vy / mag_b
    ) / mag_b)),
    p_C_zc1 = -(((
      b_vz / mag_b - C_val * a_vz / mag_a
    ) / mag_a) + ((
      a_vz / mag_a - C_val * b_vz / mag_b
    ) / mag_b))
  )]
  work_dt[, `:=`(
    p_C_xc2 = (b_vx / mag_b - C_val * a_vx / mag_a) / mag_a,
    p_C_yc2 = (b_vy / mag_b - C_val * a_vy / mag_a) / mag_a,
    p_C_zc2 = (b_vz / mag_b - C_val * a_vz / mag_a) / mag_a
  )]
  work_dt[, `:=`(
    p_C_xc3 = (a_vx / mag_a - C_val * b_vx / mag_b) / mag_b,
    p_C_yc3 = (a_vy / mag_a - C_val * b_vy / mag_b) / mag_b,
    p_C_zc3 = (a_vz / mag_a - C_val * b_vz / mag_b) / mag_b
  )]
  work_dt[, s_C_sq := (p_C_xc1^2 * s_xc1_sq) + (p_C_yc1^2 * s_yc1_sq) +
            (p_C_zc1^2 * s_zc1_sq) + (p_C_xc2^2 * s_xc2_sq) + (p_C_yc2^2 * s_yc2_sq) +
            (p_C_zc2^2 * s_zc2_sq) + (p_C_xc3^2 * s_xc3_sq) + (p_C_yc3^2 * s_yc3_sq) +
            (p_C_zc3^2 * s_zc3_sq)]
  work_dt[C_val^2 >= 1.0, s_theta_sq := 0][C_val^2 < 1.0, s_theta_sq := s_C_sq /
                                             (1 - C_val^2)]
  work_dt[, AngleError := sqrt(pmax(0, s_theta_sq)) * 180 / pi]
  return(merge(
    bond_angles,
    work_dt[, .(CentralAtom, Neighbor1, Neighbor2, AngleError)],
    by = c("CentralAtom", "Neighbor1", "Neighbor2"),
    all.x = TRUE
  ))
}

#' @title Identify Atomic Bonds using the Minimum Distance Method
#' @description Identifies bonded atoms by finding the nearest neighbor distance (d_min)
#'   for each central atom and defining a cutoff distance (d_cut) as
#'   d_cut = (1 + delta) * d_min.
#' @param distances A `data.table` of interatomic distances from `calculate_distances`.
#' @param delta The relative tolerance parameter (default 0.1 or ten percent).
#' @return A `data.table` of bonded pairs.
#' @family bonding algorithms
#' @export
#' @examples
#' cif_file <- system.file("extdata", "ICSD422.cif", package = "crysmal")
#' if (file.exists(cif_file)) {
#'   # Setup code
#'   cif_content <- read_cif_files(cif_file)[[1]]
#'   atoms <- extract_atomic_coordinates(cif_content)
#'   metrics <- extract_unit_cell_metrics(cif_content)
#'   sym_ops <- extract_symmetry_operations(cif_content)
#'   full_cell <- apply_symmetry_operations(atoms, sym_ops)
#'   super_cell <- expand_transformed_coords(full_cell)
#'   dists <- calculate_distances(atoms, super_cell, metrics)
#'
#'   # Apply Minimum Distance Method
#'   bonded <- minimum_distance(dists, delta = 0.1)
#'   print(head(bonded))
#' }
minimum_distance <- function(distances, delta = 0.1) {
  dmin <- distances[, .(dmin = min(Distance)), by = .(Atom1)]
  dmin[, dcut := (1 + delta) * dmin]
  bonded_pairs <- distances[dmin, on = .(Atom1), allow.cartesian = TRUE][Distance <= dcut, .(Atom1, Atom2, Distance, DeltaX, DeltaY, DeltaZ, dcut, dmin)]
  return(bonded_pairs)
}

#' @title Identify Atomic Bonds using Brunner's Method
#' @description (In Progress) An alternative bonding detection method using the
#'   largest reciprocal gap.
#' @param distances A `data.table` of interatomic distances.
#' @param delta A small distance offset (default 0.0001).
#' @return A `data.table` of bonded pairs.
#' @family bonding algorithms
#' @export
brunner <- function(distances, delta = 0.0001) {
  bonds <- list()
  unique_atoms <- unique(distances$Atom1)
  for (atom in unique_atoms) {
    atom_distances <- distances[Atom1 == atom][order(Distance)]
    largest_gap <- -Inf
    j_max <- NA
    for (j in 1:(nrow(atom_distances) - 1)) {
      reciprocal_gap <- 1 / atom_distances$Distance[j] - 1 / atom_distances$Distance[j + 1]
      if (reciprocal_gap > largest_gap) {
        largest_gap <- reciprocal_gap
        j_max <- j
      }
    }
    d_cut <- atom_distances$Distance[j_max] + delta
    bonds[[atom]] <- atom_distances[Distance <= d_cut]
  }
  return(rbindlist(bonds, fill = TRUE))
}

#' @title Identify Atomic Bonds using Hoppe's Method
#' @description (In Progress) An alternative bonding detection method using
#'   Effective Coordination Numbers.
#' @param distances A `data.table` of interatomic distances.
#' @param delta A bond strength cutoff (default 0.5).
#' @param tolerance The convergence tolerance for iterative calculation.
#' @return A `data.table` of bonded pairs.
#' @family bonding algorithms
#' @export
hoppe <- function(distances,
                  delta = 0.5,
                  tolerance = 0.001) {
  bonded_pairs <- list()
  unique_atoms <- unique(distances$Atom1)
  for (atom in unique_atoms) {
    atom_distances <- distances[Atom1 == atom]
    dmin <- min(atom_distances$Distance)
    davg <- sum(atom_distances$Distance * exp(1 - (atom_distances$Distance /
                                                     dmin)^6)) / sum(exp(1 - (atom_distances$Distance / dmin)^6))
    while (TRUE) {
      prev_davg <- davg
      davg <- sum(atom_distances$Distance * exp(1 - (
        atom_distances$Distance / prev_davg
      )^6)) / sum(exp(1 - (
        atom_distances$Distance / prev_davg
      )^6))
      if (abs(davg - prev_davg) <= tolerance)
        break
    }
    atom_distances[, BondStrength := exp(1 - (Distance / davg)^6)]
    bonded_pairs[[atom]] <- atom_distances[BondStrength >= delta, .(Atom1, Atom2, Distance)]
  }
  return(rbindlist(bonded_pairs, fill = TRUE))
}

#' @title Filter Data by Atom Symbol Interactively
#' @description Prompts the user to select chemical elements to keep in a data table
#'   of bonds or angles. Filtering is based on matching the base chemical symbol
#'   in a specified column (e.g., "CentralAtom").
#'
#' @details
#' The function first identifies all unique base chemical symbols from the atom
#' labels in the specified column (e.g., it extracts 'C' from 'C1', 'Si' from 'Si2_1').
#' It then presents these symbols to the user and asks for a comma-separated list
#' of the symbols they wish to retain.
#'
#' The matching logic is designed to be specific to avoid ambiguity between elements.
#' For example, if the user enters 'C', the function will match labels like 'C1',
#' 'C2', 'C_10', or a lone 'C'. However, it will *not* match labels for different
#' elements that start with C, such as 'Cr1' or 'Ca2'. This is achieved by
#' constructing a regular expression that ensures the character(s) immediately
#' following the selected symbol are not alphabetical letters.
#'
#' This function is intended for interactive use.
#'
#' @param data_table A `data.table` object containing atomic information, such as
#'   the output from `calculate_angles` or `minimum_distance`.
#' @param atom_col A character string specifying the name of the column in
#'   `data_table` that contains the atom labels to filter by. Defaults to "CentralAtom".
#' @return A `data.table` filtered to include only the rows where the atom
#'   label in `atom_col` corresponds to one of the user-selected chemical symbols.
#'   If the user provides no input, an empty `data.table` is returned.
#' @family property calculators
#' @export
#' @examples
#' # This example demonstrates how the function would work.
#' # In a real session, you would call: filter_atoms_by_symbol(my_data)
#'
#' # 1. Create a sample data.table of bond angles
#' sample_angles <- data.table::data.table(
#'   CentralAtom = c("C1", "C2", "Si1", "Cr1", "O1", "O2", "C"),
#'   Neighbor1 = c("O1", "O2", "O1", "N1", "C1", "C2", "H1"),
#'   Neighbor2 = c("H1", "H2", "O2", "N2", "H3", "H4", "H2"),
#'   Angle = c(109.5, 109.5, 120.0, 90.0, 104.5, 104.5, 120)
#' )
#'
#' # 2. In an interactive R session, the function would prompt the user.
#' # For example, if the user sees the available symbols (C, Si, Cr, O) and
#' # enters "C,Si" at the prompt, the function would return the rows for
#' # "C1", "C2", "C", and "Si1".
#'
#' if (interactive()) {
#'   filtered_data <- filter_atoms_by_symbol(sample_angles, atom_col = "CentralAtom")
#'   print(filtered_data)
#' }
#'
filter_atoms_by_symbol <- function(data_table, atom_col = "CentralAtom") {
  # Check if the column exists in the data.table
  if (!atom_col %in% names(data_table)) {
    stop(paste("Column '", atom_col, "' not found in the data table.", sep = ""))
  }

  # Extract unique base chemical symbols from the specified column to show the user
  unique_labels <- unique(data_table[[atom_col]])
  # This regex extracts the initial sequence of one or more letters
  base_symbols <- sort(unique(stringr::str_extract(unique_labels, "^[A-Za-z]+")))

  if (length(base_symbols) == 0) {
    message("No recognizable atom symbols found in the specified column.")
    return(data_table) # Return original table if no symbols are found
  }

  # Interactively prompt the user for their choice of atoms
  cat("Available base atom symbols found:",
      paste(base_symbols, collapse = ", "),
      "\n")
  cat(
    "Please enter the chemical symbols you want to filter for, separated by commas (e.g., C,Si,O):\n"
  )
  user_input <- readline()

  # Process the user's input string
  symbols_to_keep <- trimws(strsplit(user_input, ",")[[1]])
  # Remove any empty strings that might result from trailing commas or empty input
  symbols_to_keep <- symbols_to_keep[symbols_to_keep != ""]

  # If the user did not provide any valid symbols, return an empty data.table
  if (length(symbols_to_keep) == 0) {
    message("No symbols entered. Returning an empty data table.")
    return(data_table[0, ])
  }

  # Construct the regular expression pattern for grepl
  # For each symbol (e.g., "C"), the pattern is (^C$)|(^C[^A-Za-z])
  # This matches the exact symbol OR the symbol followed by a non-alphabetic character.
  # This correctly distinguishes 'C' from 'Cr', for example.
  patterns <- sapply(symbols_to_keep, function(sym) {
    paste0("(^", sym, "$)|(^", sym, "[^A-Za-z])")
  })
  full_pattern <- paste(patterns, collapse = "|")

  # Filter the data.table using the generated regex pattern on the specified column
  filtered_dt <- data_table[grepl(full_pattern, get(atom_col))]

  return(filtered_dt)
}
