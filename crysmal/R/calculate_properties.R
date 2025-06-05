#' Calculate Interatomic Distances
#'
#' Calculates distances between atoms in `atomic_coordinates` (reference set) and
#' `expanded_coords` (search set) using the unit cell metric tensor.
#'
#' @param atomic_coordinates A data.table of reference atomic coordinates (Label, x_a, y_b, z_c).
#' @param expanded_coords A data.table of search atomic coordinates (Label, x_a, y_b, z_c),
#'                        typically including original atoms and atoms in neighboring cells.
#' @param unit_cell_metrics A data.table or list with unit cell parameters:
#'                          `_cell_length_a`, `_cell_length_b`, `_cell_length_c`,
#'                          `_cell_angle_alpha`, `_cell_angle_beta`, `_cell_angle_gamma`.
#' @return A data.table with Atom1 (from `atomic_coordinates`), Atom2 (from `expanded_coords`),
#'         Distance, and fractional coordinate differences (DeltaX, DeltaY, DeltaZ), and cosines of cell angles.
#'         Returns NULL if inputs are invalid.
#' @export
#' @examples
#' # ucm <- data.table("_cell_length_a"=10, "_cell_length_b"=10, "_cell_length_c"=10,
#' #                   "_cell_angle_alpha"=90, "_cell_angle_beta"=90, "_cell_angle_gamma"=90)
#' # ac <- data.table(Label="A1", x_a=0.1, y_b=0.1, z_c=0.1)
#' # ec <- data.table(Label=c("A1_exp","A2_exp"), x_a=c(0.1,0.5), y_b=c(0.1,0.5), z_c=c(0.1,0.5))
#' # calculate_distances(ac, ec, ucm)
calculate_distances <- function(atomic_coordinates, expanded_coords, unit_cell_metrics) {
  if (is.null(atomic_coordinates) || nrow(atomic_coordinates) == 0 ||
      is.null(expanded_coords) || nrow(expanded_coords) == 0 ||
      is.null(unit_cell_metrics) || nrow(unit_cell_metrics) == 0) {
    return(NULL)
  }

  # Extract unit cell parameters
  a <- unit_cell_metrics[["_cell_length_a"]]
  b <- unit_cell_metrics[["_cell_length_b"]]
  c <- unit_cell_metrics[["_cell_length_c"]]
  alpha_deg <- unit_cell_metrics[["_cell_angle_alpha"]]
  beta_deg <- unit_cell_metrics[["_cell_angle_beta"]]
  gamma_deg <- unit_cell_metrics[["_cell_angle_gamma"]]

  if (any(is.na(c(a, b, c, alpha_deg, beta_deg, gamma_deg)))) {
    warning("One or more unit cell parameters are NA. Cannot calculate distances.")
    return(NULL)
  }

  alpha <- alpha_deg * pi / 180 # Convert degrees to radians
  beta <- beta_deg * pi / 180
  gamma <- gamma_deg * pi / 180

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

  tol <- 1e-10 # Tolerance for checking if cosine is effectively zero
  if (abs(cos_alpha) < tol) cos_alpha <- 0
  if (abs(cos_beta) < tol) cos_beta <- 0
  if (abs(cos_gamma) < tol) cos_gamma <- 0

  r2 <- (a^2 * delta_x^2) + (b^2 * delta_y^2) + (c^2 * delta_z^2) +
    (2 * b * c * cos_alpha * delta_y * delta_z) +
    (2 * c * a * cos_beta * delta_z * delta_x) +
    (2 * a * b * cos_gamma * delta_x * delta_y)

  # Handle potential negative r2 due to floating point issues if distances are very small
  r2[r2 < 0 & abs(r2) < tol] <- 0
  if (any(r2 < 0)) {
    warning("Negative squared distances encountered; check coordinates or cell parameters.")
  }
  r <- sqrt(r2)

  atom_pairs <- CJ(Atom1 = labels_atomic, Atom2 = labels_expanded)

  distances_dt <- data.table(
    Atom1 = atom_pairs$Atom1,
    Atom2 = atom_pairs$Atom2,
    Distance = as.vector(r),
    DeltaX = as.vector(delta_x), # Fractional differences
    DeltaY = as.vector(delta_y),
    DeltaZ = as.vector(delta_z),
    CosAlpha = cos_alpha, # These are constants for the cell
    CosBeta = cos_beta,
    CosGamma = cos_gamma
  )

  # Remove self-distances (Distance effectively zero)
  distances_dt <- distances_dt[Distance > tol]

  return(distances_dt)
}

#' Minimum Distance Bonding
#' @param distances A data.table.
#' @param delta A numeric tolerance.
#' @return A data.table.
#' @export
#' @examples
#' # See package vignette for examples.
minimum_distance <- function(distances, delta = 0.1) {
  # Add a simple title line for the description if Rd still says "no description"
  # #' @description Identifies bonded pairs using the Minimum Distance method.
  # The function body remains the same:
  if (is.null(distances) || nrow(distances) == 0) return(NULL)
  if (!all(c("Atom1", "Atom2", "Distance") %in% names(distances))) {
    warning("Distances table must contain Atom1, Atom2, Distance columns.")
    return(NULL)
  }

  dmin_table <- distances[, .(dmin = min(Distance)), by = Atom1]
  dmin_table[, dcut := (1 + delta) * dmin]

  setkey(dmin_table, Atom1)
  setkey(distances, Atom1)

  bonded_pairs <- dmin_table[distances, on = .(Atom1), nomatch = 0][Distance <= dcut,
                                                                    .(Atom1, Atom2, Distance, dcut, dmin)]

  if (nrow(bonded_pairs) == 0) return(NULL)
  return(bonded_pairs)
}

#' Apply Brunner's Largest Reciprocal Gap Method for Bonding
#'
#' Identifies bonded pairs using Brunner's method, which looks for the largest
#' gap in the reciprocals of sorted distances for each central atom.
#'
#' @param distances A data.table from \code{calculate_distances} with Atom1, Atom2, Distance.
#' @param delta A small tolerance added to the distance cut-off (default 0.0001).
#' @return A data.table of bonded pairs. Returns NULL if input is invalid or no bonds found.
#' @export
#' @examples
#' # dist_data <- data.table(Atom1=rep("A",5), Atom2=paste0("N",1:5),
#' #                         Distance=c(2.0, 2.1, 2.5, 3.0, 3.1))
#' # brunner(dist_data)
brunner <- function(distances, delta = 0.0001) {
  if (is.null(distances) || nrow(distances) == 0) return(NULL)
  if (!all(c("Atom1", "Atom2", "Distance") %in% names(distances))) {
    warning("Distances table must contain Atom1, Atom2, Distance columns.")
    return(NULL)
  }

  bonds_list <- list()
  unique_central_atoms <- unique(distances$Atom1)

  for (atom in unique_central_atoms) {
    atom_distances <- distances[Atom1 == atom][order(Distance)]
    if (nrow(atom_distances) == 0) next
    if (nrow(atom_distances) == 1) { # Only one neighbor, considered bonded
      bonds_list[[atom]] <- atom_distances
      next
    }

    largest_gap <- -Inf
    j_max <- 1 # Default to first atom if only one distance or no significant gap

    if (nrow(atom_distances) > 1) { # Need at least two distances to calculate a gap
      for (j in 1:(nrow(atom_distances) - 1)) {
        # Distances should be positive (filtered in calculate_distances)
        reciprocal_gap <- (1 / atom_distances$Distance[j]) - (1 / atom_distances$Distance[j + 1])
        if (reciprocal_gap > largest_gap) {
          largest_gap <- reciprocal_gap
          j_max <- j
        }
      }
    } else {
      j_max <- nrow(atom_distances) # Only one distance
    }

    d_cut <- atom_distances$Distance[j_max] + delta
    bonded_for_atom <- atom_distances[Distance <= d_cut]
    bonds_list[[atom]] <- bonded_for_atom
  }

  if (length(bonds_list) == 0) return(NULL)
  all_bonds <- rbindlist(bonds_list, fill = TRUE)
  if (nrow(all_bonds) == 0) return(NULL)
  return(all_bonds)
}

#' Apply Hoppe's Method of Effective Coordination Numbers for Bonding
#'
#' Identifies bonded pairs using Hoppe's ECON method. Bond strength is calculated
#' based on an iteratively determined average distance (d_avg).
#'
#' @param distances A data.table from \code{calculate_distances} with Atom1, Atom2, Distance.
#' @param bond_strength_threshold Delta, the minimum bond strength to consider a pair bonded (default 0.5).
#' @param tolerance Convergence tolerance for iterative d_avg calculation (default 0.001).
#' @return A data.table of bonded pairs (Atom1, Atom2, Distance, BondStrength).
#'         Returns NULL if input is invalid or no bonds found.
#' @export
#' @examples
#' # dist_data <- data.table(Atom1="A", Atom2=c("N1","N2","N3"), Distance=c(2.0, 2.2, 3.5))
#' # hoppe(dist_data)
hoppe <- function(distances, bond_strength_threshold = 0.5, tolerance = 0.001) {
  if (is.null(distances) || nrow(distances) == 0) return(NULL)
  if (!all(c("Atom1", "Atom2", "Distance") %in% names(distances))) {
    warning("Distances table must contain Atom1, Atom2, Distance columns.")
    return(NULL)
  }

  bonded_pairs_list <- list()
  unique_central_atoms <- unique(distances$Atom1)

  for (atom in unique_central_atoms) {
    atom_distances_dt <- distances[Atom1 == atom]
    if (nrow(atom_distances_dt) == 0) next

    dists_vec <- atom_distances_dt$Distance
    dmin_val <- min(dists_vec)
    if (dmin_val <= 0) { # Avoid division by zero or log of non-positive
      warning(paste("Atom", atom, "has non-positive minimum distance. Skipping for Hoppe's method."))
      next
    }

    weights_initial <- exp(1 - (dists_vec / dmin_val)^6)
    davg <- sum(dists_vec * weights_initial) / sum(weights_initial)
    if(is.na(davg) || !is.finite(davg)) {
      warning(paste("Initial d_avg calculation failed for atom", atom, ". Skipping."))
      next
    }

    iter_count <- 0
    max_iters <- 100 # Prevent potential infinite loop
    while (iter_count < max_iters) {
      prev_davg <- davg
      weights_iter <- exp(1 - (dists_vec / prev_davg)^6)
      davg <- sum(dists_vec * weights_iter) / sum(weights_iter)

      if(is.na(davg) || !is.finite(davg)) {
        warning(paste("Iterative d_avg calculation failed for atom", atom, "at iteration", iter_count, ". Using previous d_avg."))
        davg <- prev_davg # Revert if calculation fails
        break
      }
      if (abs(davg - prev_davg) <= tolerance) break
      iter_count <- iter_count + 1
    }
    if(iter_count == max_iters) {
      warning(paste("Hoppe's method d_avg did not converge for atom", atom, "within", max_iters, "iterations."))
    }

    atom_distances_dt[, BondStrength := exp(1 - (Distance / davg)^6)]
    bonded_for_atom <- atom_distances_dt[BondStrength >= bond_strength_threshold, .(Atom1, Atom2, Distance, BondStrength)]

    if (nrow(bonded_for_atom) > 0) {
      bonded_pairs_list[[atom]] <- bonded_for_atom
    }
  }

  if (length(bonded_pairs_list) == 0) return(NULL)
  all_bonded_pairs <- rbindlist(bonded_pairs_list, fill = TRUE)
  if (nrow(all_bonded_pairs) == 0) return(NULL)
  return(all_bonded_pairs)
}

#' Calculate Nearest Neighbor Counts
#'
#' Counts the number of bonded neighbors for each central atom.
#'
#' @param bonded_pairs_table A data.table of bonded pairs, typically output from
#'                           a bonding algorithm function (e.g., `minimum_distance`).
#'                           Must contain an 'Atom1' column for the central atom.
#' @return A data.table with Atom and NeighborCount.
#' @export
#' @examples
#' # bonds <- data.table(Atom1=c("A","A","B"), Atom2=c("N1","N2","N3"), Distance=c(1,1,1))
#' # calculate_neighbor_counts(bonds)
calculate_neighbor_counts <- function(bonded_pairs_table) {
  if (is.null(bonded_pairs_table) || nrow(bonded_pairs_table) == 0) {
    return(data.table(Atom = character(), NeighborCount = integer()))
  }
  if (!"Atom1" %in% names(bonded_pairs_table)) {
    warning("bonded_pairs_table must contain 'Atom1' column.")
    return(data.table(Atom = character(), NeighborCount = integer()))
  }

  neighbor_counts <- bonded_pairs_table[, .(NeighborCount = .N), by = .(Atom1)]
  setnames(neighbor_counts, "Atom1", "Atom")
  return(neighbor_counts)
}

#' Calculate Bond Angles
#'
#' Calculates bond angles (Neighbor1 - CentralAtom - Neighbor2) based on bonded pairs
#' and atomic coordinates, using the unit cell metric tensor.
#'
#' @param bonded_pairs A data.table of bonded pairs, with columns Atom1 (central) and Atom2 (neighbor).
#' @param atomic_coordinates A data.table of the original (asymmetric unit) atomic coordinates.
#'                           These are used as the central atoms for angle calculations.
#' @param expanded_coords A data.table of expanded coordinates (including those in neighboring cells).
#'                        This table is used to look up the positions of all atoms involved in angles.
#' @param unit_cell_metrics A data.table or list with unit cell parameters.
#' @return A data.table with CentralAtom, Neighbor1, Neighbor2, and Angle (in degrees).
#'         Returns NULL if inputs are invalid or no angles can be formed.
#' @export
#' @importFrom utils combn
#' @examples
#' # ucm <- data.table("_cell_length_a"=10, "_cell_length_b"=10, "_cell_length_c"=10,
#' #                   "_cell_angle_alpha"=90, "_cell_angle_beta"=90, "_cell_angle_gamma"=90)
#' # ac <- data.table(Label=c("C1", "H1", "H2"), x_a=c(0,0,0.1), y_b=c(0,0.1,0), z_c=c(0,0.05,0.05))
#' # ec <- ac # Simplified example, expanded_coords should contain all relevant atom positions
#' # bonds <- data.table(Atom1=c("C1","C1"), Atom2=c("H1","H2"))
#' # calculate_angles(bonds, ac, ec, ucm)
calculate_angles <- function(bonded_pairs, atomic_coordinates, expanded_coords, unit_cell_metrics) {
  if (is.null(bonded_pairs) || nrow(bonded_pairs) == 0 ||
      is.null(atomic_coordinates) || nrow(atomic_coordinates) == 0 ||
      is.null(expanded_coords) || nrow(expanded_coords) == 0 ||
      is.null(unit_cell_metrics) || nrow(unit_cell_metrics) == 0) {
    return(NULL)
  }

  a <- unit_cell_metrics[["_cell_length_a"]]
  b <- unit_cell_metrics[["_cell_length_b"]]
  c <- unit_cell_metrics[["_cell_length_c"]]
  alpha_rad <- unit_cell_metrics[["_cell_angle_alpha"]] * pi / 180
  beta_rad <- unit_cell_metrics[["_cell_angle_beta"]] * pi / 180
  gamma_rad <- unit_cell_metrics[["_cell_angle_gamma"]] * pi / 180

  if (any(is.na(c(a, b, c, alpha_rad, beta_rad, gamma_rad)))) {
    warning("Unit cell parameters for angle calculation are NA.")
    return(NULL)
  }

  cos_alpha <- cos(alpha_rad); cos_beta <- cos(beta_rad); cos_gamma <- cos(gamma_rad)
  tol <- 1e-10 # Tolerance for floating point comparisons
  if (abs(cos_alpha) < tol) cos_alpha <- 0
  if (abs(cos_beta) < tol) cos_beta <- 0
  if (abs(cos_gamma) < tol) cos_gamma <- 0

  cols_to_select <- c("Label", "x_a", "y_b", "z_c")
  if (!all(cols_to_select %in% names(atomic_coordinates)) ||
      !all(cols_to_select %in% names(expanded_coords))) {
    warning("Atomic coordinates and expanded_coords must have Label, x_a, y_b, z_c columns.")
    return(NULL)
  }

  all_coords_for_angles <- rbindlist(list(
    atomic_coordinates[, .SD, .SDcols = cols_to_select],
    expanded_coords[, .SD, .SDcols = cols_to_select]
  ), use.names = TRUE, fill = FALSE)

  all_coords_for_angles <- unique(all_coords_for_angles, by = "Label")
  setkey(all_coords_for_angles, Label) # Keying is still useful for efficiency if many lookups

  angle_list <- list()
  unique_central_atoms <- unique(bonded_pairs$Atom1)

  for (central_atom_label in unique_central_atoms) {
    neighbors_labels <- bonded_pairs[Atom1 == central_atom_label, Atom2]
    if (length(neighbors_labels) < 2) next

    coord_central <- all_coords_for_angles[Label == central_atom_label, .(x_a, y_b, z_c)]
    if(nrow(coord_central) == 0) {
      warning(paste("Central atom", central_atom_label, "not found in coordinate lookup table for angle calculation."))
      next
    }

    neighbor_combinations <- utils::combn(neighbors_labels, 2, simplify = FALSE)

    for (pair in neighbor_combinations) {
      neighbor1_label <- pair[1]
      neighbor2_label <- pair[2]

      coord_n1 <- all_coords_for_angles[Label == neighbor1_label, .(x_a, y_b, z_c)]
      coord_n2 <- all_coords_for_angles[Label == neighbor2_label, .(x_a, y_b, z_c)]

      if(nrow(coord_n1) == 0 || nrow(coord_n2) == 0) {
        warning(paste("One or both neighbors", neighbor1_label, ",", neighbor2_label, "for central atom", central_atom_label, "not found in coordinate lookup table."))
        next
      }

      v1_frac <- as.numeric(coord_n1 - coord_central)
      v2_frac <- as.numeric(coord_n2 - coord_central)

      xf1 <- v1_frac[1]; yf1 <- v1_frac[2]; zf1 <- v1_frac[3]
      xf2 <- v2_frac[1]; yf2 <- v2_frac[2]; zf2 <- v2_frac[3]

      dot_product <- (xf1*xf2*a^2 + yf1*yf2*b^2 + zf1*zf2*c^2 +
                        (xf1*yf2 + yf1*xf2)*a*b*cos_gamma +
                        (xf1*zf2 + zf1*xf2)*a*c*cos_beta +
                        (yf1*zf2 + zf1*yf2)*b*c*cos_alpha)

      mag_sq1 <- (xf1^2*a^2 + yf1^2*b^2 + zf1^2*c^2 +
                    2*xf1*yf1*a*b*cos_gamma + 2*xf1*zf1*a*c*cos_beta + 2*yf1*zf1*b*c*cos_alpha)
      mag_sq2 <- (xf2^2*a^2 + yf2^2*b^2 + zf2^2*c^2 +
                    2*xf2*yf2*a*b*cos_gamma + 2*xf2*zf2*a*c*cos_beta + 2*yf2*zf2*b*c*cos_alpha)

      if (mag_sq1 <= tol || mag_sq2 <= tol) {
        warning(paste("Zero magnitude vector encountered for angle:", central_atom_label, "-", neighbor1_label, "-", neighbor2_label))
        next
      }

      cos_theta_val <- dot_product / (sqrt(mag_sq1) * sqrt(mag_sq2))
      cos_theta_val <- min(max(cos_theta_val, -1.0), 1.0)

      angle_degrees <- acos(cos_theta_val) * 180 / pi

      if (!is.na(angle_degrees)) {
        angle_list[[length(angle_list) + 1]] <- data.table(
          CentralAtom = central_atom_label,
          Neighbor1 = neighbor1_label,
          Neighbor2 = neighbor2_label,
          Angle = angle_degrees
        )
      }
    }
  }

  if (length(angle_list) > 0) {
    angle_results <- rbindlist(angle_list)
    return(angle_results[order(CentralAtom, Neighbor1, Neighbor2)])
  } else {
    return(data.table(CentralAtom=character(), Neighbor1=character(), Neighbor2=character(), Angle=numeric()))
  }
}
