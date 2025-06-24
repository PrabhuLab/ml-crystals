#' Calculate Geometric Properties and Bonding
#'
#' A suite of functions to calculate interatomic distances, identify bonds
#' using various algorithms, and compute derivative properties like neighbor counts
#' and bond angles.
#' @name calculate_properties
NULL

#' @rdname calculate_properties
#' @param atomic_coordinates A data.table of reference atomic coordinates (Label, x_a, y_b, z_c).
#' @param expanded_coords A data.table of search atomic coordinates (Label, x_a, y_b, z_c).
#' @param unit_cell_metrics A data.table with unit cell parameters.
#' @return A data.table with Atom1, Atom2, Distance, and deltas.
#' @export
calculate_distances <- function(atomic_coordinates, expanded_coords, unit_cell_metrics) {
  if (is.null(atomic_coordinates) || nrow(atomic_coordinates) == 0 ||
      is.null(expanded_coords) || nrow(expanded_coords) == 0 ||
      is.null(unit_cell_metrics) || nrow(unit_cell_metrics) == 0) {
    return(NULL)
  }
  a <- unit_cell_metrics[["_cell_length_a"]]; b <- unit_cell_metrics[["_cell_length_b"]]; c <- unit_cell_metrics[["_cell_length_c"]]
  alpha <- unit_cell_metrics[["_cell_angle_alpha"]] * pi / 180
  beta <- unit_cell_metrics[["_cell_angle_beta"]] * pi / 180
  gamma <- unit_cell_metrics[["_cell_angle_gamma"]] * pi / 180
  if (any(is.na(c(a, b, c, alpha, beta, gamma)))) return(NULL)
  coords_atomic <- as.matrix(atomic_coordinates[, .(x_a, y_b, z_c)]); coords_expanded <- as.matrix(expanded_coords[, .(x_a, y_b, z_c)])
  labels_atomic <- atomic_coordinates$Label; labels_expanded <- expanded_coords$Label
  delta_x <- outer(coords_atomic[, 1], coords_expanded[, 1], "-"); delta_y <- outer(coords_atomic[, 2], coords_expanded[, 2], "-"); delta_z <- outer(coords_atomic[, 3], coords_expanded[, 3], "-")
  cos_alpha <- cos(alpha); cos_beta <- cos(beta); cos_gamma <- cos(gamma)
  tol <- 1e-10
  if (abs(cos_alpha) < tol) cos_alpha <- 0; if (abs(cos_beta) < tol) cos_beta <- 0; if (abs(cos_gamma) < tol) cos_gamma <- 0
  r2 <- (a^2*delta_x^2) + (b^2*delta_y^2) + (c^2*delta_z^2) + (2*b*c*cos_alpha*delta_y*delta_z) + (2*c*a*cos_beta*delta_z*delta_x) + (2*a*b*cos_gamma*delta_x*delta_y)
  r2[r2 < 0 & abs(r2) < tol] <- 0
  r <- sqrt(r2)
  atom_pairs <- CJ(Atom1 = labels_atomic, Atom2 = labels_expanded)
  distances_dt <- data.table(Atom1 = atom_pairs$Atom1, Atom2 = atom_pairs$Atom2, Distance = as.vector(r),
                             DeltaX = as.vector(delta_x), DeltaY = as.vector(delta_y), DeltaZ = as.vector(delta_z),
                             CosAlpha = cos_alpha, CosBeta = cos_beta, CosGamma = cos_gamma)
  return(distances_dt[Distance > 0])
}

#' Identify Bonded Pairs
#'
#' Functions to identify bonded pairs from a table of interatomic distances using
#' different algorithms.
#' @param distances A data.table from \code{calculate_distances}.
#' @param delta A numeric tolerance or threshold parameter specific to each method.
#' @return A data.table of bonded pairs.
#' @name bonding
NULL

#' @rdname bonding
#' @export
minimum_distance <- function(distances, delta = 0.1) {
  # Input validation: Return an empty data.table with a matching schema if input is empty.
  # This is more robust than returning NULL.
  if (is.null(distances) || nrow(distances) == 0) {
    return(
      data.table::data.table(
        Atom1 = character(),
        Atom2 = character(),
        Distance = numeric(),
        dcut = numeric(),
        dmin = numeric()
      )
    )
  }

  # Step 1: Calculate dmin for each central atom (Atom1).
  # This creates a new, small summary table.
  dmin_table <- distances[, .(dmin = min(Distance)), by = Atom1]

  # Step 2: Calculate the cutoff distance (dcut) for each Atom1 in the summary table.
  dmin_table[, dcut := (1 + delta) * dmin]

  # Step 3: Perform a join and filter to find bonded pairs.
  # This is the idiomatic `data.table` way that avoids modifying the original `distances` table.
  # - It joins the dcut/dmin info from dmin_table to the distances table.
  # - `allow.cartesian=TRUE` is needed because each Atom1 in dmin_table matches multiple rows in distances.
  # - It then immediately filters rows where the distance is within the calculated cutoff.
  bonded_pairs <- distances[dmin_table, on = .(Atom1), allow.cartesian = TRUE][
    Distance <= dcut,
    .(Atom1, Atom2, Distance, dcut, dmin)
  ]

  return(bonded_pairs)
}

#' @rdname bonding
#' @export
brunner <- function(distances, delta = 0.0001) {
  if (is.null(distances) || nrow(distances) == 0) return(NULL)
  bonds_list <- list()
  for (atom in unique(distances$Atom1)) {
    atom_distances <- distances[Atom1 == atom][order(Distance)]
    if (nrow(atom_distances) <= 1) {
      if(nrow(atom_distances) == 1) bonds_list[[atom]] <- atom_distances
      next
    }
    largest_gap <- -Inf; j_max <- 1
    if(nrow(atom_distances) > 1){
      for (j in 1:(nrow(atom_distances) - 1)) {
        reciprocal_gap <- (1 / atom_distances$Distance[j]) - (1 / atom_distances$Distance[j + 1])
        if (reciprocal_gap > largest_gap) { largest_gap <- reciprocal_gap; j_max <- j }
      }
    }
    d_cut <- atom_distances$Distance[j_max] + delta
    bonds_list[[atom]] <- atom_distances[Distance <= d_cut]
  }
  if (length(bonds_list) == 0) return(NULL)
  return(rbindlist(bonds_list, fill = TRUE))
}

#' @rdname bonding
#' @param bond_strength_threshold For `hoppe`, the minimum bond strength to consider a pair bonded.
#' @param tolerance For `hoppe`, the convergence tolerance for iterative d_avg calculation.
#' @export
hoppe <- function(distances, bond_strength_threshold = 0.5, tolerance = 0.001) {
  if (is.null(distances) || nrow(distances) == 0) return(NULL)
  bonded_pairs_list <- list()
  for (atom in unique(distances$Atom1)) {
    atom_distances_dt <- distances[Atom1 == atom]
    if (nrow(atom_distances_dt) == 0) next
    dists_vec <- atom_distances_dt$Distance; dmin_val <- min(dists_vec)
    if (dmin_val <= 1e-6) next
    weights_initial <- exp(1 - (dists_vec / dmin_val)^6)
    davg <- sum(dists_vec * weights_initial) / sum(weights_initial)
    if(is.na(davg) || !is.finite(davg)) next
    iter_count <- 0
    while (iter_count < 100) {
      prev_davg <- davg
      weights_iter <- exp(1 - (dists_vec / prev_davg)^6)
      davg <- sum(dists_vec * weights_iter) / sum(weights_iter)
      if(is.na(davg) || !is.finite(davg) || abs(davg - prev_davg) <= tolerance) break
      iter_count <- iter_count + 1
    }
    atom_distances_dt[, BondStrength := exp(1 - (Distance / davg)^6)]
    bonded_for_atom <- atom_distances_dt[BondStrength >= bond_strength_threshold, .(Atom1, Atom2, Distance, BondStrength)]
    if (nrow(bonded_for_atom) > 0) bonded_pairs_list[[atom]] <- bonded_for_atom
  }
  if (length(bonded_pairs_list) == 0) return(NULL)
  return(rbindlist(bonded_pairs_list, fill = TRUE))
}

#' @rdname calculate_properties
#' @param bonded_pairs_table A data.table of bonded pairs.
#' @return A data.table with Atom and NeighborCount.
#' @export
calculate_neighbor_counts <- function(bonded_pairs_table) {
  if (is.null(bonded_pairs_table) || nrow(bonded_pairs_table) == 0) {
    return(data.table(Atom = character(), NeighborCount = integer()))
  }
  neighbor_counts <- bonded_pairs_table[, .(NeighborCount = .N), by = .(Atom1)]
  setnames(neighbor_counts, "Atom1", "Atom")
  return(neighbor_counts)
}

#' @rdname calculate_properties
#' @param bonded_pairs A data.table of bonded pairs.
#' @return A data.table with CentralAtom, Neighbor1, Neighbor2, and Angle.
#' @export
#' @importFrom utils combn
calculate_angles <- function(bonded_pairs, atomic_coordinates, expanded_coords, unit_cell_metrics) {
  if (is.null(bonded_pairs) || nrow(bonded_pairs) == 0) return(NULL)
  if (is.null(atomic_coordinates) || is.null(expanded_coords) || is.null(unit_cell_metrics)) return(NULL)
  a <- unit_cell_metrics[["_cell_length_a"]]; b <- unit_cell_metrics[["_cell_length_b"]]; c <- unit_cell_metrics[["_cell_length_c"]]
  alpha_rad <- unit_cell_metrics[["_cell_angle_alpha"]] * pi / 180; beta_rad <- unit_cell_metrics[["_cell_angle_beta"]] * pi / 180; gamma_rad <- unit_cell_metrics[["_cell_angle_gamma"]] * pi / 180
  cos_alpha <- cos(alpha_rad); cos_beta <- cos(beta_rad); cos_gamma <- cos(gamma_rad)
  all_coords <- rbind(atomic_coordinates[, .(Label, x_a, y_b, z_c)], expanded_coords[, .(Label, x_a, y_b, z_c)], fill = TRUE)
  all_coords <- unique(all_coords, by = "Label"); setkey(all_coords, Label)
  calculate_angle_metric <- function(central_atom, atom2, atom3) {
    coord1 <- all_coords[central_atom, .(x_a, y_b, z_c)]; coord2 <- all_coords[atom2, .(x_a, y_b, z_c)]; coord3 <- all_coords[atom3, .(x_a, y_b, z_c)]
    if (anyNA(coord1) || anyNA(coord2) || anyNA(coord3)) return(NA_real_)
    v1_frac <- as.numeric(coord2 - coord1); v2_frac <- as.numeric(coord3 - coord1)
    xf1 <- v1_frac[1]; yf1 <- v1_frac[2]; zf1 <- v1_frac[3]; xf2 <- v2_frac[1]; yf2 <- v2_frac[2]; zf2 <- v2_frac[3]
    dot_product <- (xf1*xf2*a^2 + yf1*yf2*b^2 + zf1*zf2*c^2 + (xf1*yf2 + yf1*xf2)*a*b*cos_gamma + (xf1*zf2 + zf1*xf2)*a*c*cos_beta + (yf1*zf2 + zf1*yf2)*b*c*cos_alpha)
    mag_sq1 <- (xf1^2*a^2 + yf1^2*b^2 + zf1^2*c^2 + 2*xf1*yf1*a*b*cos_gamma + 2*xf1*zf1*a*c*cos_beta + 2*yf1*zf1*b*c*cos_alpha)
    mag_sq2 <- (xf2^2*a^2 + yf2^2*b^2 + zf2^2*c^2 + 2*xf2*yf2*a*b*cos_gamma + 2*xf2*zf2*a*c*cos_beta + 2*yf2*zf2*b*c*cos_alpha)
    if (mag_sq1 <= 1e-10 || mag_sq2 <= 1e-10) return(NA_real_)
    cos_theta <- min(max(dot_product / (sqrt(mag_sq1) * sqrt(mag_sq2)), -1.0), 1.0)
    return(acos(cos_theta) * 180 / pi)
  }
  angle_list <- list()
  for (central_atom_label in unique(bonded_pairs$Atom1)) {
    neighbors_labels <- bonded_pairs[Atom1 == central_atom_label, Atom2]
    if (length(neighbors_labels) >= 2) {
      for (pair in combn(neighbors_labels, 2, simplify = FALSE)) {
        angle <- calculate_angle_metric(central_atom_label, pair[1], pair[2])
        if (!is.na(angle)) {
          angle_list[[length(angle_list) + 1]] <- data.table(CentralAtom = central_atom_label, Neighbor1 = pair[1], Neighbor2 = pair[2], Angle = angle)
        }
      }
    }
  }
  if (length(angle_list) > 0) return(rbindlist(angle_list)[order(CentralAtom, Neighbor1, Neighbor2)])
  return(data.table(CentralAtom=character(), Neighbor1=character(), Neighbor2=character(), Angle=numeric()))
}
