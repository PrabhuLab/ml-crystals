#' @title Identify Atomic Bonds using the Minimum Distance Method
#' @description Identifies bonded atoms by finding the nearest neighbor distance (d_min)
#'   for each central atom and defining a cutoff distance (d_cut) as
#'   d_cut = (1 + delta) * d_min.
#' @param distances A `data.table` of interatomic distances from `calculate_distances`.
#' @param delta The relative tolerance parameter (default 0.1).
#' @return A `data.table` of bonded pairs.
#' @family bonding algorithms
#' @export
minimum_distance <- function(distances, delta = 0.1) {
  if (nrow(distances) == 0)
    return(distances)

  dmin <- distances[, .(dmin = min(Distance)), by = .(Atom1)]
  dmin[, dcut := dmin * (1 + delta)]

  bonded_pairs <- distances[dmin, on = .(Atom1)]
  bonded_pairs <- bonded_pairs[Distance <= dcut, .(Atom1,
                                                   Atom2,
                                                   Distance,
                                                   DeltaX,
                                                   DeltaY,
                                                   DeltaZ,
                                                   Weight = dmin / Distance)]
  return(bonded_pairs)
}

#' @title Identify Atomic Bonds using Brunner's Method (Reciprocal)
#' @description An alternative bonding detection method using the
#'   largest reciprocal gap. Matches `pymatgen.analysis.local_env.BrunnerNNReciprocal`.
#' @param distances A `data.table` of interatomic distances.
#' @param tol A small distance offset (default 0.0) to add to the cutoff.
#' @return A `data.table` of bonded pairs.
#' @family bonding algorithms
#' @export
brunner_nn_reciprocal <- function(distances, tol = 0.0) {
  if (nrow(distances) == 0)
    return(distances)

  dt <- copy(distances)[order(Atom1, Distance)]

  # Calculate Reciprocal Gaps
  # Use fill=NA to prevent misalignment at group boundaries
  dt[, RecipGap := (1 / Distance) - (1 / data.table::shift(Distance, type = "lead", fill = NA)), by = Atom1]

  # Handle single neighbors where gap is NA
  dt[is.na(RecipGap), RecipGap := -Inf]

  dt[, IsMaxGap := (RecipGap == max(RecipGap, na.rm = TRUE)), by = Atom1]

  # Extract cutoff distance (distance of the atom BEFORE the gap)
  cutoff_map <- dt[IsMaxGap == TRUE, .(CutoffDist = Distance[1]), by = Atom1]

  result <- distances[cutoff_map, on = "Atom1"]
  result <- result[Distance <= CutoffDist + tol]

  result[, MinDist := min(Distance), by = Atom1]
  result[, Weight := MinDist / Distance]

  return(result[, .(Atom1, Atom2, Distance, DeltaX, DeltaY, DeltaZ, Weight)])
}

#' @title Identify Atomic Bonds using Hoppe's EconNN Method
#' @description Uses Effective Coordination Numbers (ECoN) and Mean Fictive Ionic Radii (MEFIR).
#' Matches pymatgen logic regarding MEFIR initialization and iteration.
#' @param distances A `data.table` of interatomic distances.
#' @param atomic_coordinates A `data.table` of atomic coordinates used to map species to radii.
#' @param tol A bond strength cutoff (default 0.2).
#' @param use_fictive_radius Logical. If `TRUE`, calculates Hoppe's fictive ionic radius based
#' on the ratio of central/neighbor radii. If `FALSE` (default), uses the raw bond distance,
#' matching pymatgen's default behavior.
#' @return A `data.table` of bonded pairs.
#' @family bonding algorithms
#' @export
econ_nn <- function(distances,
                    atomic_coordinates,
                    tol = 0.2,
                    use_fictive_radius = FALSE) {
  if (nrow(distances) == 0)
    return(distances)

  dt <- copy(distances)

  if (use_fictive_radius) {
    # Check if OxidationState exists, else create it as numeric NA
    atom_info <- copy(atomic_coordinates)
    if (!"OxidationState" %in% names(atom_info))
      atom_info[, OxidationState := NA_real_]

    # Helper to vectorise get_radius_pymatgen_style
    get_r_vec <- Vectorize(get_radius_decision_tree,
                           vectorize.args = c("symbol", "oxidation_state"))

    dt[, `:=`(
      Sym1 = sub("[0-9].*", "", sub("_.*", "", Atom1)),
      Sym2 = sub("[0-9].*", "", sub("_.*", "", Atom2)),
      Parent1 = sub("_.*", "", Atom1),
      Parent2 = sub("_.*", "", Atom2)
    )]

    # Merge Oxidation States
    dt[atom_info, on = c("Parent1" = "Label"), OS1 := i.OxidationState]
    dt[atom_info, on = c("Parent2" = "Label"), OS2 := i.OxidationState]

    # Get Radii using standard Pymatgen logic
    dt[, R1 := get_r_vec(Sym1, OS1, cn = 6)]
    dt[, R2 := get_r_vec(Sym2, OS2, cn = 6)]

    # Fictive Radii (Hoppe Eq 1)
    dt[, RadiusSum := R1 + R2]
    dt[, FictiveR := ifelse(RadiusSum > 0, Distance * (R1 / RadiusSum), 0)]
  } else {
    # If not using fictive radius, MEFIR is calculated on raw distances
    dt[, FictiveR := Distance]
  }

  # Iterative MEFIR Solver (matches pymatgen _get_mean_fictive_ionic_radius)
  calc_mefir_group <- function(fictive_r) {
    min_v <- min(fictive_r, na.rm = TRUE)
    if (min_v <= 0 || !is.finite(min_v))
      return(0)

    mean_fir <- min_v

    for (i in 1:100) {
      prev_mean_fir <- mean_fir
      # Hoppe Eq 2: w = exp(1 - (fir / min_fir)^6) where min_fir updates to mean
      ratio <- (fictive_r / prev_mean_fir)^6
      w <- exp(1 - ratio)
      sum_w <- sum(w, na.rm = TRUE)
      if (sum_w == 0)
        return(prev_mean_fir)
      mean_fir <- sum(fictive_r * w, na.rm = TRUE) / sum_w
      if (abs(mean_fir - prev_mean_fir) < 1e-4)
        break
    }
    return(mean_fir)
  }

  # Calculate MEFIR per Central Atom
  mefir_vals <- dt[, .(MEFIR = calc_mefir_group(FictiveR)), by = Atom1]
  dt[mefir_vals, on = "Atom1", MEFIR := i.MEFIR]

  # Calculate Final Weight
  dt[, Weight := ifelse(MEFIR > 0, exp(1 - (FictiveR / MEFIR)^6), 0)]

  return(dt[Weight > tol, .(Atom1, Atom2, Distance, DeltaX, DeltaY, DeltaZ, Weight)])
}

#' @title Identify Atomic Bonds using Voronoi Tessellation
#' @description Performs 3D Voronoi analysis on the supercell.
#' strictly mapping Asymmetric Unit atoms to their Supercell equivalents for processing.
#' Ensures potential neighbors are within cutoff to avoid edge artifacts.
#' Uses QJ (Joggled input) to handle degenerate structures (FCC/BCC) and filters/merges
#' the resulting micro-faces to produce accurate Solid Angles.
#' @param atomic_coordinates A `data.table` of the primary (asymmetric) atom set.
#' @param expanded_coords A `data.table` of atoms in the expanded supercell.
#' @param unit_cell_metrics A `data.table` with cell parameters.
#' @param cutoff Distance cutoff (default 13.0).
#' @param tol Tolerance for solid angle weights (default 0).
#' @return A `data.table` of bonded pairs.
#' @family bonding algorithms
#' @export
voronoi_nn <- function(atomic_coordinates,
                       expanded_coords,
                       unit_cell_metrics,
                       cutoff = 13.0,
                       tol = 0) {
  if (!requireNamespace("geometry", quietly = TRUE))
    stop("The 'geometry' package is required.")

  # 1. Convert Supercell to Cartesian Coords
  a <- as.numeric(unit_cell_metrics$`_cell_length_a`)
  b <- as.numeric(unit_cell_metrics$`_cell_length_b`)
  c <- as.numeric(unit_cell_metrics$`_cell_length_c`)
  alpha <- as.numeric(unit_cell_metrics$`_cell_angle_alpha`) * pi / 180
  beta <- as.numeric(unit_cell_metrics$`_cell_angle_beta`) * pi / 180
  gamma <- as.numeric(unit_cell_metrics$`_cell_angle_gamma`) * pi / 180
  v_sq <- 1 - cos(alpha)^2 - cos(beta)^2 - cos(gamma)^2 + 2 * cos(alpha) * cos(beta) * cos(gamma)
  v <- sqrt(max(0, v_sq))

  M <- matrix(
    c(
      a,
      b * cos(gamma),
      c * cos(beta),
      0,
      b * sin(gamma),
      c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma),
      0,
      0,
      c * v / sin(gamma)
    ),
    nrow = 3,
    byrow = TRUE
  )

  coords_mat <- as.matrix(expanded_coords[, .(x_a, y_b, z_c)]) %*% t(M)
  if (anyNA(coords_mat))
    return(NULL)

  # 2. Delaunay Tessellation ("QJ" is crucial)
  # QJ handles degeneracy (cospherical points in FCC/BCC) by perturbing input.
  # This creates micro-faces which must be merged/filtered later.
  tess <- tryCatch({
    geometry::delaunayn(coords_mat, options = "QJ Pp")
  }, error = function(e)
    NULL)
  if (is.null(tess))
    return(NULL)

  # 3. Calculate Circumcenters
  n_tets <- nrow(tess)
  circumcenters <- matrix(NA_real_, nrow = n_tets, ncol = 3)
  for (i in 1:n_tets) {
    circumcenters[i, ] <- calculate_circumcenter(coords_mat[tess[i, ], ])
  }

  atom_to_tets <- split(rep(1:n_tets, 4), as.vector(tess))

  # 4. Identify Asymmetric Unit Atoms
  asy_frac <- as.matrix(atomic_coordinates[, .(x_a, y_b, z_c)])
  exp_frac <- as.matrix(expanded_coords[, .(x_a, y_b, z_c)])
  target_indices <- integer(nrow(asy_frac))

  for (i in 1:nrow(asy_frac)) {
    diffs <- t(t(exp_frac) - asy_frac[i, ])
    dists_sq <- rowSums(diffs^2)
    match_idx <- which.min(dists_sq)
    # Match tolerance aligned with 1e-4 deduplication
    if (dists_sq[match_idx] < 1e-5)
      target_indices[i] <- match_idx
    else
      target_indices[i] <- NA
  }
  target_indices <- na.omit(target_indices)

  results_list <- list()

  # 5. Loop over Asymmetric Unit Atoms
  for (center_idx in target_indices) {
    center_coords <- coords_mat[center_idx, ]
    full_label <- expanded_coords$Label[center_idx]
    prim_label <- sub("_[0-9-]+_[0-9-]+_[0-9-]+$", "", full_label)

    my_tets <- atom_to_tets[[as.character(center_idx)]]
    if (is.null(my_tets))
      next

    potential_neighbors <- unique(as.vector(tess[my_tets, ]))
    potential_neighbors <- potential_neighbors[potential_neighbors != center_idx]
    if (length(potential_neighbors) == 0)
      next

    # CRITICAL FILTERING STEP:
    # Filter neighbors by cutoff (removes long connections to surface atoms)
    dists <- sqrt(rowSums((t(
      t(coords_mat[potential_neighbors, , drop = FALSE]) - center_coords
    ))^2))

    valid_mask <- dists <= cutoff
    neighbor_indices <- potential_neighbors[valid_mask]

    for (neigh_idx in neighbor_indices) {
      shared_tets <- intersect(my_tets, atom_to_tets[[as.character(neigh_idx)]])
      if (length(shared_tets) < 3)
        next

      raw_verts <- circumcenters[shared_tets, , drop = FALSE]
      raw_verts <- na.omit(raw_verts)
      if (nrow(raw_verts) < 3)
        next

      normal <- coords_mat[neigh_idx, ] - center_coords
      dist_sq <- sum(normal^2)
      if (is.na(dist_sq) || dist_sq < 1e-12)
        next

      # --- FIX 1: DEDUPLICATE VERTICES (Handle QJ Noise) ---
      unique_verts <- raw_verts[1, , drop = FALSE]
      if (nrow(raw_verts) > 1) {
        for (k in 2:nrow(raw_verts)) {
          pt <- raw_verts[k, ]
          dists <- sqrt(rowSums(t(t(
            unique_verts
          ) - pt)^2))
          if (all(dists > 1e-5)) {
            unique_verts <- rbind(unique_verts, pt)
          } else {
            # Average with existing to smooth QJ noise
            match_idx <- which.min(dists)
            unique_verts[match_idx, ] <- (unique_verts[match_idx, ] + pt) / 2
          }
        }
      }

      if (nrow(unique_verts) < 3)
        next

      # --- FIX 2: ORDER VERTICES ---
      # Define plane basis vectors u, v
      face_centroid <- colMeans(unique_verts)
      u <- if (abs(normal[1]) < 0.9)
        c(1, 0, 0)
      else
        c(0, 1, 0)
      u <- cross_product(normal, u)
      u <- u / sqrt(sum(u^2))
      v <- cross_product(normal, u)
      v <- v / sqrt(sum(v^2))

      # Calculate angles for ordering
      rel <- t(t(unique_verts) - face_centroid)
      angles <- atan2(rel %*% v, rel %*% u)

      # Order vertices BEFORE calculating Area OR Solid Angle
      face_verts_ordered <- unique_verts[order(angles), , drop = FALSE]

      # --- FIX 3: CALCULATE METRICS ON ORDERED POLYGON ---
      sa <- calculate_solid_angle(center_coords, face_verts_ordered)

      area <- 0
      n_fv <- nrow(face_verts_ordered)

      for (k in 1:n_fv) {
        p1 <- face_centroid
        p2 <- face_verts_ordered[k, ]
        p3 <- face_verts_ordered[if (k == n_fv)
          1
          else
            k + 1, ]
        cp <- cross_product(p2 - p1, p3 - p1)
        area <- area + 0.5 * sqrt(sum(cp^2))
      }

      # Filter micro-faces created by QJ option (Joggle)
      if (area < 1e-5 || sa < 1e-5)
        next

      results_list[[length(results_list) + 1]] <- list(
        Atom1 = prim_label,
        Atom2 = expanded_coords$Label[neigh_idx],
        Distance = sqrt(dist_sq),
        SolidAngle = sa,
        Area = area,
        DeltaX = normal[1],
        DeltaY = normal[2],
        DeltaZ = normal[3]
      )
    }
  }

  if (length(results_list) == 0)
    return(NULL)
  dt <- rbindlist(results_list)

  dt[, MaxSA := max(SolidAngle, na.rm = TRUE), by = Atom1]
  dt[, Weight := ifelse(MaxSA > 0, SolidAngle / MaxSA, 0)]

  return(dt[Weight > tol])
}

#' @title Identify Atomic Bonds using CrystalNN
#' @description Poer of Pymatgen's `CrystalNN` algorithm from the ground up.
#' It uses a Voronoi-based algorithm and solid angle weights to determine the
#' probability of various coordination environments. It modifies probability using
#' smooth distance cutoffs and Pauling electronegativity differences. The output
#' is either the most probable coordination environment (`weighted_cn = FALSE`)
#' or a weighted list of coordination environments (`weighted_cn = TRUE`).
#' @param distances Ignored (CrystalNN generates its own Voronoi basis).
#' @param atomic_coordinates Primary atom set.
#' @param expanded_coords Expanded supercell.
#' @param unit_cell_metrics Cell parameters.
#' @param cutoff_length Numeric. Cutoff in Angstroms for initial neighbor search
#' (default 7.0).
#' @param x_diff_weight Numeric. Electronegativity difference weight (default 3.0).
#' @param porous_adjustment Logical. If TRUE, adjusts Voronoi weights to
#' better describe layered / porous structures (default TRUE).
#' @param distance_cutoffs Numeric vector. Penalizes neighbor distances greater than
#' sum of covalent radii plus these cutoffs. Set to NULL to disable.
#' @param cation_anion Logical. If TRUE, restricts bonding targets to sites
#' with opposite or zero charge (requires oxidation states).
#' @param weighted_cn Logical. If FALSE (default), returns neighbors for the
#' most probable coordination environment with weight 1.0. If TRUE, returns all
#' potential neighbors with fractional probabilities.
#' @return A `data.table` of bonded pairs.
#' @family bonding algorithms
#' @export
crystal_nn <- function(distances,
                       atomic_coordinates,
                       expanded_coords,
                       unit_cell_metrics,
                       cutoff_length = 7.0,
                       x_diff_weight = 3.0,
                       porous_adjustment = TRUE,
                       distance_cutoffs = c(0.5, 1.0),
                       cation_anion = FALSE,
                       weighted_cn = FALSE) {

  # 1. Get Base Voronoi Neighbors
  base <- voronoi_nn(
    atomic_coordinates,
    expanded_coords,
    unit_cell_metrics,
    cutoff = cutoff_length,
    tol = 0
  )
  if (is.null(base) || nrow(base) == 0)
    return(NULL)

  dt <- copy(base)
  dt[, `:=`(Sym1 = sub("[0-9].*", "", sub("_.*", "", Atom1)),
            Sym2 = sub("[0-9].*", "", sub("_.*", "", Atom2)))]

  # Merge Oxidation States if available (for radii/EN logic)
  atom_info <- copy(atomic_coordinates)
  if (!"OxidationState" %in% names(atom_info))
    atom_info[, OxidationState := NA_real_]

  dt[, Parent1 := sub("_.*", "", Atom1)]
  dt[, Parent2 := sub("_.*", "", Atom2)]
  dt[atom_info, on = c("Parent1" = "Label"), OS1 := i.OxidationState]
  dt[atom_info, on = c("Parent2" = "Label"), OS2 := i.OxidationState]

  # --- Optional: Cation/Anion constraint ---
  if (cation_anion) {
    if (all(is.na(dt$OS1)) || all(is.na(dt$OS2))) {
      warning("CrystalNN: cation_anion=TRUE but oxidation states are missing. Ignoring constraint.")
    } else {
      # Keep if opposite charge, if one is neutral (0), or if NA to fail gracefully
      dt <- dt[is.na(OS1) | is.na(OS2) | (OS1 * OS2 <= 0)]
      if (nrow(dt) == 0) return(NULL)
    }
  }

  # --- Step 1: Porous Adjustment ---
  if (porous_adjustment) {
    dt[, RawScore := ifelse(Area > 1e-8, (SolidAngle * SolidAngle) / Area, SolidAngle)]
  } else {
    dt[, RawScore := SolidAngle]
  }

  # --- Step 2: Chemical Adjustment (Electronegativity) ---
  if (!is.null(x_diff_weight) && x_diff_weight > 0) {
    dt[, EN1 := get_electronegativity(Sym1)]
    dt[, EN2 := get_electronegativity(Sym2)]

    dt[, ChemMod := 1 + x_diff_weight * sqrt(abs(EN1 - EN2) / 3.3)]
    dt[is.na(ChemMod), ChemMod := 1.0]

    dt[, RawScore := RawScore * ChemMod]
  }

  # --- Step 3: Normalize (Highest weight = 1.0) ---
  dt[, MaxScore := max(RawScore, na.rm = TRUE), by = Atom1]
  dt[, Weight := ifelse(MaxScore > 0, RawScore / MaxScore, 0)]

  # --- Step 4: Distance Cutoffs ---
  if (!is.null(distance_cutoffs)) {
    atomic_rads <- get_radii_data()

    get_r_vec <- Vectorize(get_radius_decision_tree,
                           vectorize.args = c("symbol", "oxidation_state"))

    dt[, R1 := get_r_vec(Sym1, OS1, cn = 6)]
    dt[, R2 := get_r_vec(Sym2, OS2, cn = 6)]

    get_default_r <- Vectorize(function(sym) {
      r <- atomic_rads[Symbol == sym & Type == "covalent"]$Radius
      if (length(r) == 0)
        r <- atomic_rads[Symbol == sym & Type == "atomic"]$Radius
      if (length(r) == 0)
        return(0)
      return(r[1])
    }, vectorize.args = "sym")

    dt[R1 == 0, R1 := get_default_r(Sym1)]
    dt[R2 == 0, R2 := get_default_r(Sym2)]

    dt[, Diameter := R1 + R2]
    dt[, CutLow := Diameter + distance_cutoffs[1]]
    dt[, CutHigh := Diameter + distance_cutoffs[2]]

    dt[, DistMod := data.table::fcase(
      Distance <= CutLow, 1.0,
      Distance >= CutHigh, 0.0,
      default = 0.5 * (cos((Distance - CutLow) / (CutHigh - CutLow) * pi) + 1)
    )]
    dt[is.na(DistMod), DistMod := 0]

    dt[, Weight := Weight * DistMod]
  }

  # --- Step 5: Semicircle Integral and CN Probability ---
  dt[, Weight := round(Weight, 3)]
  dt <- dt[Weight > 0]
  if (nrow(dt) == 0) return(NULL)

  # Internal function to get an integral between two bounds of a unit semicircle.
  # Used in algorithm to determine bond probabilities.
  semicircle_integral <- function(x1, x2) {
    radius <- 1.0
    calc_area <- function(x) {
      res <- numeric(length(x))
      res[x >= 1.0] <- 0.25 * pi * radius^2
      res[x <= 0.0] <- 0.0
      mid <- x > 0.0 & x < 1.0
      x_m <- x[mid]
      res[mid] <- 0.5 * (x_m * sqrt(radius^2 - x_m^2) + radius^2 * atan(x_m / sqrt(radius^2 - x_m^2)))
      return(res)
    }
    return((calc_area(x1) - calc_area(x2)) / (0.25 * pi * radius^2))
  }

  result_list <- lapply(split(dt, dt$Atom1), function(sub_dt) {
    sub_dt <- sub_dt[order(-Weight)]

    if (weighted_cn) {
      # Mathematical simplification of pymatgen algorithm: The fractional weight of
      # a neighbor is the sum of the probabilities of all coordination environments
      # it participates in. Because the discrete CN probability bins partition the
      # interval [0, 1], this exactly simplifies to the integral of the unit
      # semicircle starting from the neighbor's weight down to 0.0.
      sub_dt[, FinalWeight := semicircle_integral(Weight, 0.0)]
      return(sub_dt)

    } else {
      # Map distinct weights into discrete Coordination Environments (CN probabilities)
      w_vals <- sub_dt$Weight
      dist_bins <- unique(w_vals)
      dist_bins <- c(dist_bins, 0.0)

      cn_weights <- numeric(nrow(sub_dt) + 1) # index 1 = CN0, index 2 = CN1, etc.

      for (i in seq_len(length(dist_bins) - 1)) {
        x1 <- dist_bins[i]
        x2 <- dist_bins[i + 1]
        prob <- semicircle_integral(x1, x2)
        cn <- sum(w_vals >= x1)
        cn_weights[cn + 1] <- cn_weights[cn + 1] + prob
      }

      cn_weights[1] <- max(0, 1.0 - sum(cn_weights)) # CN0 (No coordinated neighbors)
      best_cn <- which.max(cn_weights) - 1

      if (best_cn == 0) {
        return(NULL)
      } else {
        # Standard approach: Fix CN rigidly at the most probable value and assign a weight of 1
        res <- head(sub_dt, best_cn)
        res[, FinalWeight := 1.0]
        return(res)
      }
    }
  })

  final_dt <- rbindlist(result_list)
  if (nrow(final_dt) == 0) return(NULL)

  final_dt[, Weight := FinalWeight]
  return(final_dt[, .(Atom1, Atom2, Distance, DeltaX, DeltaY, DeltaZ, Weight)])
}
