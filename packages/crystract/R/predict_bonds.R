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
  dt[, RecipGap := (1 / Distance) - (1 / data.table::shift(Distance, type = "lead", fill = NA)), by = Atom1]
  dt[is.na(RecipGap), RecipGap := -Inf]

  dt[, IsMaxGap := (RecipGap == max(RecipGap, na.rm = TRUE)), by = Atom1]
  cutoff_map <- dt[IsMaxGap == TRUE, .(CutoffDist = Distance[1]), by = Atom1]

  result <- distances[cutoff_map, on = "Atom1"]
  result <- result[Distance <= CutoffDist + tol]

  result[, MinDist := min(Distance), by = Atom1]
  result[, Weight := MinDist / Distance]

  return(result[, .(Atom1, Atom2, Distance, DeltaX, DeltaY, DeltaZ, Weight)])
}

#' @title Identify Atomic Bonds using Hoppe's EconNN Method
#' @description Uses Effective Coordination Numbers (ECoN) and Mean Fictive Ionic Radii (MEFIR).
#' @param distances A `data.table` of interatomic distances.
#' @param atomic_coordinates A `data.table` of atomic coordinates used to map species to radii.
#' @param tol A bond strength cutoff (default 0.2).
#' @param cation_anion Logical. If TRUE, restricts bonding to sites with opposite or zero charge.
#' @return A `data.table` of bonded pairs.
#' @family bonding algorithms
#' @export
econ_nn <- function(distances,
                    atomic_coordinates,
                    tol = 0.2,
                    cation_anion = FALSE) {
  if (nrow(distances) == 0)
    return(distances)

  atom_info <- copy(atomic_coordinates)
  if (!"OxidationState" %in% names(atom_info))
    atom_info[, OxidationState := NA_real_]

  # Vectorized radius lookup helper
  get_r_vec <- Vectorize(get_radius, vectorize.args = c("symbol", "oxidation_state"))

  dt <- copy(distances)
  dt[, `:=`(
    Sym1 = sub("[0-9].*", "", sub("_.*", "", Atom1)),
    Sym2 = sub("[0-9].*", "", sub("_.*", "", Atom2)),
    Parent1 = sub("_.*", "", Atom1),
    Parent2 = sub("_.*", "", Atom2)
  )]

  # Merge Oxidation States
  dt[atom_info, on = c("Parent1" = "Label"), OS1 := i.OxidationState]
  dt[atom_info, on = c("Parent2" = "Label"), OS2 := i.OxidationState]

  # 1. Cation-Anion Filter (Pymatgen Logic)
  # If enabled, filter out like-charged neighbors (excluding neutrals)
  if (cation_anion) {
    # Logic: Keep if any OS is NA, OR signs are opposite/zero
    dt <- dt[is.na(OS1) |
               is.na(OS2) |
               (OS1 >= 0 & OS2 <= 0) | (OS1 < 0 & OS2 >= 0)]
  }

  if (nrow(dt) == 0)
    return(NULL)

  # 2. Get Radii (Using new interpolated get_radius logic)
  dt[, R1 := get_r_vec(Sym1, OS1)]
  dt[, R2 := get_r_vec(Sym2, OS2)]

  # 3. Fictive Radii (Hoppe Eq 1)
  dt[, RadiusSum := R1 + R2]
  # Fallback to bond distance if radii sum is invalid
  dt[, FictiveR := ifelse(RadiusSum > 0, Distance * (R1 / RadiusSum), Distance)]

  # 4. Iterative MEFIR Solver
  # Updates the "minimum_fir" reference to the current mean in each loop
  calc_mefir_group <- function(fictive_r) {
    min_v <- min(fictive_r, na.rm = TRUE)
    if (min_v <= 0 || !is.finite(min_v))
      return(0)

    mean_fir <- min_v

    for (i in 1:100) {
      prev_mean_fir <- mean_fir
      # Weighting based on previous MEFIR (pymatgen logic: minimum_fir = mean_fictive_ionic_radius)
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

  mefir_vals <- dt[, .(MEFIR = calc_mefir_group(FictiveR)), by = Atom1]
  dt[mefir_vals, on = "Atom1", MEFIR := i.MEFIR]

  # 5. Final Weight
  dt[, Weight := ifelse(MEFIR > 0, exp(1 - (FictiveR / MEFIR)^6), 0)]

  return(dt[Weight > tol, .(Atom1, Atom2, Distance, DeltaX, DeltaY, DeltaZ, Weight)])
}

#' @title Identify Atomic Bonds using Voronoi Tessellation
#' @description Performs 3D Voronoi analysis on the supercell.
#' strictly mapping Asymmetric Unit atoms to their Supercell equivalents.
#' @param atomic_coordinates A `data.table` of the primary (asymmetric) atom set.
#' @param expanded_coords A `data.table` of atoms in the expanded supercell.
#' @param unit_cell_metrics A `data.table` with cell parameters.
#' @param cutoff Distance cutoff (default 13.0).
#' @param tol Tolerance for solid angle weights (default 0).
#' @param allow_pathological Logical. If FALSE, raises error for infinite vertices.
#' @return A `data.table` of bonded pairs.
#' @family bonding algorithms
#' @export
voronoi_nn <- function(atomic_coordinates,
                       expanded_coords,
                       unit_cell_metrics,
                       cutoff = 13.0,
                       tol = 0,
                       allow_pathological = FALSE) {
  if (!requireNamespace("geometry", quietly = TRUE))
    stop("The 'geometry' package is required.")

  # 1. Convert Supercell to Cartesian Coords
  a <- as.numeric(unit_cell_metrics$`_cell_length_a`)
  b <- as.numeric(unit_cell_metrics$`_cell_length_b`)
  c <- as.numeric(unit_cell_metrics$`_cell_length_c`)
  alpha_r <- as.numeric(unit_cell_metrics$`_cell_angle_alpha`) * pi / 180
  beta_r <- as.numeric(unit_cell_metrics$`_cell_angle_beta`) * pi / 180
  gamma_r <- as.numeric(unit_cell_metrics$`_cell_angle_gamma`) * pi / 180

  v_sq <- 1 - cos(alpha_r)^2 - cos(beta_r)^2 - cos(gamma_r)^2 + 2 * cos(alpha_r) * cos(beta_r) * cos(gamma_r)
  v <- sqrt(max(0, v_sq))

  M <- matrix(
    c(
      a,
      b * cos(gamma_r),
      c * cos(beta_r),
      0,
      b * sin(gamma_r),
      c * (cos(alpha_r) - cos(beta_r) * cos(gamma_r)) / sin(gamma_r),
      0,
      0,
      c * v / sin(gamma_r)
    ),
    nrow = 3,
    byrow = TRUE
  )

  coords_mat <- as.matrix(expanded_coords[, .(x_a, y_b, z_c)]) %*% t(M)

  if (anyNA(coords_mat))
    return(NULL)

  # 2. Delaunay Tessellation
  tess <- tryCatch({
    geometry::delaunayn(coords_mat, options = "QJ Pp")
  }, error = function(e)
    NULL)

  if (is.null(tess))
    return(NULL)

  # 3. Calculate Circumcenters (Voronoi Vertices)
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
    dists_sq <- rowSums(t(t(exp_frac) - asy_frac[i, ])^2)
    target_indices[i] <- which.min(dists_sq)
  }
  target_indices <- na.omit(unique(target_indices))

  # Pathological Check: If atoms are on the convex hull, cells are infinite
  if (!allow_pathological) {
    hull_inds <- tryCatch(
      unique(as.vector(geometry::convhulln(coords_mat))),
      error = function(e)
        NULL
    )
    if (!is.null(hull_inds) && any(target_indices %in% hull_inds)) {
      stop(
        "Structure is pathological: Voronoi cell contains infinite vertices. Increase supercell size or set allow_pathological=TRUE."
      )
    }
  }

  results_list <- list()

  # 5. Loop ONLY over Asymmetric Unit Atoms
  for (center_idx in target_indices) {
    center_coords <- coords_mat[center_idx, ]
    prim_label <- sub("_[0-9-]+_[0-9-]+_[0-9-]+$", "", expanded_coords$Label[center_idx])

    my_tets <- atom_to_tets[[as.character(center_idx)]]
    if (is.null(my_tets))
      next

    potential_neighbors <- unique(as.vector(tess[my_tets, ]))
    potential_neighbors <- potential_neighbors[potential_neighbors != center_idx]

    # Pre-filter by cutoff
    dists <- sqrt(rowSums((t(
      t(coords_mat[potential_neighbors, , drop = FALSE]) - center_coords
    ))^2))
    valid_mask <- dists <= cutoff
    neighbor_indices <- potential_neighbors[valid_mask]

    for (neigh_idx in neighbor_indices) {
      shared_tets <- intersect(my_tets, atom_to_tets[[as.character(neigh_idx)]])
      if (length(shared_tets) < 3)
        next

      # Face vertices = circumcenters of shared tetrahedra
      face_verts <- circumcenters[shared_tets, , drop = FALSE]
      face_verts <- na.omit(face_verts)
      if (nrow(face_verts) < 3)
        next

      normal <- coords_mat[neigh_idx, ] - center_coords
      face_dist <- sqrt(sum(normal^2)) / 2 # Dist to face = half bond length

      # Solid Angle (Oosterom & Strackee)
      sa <- calculate_solid_angle(center_coords, face_verts)

      # Area Calculation: A = 3 * Volume / distance
      # We must sum volumes of subtetrahedra formed by center + face triangles
      # To form triangles, order the face vertices around the normal

      # Basis construction
      u <- if (abs(normal[1]) < 0.9)
        c(1, 0, 0)
      else
        c(0, 1, 0)
      u <- cross_product(normal, u)
      u <- u / sqrt(sum(u^2))
      v <- cross_product(normal, u)
      v <- v / sqrt(sum(v^2))

      # Project vertices to 2D for ordering
      fc <- colMeans(face_verts)
      rel <- t(t(face_verts) - fc)
      angles <- atan2(rel %*% v, rel %*% u)
      ordered_verts <- face_verts[order(angles), , drop = FALSE]

      # Sum volumes
      total_vol <- 0
      n_v <- nrow(ordered_verts)
      for (k in 1:n_v) {
        p1 <- ordered_verts[1, ] # Fan center (on face)
        p2 <- ordered_verts[k, ]
        p3 <- ordered_verts[if (k == n_v)
          1
          else
            k + 1, ]
        # Tetrahedron formed by AtomCenter, and the 3 face points
        total_vol <- total_vol + vol_tetra(center_coords, p1, p2, p3)
      }

      area <- 3 * total_vol / face_dist

      results_list[[length(results_list) + 1]] <- list(
        Atom1 = prim_label,
        Atom2 = expanded_coords$Label[neigh_idx],
        Distance = face_dist * 2,
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

  # Normalize Weights
  dt[, MaxSA := max(SolidAngle, na.rm = TRUE), by = Atom1]
  dt[, Weight := ifelse(MaxSA > 0, SolidAngle / MaxSA, 0)]

  return(dt[Weight > tol])
}

#' @title Identify Atomic Bonds using CrystalNN
#' @description Port of Pymatgen's `CrystalNN` weighting logic.
#' @param distances Ignored.
#' @param atomic_coordinates Primary atom set.
#' @param expanded_coords Expanded supercell.
#' @param unit_cell_metrics Cell parameters.
#' @param cutoff_length Voronoi cutoff.
#' @param x_diff_weight Electronegativity weight (default 3.0).
#' @param porosity_adjustment Logical (default TRUE).
#' @param distance_cutoffs Vector c(0.5, 1.0).
#' @return A `data.table` of bonded pairs.
#' @family bonding algorithms
#' @export
crystal_nn <- function(distances,
                       atomic_coordinates,
                       expanded_coords,
                       unit_cell_metrics,
                       cutoff_length = 13.0,
                       x_diff_weight = 3.0,
                       porosity_adjustment = TRUE,
                       distance_cutoffs = c(0.5, 1.0)) {
  # 1. Get Base Voronoi Neighbors
  base <- voronoi_nn(
    atomic_coordinates,
    expanded_coords,
    unit_cell_metrics,
    cutoff = cutoff_length,
    tol = 0,
    allow_pathological = TRUE # Handled loosely here
  )
  if (is.null(base) || nrow(base) == 0)
    return(NULL)

  dt <- copy(base)
  dt[, `:=`(Sym1 = sub("[0-9].*", "", sub("_.*", "", Atom1)),
            Sym2 = sub("[0-9].*", "", sub("_.*", "", Atom2)))]

  # 2. Porosity Adjustment
  # Pymatgen: weight = solid_angle / area (then normalize)
  # R VoronoiNN returns Weight=SA/MaxSA.
  # We construct RawScore ~ SA * (SA/Area) to match logic
  if (porosity_adjustment) {
    dt[, RawScore := ifelse(Area > 1e-8, SolidAngle * (SolidAngle / Area), SolidAngle)]
  } else {
    dt[, RawScore := SolidAngle]
  }

  # 3. Chemical Adjustment (Electronegativity)
  if (x_diff_weight > 0) {
    dt[, EN1 := get_electronegativity(Sym1)]
    dt[, EN2 := get_electronegativity(Sym2)]

    dt[, ChemMod := 1 + x_diff_weight * sqrt(abs(EN1 - EN2) / 3.3)]
    dt[is.na(ChemMod), ChemMod := 1.0]

    dt[, RawScore := RawScore * ChemMod]
  }

  # 4. Normalize (Highest weight = 1.0)
  dt[, MaxScore := max(RawScore, na.rm = TRUE), by = Atom1]
  dt[, Weight := ifelse(MaxScore > 0, RawScore / MaxScore, 0)]

  # 5. Distance Cutoffs
  if (!is.null(distance_cutoffs)) {
    atom_info <- copy(atomic_coordinates)
    if (!"OxidationState" %in% names(atom_info))
      atom_info[, OxidationState := NA_real_]
    get_r_vec <- Vectorize(get_radius, vectorize.args = c("symbol", "oxidation_state"))

    dt[, Parent1 := sub("_.*", "", Atom1)]
    dt[, Parent2 := sub("_.*", "", Atom2)]
    dt[atom_info, on = c("Parent1" = "Label"), OS1 := i.OxidationState]
    dt[atom_info, on = c("Parent2" = "Label"), OS2 := i.OxidationState]

    dt[, R1 := get_r_vec(Sym1, OS1)]
    dt[, R2 := get_r_vec(Sym2, OS2)]

    dt[, Diameter := R1 + R2]
    # Fallback to covalent sum if 0
    dt[Diameter == 0, Diameter := get_radius(Sym1, NA) + get_radius(Sym2, NA)]

    dt[, CutLow := Diameter + distance_cutoffs[1]]
    dt[, CutHigh := Diameter + distance_cutoffs[2]]

    dt[, DistMod := data.table::fcase(Distance <= CutLow,
                                      1.0,
                                      Distance >= CutHigh,
                                      0.0,
                                      default = 0.5 * (cos((Distance - CutLow) / (CutHigh - CutLow) * pi
                                      ) + 1))]
    dt[is.na(DistMod), DistMod := 0]
    dt[, Weight := Weight * DistMod]
  }

  return(dt[Weight > 0, .(Atom1, Atom2, Distance, DeltaX, DeltaY, DeltaZ, Weight)])
}
