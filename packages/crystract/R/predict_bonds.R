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
#' @return A `data.table` of bonded pairs.
#' @family bonding algorithms
#' @export
econ_nn <- function(distances, atomic_coordinates, tol = 0.2) {
  if (nrow(distances) == 0)
    return(distances)

  # Check if OxidationState exists, else create it as numeric NA
  atom_info <- copy(atomic_coordinates)
  if (!"OxidationState" %in% names(atom_info))
    atom_info[, OxidationState := NA_real_]

  # Helper to vectorise get_radius
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

  # Get Radii
  dt[, R1 := get_r_vec(Sym1, OS1)]
  dt[, R2 := get_r_vec(Sym2, OS2)]

  # Fictive Radii (Hoppe Eq 1)
  # Handle cases where R1+R2 is 0 (though rare with defaults)
  dt[, RadiusSum := R1 + R2]
  dt[, FictiveR := ifelse(RadiusSum > 0, Distance * (R1 / RadiusSum), 0)]

  # Iterative MEFIR Solver (matches pymatgen _get_mean_fictive_ionic_radius)
  # Returns vector of MEFIRs corresponding to input groups
  calc_mefir_group <- function(fictive_r) {
    # Pymatgen Init: minimum_fir = min(fictive_ionic_radii)
    min_v <- min(fictive_r, na.rm = TRUE)
    if (min_v <= 0 || !is.finite(min_v))
      return(0)

    # First iteration uses min_v as the "previous" MEFIR
    mean_fir <- min_v

    # Iteration loop
    for (i in 1:100) {
      prev_mean_fir <- mean_fir

      # Hoppe Eq 2: Weighting based on previous MEFIR
      # Pymatgen: weighted_sum += fir * exp(1 - (fir/prev)^6)
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

  # Transformation matrix (Fractional -> Cartesian)
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

  # 2. Delaunay Tessellation
  # "QJ" (Joggle) is crucial for symmetric crystals to avoid coplanar/cospherical failures
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

  # Map atoms to tetrahedra
  atom_to_tets <- split(rep(1:n_tets, 4), as.vector(tess))

  # 4. Identify Asymmetric Unit Atoms within Supercell
  # We identify indices in expanded_coords that match atomic_coordinates
  # by distance (handling floating point tolerance).
  # Note: atomic_coordinates are fractional.

  asy_frac <- as.matrix(atomic_coordinates[, .(x_a, y_b, z_c)])
  exp_frac <- as.matrix(expanded_coords[, .(x_a, y_b, z_c)])

  # Find indices in expanded_coords that correspond to atomic_coordinates
  # We assume the asymmetric unit is present in the supercell (usually near origin)
  target_indices <- integer(nrow(asy_frac))

  for (i in 1:nrow(asy_frac)) {
    diffs <- t(t(exp_frac) - asy_frac[i, ])
    dists_sq <- rowSums(diffs^2)
    # Find closest match
    match_idx <- which.min(dists_sq)
    if (dists_sq[match_idx] < 1e-5) {
      target_indices[i] <- match_idx
    } else {
      # Fallback: if not found (rare), skip
      target_indices[i] <- NA
    }
  }
  target_indices <- na.omit(target_indices)

  results_list <- list()

  # 5. Loop ONLY over Asymmetric Unit Atoms
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

    # Calculate Distances
    dists <- sqrt(rowSums((t(
      t(coords_mat[potential_neighbors, , drop = FALSE]) - center_coords
    ))^2))

    # Filter by cutoff
    valid_mask <- dists <= cutoff
    neighbor_indices <- potential_neighbors[valid_mask]

    for (neigh_idx in neighbor_indices) {
      shared_tets <- intersect(my_tets, atom_to_tets[[as.character(neigh_idx)]])
      # A face must share at least 3 tetrahedra in Delaunay dual
      if (length(shared_tets) < 3)
        next

      face_verts <- circumcenters[shared_tets, , drop = FALSE]
      face_verts <- na.omit(face_verts)
      if (nrow(face_verts) < 3)
        next

      normal <- coords_mat[neigh_idx, ] - center_coords
      dist_sq <- sum(normal^2)
      if (is.na(dist_sq) || dist_sq < 1e-12)
        next

      # Solid Angle Calculation
      sa <- calculate_solid_angle(center_coords, face_verts)

      # Area Calculation
      # Project vertices to 2D plane perpendicular to normal for area,
      # or sum triangle areas from centroid.
      face_centroid <- colMeans(face_verts)
      area <- 0
      n_fv <- nrow(face_verts)

      # Order vertices for correct area calculation
      # 1. Create Basis u, v
      u <- if (abs(normal[1]) < 0.9)
        c(1, 0, 0)
      else
        c(0, 1, 0)
      u <- cross_product(normal, u)
      u <- u / sqrt(sum(u^2))
      v <- cross_product(normal, u)
      v <- v / sqrt(sum(v^2))

      # 2. Project relative vectors
      rel <- t(t(face_verts) - face_centroid)
      angles <- atan2(rel %*% v, rel %*% u)
      face_verts_ordered <- face_verts[order(angles), , drop = FALSE]

      # 3. Sum Triangles
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

  # Normalize Weights (SolidAngle / Max SolidAngle)
  # CrystalNN needs absolute SolidAngle, so we keep that column distinct
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
    tol = 0
  )
  if (is.null(base) || nrow(base) == 0)
    return(NULL)

  dt <- copy(base)
  dt[, `:=`(Sym1 = sub("[0-9].*", "", sub("_.*", "", Atom1)),
            Sym2 = sub("[0-9].*", "", sub("_.*", "", Atom2)))]

  # 2. Porosity Adjustment (Apply BEFORE normalization)
  # Pymatgen: weight = solid_angle / area
  # Note: Pymatgen calculates this, then normalizes by the Max(solid_angle/area) later
  if (porosity_adjustment) {
    dt[, RawScore := ifelse(Area > 1e-8, SolidAngle / Area, SolidAngle)]
  } else {
    dt[, RawScore := SolidAngle]
  }

  # 3. Chemical Adjustment (Electronegativity)
  if (x_diff_weight > 0) {
    dt[, EN1 := get_electronegativity(Sym1)]
    dt[, EN2 := get_electronegativity(Sym2)]

    dt[, ChemMod := 1 + x_diff_weight * sqrt(abs(EN1 - EN2) / 3.3)]
    dt[is.na(ChemMod), ChemMod := 1.0] # Handle missing EN

    dt[, RawScore := RawScore * ChemMod]
  }

  # 4. Normalize (Highest weight = 1.0)
  dt[, MaxScore := max(RawScore, na.rm = TRUE), by = Atom1]
  dt[, Weight := ifelse(MaxScore > 0, RawScore / MaxScore, 0)]

  # 5. Distance Cutoffs (Apply penalty AFTER normalization)
  if (!is.null(distance_cutoffs)) {
    # Get Radii
    atom_info <- copy(atomic_coordinates)
    if (!"OxidationState" %in% names(atom_info))
      atom_info[, OxidationState := NA_real_]
    get_r_vec <- Vectorize(get_radius, vectorize.args = c("symbol", "oxidation_state"))

    dt[, Parent2 := sub("_.*", "", Atom2)]

    # Safe lookup map
    rad_map <- data.table(Label = atom_info$Label,
                          R = get_r_vec(
                            sub("[0-9]+$", "", atom_info$Label),
                            atom_info$OxidationState
                          ))

    # Simplified Radius Lookup:
    dt[, R1 := get_r_vec(Sym1, NA)] # Fallback to covalent/atomic if matching fails
    dt[, R2 := get_r_vec(Sym2, NA)]

    dt[, Diameter := R1 + R2]
    dt[, CutLow := Diameter + distance_cutoffs[1]]
    dt[, CutHigh := Diameter + distance_cutoffs[2]]
    dt[is.na(Diameter), Diameter := 0]

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
