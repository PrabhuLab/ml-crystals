#' @title Identify Atomic Bonds using the Minimum Distance Method
#' @description Identifies bonded atoms by finding the nearest neighbor distance (d_min)
#'   for each central atom and defining a cutoff distance (d_cut) as
#'   d_cut = (1 + delta) * d_min.
#' @param distances A `data.table` of interatomic distances from `calculate_distances`.
#' @param delta The relative tolerance parameter (default 0.1).
#' @return A `data.table` of bonded pairs.
#' @family bonding algorithms
#' @export
#' @examples
#' dists <- data.table::data.table(Atom1 = c("A", "A"), Atom2 = c("B", "C"),
#'                                 Distance = c(1.5, 2.5),
#'                                 DeltaX = c(1, 0), DeltaY = c(0, 1), DeltaZ = c(0, 0))
#' minimum_distance(dists, delta = 0.1)
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
#' @references
#' Brunner, G. O. (1977). "A Definition of Coordination and Its Relevance in the
#' Structure Types AlB2 and NiAs." \emph{Acta Crystallographica Section A}, 33(1),
#' 226–227. \doi{10.1107/S0567739477000461}
#' @family bonding algorithms
#' @export
#' @examples
#' dists <- data.table::data.table(Atom1 = c("A", "A"), Atom2 = c("B", "C"),
#'                                 Distance = c(1.5, 2.5),
#'                                 DeltaX = c(1, 0), DeltaY = c(0, 1), DeltaZ = c(0, 0))
#' brunner_nn_reciprocal(dists)
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
#' @param distances A `data.table` of interatomic distances.
#' @param atomic_coordinates A `data.table` of atomic coordinates used to map species to radii.
#' @param tol A bond strength cutoff (default 0.2).
#' @param use_fictive_radius Logical. If `TRUE`, calculates Hoppe's fictive ionic radius.
#' @return A `data.table` of bonded pairs.
#' @references
#' Hoppe, R. (1979). "Effective Coordination Numbers (ECoN) and Mean Fictive Ionic
#' Radii (MEFIR)." \emph{Zeitschrift Für Kristallographie}, 150(1–4), 23–52.
#' \doi{10.1524/zkri.1979.150.14.23}
#'
#' Shannon, R. D. (1976). "Revised Effective Ionic Radii and Systematic Studies of
#' Interatomic Distances in Halides and Chalcogenides." \emph{Acta Crystallographica
#' Section A}, 32(5), 751–767. \doi{10.1107/S0567739476001551}
#' @family bonding algorithms
#' @export
#' @examples
#' dists <- data.table::data.table(Atom1 = c("Si1", "Si1"), Atom2 = c("O1", "O2"),
#'                                 Distance = c(1.6, 2.0),
#'                                 DeltaX = c(1, 0), DeltaY = c(0, 1), DeltaZ = c(0, 0))
#' ac <- data.table::data.table(Label = c("Si1", "O1", "O2"),
#'                              OxidationState = c(4, -2, -2))
#' econ_nn(dists, ac)
econ_nn <- function(distances,
                    atomic_coordinates,
                    tol = 0.2,
                    use_fictive_radius = FALSE) {
  if (nrow(distances) == 0)
    return(distances)

  dt <- copy(distances)

  if (use_fictive_radius) {
    atom_info <- copy(atomic_coordinates)
    if (!"OxidationState" %in% names(atom_info))
      atom_info[, OxidationState := NA_real_]

    get_r_ionic <- Vectorize(get_radius_pymatgen_style,
                             vectorize.args = c("symbol", "oxidation_state"))
    get_r_def <- Vectorize(function(sym)
      get_radius_pymatgen_style(sym, NA_real_, 6),
      vectorize.args = "sym")

    dt[, `:=`(
      Sym1 = sub("[0-9].*", "", sub("_.*", "", Atom1)),
      Sym2 = sub("[0-9].*", "", sub("_.*", "", Atom2)),
      Parent1 = sub("_.*", "", Atom1),
      Parent2 = sub("_.*", "", Atom2)
    )]

    dt[atom_info, on = c("Parent1" = "Label"), OS1 := i.OxidationState]
    dt[atom_info, on = c("Parent2" = "Label"), OS2 := i.OxidationState]

    dt[, R1 := get_r_ionic(Sym1, OS1, cn = 6)]
    dt[, R2 := get_r_ionic(Sym2, OS2, cn = 6)]

    # Fallbacks for fictive radius math
    dt[R1 == 0, R1 := get_r_def(Sym1)]
    dt[R2 == 0, R2 := get_r_def(Sym2)]

    dt[, RadiusSum := R1 + R2]
    dt[, FictiveR := ifelse(RadiusSum > 0, Distance * (R1 / RadiusSum), 0)]
  } else {
    dt[, FictiveR := Distance]
  }

  calc_mefir_group <- function(fictive_r) {
    min_v <- min(fictive_r, na.rm = TRUE)
    if (min_v <= 0 || !is.finite(min_v))
      return(0)

    mean_fir <- min_v
    for (i in 1:100) {
      prev_mean_fir <- mean_fir
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
  dt[, Weight := ifelse(MEFIR > 0, exp(1 - (FictiveR / MEFIR)^6), 0)]

  return(dt[Weight > tol, .(Atom1, Atom2, Distance, DeltaX, DeltaY, DeltaZ, Weight)])
}

#' @title Identify Atomic Bonds using Voronoi Tessellation
#' @description Performs 3D Voronoi analysis on the supercell.
#' @param atomic_coordinates A `data.table` of the primary (asymmetric) atom set.
#' @param expanded_coords A `data.table` of atoms in the expanded supercell.
#' @param unit_cell_metrics A `data.table` with cell parameters.
#' @param cutoff Distance cutoff (default 13.0).
#' @param tol Tolerance for solid angle weights (default 0).
#' @return A `data.table` of bonded pairs.
#' @references
#' O'Keeffe, M. (1979). "A Proposed Rigorous Definition of Coordination Number."
#' \emph{Acta Crystallographica Section A}, 35(5), 772–775. \doi{10.1107/S0567739479001765}
#'
#' Aurenhammer, F., Klein, R., & Lee, D.-T. (2013). \emph{Voronoi Diagrams and
#' Delaunay Triangulations}. World Scientific Publishing Company.
#' @family bonding algorithms
#' @export
#' @examples
#' \donttest{
#' cif_path <- system.file("extdata", "1590946.cif", package = "crystract")
#' if (file.exists(cif_path)) {
#'   old_threads <- data.table::setDTthreads(2)
#'   res <- analyze_single_cif(cif_path, bonding_algorithms = "voronoi")
#'   print(head(res$bonds_voronoi[[1]]))
#'   data.table::setDTthreads(old_threads)
#' }
#' }
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

  # 2. Delaunay Tessellation ("QJ" is crucial to prevent crashes on symmetric grids)
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

      unique_verts <- raw_verts[1, , drop = FALSE]
      if (nrow(raw_verts) > 1) {
        for (k in 2:nrow(raw_verts)) {
          pt <- raw_verts[k, ]
          dists_v <- sqrt(rowSums(t(t(
            unique_verts
          ) - pt)^2))
          if (all(dists_v > 1e-5)) {
            unique_verts <- rbind(unique_verts, pt)
          } else {
            match_idx <- which.min(dists_v)
            unique_verts[match_idx, ] <- (unique_verts[match_idx, ] + pt) / 2
          }
        }
      }
      if (nrow(unique_verts) < 3)
        next

      face_centroid <- colMeans(unique_verts)
      u <- if (abs(normal[1]) < 0.9)
        c(1, 0, 0)
      else
        c(0, 1, 0)
      u <- cross_product(normal, u)
      u <- u / sqrt(sum(u^2))
      v <- cross_product(normal, u)
      v <- v / sqrt(sum(v^2))

      rel <- t(t(unique_verts) - face_centroid)
      angles <- atan2(rel %*% v, rel %*% u)
      face_verts_ordered <- unique_verts[order(angles), , drop = FALSE]

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

      if (area < 1e-10 || sa < 1e-10)
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
  raw_dt <- rbindlist(results_list)

  # Reconstruct true mathematical faces by summing the shattered shards belonging to the same atom pair.
  dt <- raw_dt[, .(
    Distance = min(Distance),
    SolidAngle = sum(SolidAngle),
    Area = sum(Area),
    DeltaX = mean(DeltaX),
    DeltaY = mean(DeltaY),
    DeltaZ = mean(DeltaZ)
  ), by = .(Atom1, Atom2)]

  # Apply the physical filter now that the faces are reconstructed
  dt <- dt[Area >= 1e-5 & SolidAngle >= 1e-5]
  if (nrow(dt) == 0)
    return(NULL)

  dt[, MaxSA := max(SolidAngle, na.rm = TRUE), by = Atom1]
  dt[, Weight := ifelse(MaxSA > 0, SolidAngle / MaxSA, 0)]

  return(dt[Weight > tol])
}

#' @title Identify Atomic Bonds using CrystalNN
#' @description Port of Pymatgen's `CrystalNN` weighting logic. Includes Porosity,
#' Electronegativity, and Distance penalties.
#' Uses specific Decision Tree logic for radii (Shannon -> Covalent -> Atomic) per
#' pymatgen's fallback logic.
#' @param distances Ignored.
#' @param atomic_coordinates Primary atom set.
#' @param expanded_coords Expanded supercell.
#' @param unit_cell_metrics Cell parameters.
#' @param cutoff_length Voronoi cutoff. Default 7.0 matches Pymatgen CrystalNN default.
#' @param x_diff_weight Electronegativity weight (default 3.0).
#' @param porosity_adjustment Logical (default TRUE).
#' @param distance_cutoffs Vector c(0.5, 1.0).
#' @return A `data.table` of bonded pairs.
#' @references
#' Shannon, R. D. (1976). "Revised effective ionic radii and systematic studies of
#' interatomic distances in halides and chalcogenides." \emph{Acta Crystallographica
#' Section A}, 32(5), 751-767.
#'
#' Pauling, L. (1960). \emph{The Nature of the Chemical Bond and the Structure of
#' Molecules and Crystals: An Introduction to Modern Structural Chemistry}.
#' Cornell University Press.
#'
#' Pan, H., Ganose, A. M., Horton, M., et al. (2021). "Benchmarking Coordination Number
#' Prediction Algorithms on Inorganic Crystal Structures." \emph{Inorganic Chemistry},
#' 60(3), 1590–1603. \doi{10.1021/acs.inorgchem.0c02996}
#' @family bonding algorithms
#' @export
#' @examples
#' \donttest{
#' cif_path <- system.file("extdata", "1590946.cif", package = "crystract")
#' if (file.exists(cif_path)) {
#'   old_threads <- data.table::setDTthreads(2)
#'   res <- analyze_single_cif(cif_path, bonding_algorithms = "crystal_nn")
#'   print(head(res$bonds_crystal_nn[[1]]))
#'   data.table::setDTthreads(old_threads)
#' }
#' }
crystal_nn <- function(distances,
                       atomic_coordinates,
                       expanded_coords,
                       unit_cell_metrics,
                       cutoff_length = 7.0,
                       # Pymatgen default for CrystalNN is 7.0 (search_cutoff)
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

  # Merge Oxidation States if available (for radii/EN logic)
  atom_info <- copy(atomic_coordinates)
  if (!"OxidationState" %in% names(atom_info))
    atom_info[, OxidationState := NA_real_]

  dt[, Parent1 := sub("_.*", "", Atom1)]
  dt[, Parent2 := sub("_.*", "", Atom2)]
  dt[atom_info, on = c("Parent1" = "Label"), OS1 := i.OxidationState]
  dt[atom_info, on = c("Parent2" = "Label"), OS2 := i.OxidationState]

  # --- Step 1: Porosity Adjustment ---
  # Pymatgen Logic: x['weight'] *= x['poly_info']['solid_angle'] / x['poly_info']['area']
  # Since x['weight'] IS the SolidAngle (from VoronoiNN),
  # result = SolidAngle * (SolidAngle / Area) = SolidAngle^2 / Area.
  if (porosity_adjustment) {
    dt[, RawScore := ifelse(Area > 1e-8, (SolidAngle * SolidAngle) / Area, SolidAngle)]
  } else {
    dt[, RawScore := SolidAngle]
  }

  # --- Step 2: Chemical Adjustment (Electronegativity) ---
  if (x_diff_weight > 0) {
    dt[, EN1 := get_electronegativity(Sym1)]
    dt[, EN2 := get_electronegativity(Sym2)]

    # 3.3 is max deltaX constant from pymatgen
    dt[, ChemMod := 1 + x_diff_weight * sqrt(abs(EN1 - EN2) / 3.3)]
    dt[is.na(ChemMod), ChemMod := 1.0]

    dt[, RawScore := RawScore * ChemMod]
  }

  # --- Step 3: Normalize (Highest weight = 1.0) ---
  dt[, MaxScore := max(RawScore, na.rm = TRUE), by = Atom1]
  dt[, Weight := ifelse(MaxScore > 0, RawScore / MaxScore, 0)]

  # --- Step 4: Distance Cutoffs (Penalty applied AFTER normalization) ---
  if (!is.null(distance_cutoffs)) {
    # DEFINE atomic_rads IN LOCAL SCOPE FOR Vectorize TO FIND
    atomic_rads <- get_radii_data()

    get_r_vec <- Vectorize(get_radius_decision_tree,
                           vectorize.args = c("symbol", "oxidation_state"))

    # Use Decision Tree Logic: Shannon -> Covalent -> Atomic
    dt[, R1 := get_r_vec(Sym1, OS1, cn = 6)]
    dt[, R2 := get_r_vec(Sym2, OS2, cn = 6)]

    # Pymatgen Logic: If specific radii failed (returned 0 by fallback),
    # fallback to default Covalent/Atomic
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

    # Calculate penalty
    dt[, DistMod := data.table::fcase(Distance <= CutLow,
                                      1.0,
                                      Distance >= CutHigh,
                                      0.0,
                                      default = 0.5 * (cos((Distance - CutLow) / (CutHigh - CutLow) * pi
                                      ) + 1))]
    dt[is.na(DistMod), DistMod := 0]

    dt[, Weight := Weight * DistMod]
  }

  # --- Step 5: Pymatgen Rounding Logic ---
  # Pymatgen explicitly rounds weights to 3 decimal places.
  # This acts as a filter for very small weights (e.g., 0.0004 -> 0).
  dt[, Weight := round(Weight, 3)]

  return(dt[Weight > 0, .(Atom1, Atom2, Distance, DeltaX, DeltaY, DeltaZ, Weight)])
}
