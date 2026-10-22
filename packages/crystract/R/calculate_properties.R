#' @title Calculate Interatomic Distances
#' @description Computes the distances between a central set of atoms and an
#'   expanded set, using the metric tensor for accuracy.
#' @details
#' With \eqn{\Delta x=x_{f1}-x_{f2}}, \eqn{\Delta y=y_{f1}-y_{f2}} and
#' \eqn{\Delta z=z_{f1}-z_{f2}}, the interatomic distance is
#' \deqn{d=\left[a^2\Delta x^2+b^2\Delta y^2+c^2\Delta z^2+
#' 2ab\Delta x\Delta y\cos\gamma+2ac\Delta x\Delta z\cos\beta+
#' 2bc\Delta y\Delta z\cos\alpha\right]^{1/2}.}
#' The columns `x_a`, `y_b`, `z_c` contain fractional coordinates
#' \eqn{x_f,y_f,z_f}; `z_c` is not a Cartesian coordinate here.
#' @param atomic_coordinates A `data.table` of the primary (asymmetric) atom set.
#' @param expanded_coords A `data.table` of atoms in the expanded supercell.
#' @param unit_cell_metrics A `data.table` with cell parameters.
#' @param tolerance A numeric cutoff (default 1e-6). Distances smaller than this
#'   value are considered floating-point noise (overlapping atoms) and are
#'   filtered out.
#' @return A `data.table` of all non-zero distances.
#' @family property calculators
#' @export
#' @examples
#' ac <- data.table::data.table(Label = "Si1", x_a = 0, y_b = 0, z_c = 0)
#' ec <- data.table::data.table(Label = "O1_1", x_a = 0.5, y_b = 0.5, z_c = 0.5)
#' uc <- data.table::data.table(`_cell_length_a` = 10, `_cell_length_b` = 10,
#'                              `_cell_length_c` = 10, `_cell_angle_alpha` = 90,
#'                              `_cell_angle_beta` = 90, `_cell_angle_gamma` = 90)
#' calculate_distances(ac, ec, uc)
calculate_distances <- function(atomic_coordinates,
                                expanded_coords,
                                unit_cell_metrics,
                                tolerance = 1e-6) {
  # Safe extraction of cell metrics
  a <- as.numeric(unit_cell_metrics$`_cell_length_a`)
  b <- as.numeric(unit_cell_metrics$`_cell_length_b`)
  c <- as.numeric(unit_cell_metrics$`_cell_length_c`)
  alpha <- as.numeric(unit_cell_metrics$`_cell_angle_alpha`) * pi / 180
  beta <- as.numeric(unit_cell_metrics$`_cell_angle_beta`) * pi / 180
  gamma <- as.numeric(unit_cell_metrics$`_cell_angle_gamma`) * pi / 180

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

  # FIX: stringsAsFactors = FALSE prevents "Ops.factor" warnings downstream
  atom_pairs <- expand.grid(
    Atom1 = labels_atomic,
    Atom2 = labels_expanded,
    KEEP.OUT.ATTRS = TRUE,
    stringsAsFactors = FALSE
  )

  distances <- data.table(
    Atom1 = atom_pairs$Atom1,
    Atom2 = atom_pairs$Atom2,
    Distance = as.vector(r),
    DeltaX = as.vector(delta_x),
    DeltaY = as.vector(delta_y),
    DeltaZ = as.vector(delta_z)
  )

  # Filter out self-interactions or near-zero distances based on tolerance
  distances <- distances[Distance > tolerance]

  return(distances)
}

#' @title Calculate Coordination Numbers
#' @description Counts the number of nearest neighbors for each central atom based on
#'   a table of bonded pairs.
#' @param bonded_pairs_table A `data.table` of bonded pairs.
#' @return A `data.table` with columns 'Atom' and 'CoordinationNumber'.
#' @family property calculators
#' @export
#' @examples
#' bp <- data.table::data.table(Atom1 = c("Si1", "Si1", "O1"),
#'                              Atom2 = c("O1", "O2", "Si1"))
#' calculate_neighbor_counts(bp)
calculate_neighbor_counts <- function(bonded_pairs_table) {
  if (is.null(bonded_pairs_table) || nrow(bonded_pairs_table) == 0) {
    return(data.table(Atom = character(), CoordinationNumber = integer()))
  }

  neighbor_counts <- bonded_pairs_table[, .(CoordinationNumber = .N), by = .(Atom1)]

  setnames(neighbor_counts, "Atom1", "Atom")

  return(neighbor_counts)
}

#' @title Calculate Weighted Coordination Numbers
#' @description Counts the weighted number of nearest neighbors for each central atom based on
#'   a table of bonded pairs, accounting for the fractional occupancy of the neighbor sites.
#'   This is particularly useful for analyzing disordered crystal structures.
#'
#'   If all occupancies are 1.0, this efficiently returns the standard coordination number.
#'   If the fractional occupancies on any single crystallographic site sum to greater than 1.0,
#'   this indicates an invalid or physically impossible structure; the function will throw a
#'   warning and return `NULL` so the file can be discarded by the analysis pipeline.
#'
#' @details
#' The occupancy-corrected coordination number is
#' \deqn{\mathrm{CN}_{\mathrm{occ}}=\sum_{i\in\mathrm{neighbour\ sites}}m_i\,\mathrm{occ}_i.}
#' Here \eqn{m_i} is the local number of appearances of asymmetric neighbour
#' site \eqn{i} in the central atom's coordination sphere, not its unit-cell
#' Wyckoff multiplicity, and \eqn{\mathrm{occ}_i} is its site occupancy.
#' @param bonded_pairs_table A `data.table` of bonded pairs (e.g., from `minimum_distance`).
#' @param atomic_coordinates A `data.table` of asymmetric atoms from `extract_atomic_coordinates`.
#' @return A `data.table` with columns 'Atom', 'CoordinationNumber', and 'WeightedCoordinationNumber'.
#'   Returns `NULL` if the occupancies sum to > 1 on any shared site.
#' @family property calculators
#' @export
#' @examples
#' bp <- data.table::data.table(Atom1 = c("Si1", "Si1"),
#'                              Atom2 = c("O1_1_0_0", "O2_1_0_0"))
#' ac <- data.table::data.table(Label = c("Si1", "O1", "O2"),
#'                              Occupancy = c(1.0, 0.5, 0.5),
#'                              x_a=c(0,0,0), y_b=c(0,0,0), z_c=c(0,0,0))
#' calculate_weighted_neighbor_counts(bp, ac)
calculate_weighted_neighbor_counts <- function(bonded_pairs_table, atomic_coordinates) {
  if (is.null(bonded_pairs_table) || nrow(bonded_pairs_table) == 0) {
    return(
      data.table(
        Atom = character(),
        CoordinationNumber = integer(),
        WeightedCoordinationNumber = numeric()
      )
    )
  }

  coords <- copy(atomic_coordinates)

  # Ensure Occupancy column exists and handle NAs
  if (!"Occupancy" %in% names(coords)) {
    coords[, Occupancy := 1.0]
  } else {
    coords[is.na(Occupancy), Occupancy := 1.0]
  }

  # --- 1. VALIDATION ---
  # Check if fractional occupancies sum to > 1 at any physical crystallographic site.
  # Grouping by exact snapped coordinates identifies atoms sharing the same site.
  site_occ <- coords[, .(SumOccupancy = sum(Occupancy)), by = .(x_a, y_b, z_c)]

  # We use 1.01 to allow for minor floating-point rounding errors (e.g., 1/3 + 2/3)
  if (any(site_occ$SumOccupancy > 1.01)) {
    warning(
      "Fractional occupancies add up to greater than 1.0 on one or more sites. Discarding file."
    )
    return(NULL)
  }

  # --- 2. OPTIMIZATION ---
  # If all occupancies are fully occupied (1.0), skip the fractional math (WCN == CN)
  if (all(abs(coords$Occupancy - 1.0) < 1e-5)) {
    neighbor_counts <- bonded_pairs_table[, .(CoordinationNumber = .N,
                                              WeightedCoordinationNumber = as.numeric(.N)), by = .(Atom1)]

    setnames(neighbor_counts, "Atom1", "Atom")
    return(neighbor_counts)
  }

  # --- 3. CALCULATE WEIGHTED CN ---
  work_dt <- copy(bonded_pairs_table)

  # Extract the parent site label (e.g., "Fe1" from "Fe1_1_0_0")
  work_dt[, NeighborSite := sub("_.*", "", Atom2)]

  # Count how many times a specific asymmetric site appears as a neighbor (multiplicity m_i)
  site_counts <- work_dt[, .(m_i = .N), by = .(Atom1, NeighborSite)]

  # Merge occupancy info from the asymmetric unit coordinates
  site_counts <- merge(
    site_counts,
    coords[, .(Label, Occupancy)],
    by.x = "NeighborSite",
    by.y = "Label",
    all.x = TRUE
  )

  # Handle any missing occupancies resulting from the merge just in case
  site_counts[is.na(Occupancy), Occupancy := 1.0]

  # Weighted Sum: Sum(Multiplicity * Occupancy)
  w_cn <- site_counts[, .(
    CoordinationNumber = sum(m_i),
    WeightedCoordinationNumber = sum(m_i * Occupancy)
  ), by = .(Atom1)]

  setnames(w_cn, "Atom1", "Atom")

  return(w_cn)
}

#' @title Calculate Bond Angles
#' @description Calculates all bond angles centered on each atom, formed by pairs
#'   of its bonded neighbors.
#' @details
#' Atom \eqn{A_1} is central and
#' \eqn{\mathbf{c}_k=(x_{ck},y_{ck},z_{ck})^{\mathsf T}=\mathbf{M}\mathbf{f}_k}.
#' Set \eqn{\mathbf{a}=\mathbf{c}_2-\mathbf{c}_1} and
#' \eqn{\mathbf{b}=\mathbf{c}_3-\mathbf{c}_1}. Then
#' \deqn{C=\frac{\mathbf{a}\cdot\mathbf{b}}{\lVert\mathbf{a}\rVert\lVert\mathbf{b}\rVert},
#' \qquad \theta=\arccos C.}
#' Bold \eqn{\mathbf{a},\mathbf{b}} are bond vectors; scalar \eqn{a,b,c}
#' are cell lengths. The implementation uses the equivalent fractional-coordinate
#' metric product, with \eqn{\mathbf{G}=\mathbf{M}^{\mathsf T}\mathbf{M}},
#' and returns \eqn{\theta} in degrees.
#' @param bonded_pairs Data.table of bonded atoms.
#' @param atomic_coordinates Data.table of asymmetric atom coordinates.
#' @param expanded_coords Data.table of supercell atom coordinates.
#' @param unit_cell_metrics Data.table with unit cell parameters.
#' @return A `data.table` of all unique bond angles.
#' @family property calculators
#' @export
#' @examples
#' bp <- data.table::data.table(Atom1 = c("Si1", "Si1"), Atom2 = c("O1", "O2"))
#' ac <- data.table::data.table(Label = "Si1", x_a = 0, y_b = 0, z_c = 0)
#' ec <- data.table::data.table(Label = c("O1", "O2"),
#'                              x_a = c(1, 0), y_b = c(0, 1), z_c = c(0, 0))
#' uc <- data.table::data.table(`_cell_length_a` = 10, `_cell_length_b` = 10,
#'                              `_cell_length_c` = 10, `_cell_angle_alpha` = 90,
#'                              `_cell_angle_beta` = 90, `_cell_angle_gamma` = 90)
#' calculate_angles(bp, ac, ec, uc)
calculate_angles <- function(bonded_pairs,
                             atomic_coordinates,
                             expanded_coords,
                             unit_cell_metrics) {
  if (is.null(bonded_pairs) || nrow(bonded_pairs) == 0) {
    return(
      data.table(
        CentralAtom = character(),
        Neighbor1 = character(),
        Neighbor2 = character(),
        Angle = numeric()
      )
    )
  }

  # Extract unit cell parameters
  a <- unit_cell_metrics$`_cell_length_a`
  b <- unit_cell_metrics$`_cell_length_b`
  c <- unit_cell_metrics$`_cell_length_c`
  alpha_rad <- unit_cell_metrics$`_cell_angle_alpha` * pi / 180
  beta_rad <- unit_cell_metrics$`_cell_angle_beta` * pi / 180
  gamma_rad <- unit_cell_metrics$`_cell_angle_gamma` * pi / 180

  cos_alpha <- cos(alpha_rad)
  cos_beta <- cos(beta_rad)
  cos_gamma <- cos(gamma_rad)

  # Combine all unique atomic coordinates and set key for fast lookup
  # (Including expanded_coords ensures all possible neighbors are available)
  all_coords <- unique(rbind(atomic_coordinates[, .(Label, x_a, y_b, z_c)], expanded_coords[, .(Label, x_a, y_b, z_c)], fill = TRUE),
                       by = "Label")
  setkey(all_coords, Label)

  # Generate all unique (CentralAtom, Neighbor1, Neighbor2) triplets
  # Start with central atoms and their first neighbors
  angle_triplets <- bonded_pairs[, .(CentralAtom = Atom1, Neighbor1 = Atom2)]

  # Self-join to find all pairs of neighbors for each CentralAtom
  angle_triplets <- merge(
    angle_triplets,
    bonded_pairs[, .(CentralAtom = Atom1, Neighbor2 = Atom2)],
    by = "CentralAtom",
    allow.cartesian = TRUE # Needed to get all combinations
  )

  # Filter out cases where Neighbor1 is the same as Neighbor2
  angle_triplets <- angle_triplets[Neighbor1 != Neighbor2]

  # Ensure consistent ordering of Neighbor1 and Neighbor2 to avoid duplicate angles
  # (e.g., A-B-C is the same as C-B-A when B is central)
  angle_triplets[, `:=`(
    Neighbor1_sorted = pmin(Neighbor1, Neighbor2),
    Neighbor2_sorted = pmax(Neighbor1, Neighbor2)
  )]
  angle_triplets <- unique(angle_triplets,
                           by = c("CentralAtom", "Neighbor1_sorted", "Neighbor2_sorted"))

  # Merge fractional coordinates for all three atoms involved in each angle
  angle_triplets[all_coords, on = .(CentralAtom = Label), `:=`(cx = i.x_a, cy = i.y_b, cz = i.z_c)]
  angle_triplets[all_coords, on = .(Neighbor1_sorted = Label), `:=`(n1x = i.x_a, n1y = i.y_b, n1z = i.z_c)]
  angle_triplets[all_coords, on = .(Neighbor2_sorted = Label), `:=`(n2x = i.x_a, n2y = i.y_b, n2z = i.z_c)]

  # Calculate vectors v1 (CentralAtom -> Neighbor1) and v2 (CentralAtom -> Neighbor2)
  angle_triplets[, `:=`(
    xf1 = n1x - cx,
    yf1 = n1y - cy,
    zf1 = n1z - cz,
    xf2 = n2x - cx,
    yf2 = n2y - cy,
    zf2 = n2z - cz
  )]

  # Vectorized calculation of dot product using metric tensor
  angle_triplets[, dot_product := (
    xf1 * xf2 * a^2 + yf1 * yf2 * b^2 + zf1 * zf2 * c^2 +
      (xf1 * yf2 + yf1 * xf2) * a * b * cos_gamma +
      (xf1 * zf2 + zf1 * xf2) * a * c * cos_beta +
      (yf1 * zf2 + zf1 * yf2) * b * c * cos_alpha
  )]

  # Vectorized calculation of squared magnitudes for both vectors
  angle_triplets[, mag_sq1 := (
    xf1^2 * a^2 + yf1^2 * b^2 + zf1^2 * c^2 +
      2 * xf1 * yf1 * a * b * cos_gamma +
      2 * xf1 * zf1 * a * c * cos_beta +
      2 * yf1 * zf1 * b * c * cos_alpha
  )]
  angle_triplets[, mag_sq2 := (
    xf2^2 * a^2 + yf2^2 * b^2 + zf2^2 * c^2 +
      2 * xf2 * yf2 * a * b * cos_gamma +
      2 * xf2 * zf2 * a * c * cos_beta +
      2 * yf2 * zf2 * b * c * cos_alpha
  )]

  # Calculate magnitudes, guarding against division by zero for very small lengths
  angle_triplets[, `:=`(
    mag1 = sqrt(pmax(mag_sq1, 0)),
    # pmax(...,0) to handle potential negative due to FP error
    mag2 = sqrt(pmax(mag_sq2, 0))
  )]

  # Calculate cos_theta, guarding against division by zero or invalid acos input range [-1, 1]
  angle_triplets[, cos_theta := pmin(pmax(dot_product / (mag1 * mag2), -1.0), 1.0)]
  angle_triplets[mag1 <= 1e-10 |
                   mag2 <= 1e-10, cos_theta := NA_real_] # Set to NA if magnitudes are too small

  # Calculate angle in degrees
  angle_triplets[, Angle := acos(cos_theta) * 180 / pi]

  # Filter out NA angles (if any, due to impossible geometry) and select final columns
  result_angles <- angle_triplets[!is.na(Angle), .(
    CentralAtom = CentralAtom,
    Neighbor1 = Neighbor1_sorted,
    Neighbor2 = Neighbor2_sorted,
    Angle
  )]

  return(result_angles[order(CentralAtom, Neighbor1, Neighbor2)])
}

# Internal geometry for both uncertainty propagators. The cell angles are radians.
# Sources: IUCr Cartesian coordinates; Angel et al.
# (2025), doi:10.1107/S2053273325002682; JCGM 100:2008.
uncertainty_cell <- function(unit_cell_metrics) {
  cell <- c(a = unit_cell_metrics$`_cell_length_a`,
            b = unit_cell_metrics$`_cell_length_b`,
            c = unit_cell_metrics$`_cell_length_c`,
            alpha = unit_cell_metrics$`_cell_angle_alpha` * pi / 180,
            beta = unit_cell_metrics$`_cell_angle_beta` * pi / 180,
            gamma = unit_cell_metrics$`_cell_angle_gamma` * pi / 180)
  with(as.list(cell), {
    ca <- cos(alpha); cb <- cos(beta); cg <- cos(gamma)
    sa <- sin(alpha); sb <- sin(beta); sg <- sin(gamma)
    D <- 1 - ca^2 - cb^2 - cg^2 + 2 * ca * cb * cg
    if (any(!is.finite(cell)) || D <= 0 || abs(sg) < 1e-12)
      stop("Invalid unit cell for uncertainty propagation")
    v <- sqrt(D)
    # Stage 1: convert fractional position f to Cartesian c = M f.
    M <- rbind(c(a, b * cg, c * cb),
               c(0, b * sg, c * (ca - cb * cg) / sg),
               c(0, 0, c * v / sg))
    # Stage 2: differentiate M in each independent cell parameter. In
    # particular, dz_c/dalpha and dz_c/dbeta include their sine factors.
    dM <- list(
      a = rbind(c(1, 0, 0), c(0, 0, 0), c(0, 0, 0)),
      b = rbind(c(0, cg, 0), c(0, sg, 0), c(0, 0, 0)),
      c = rbind(c(0, 0, cb), c(0, 0, (ca - cb * cg) / sg), c(0, 0, v / sg)),
      alpha = rbind(c(0, 0, 0), c(0, 0, -c * sa / sg),
                    c(0, 0, c * sa * (ca - cb * cg) / (v * sg))),
      beta = rbind(c(0, 0, -c * sb), c(0, 0, c * sb * cg / sg),
                   c(0, 0, c * sb * (cb - ca * cg) / (v * sg))),
      gamma = rbind(c(0, -b * sg, 0),
                    c(0, b * cg, c * (cb - ca * cg) / sg^2),
                    c(0, 0, c * (sg^2 * (cg - ca * cb) - D * cg) / (v * sg^2)))
    )
    list(M = M, dM = dM)
  })
}

# A site stores df/dA_j for every original asymmetric-site parent A_j.
uncertainty_site <- function(label, atomic_coordinates, expanded_coords = NULL) {
  i <- match(label, atomic_coordinates$Label)
  source <- atomic_coordinates
  if (is.na(i) && !is.null(expanded_coords)) {
    i <- match(label, expanded_coords$Label)
    source <- expanded_coords
  }
  if (is.na(i)) return(NULL)
  map <- if ("ParentJacobian" %in% names(source)) source$ParentJacobian[[i]] else NULL
  if (is.null(map) && identical(source, atomic_coordinates))
    map <- stats::setNames(list(diag(3)), label)
  if (is.null(map)) return(NULL) # A generated label alone cannot reveal W.
  list(f = as.numeric(source[i, .(x_a, y_b, z_c)]), map = map)
}

uncertainty_add_parent <- function(gradient, site, fractional_gradient) {
  # Stage 3: apply the chain rule and add shared-parent terms before squaring.
  for (parent in names(site$map)) {
    contribution <- as.numeric(t(site$map[[parent]]) %*% fractional_gradient)
    for (axis in seq_along(contribution)) {
      key <- paste0(parent, ":", c("x", "y", "z")[axis])
      gradient[key] <- if (key %in% names(gradient)) gradient[key] + contribution[axis] else contribution[axis]
    }
  }
  gradient
}

uncertainty_variance <- function(gradient, unit_cell_metrics,
                                 independent_coordinates, covariance = NULL) {
  # Stage 4: keep one entry per independent input p. A full named Sigma_p
  # may be supplied; otherwise use reported CIF variances on the diagonal.
  if (any(!is.finite(gradient))) return(NA_real_)
  gradient <- gradient[abs(gradient) > 1e-12]
  if (!length(gradient)) return(0)
  if (!is.null(covariance)) {
    if (is.null(rownames(covariance)) || is.null(colnames(covariance)) ||
        !all(names(gradient) %in% rownames(covariance)) ||
        !all(names(gradient) %in% colnames(covariance)))
      stop("covariance must name every independent input with a nonzero gradient")
    Sigma <- covariance[names(gradient), names(gradient), drop = FALSE]
    if (any(!is.finite(Sigma))) return(NA_real_)
    return(as.numeric(t(gradient) %*% Sigma %*% gradient))
  }
  cell_su <- vapply(c("a", "b", "c", "alpha", "beta", "gamma"), function(name) {
    field <- paste0(if (name %in% c("a", "b", "c"))
      "_cell_length_" else "_cell_angle_", name, "_error")
    value <- unit_cell_metrics[[field]]
    if (!length(value)) return(NA_real_)
    as.numeric(value[[1]]) * if (name %in% c("a", "b", "c")) 1 else pi / 180
  }, numeric(1))
  names(cell_su) <- paste0("cell_", names(cell_su))
  atom_su <- unlist(lapply(seq_len(nrow(independent_coordinates)), function(i) {
    x <- unlist(independent_coordinates[i, .(x_error, y_error, z_error)], use.names = FALSE)
    stats::setNames(x, paste0(independent_coordinates$Label[i], ":", c("x", "y", "z")))
  }))
  su <- c(cell_su, atom_su)[names(gradient)]
  if (anyNA(su) || any(!is.finite(su)) || any(su < 0)) return(NA_real_)
  sum((gradient * su)^2)
}

#' @title Propagate Distance Error
#' @description Calculates the standard uncertainty for each interatomic distance.
#' @details
#' The uncertainty is calculated using each independently refined
#' input once in \eqn{\mathbf{p}}:
#' \deqn{\sigma_d^2=(\nabla_{\mathbf{p}}d)^{\mathsf T}
#' \boldsymbol{\Sigma}_{\mathbf{p}}(\nabla_{\mathbf{p}}d).}
#' If only individual standard uncertainties are available, unavailable input
#' covariances are approximated as zero, giving
#' \eqn{\sigma_d^2=\sum_i(\partial d/\partial p_i)^2\sigma_{p_i}^2}
#' only after known constraints have been mapped to independent inputs.
#' For \eqn{\mathbf{f}_k=\mathbf{W}_k\mathbf{A}_{j(k)}+\mathbf{t}_k},
#' \deqn{\nabla_{\mathbf{A}_j}d=\sum_{k:j(k)=j}\mathbf{W}_k^{\mathsf T}
#' \nabla_{\mathbf{f}_k}d.}
#' Here \eqn{\mathbf{A}_j} is the parent fractional-coordinate vector,
#' \eqn{\mathbf{W}_k} is the symmetry matrix and \eqn{\mathbf{t}_k} is fixed.
#' With the fractional-coordinate differences defined above, the cell derivative is
#' \deqn{\frac{\partial d}{\partial c}=
#' \frac{c\Delta z^2+a\Delta x\Delta z\cos\beta+b\Delta y\Delta z\cos\alpha}{d}.}
#' Cell angles and their standard uncertainties are converted to radians.
#'
#' Shared generated-site parents are combined before evaluating this variance.
#' Missing standard uncertainties yield `NA` unless a named covariance is supplied.
#' @references Angel et al. (2025), doi:10.1107/S2053273325002682;
#' JCGM 100:2008, doi:10.59161/JCGM100-2008E; IUCr, "Matrices, Mappings
#' and Crystallographic Symmetry" (2002). The package vignette, "Error
#' Propagation", gives the derivative equations.
#' @param bonded_pairs Data.table of bonded atoms with their distances.
#' @param atomic_coordinates Data.table with fractional coordinates and errors.
#' @param unit_cell_metrics Data.table with unit cell parameters and errors.
#' @param expanded_coords Generated sites with `ParentJacobian` provenance.
#' @param independent_coordinates Original asymmetric sites and their uncertainties.
#' @param covariance Optional named covariance matrix for independent inputs.
#' @return The input `bonded_pairs` data.table with a new 'DistanceError' column.
#' @family error propagators
#' @export
#' @examples
#' bp <- data.table::data.table(Atom1 = "Si1", Atom2 = "O1", Distance = sqrt(2),
#'                              DeltaX = 0.1, DeltaY = 0.1, DeltaZ = 0)
#' ac <- data.table::data.table(Label = c("Si1", "O1"),
#'                              x_a = c(0, 0.1), y_b = c(0, 0.1), z_c = c(0, 0),
#'                              x_error = c(0.01, 0.01),
#'                              y_error = c(0.01, 0.01),
#'                              z_error = c(0.01, 0.01))
#' uc <- data.table::data.table(`_cell_length_a` = 10, `_cell_length_a_error` = 0.1,
#'                              `_cell_length_b` = 10, `_cell_length_b_error` = 0.1,
#'                              `_cell_length_c` = 10, `_cell_length_c_error` = 0.1,
#'                              `_cell_angle_alpha` = 90, `_cell_angle_alpha_error` = 0,
#'                              `_cell_angle_beta` = 90, `_cell_angle_beta_error` = 0,
#'                              `_cell_angle_gamma` = 90, `_cell_angle_gamma_error` = 0)
#' propagate_distance_error(bp, ac, uc)
propagate_distance_error <- function(bonded_pairs,
                                     atomic_coordinates,
                                     unit_cell_metrics,
                                     expanded_coords = NULL,
                                     independent_coordinates = atomic_coordinates,
                                     covariance = NULL) {
  if (is.null(bonded_pairs)) return(NULL)
  result <- copy(bonded_pairs)
  result[, DistanceError := NA_real_]
  if (!nrow(result)) return(result)
  geometry <- uncertainty_cell(unit_cell_metrics)
  M <- geometry$M
  G <- t(M) %*% M
  errors <- vapply(seq_len(nrow(result)), function(i) {
    first <- uncertainty_site(result$Atom1[i], atomic_coordinates, expanded_coords)
    second <- uncertainty_site(result$Atom2[i], atomic_coordinates, expanded_coords)
    if (is.null(first) || is.null(second)) return(NA_real_)
    # Stage 1: fractional displacement gives d^2 = Delta_f^T G Delta_f.
    delta <- first$f - second$f
    d <- sqrt(as.numeric(t(delta) %*% G %*% delta))
    if (!is.finite(d) || d <= 1e-12) return(NA_real_)
    # Stage 2: differentiate the metric distance in all six cell inputs.
    # This is the vignette's explicit dd/da, dd/db, dd/dc and angle formulas.
    cart_delta <- as.numeric(M %*% delta)
    gradient <- stats::setNames(vapply(geometry$dM, function(dMs)
      sum(cart_delta * as.numeric(dMs %*% delta)) / d, numeric(1)),
      paste0("cell_", names(geometry$dM)))
    # Stage 3: h_1 = G Delta_f/d, h_2 = -h_1. Shared parents cancel or
    # reinforce through W^T h before calculating their variance.
    h <- as.numeric(G %*% delta) / d
    gradient <- uncertainty_add_parent(gradient, first, h)
    gradient <- uncertainty_add_parent(gradient, second, -h)
    variance <- uncertainty_variance(gradient, unit_cell_metrics,
                                     independent_coordinates, covariance)
    if (is.na(variance) || variance < -1e-10) NA_real_ else sqrt(max(0, variance))
  }, numeric(1))
  result[, DistanceError := errors]
  if (anyNA(errors)) warning("Distance uncertainty unavailable for rows with missing input uncertainty or provenance")
  result
}

#' @title Propagate Angle Error
#' @description Calculates the standard uncertainty for each bond angle.
#' @details
#' Use the notation \eqn{\mathbf{c}_k=\mathbf{M}\mathbf{f}_k},
#' \eqn{\mathbf{a}=\mathbf{c}_2-\mathbf{c}_1},
#' \eqn{\mathbf{b}=\mathbf{c}_3-\mathbf{c}_1} and
#' \eqn{C=\mathbf{a}\cdot\mathbf{b}/(\lVert\mathbf{a}\rVert\lVert\mathbf{b}\rVert)}.
#' The Cartesian gradients are
#' \deqn{\nabla_{\mathbf{c}_2}C=\frac{\mathbf{b}}{\lVert\mathbf{a}\rVert\lVert\mathbf{b}\rVert}
#' -\frac{C\mathbf{a}}{\lVert\mathbf{a}\rVert^2},\qquad
#' \nabla_{\mathbf{c}_3}C=\frac{\mathbf{a}}{\lVert\mathbf{a}\rVert\lVert\mathbf{b}\rVert}
#' -\frac{C\mathbf{b}}{\lVert\mathbf{b}\rVert^2}.}
#' \deqn{\nabla_{\mathbf{c}_1}C=-\nabla_{\mathbf{c}_2}C-\nabla_{\mathbf{c}_3}C,
#' \qquad \mathbf{g}_k=\nabla_{\mathbf{c}_k}\theta=-\frac{\nabla_{\mathbf{c}_k}C}{\sqrt{1-C^2}}.}
#' For each cell parameter \eqn{s\in\{a,b,c,\alpha,\beta,\gamma\}} and
#' parent \eqn{\mathbf{A}_j}, accumulate contributions before forming the variance:
#' \deqn{\frac{\partial\theta}{\partial s}=\sum_{k=1}^3\mathbf{g}_k^{\mathsf T}
#' \frac{\partial\mathbf{M}}{\partial s}\mathbf{f}_k,\qquad
#' \nabla_{\mathbf{A}_j}\theta=\sum_{k:j(k)=j}\mathbf{W}_k^{\mathsf T}
#' \mathbf{M}^{\mathsf T}\mathbf{g}_k.}
#' After mapping any further known constraints to independent inputs \eqn{\mathbf{p}},
#' \deqn{\sigma_\theta^2=(\nabla_{\mathbf{p}}\theta)^{\mathsf T}
#' \boldsymbol{\Sigma}_{\mathbf{p}}(\nabla_{\mathbf{p}}\theta).}
#' A diagonal approximation applies only to independent inputs. Equivalently,
#' retain \eqn{\boldsymbol{\Sigma}_{\mathbf{q}}=\mathbf{J}\boldsymbol{\Sigma}_{\mathbf{p}}\mathbf{J}^{\mathsf T}}
#' for all nine Cartesian coordinates, including off-diagonal entries.
#' With \eqn{D=1-\cos^2\alpha-\cos^2\beta-\cos^2\gamma+2\cos\alpha\cos\beta\cos\gamma},
#' the corrected transformation derivatives include
#' \deqn{\frac{\partial z_c}{\partial\alpha}=
#' \frac{cz_f\sin\alpha(\cos\alpha-\cos\beta\cos\gamma)}{\sqrt{D}\sin\gamma},\qquad
#' \frac{\partial z_c}{\partial\beta}=
#' \frac{cz_f\sin\beta(\cos\beta-\cos\alpha\cos\gamma)}{\sqrt{D}\sin\gamma}.}
#' Differentiate in radians and return \eqn{\sigma_{\theta,\mathrm{deg}}=
#' \sigma_{\theta,\mathrm{rad}}180/\pi}. These first-order gradients require
#' nonzero bond lengths and \eqn{|C|<1}; they are singular at 0 and 180 degrees.
#'
#' Shared generated-site parents are combined before evaluating this variance.
#' Singular angles and missing input uncertainty yield `NA`.
#' @references Angel et al. (2025), doi:10.1107/S2053273325002682;
#' JCGM 100:2008, doi:10.59161/JCGM100-2008E; IUCr recommendations on
#' uncertainty of derived angles. The package vignette, "Error Propagation",
#' gives the angle gradients.
#' @param bond_angles Data.table of calculated bond angles.
#' @param atomic_coordinates Data.table with fractional coordinates and errors.
#' @param expanded_coords Data.table of supercell atom coordinates.
#' @param unit_cell_metrics Data.table with unit cell parameters and errors.
#' @param independent_coordinates Original asymmetric sites and their uncertainties.
#' @param covariance Optional named covariance matrix for independent inputs.
#' @return The input `bond_angles` data.table with a new 'AngleError' column.
#' @family error propagators
#' @export
#' @examples
#' # 1. Create dummy bond angles
#' ba <- data.table::data.table(
#'   CentralAtom = "Si1", Neighbor1 = "O1", Neighbor2 = "O2", Angle = 109.5
#' )
#'
#' # 2. Create dummy atomic coordinates with errors
#' ac <- data.table::data.table(
#'   Label = c("Si1", "O1", "O2"),
#'   x_a = c(0, 0.1, 0), y_b = c(0, 0, 0.1), z_c = c(0, 0, 0),
#'   x_error = c(0.01, 0.01, 0.01),
#'   y_error = c(0.01, 0.01, 0.01),
#'   z_error = c(0.01, 0.01, 0.01)
#' )
#'
#' # 3. Create dummy expanded coordinates
#' ec <- data.table::data.table(
#'   Label = c("O1", "O2"),
#'   x_a = c(0.1, 0), y_b = c(0, 0.1), z_c = c(0, 0)
#' )
#'
#' # 4. Create dummy unit cell metrics
#' uc <- data.table::data.table(
#'   `_cell_length_a` = 10, `_cell_length_a_error` = 0.1,
#'   `_cell_length_b` = 10, `_cell_length_b_error` = 0.1,
#'   `_cell_length_c` = 10, `_cell_length_c_error` = 0.1,
#'   `_cell_angle_alpha` = 90, `_cell_angle_alpha_error` = 0,
#'   `_cell_angle_beta` = 90, `_cell_angle_beta_error` = 0,
#'   `_cell_angle_gamma` = 90, `_cell_angle_gamma_error` = 0
#' )
#'
#' # 5. Run the error propagation
#' propagate_angle_error(ba, ac, ec, uc)
propagate_angle_error <- function(bond_angles,
                                  atomic_coordinates,
                                  expanded_coords,
                                  unit_cell_metrics,
                                  independent_coordinates = atomic_coordinates,
                                  covariance = NULL) {
  if (is.null(bond_angles)) return(NULL)
  result <- copy(bond_angles)
  result[, AngleError := NA_real_]
  if (!nrow(result)) return(result)
  geometry <- uncertainty_cell(unit_cell_metrics)
  M <- geometry$M
  errors <- vapply(seq_len(nrow(result)), function(i) {
    labels <- c(result$CentralAtom[i], result$Neighbor1[i], result$Neighbor2[i])
    sites <- lapply(labels, uncertainty_site, atomic_coordinates = atomic_coordinates,
                    expanded_coords = expanded_coords)
    if (any(vapply(sites, is.null, logical(1)))) return(NA_real_)
    # Stage 1: convert all three generated positions c_k = M f_k.
    positions <- lapply(sites, function(site) as.numeric(M %*% site$f))
    a_vec <- positions[[2]] - positions[[1]]
    b_vec <- positions[[3]] - positions[[1]]
    na <- sqrt(sum(a_vec^2)); nb <- sqrt(sum(b_vec^2))
    if (!is.finite(na * nb) || na <= 1e-12 || nb <= 1e-12) return(NA_real_)
    C <- sum(a_vec * b_vec) / (na * nb)
    if (!is.finite(C) || abs(C) >= 1 - 1e-12) return(NA_real_)
    # Stage 2: differentiate C = (a.b)/(|a||b|), then theta = acos(C).
    # The central gradient is minus the sum of the two neighbor gradients.
    dC2 <- b_vec / (na * nb) - C * a_vec / na^2
    dC3 <- a_vec / (na * nb) - C * b_vec / nb^2
    g2 <- -dC2 / sqrt(1 - C^2)
    g3 <- -dC3 / sqrt(1 - C^2)
    g <- list(-g2 - g3, g2, g3)
    # Stage 3: each cell input moves all three positions through dM/ds.
    gradient <- stats::setNames(vapply(geometry$dM, function(dMs)
      sum(vapply(seq_along(sites), function(k)
        sum(g[[k]] * as.numeric(dMs %*% sites[[k]]$f)), numeric(1))), numeric(1)),
      paste0("cell_", names(geometry$dM)))
    # Stage 4: M^T g_k returns each Cartesian gradient to fractional space;
    # the stored W mappings combine contributions from repeated parents.
    for (k in seq_along(sites))
      gradient <- uncertainty_add_parent(gradient, sites[[k]],
                                         as.numeric(t(M) %*% g[[k]]))
    variance <- uncertainty_variance(gradient, unit_cell_metrics,
                                     independent_coordinates, covariance)
    # Stage 5: the gradient is in radians; report the angle SU in degrees.
    if (is.na(variance) || variance < -1e-10) NA_real_ else sqrt(max(0, variance)) * 180 / pi
  }, numeric(1))
  result[, AngleError := errors]
  if (anyNA(errors)) warning("Angle uncertainty unavailable for singular geometry, missing input uncertainty, or provenance")
  result
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
#' # 1. Create a sample data.table of bond angles
#' sample_angles <- data.table::data.table(
#'   CentralAtom = c("C1", "C2", "Si1", "Cr1", "O1", "O2", "C"),
#'   Neighbor1 = c("O1", "O2", "O1", "N1", "C1", "C2", "H1"),
#'   Neighbor2 = c("H1", "H2", "O2", "N2", "H3", "H4", "H2"),
#'   Angle = c(109.5, 109.5, 120.0, 90.0, 104.5, 104.5, 120)
#' )
#'
#' # 2. In an interactive R session, the function would prompt the user.
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
  message("Available base atom symbols found: ",
          paste(base_symbols, collapse = ", "))

  user_input <- readline(prompt = "Please enter the chemical symbols you want to filter for, separated by commas (e.g., C,Si,O):\n")

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

#' @title Filter Data by Wyckoff Symbol
#' @description Filters a data table (e.g., bonds or angles) to include only
#'   entries where a specified atom occupies one of the given Wyckoff sites.
#'
#' @details
#' This function is designed to facilitate analysis based on crystallographic
#' site symmetry. It is particularly useful for analyzing complex structures like
#' clathrates where different atoms perform distinct structural roles based
#' on their site symmetry.
#'
#' The function works by:
#' 1.  Creating a full Wyckoff label (e.g., "4c", "24k") by combining the
#'     `WyckoffMultiplicity` and `WyckoffSymbol` columns from the
#'     `atomic_coordinates` table.
#' 2.  Identifying the parent atom for each entry in the `data_table`.
#' 3.  Merging this Wyckoff information into the results table.
#' 4.  Filtering to keep only rows where the atom's full Wyckoff label matches
#'     one of the symbols provided in the `wyckoff_symbols` vector.
#'
#' @param data_table A `data.table` object, such as one produced by
#'   `minimum_distance` or
#'   `calculate_angles`.
#' @param atomic_coordinates A `data.table` from `extract_atomic_coordinates`
#'   containing `WyckoffMultiplicity` and `WyckoffSymbol` columns.
#' @param atom_col A character string specifying the column in `data_table`
#'   that contains the atom labels to filter by (e.g., "Atom1", "CentralAtom").
#' @param wyckoff_symbols A character vector of the full Wyckoff symbols to
#'   keep (e.g., `c("4c")` or `c("6c", "16i", "24k")`).
#' @return A `data.table` filtered to include only rows where the specified
#'   atom occupies one of the desired Wyckoff sites.
#' @family property calculators
#' @export
#' @examples
#' cif_file <- system.file("extdata", "1590946.cif", package = "crystract")
#' if (file.exists(cif_file)) {
#'   old_threads <- data.table::setDTthreads(2)
#'   # 1. Perform a standard analysis to get bond and coordinate tables
#'   cif_content <- read_cif_files(cif_file)[[1]]
#'   atoms <- extract_atomic_coordinates(cif_content)
#'   metrics <- extract_unit_cell_metrics(cif_content)
#'   sym_ops <- extract_symmetry_operations(cif_content)
#'   full_cell <- apply_symmetry_operations(atoms, sym_ops, metrics)
#'   super_cell <- expand_transformed_coords(full_cell)
#'   dists <- calculate_distances(atoms, super_cell, metrics)
#'   bonds <- minimum_distance(dists)
#'
#'   # 2. Mock Wyckoff sites since 1590946 doesn't explicitly declare them
#'   atoms[, WyckoffSymbol := "c"]
#'   atoms[, WyckoffMultiplicity := 4]
#'
#'   print("Original atomic coordinates showing Wyckoff sites:")
#'   print(atoms[, .(Label, WyckoffSymbol, WyckoffMultiplicity)])
#'
#'   filtered_bonds <- filter_by_wyckoff_symbol(
#'     data_table = bonds,
#'     atomic_coordinates = atoms,
#'     atom_col = "Atom1",
#'     wyckoff_symbols = "4c"
#'   )
#'
#'   cat("\nNumber of bonds in original table:", nrow(bonds), "\n")
#'   cat("Number of bonds after filtering for '4c' site:", nrow(filtered_bonds), "\n")
#'   data.table::setDTthreads(old_threads)
#' }
filter_by_wyckoff_symbol <- function(data_table,
                                     atomic_coordinates,
                                     atom_col,
                                     wyckoff_symbols) {
  # --- Input Validation ---
  if (!atom_col %in% names(data_table)) {
    stop(paste("Column '", atom_col, "' not found in the data table.", sep = ""))
  }
  required_cols <- c("WyckoffSymbol", "WyckoffMultiplicity", "Label")
  if (!all(required_cols %in% names(atomic_coordinates))) {
    stop(paste(
      "`atomic_coordinates` must contain the columns:",
      paste(required_cols, collapse = ", ")
    ))
  }
  if (!is.character(wyckoff_symbols)) {
    stop("`wyckoff_symbols` must be a character vector.")
  }

  # --- Data Preparation ---
  work_dt <- copy(data_table)

  # Extract parent labels to merge with atomic_coordinates
  work_dt[, ParentLabel := sub("_.*", "", get(atom_col))]

  # Prepare the Wyckoff info lookup table
  wyckoff_info <- atomic_coordinates[, .(Label, FullWyckoff = paste0(WyckoffMultiplicity, WyckoffSymbol))]

  # --- Merging and Filtering ---
  merged_dt <- merge(
    work_dt,
    wyckoff_info,
    by.x = "ParentLabel",
    by.y = "Label",
    all.x = TRUE
  )

  # Filter based on the provided list of full Wyckoff symbols
  filtered_dt <- merged_dt[FullWyckoff %in% wyckoff_symbols, ]

  # --- Clean and Return ---
  # Remove the temporary columns used for the merge
  filtered_dt[, `:=`(ParentLabel = NULL, FullWyckoff = NULL)]

  return(filtered_dt)
}

#' @title Filter Distances by Element Symbols
#' @description Removes any distance pair where at least one of the atoms
#'   corresponds to a specified list of element symbols.
#' @param distances A `data.table` of interatomic distances.
#' @param atomic_coordinates A `data.table` of atomic coordinates.
#' @param elements_to_exclude A character vector of element symbols to exclude.
#' @return A `data.table` of distances with the specified elements removed.
#' @family property calculators
#' @export
#' @examples
#' dists <- data.table::data.table(Atom1 = "Si1_1_0_0", Atom2 = "O1_1_0_0", Distance = 1.6)
#' ac <- data.table::data.table(Label = c("Si1", "O1"))
#' filter_by_elements(dists, ac, elements_to_exclude = "O")
filter_by_elements <- function(distances,
                               atomic_coordinates,
                               elements_to_exclude) {
  if (is.null(elements_to_exclude) ||
      length(elements_to_exclude) == 0) {
    return(distances)
  }
  element_info <- atomic_coordinates[, .(ParentLabel = Label,
                                         Element = sub("[0-9].*", "", Label))]
  dist_with_elements <- copy(distances)
  dist_with_elements[, Parent1 := sub("_.*", "", Atom1)]
  dist_with_elements[, Parent2 := sub("_.*", "", Atom2)]
  merged <- merge(
    dist_with_elements,
    element_info,
    by.x = "Parent1",
    by.y = "ParentLabel",
    all.x = TRUE
  )
  setnames(merged, "Element", "Element1")
  merged <- merge(
    merged,
    element_info,
    by.x = "Parent2",
    by.y = "ParentLabel",
    all.x = TRUE
  )
  setnames(merged, "Element", "Element2")
  filtered <- merged[!(Element1 %in% elements_to_exclude) &
                       !(Element2 %in% elements_to_exclude)]
  return(filtered[, names(distances), with = FALSE])
}

#' @title Filter Ghost Distances Using Atomic Radii
#' @description Cleans a distance table by removing physically implausible distances.
#' It uses a table of atomic radii to establish a plausible bond length
#' range for each atom pair. Any calculated distance falling outside this range
#' (defined by a `margin`) is considered a "ghost" distance and is removed. This
#' is particularly useful for cleaning data from disordered crystal structures.
#' @param distances A `data.table` of interatomic distances, typically from
#' `calculate_distances`. Must contain `Atom1`, `Atom2`, and `Distance`.
#' @param atomic_coordinates A `data.table` of asymmetric atoms from
#' `extract_atomic_coordinates`. Used to link atom labels to element types.
#' @param margin A numeric value (default 0.1) specifying the tolerance. A
#' distance `d` between atoms with radii `r1` and `r2` is kept if
#' `(r1+r2)*(1-margin) <= d <= (r1+r2)*(1+margin)`.
#' @param radii_type A character string specifying the type of radius to use for
#' the calculation (e.g., `"covalent"`, `"ionic"`). This value must correspond
#' to an entry in the `Type` column of the active radii table. Defaults to `"covalent"`.
#' The radii table can be customized for the session using `set_radii_data()`.
#' @return A list containing two `data.table`s:
#' \item{kept}{The distances considered physically plausible.}
#' \item{removed}{The "ghost" distances that were filtered out, with columns
#' explaining the reason for removal.}
#' @family post-processing
#' @export
#' @examples
#' # Create minimal dummy data for demonstration
#' distances <- data.table::data.table(
#'   Atom1 = c("Si1_1_0_0", "O1_1_0_0"),
#'   Atom2 = c("O1_1_0_0", "Si1_1_0_0"),
#'   Distance = c(1.6, 0.5) # 0.5 is implausibly short
#' )
#'
#' atomic_coords <- data.table::data.table(
#'   Label = c("Si1", "O1")
#' )
#'
#' # Run the filter
#' result <- filter_ghost_distances(distances, atomic_coords, margin = 0.1)
#'
#' print(result$kept)
#' print(result$removed)
filter_ghost_distances <- function(distances,
                                   atomic_coordinates,
                                   margin = 0.1,
                                   radii_type = "covalent") {
  # --- 1. Augment data for calculation ---
  # Create a lookup for element type from the parent atom label
  element_info <- atomic_coordinates[, .(ParentLabel = Label,
                                         Element = sub("[0-9].*", "", Label))]

  # Use the new helper to get the active radii table (user-defined or default)
  all_radii_data <- get_radii_data()

  # Filter the radii table for the specific type requested by the user
  radii_lookup <- all_radii_data[Type == radii_type]

  if (nrow(radii_lookup) == 0) {
    stop(
      paste0(
        "Radii type '",
        radii_type,
        "' not found in the active radii table. ",
        "Available types are: ",
        paste(unique(all_radii_data$Type), collapse = ", ")
      )
    )
  }

  # Keep only Symbol and Radius for the merge
  radii_lookup <- radii_lookup[, .(Symbol, Radius)]

  merged <- copy(distances)
  # Create Parent labels in the distances table to match element_info
  merged[, `:=`(Parent1 = sub("_.*", "", Atom1),
                Parent2 = sub("_.*", "", Atom2))]
  # Merge to get Element and Radius for Atom1
  merged <- merge(
    merged,
    element_info,
    by.x = "Parent1",
    by.y = "ParentLabel",
    all.x = TRUE
  )
  setnames(merged, "Element", "Element1")
  merged <- merge(
    merged,
    radii_lookup,
    by.x = "Element1",
    by.y = "Symbol",
    all.x = TRUE
  )
  setnames(merged, "Radius", "Radius1")
  # Merge to get Element and Radius for Atom2
  merged <- merge(
    merged,
    element_info,
    by.x = "Parent2",
    by.y = "ParentLabel",
    all.x = TRUE
  )
  setnames(merged, "Element", "Element2")
  merged <- merge(
    merged,
    radii_lookup,
    by.x = "Element2",
    by.y = "Symbol",
    all.x = TRUE
  )
  setnames(merged, "Radius", "Radius2")
  # Handle cases where radii are not found (e.g., vacancies) by setting them to 0
  merged[is.na(Radius1), Radius1 := 0][is.na(Radius2), Radius2 := 0]

  # --- 2. Calculate the plausible distance range ---
  merged[, expected_dist := Radius1 + Radius2]
  merged[, lower_bound := expected_dist * (1 - margin)]
  merged[, upper_bound := expected_dist * (1 + margin)]

  # --- 3. Apply the filter ---
  kept_mask <- merged$Distance >= merged$lower_bound &
    merged$Distance <= merged$upper_bound &
    merged$expected_dist > 0
  kept_dt <- merged[kept_mask, names(distances), with = FALSE]
  # Create a detailed table of removed rows for quality control
  removed_dt <- merged[!kept_mask, .(
    Atom1,
    Atom2,
    Distance,
    expected_dist,
    lower_bound,
    upper_bound,
    Reason = ifelse(
      Distance < lower_bound,
      "Distance is TOO SHORT",
      "Distance is TOO LONG"
    )
  )]
  return(list(kept = kept_dt, removed = removed_dt))
}

#' @title Calculate Weighted Average Network Bond Distance
#' @description Computes a single, representative bond length for a specified atomic
#'   network. This function precisely implements the validated logic that accounts
#'   for site multiplicity and occupancy of the central atoms.
#' @details
#' The occupancy-corrected average interatomic distance is
#' \deqn{\bar d=\frac{\sum_j m_j\,\mathrm{occ}_j\sum_{i=1}^{n_j}d_{ij}}
#' {\sum_j m_j\,\mathrm{occ}_j\,n_j}.}
#' The index \eqn{j} runs over symmetry-inequivalent central atoms in the selected
#' network; \eqn{m_j} is the Wyckoff multiplicity, \eqn{\mathrm{occ}_j} the central
#' site occupancy, \eqn{n_j} the retained distance count and \eqn{d_{ij}} a distance
#' to a neighbour. The inner sum is not weighted by neighbour occupancy.
#' @param distances A `data.table` of interatomic distances filtered to include
#'   **only bonded pairs** (e.g., from `minimum_distance`).
#' @param atomic_coordinates A `data.table` of asymmetric atoms from
#'   `extract_atomic_coordinates`.
#' @param wyckoff_symbols A character vector of Wyckoff symbols defining the
#'   atomic network (e.g., `c("6c", "16i", "24k")`). Must be the full symbol.
#' @return A single numeric value representing the weighted average bond distance.
#' @family post-processing
#' @export
#' @examples
#' dists <- data.table::data.table(Atom1 = "Si1", Atom2 = "O1", Distance = 1.6)
#' ac <- data.table::data.table(Label = c("Si1", "O1"),
#'                              WyckoffMultiplicity = c(4, 4),
#'                              WyckoffSymbol = c("c", "c"),
#'                              Occupancy = c(1, 1))
#' calculate_weighted_average_network_distance(dists, ac, wyckoff_symbols = "4c")
calculate_weighted_average_network_distance <- function(distances,
                                                        atomic_coordinates,
                                                        wyckoff_symbols) {
  # --- 1. Prepare atomic info and identify network atoms ---
  atom_info <- copy(atomic_coordinates)
  atom_info[, FullWyckoff := paste0(WyckoffMultiplicity, WyckoffSymbol)]
  network_atom_labels <- atom_info[FullWyckoff %in% wyckoff_symbols, Label]

  if (length(network_atom_labels) == 0) {
    return(NA_real_)
  }

  # --- 2. Filter bonds originating from network atoms ---
  network_distances <- distances[sub("_.*", "", Atom1) %in% network_atom_labels]

  if (nrow(network_distances) == 0) {
    return(NA_real_)
  }

  # --- 3. Sum lengths and count neighbors for each central site. Neighbor
  # occupancy does not enter the network-distance formula. ---
  atom_dist_summary <- network_distances[, .(
    sum_d = sum(Distance),
    # Simple, unweighted sum of bond lengths
    n = .N               # Simple, unweighted coordination number
  ), by = Atom1]

  atom_dist_summary[, ParentLabel := sub("_.*", "", Atom1)]

  # --- 4. Merge summary with main atomic info ---
  unique_summary <- unique(atom_dist_summary, by = "ParentLabel")
  merged_data <- merge(
    atom_info[Label %in% network_atom_labels],
    unique_summary[, .(ParentLabel, sum_d, n)],
    by.x = "Label",
    by.y = "ParentLabel",
    all.x = TRUE
  )

  merged_data <- merged_data[!is.na(sum_d)]
  if (nrow(merged_data) == 0)
    return(NA_real_)

  # --- 5. Calculate the final weighted average ---
  numerator <- sum(
    merged_data$WyckoffMultiplicity * merged_data$Occupancy * merged_data$sum_d,
    na.rm = TRUE
  )
  denominator <- sum(merged_data$WyckoffMultiplicity * merged_data$Occupancy * merged_data$n,
                     na.rm = TRUE)

  if (is.na(denominator) || denominator == 0) {
    return(NA_real_)
  }

  final_result <- numerator / denominator
  return(final_result)
}
