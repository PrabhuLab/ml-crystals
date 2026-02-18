#' @title Merge Close Atoms (Pymatgen Style)
#' @description Merges atomic sites that are within a specified distance tolerance,
#' accounting for Periodic Boundary Conditions (PBC). This mimics pymatgen's
#' `Structure.merge_sites` logic: it performs hierarchical clustering and then
#' averages the coordinates of the merged sites, taking care of PBC wrapping.
#'
#' @param atomic_coordinates A data.table with columns `Label`, `x_a`, `y_b`, `z_c`.
#' @param unit_cell_metrics A data.table containing cell lengths and angles.
#' @param tol Numeric. Distance tolerance in Angstroms. Sites closer than this will be merged.
#' Default is 1e-4.
#' @return A deduplicated data.table with coordinates averaged.
#' @export
merge_sites_pbc <- function(atomic_coordinates,
                            unit_cell_metrics,
                            tol = 1e-4) {
  if (nrow(atomic_coordinates) < 2)
    return(atomic_coordinates)

  # 1. Extract Cell Parameters
  a <- as.numeric(unit_cell_metrics$`_cell_length_a`)
  b <- as.numeric(unit_cell_metrics$`_cell_length_b`)
  c <- as.numeric(unit_cell_metrics$`_cell_length_c`)
  alpha <- as.numeric(unit_cell_metrics$`_cell_angle_alpha`) * pi / 180
  beta  <- as.numeric(unit_cell_metrics$`_cell_angle_beta`) * pi / 180
  gamma <- as.numeric(unit_cell_metrics$`_cell_angle_gamma`) * pi / 180

  # 2. Metric Tensor components
  g11 <- a^2
  g22 <- b^2
  g33 <- c^2
  g12 <- a * b * cos(gamma)
  g13 <- a * c * cos(beta)
  g23 <- b * c * cos(alpha)

  # 3. Calculate Distance Matrix with PBC manually
  coords <- as.matrix(atomic_coordinates[, .(x_a, y_b, z_c)])
  n <- nrow(coords)
  dist_mat <- matrix(0, nrow = n, ncol = n)

  # Store displacements for use in averaging later
  # dx_mat[i, j] = dx from i to j
  dx_mat <- matrix(0, nrow = n, ncol = n)
  dy_mat <- matrix(0, nrow = n, ncol = n)
  dz_mat <- matrix(0, nrow = n, ncol = n)

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      dx <- coords[i, 1] - coords[j, 1]
      dy <- coords[i, 2] - coords[j, 2]
      dz <- coords[i, 3] - coords[j, 3]

      # Apply Minimum Image Convention
      rx <- round(dx)
      ry <- round(dy)
      rz <- round(dz)

      dx <- dx - rx
      dy <- dy - ry
      dz <- dz - rz

      d2 <- g11 * dx^2 + g22 * dy^2 + g33 * dz^2 +
        2 * g12 * dx * dy + 2 * g13 * dx * dz + 2 * g23 * dy * dz

      d <- sqrt(d2)
      dist_mat[i, j] <- d
      dist_mat[j, i] <- d

      # Store wrapped displacements
      # If dist is small, this helps us average correctly
      if (d < tol + 1e-6) {
        dx_mat[i, j] <- dx
        dy_mat[i, j] <- dy
        dz_mat[i, j] <- dz
        dx_mat[j, i] <- -dx
        dy_mat[j, i] <- -dy
        dz_mat[j, i] <- -dz
      }
    }
  }

  # 4. Hierarchical Clustering
  dist_obj <- as.dist(dist_mat)
  # Single linkage groups chain-linked neighbors (standard for merging)
  hc <- hclust(dist_obj, method = "single")
  clusters <- cutree(hc, h = tol)

  # 5. Merge Sites with Averaging (matching pymatgen structure.py logic)
  # We iterate through unique clusters. If a cluster has size > 1, we average.

  # Add cluster ID to data
  work_dt <- copy(atomic_coordinates)
  work_dt[, ClusterID := clusters]

  # Identify clusters with duplicates
  dup_clusters <- work_dt[, .N, by = ClusterID][N > 1]$ClusterID

  if (length(dup_clusters) > 0) {
    for (cid in dup_clusters) {
      # Get indices of atoms in this cluster
      indices <- which(clusters == cid)

      # Base coordinate is the first atom
      base_idx <- indices[1]
      avg_x <- coords[base_idx, 1]
      avg_y <- coords[base_idx, 2]
      avg_z <- coords[base_idx, 3]

      # For subsequent atoms, add the Displacement from base to current
      # This effectively unwraps them relative to the base atom
      # Matches pymatgen: coords += ((offset - np.round(offset)) / (site_idx + 2))
      # Here we just sum them up then divide by N

      for (k in 2:length(indices)) {
        curr_idx <- indices[k]
        # dx_mat[base, curr] represents (base - curr)
        # We want curr relative to base, so we subtract (base-curr) from base
        # No, wait.
        # base + (curr - base)_wrapped
        # (curr - base)_wrapped is -dx_mat[base, curr]

        avg_x <- avg_x + (coords[base_idx, 1] - dx_mat[base_idx, curr_idx])
        avg_y <- avg_y + (coords[base_idx, 2] - dy_mat[base_idx, curr_idx])
        avg_z <- avg_z + (coords[base_idx, 3] - dz_mat[base_idx, curr_idx])
      }

      # Final average
      n_clus <- length(indices)
      work_dt[ClusterID == cid, `:=`(
        x_a = avg_x / n_clus,
        y_b = avg_y / n_clus,
        z_c = avg_z / n_clus
      )]
    }
  }

  # Deduplicate: take first entry of each cluster, but with the averaged coordinates
  merged_dt <- unique(work_dt, by = "ClusterID")

  # Wrap back to unit cell [0, 1)
  merged_dt[, `:=`(x_a = x_a %% 1,
                   y_b = y_b %% 1,
                   z_c = z_c %% 1)]
  merged_dt[, ClusterID := NULL]

  return(merged_dt)
}

#' @title Apply Symmetry Operations to Generate a Full Unit Cell
#' @description Generates all symmetry-equivalent atomic positions within the unit cell.
#' @param atomic_coordinates A `data.table` of asymmetric atoms.
#' @param symmetry_operations A `data.table` of operations.
#' @param unit_cell_metrics A `data.table` of cell parameters (needed for distance merging).
#' @param tolerance Numeric. The distance tolerance in Angstroms for merging close atoms.
#' Default is 1e-4.
#' @return A `data.table` containing all unique atomic positions in the unit cell.
#' @family coordinate processors
#' @export
apply_symmetry_operations <- function(atomic_coordinates,
                                      symmetry_operations,
                                      unit_cell_metrics,
                                      tolerance = 1e-4) {
  if (is.null(atomic_coordinates) || nrow(atomic_coordinates) == 0) {
    return(
      data.table(
        Label = character(),
        SymOp = integer(),
        x_a = numeric(),
        y_b = numeric(),
        z_c = numeric()
      )
    )
  }

  safe_eval <- function(op_str, x, y, z) {
    expr_str <- gsub("x", sprintf("(%.16f)", x), op_str, fixed = TRUE)
    expr_str <- gsub("y", sprintf("(%.16f)", y), expr_str, fixed = TRUE)
    expr_str <- gsub("z", sprintf("(%.16f)", z), expr_str, fixed = TRUE)
    tryCatch(
      eval(parse(text = expr_str)),
      error = function(e)
        stop(paste("SymOp parse error:", op_str))
    )
  }

  n_atoms <- nrow(atomic_coordinates)
  n_ops <- nrow(symmetry_operations)

  expanded_list <- lapply(seq_len(n_atoms), function(i) {
    row <- atomic_coordinates[i, ]
    lbl <- row$Label
    val_x <- as.numeric(row$x_a)
    val_y <- as.numeric(row$y_b)
    val_z <- as.numeric(row$z_c)

    new_coords <- lapply(seq_len(n_ops), function(j) {
      op_x <- symmetry_operations$x[j]
      op_y <- symmetry_operations$y[j]
      op_z <- symmetry_operations$z[j]

      data.table(
        Label = paste(lbl, j, sep = "_"),
        SymOp = j,
        x_a = safe_eval(op_x, val_x, val_y, val_z) %% 1,
        y_b = safe_eval(op_y, val_x, val_y, val_z) %% 1,
        z_c = safe_eval(op_z, val_x, val_y, val_z) %% 1
      )
    })
    rbindlist(new_coords)
  })

  transformed_coords <- rbindlist(expanded_list)
  if (nrow(transformed_coords) == 0)
    return(transformed_coords)

  # --- DEDUPLICATION LOGIC ---
  # Pymatgen-style: Clean up using distance tolerance and averaging
  transformed_coords <- merge_sites_pbc(transformed_coords, unit_cell_metrics, tol = tolerance)

  return(transformed_coords)
}

#' @title Calculate Supercell Expansion Factors
#' @description Computes the number of unit cell repetitions required in each
#'   direction (a, b, c) to encompass a sphere of a given radius.
#' @param unit_cell_metrics A `data.table` with cell parameters.
#' @param radius Numeric. The cutoff radius (e.g., 13.0 for Voronoi).
#' @return A named numeric vector `c(a=..., b=..., c=...)`.
#' @family coordinate processors
#' @export
calculate_expansion_factors <- function(unit_cell_metrics, radius) {
  a <- as.numeric(unit_cell_metrics$`_cell_length_a`)
  b <- as.numeric(unit_cell_metrics$`_cell_length_b`)
  c <- as.numeric(unit_cell_metrics$`_cell_length_c`)
  alpha <- as.numeric(unit_cell_metrics$`_cell_angle_alpha`) * pi / 180
  beta <- as.numeric(unit_cell_metrics$`_cell_angle_beta`) * pi / 180
  gamma <- as.numeric(unit_cell_metrics$`_cell_angle_gamma`) * pi / 180

  # Volume calculation
  cos_a <- cos(alpha)
  cos_b <- cos(beta)
  cos_g <- cos(gamma)

  term <- 1 - cos_a^2 - cos_b^2 - cos_g^2 + 2 * cos_a * cos_b * cos_g
  V <- a * b * c * sqrt(max(0, term))

  # Perpendicular widths (height of the parallelepiped along each vector)
  h_a <- V / (b * c * sin(alpha))
  h_b <- V / (a * c * sin(beta))
  h_c <- V / (a * b * sin(gamma))

  # We need radius to be covered by n * height.
  # Note: Pymatgen adds a small buffer (0.15). We match that.
  r_eff <- radius + 0.15

  n_a <- ceiling(r_eff / h_a)
  n_b <- ceiling(r_eff / h_b)
  n_c <- ceiling(r_eff / h_c)

  # Ensure at least 1 (though 3x3x3 means index -1 to 1, so minimum factor is 1)
  # If factors are c(1,1,1), ranges are -1:1.
  return(c(
    a = max(1, n_a),
    b = max(1, n_b),
    c = max(1, n_c)
  ))
}

#' @title Expand Coordinates into a Supercell
#' @description Takes a set of atomic coordinates within a single unit cell and
#'   replicates them into a supercell grid defined by expansion factors.
#' @param transformed_coords A `data.table` of atom positions within one unit cell.
#' @param expansion_factors A numeric vector `c(a, b, c)` specifying the number of
#'   images in each direction (+/-). Default `c(1,1,1)` creates a 3x3x3 grid.
#' @return A `data.table` containing all atomic positions in the expanded supercell.
#' @family coordinate processors
#' @export
expand_transformed_coords <- function(transformed_coords,
                                      expansion_factors = c(1, 1, 1)) {
  if (is.null(transformed_coords) ||
      nrow(transformed_coords) == 0 ||
      !"x_a" %in% names(transformed_coords)) {
    return(
      data.table(
        Label = character(),
        SymOp = integer(),
        Tx = integer(),
        Ty = integer(),
        Tz = integer(),
        x_a = numeric(),
        y_b = numeric(),
        z_c = numeric()
      )
    )
  }

  na <- expansion_factors[1]
  nb <- expansion_factors[2]
  nc <- expansion_factors[3]

  # Generate indices from -n to +n
  cell_indices <- as.data.table(expand.grid(
    x = -na:na,
    y = -nb:nb,
    z = -nc:nc
  ))

  expanded_coords <- rbindlist(lapply(1:nrow(cell_indices), function(i) {
    cell_shift <- cell_indices[i, ]
    dt <- copy(transformed_coords)
    dt[, `:=`(Tx = cell_shift$x,
              Ty = cell_shift$y,
              Tz = cell_shift$z)]
    # Make labels unique for the supercell images
    dt[, Label := paste(Label, Tx, Ty, Tz, sep = "_")]
    dt[, `:=`(x_a = x_a + Tx,
              y_b = y_b + Ty,
              z_c = z_c + Tz)]
    return(dt)
  }))

  # Simple unique check
  expanded_coords <- unique(expanded_coords, by = c("x_a", "y_b", "z_c"))
  return(expanded_coords)
}
