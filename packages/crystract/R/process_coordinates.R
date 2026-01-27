#' @title Merge Close Atoms (Pymatgen Style)
#' @description Merges atomic sites that are within a specified distance tolerance,
#' accounting for Periodic Boundary Conditions (PBC). This mimics pymatgen's
#' `Structure.merge_sites` logic using hierarchical clustering.
#'
#' @param atomic_coordinates A data.table with columns `Label`, `x_a`, `y_b`, `z_c`.
#' @param unit_cell_metrics A data.table containing cell lengths and angles.
#' @param tol Numeric. Distance tolerance in Angstroms. Sites closer than this will be merged.
#' @return A deduplicated data.table with coordinates averaged.
#' @export
merge_sites_pbc <- function(atomic_coordinates,
                            unit_cell_metrics,
                            tol = 0.1) {
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

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      dx <- coords[i, 1] - coords[j, 1]
      dy <- coords[i, 2] - coords[j, 2]
      dz <- coords[i, 3] - coords[j, 3]

      # Apply Minimum Image Convention
      dx <- dx - round(dx)
      dy <- dy - round(dy)
      dz <- dz - round(dz)

      d2 <- g11 * dx^2 + g22 * dy^2 + g33 * dz^2 +
        2 * g12 * dx * dy + 2 * g13 * dx * dz + 2 * g23 * dy * dz

      d <- sqrt(d2)
      dist_mat[i, j] <- d
      dist_mat[j, i] <- d
    }
  }

  # 4. Hierarchical Clustering
  dist_obj <- as.dist(dist_mat)
  # Single linkage groups chain-linked neighbors (standard for merging)
  hc <- hclust(dist_obj, method = "single")
  clusters <- cutree(hc, h = tol)

  # 5. Merge Sites (Take first entry per cluster to preserve labels)
  merged_dt <- atomic_coordinates[, .SD[1], by = clusters]
  merged_dt[, clusters := NULL]

  return(merged_dt)
}

#' @title Apply Symmetry Operations to Generate a Full Unit Cell
#' @description Generates all symmetry-equivalent atomic positions within the unit cell.
#' @param atomic_coordinates A `data.table` of asymmetric atoms.
#' @param symmetry_operations A `data.table` of operations.
#' @param unit_cell_metrics A `data.table` of cell parameters (needed for distance merging).
#' @param tolerance Numeric. The distance tolerance in Angstroms for merging close atoms.
#' @return A `data.table` containing all unique atomic positions in the unit cell.
#' @family coordinate processors
#' @export
apply_symmetry_operations <- function(atomic_coordinates,
                                      symmetry_operations,
                                      unit_cell_metrics,
                                      tolerance = 0.1) {
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
  # Pymatgen-style: Clean up using distance tolerance
  transformed_coords <- merge_sites_pbc(transformed_coords, unit_cell_metrics, tol = tolerance)

  return(transformed_coords)
}


#' @title Expand Coordinates into a Supercell
#' @description Takes a set of atomic coordinates within a single unit cell and
#'   replicates them into a 3x3x3 supercell grid.
#' @param transformed_coords A `data.table` of atom positions within one unit cell.
#' @return A `data.table` containing all atomic positions in the expanded supercell.
#' @family coordinate processors
#' @export
expand_transformed_coords <- function(transformed_coords) {
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

  cell_indices <- as.data.table(expand.grid(
    x = -1:1,
    y = -1:1,
    z = -1:1
  ))

  expanded_coords <- rbindlist(lapply(1:nrow(cell_indices), function(i) {
    cell_shift <- cell_indices[i, ]
    dt <- copy(transformed_coords)
    dt[, `:=`(Tx = cell_shift$x,
              Ty = cell_shift$y,
              Tz = cell_shift$z)]
    dt[, Label := paste(Label, Tx, Ty, Tz, sep = "_")]
    dt[, `:=`(x_a = x_a + Tx,
              y_b = y_b + Ty,
              z_c = z_c + Tz)]
    return(dt)
  }))

  # Use simple unique here; proximity issues should have been resolved
  # in the unit cell via apply_symmetry_operations -> merge_sites_pbc
  expanded_coords <- unique(expanded_coords, by = c("x_a", "y_b", "z_c"))
  return(expanded_coords)
}
