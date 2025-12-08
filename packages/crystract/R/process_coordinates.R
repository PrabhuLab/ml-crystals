#' @title Apply Symmetry Operations to Generate a Full Unit Cell
#' @description Generates all symmetry-equivalent atomic positions within the unit cell.
#' @param atomic_coordinates A `data.table` of asymmetric atoms.
#' @param symmetry_operations A `data.table` of operations.
#' @return A `data.table` containing all unique atomic positions in the unit cell.
#' @family coordinate processors
#' @export
apply_symmetry_operations <- function(atomic_coordinates, symmetry_operations) {
  if (is.null(atomic_coordinates) || nrow(atomic_coordinates) == 0) {
    return(data.table(
      Label = character(),
      x_a = numeric(),
      y_b = numeric(),
      z_c = numeric()
    ))
  }

  # Helper to safely evaluate the string math (e.g. "x+1/2")
  # Uses strict formatting to prevent parsing errors
  safe_eval <- function(op_str, x, y, z) {
    # Replace coordinates with their values wrapped in parens to preserve order of operations
    # e.g., "-x" where x is -0.1 becomes "-(-0.1)" instead of "--0.1"
    expr_str <- gsub("x", sprintf("(%.10f)", x), op_str, fixed = TRUE)
    expr_str <- gsub("y", sprintf("(%.10f)", y), expr_str, fixed = TRUE)
    expr_str <- gsub("z", sprintf("(%.10f)", z), expr_str, fixed = TRUE)

    tryCatch({
      eval(parse(text = expr_str))
    }, error = function(e) {
      stop(sprintf(
        "Failed to parse symmetry operation: '%s' -> '%s'",
        op_str,
        expr_str
      ))
    })
  }

  n_atoms <- nrow(atomic_coordinates)
  n_ops <- nrow(symmetry_operations)

  # Expand structure using optimized lapply (Vectorized approach, no set() loops)
  expanded_list <- lapply(seq_len(n_atoms), function(i) {
    row <- atomic_coordinates[i, ]
    lbl <- row$Label
    val_x <- row$x_a
    val_y <- row$y_b
    val_z <- row$z_c

    # Apply all symmetry operations to this specific atom
    new_coords <- lapply(seq_len(n_ops), function(j) {
      op_x <- symmetry_operations$x[j]
      op_y <- symmetry_operations$y[j]
      op_z <- symmetry_operations$z[j]

      data.table(
        Label = paste(lbl, j, sep = "_"),
        x_a = safe_eval(op_x, val_x, val_y, val_z) %% 1,
        y_b = safe_eval(op_y, val_x, val_y, val_z) %% 1,
        z_c = safe_eval(op_z, val_x, val_y, val_z) %% 1
      )
    })

    rbindlist(new_coords)
  })

  transformed_coords <- rbindlist(expanded_list)

  # Safety check: Ensure columns exist before processing
  if (nrow(transformed_coords) == 0) {
    return(data.table(
      Label = character(),
      x_a = numeric(),
      y_b = numeric(),
      z_c = numeric()
    ))
  }

  # Round to handle floating point errors before finding unique rows
  precision <- 6
  transformed_coords[, `:=`(
    x_a = round(x_a, precision),
    y_b = round(y_b, precision),
    z_c = round(z_c, precision)
  )]

  transformed_coords <- unique(transformed_coords, by = c("x_a", "y_b", "z_c"))
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
  # Safety check: if input is empty or missing columns, return empty
  if (is.null(transformed_coords) ||
      nrow(transformed_coords) == 0 ||
      !"x_a" %in% names(transformed_coords)) {
    return(data.table(
      Label = character(),
      x_a = numeric(),
      y_b = numeric(),
      z_c = numeric()
    ))
  }

  n_cells <- 1
  cell_indices <- as.data.table(expand.grid(
    x = -n_cells:n_cells,
    y = -n_cells:n_cells,
    z = -n_cells:n_cells
  ))

  expanded_coords <- rbindlist(lapply(1:nrow(cell_indices), function(i) {
    cell_shift <- cell_indices[i, ]

    # We copy the original table and shift the coordinates
    dt <- copy(transformed_coords)
    dt[, `:=`(
      Label = paste(
        Label,
        paste(cell_shift$x, cell_shift$y, cell_shift$z, sep = "_"),
        sep = "_"
      ),
      x_a = x_a + cell_shift$x,
      y_b = y_b + cell_shift$y,
      z_c = z_c + cell_shift$z
    )]
    return(dt)
  }))

  # Round to handle floating point errors
  precision <- 6
  expanded_coords[, `:=`(
    x_a = round(x_a, precision),
    y_b = round(y_b, precision),
    z_c = round(z_c, precision)
  )]

  expanded_coords <- unique(expanded_coords, by = c("x_a", "y_b", "z_c"))
  return(expanded_coords)
}
