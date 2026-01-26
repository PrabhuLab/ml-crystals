#' @title Apply Symmetry Operations to Generate a Full Unit Cell
#' @description Generates all symmetry-equivalent atomic positions within the unit cell.
#' @param atomic_coordinates A `data.table` of asymmetric atoms.
#' @param symmetry_operations A `data.table` of operations.
#' @param precision An integer specifying the number of decimal places used to
#'   round coordinates when identifying unique atoms. Defaults to 4.
#' @return A `data.table` containing all unique atomic positions in the unit cell.
#' @family coordinate processors
#' @export
apply_symmetry_operations <- function(atomic_coordinates,
                                      symmetry_operations,
                                      precision = 4) {
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
    # Increase internal string precision to 16 to prevent truncation errors
    # with fractional coordinates (e.g., 1/3) before evaluation.
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
        # Apply modulo initially to handle operations like x+1/2
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

  # Round and modulo again to enforce unit cell boundaries and merge close atoms
  transformed_coords[, `:=`(
    x_a = round(x_a, precision) %% 1,
    y_b = round(y_b, precision) %% 1,
    z_c = round(z_c, precision) %% 1
  )]

  # Remove duplicates
  transformed_coords <- unique(transformed_coords, by = c("x_a", "y_b", "z_c"))
  return(transformed_coords)
}


#' @title Expand Coordinates into a Supercell
#' @description Takes a set of atomic coordinates within a single unit cell and
#'   replicates them into a 3x3x3 supercell grid.
#' @param transformed_coords A `data.table` of atom positions within one unit cell.
#' @param precision An integer specifying the number of decimal places used to
#'   round coordinates. Defaults to 4.
#' @return A `data.table` containing all atomic positions in the expanded supercell.
#' @family coordinate processors
#' @export
expand_transformed_coords <- function(transformed_coords, precision = 4) {
  if (is.null(transformed_coords) || nrow(transformed_coords) == 0 || !"x_a" %in% names(transformed_coords)) {
    return(data.table(Label=character(), SymOp=integer(), Tx=integer(), Ty=integer(), Tz=integer(), x_a=numeric(), y_b=numeric(), z_c=numeric()))
  }

  cell_indices <- as.data.table(expand.grid(x = -1:1, y = -1:1, z = -1:1))

  expanded_coords <- rbindlist(lapply(1:nrow(cell_indices), function(i) {
    cell_shift <- cell_indices[i, ]
    dt <- copy(transformed_coords)
    dt[, `:=`(Tx = cell_shift$x, Ty = cell_shift$y, Tz = cell_shift$z)]
    dt[, Label := paste(Label, Tx, Ty, Tz, sep = "_")]
    dt[, `:=`(x_a = x_a + Tx, y_b = y_b + Ty, z_c = z_c + Tz)]
    return(dt)
  }))

  # Round to handle floating point noise
  expanded_coords[, `:=`(
    x_a = round(x_a, precision),
    y_b = round(y_b, precision),
    z_c = round(z_c, precision)
  )]

  expanded_coords <- unique(expanded_coords, by = c("x_a", "y_b", "z_c"))
  return(expanded_coords)
}
