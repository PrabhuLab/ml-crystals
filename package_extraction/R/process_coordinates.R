#' @title Apply Symmetry Operations to Generate a Full Unit Cell
#' @description Generates all symmetry-equivalent atomic positions within the unit cell
#'   by applying the list of symmetry operations to the initial set of
#'   asymmetric atoms.
#' @param atomic_coordinates A `data.table` of asymmetric atoms from
#'   `extract_atomic_coordinates`.
#' @param symmetry_operations A `data.table` of operations from
#'   `extract_symmetry_operations`.
#' @return A `data.table` containing all unique atomic positions in the unit cell.
#' @family coordinate processors
#' @export
#' @examples
#' cif_file <- system.file("extdata", "ICSD422.cif", package = "crysmal")
#' if (file.exists(cif_file)) {
#'   cif_content <- read_cif_files(cif_file)[[1]]
#'   atoms <- extract_atomic_coordinates(cif_content)
#'   sym_ops <- extract_symmetry_operations(cif_content)
#'   full_cell_coords <- apply_symmetry_operations(atoms, sym_ops)
#'   print(full_cell_coords)
#' }
apply_symmetry_operations <- function(atomic_coordinates, symmetry_operations) {
  apply_operation <- function(operation, x, y, z) {
    operation <- gsub("x", sprintf("(%f)", x), operation)
    operation <- gsub("y", sprintf("(%f)", y), operation)
    operation <- gsub("z", sprintf("(%f)", z), operation)
    eval(parse(text = operation))
  }
  expand_coordinates <- function(row) {
    x_val <- row$x_a; y_val <- row$y_b; z_val <- row$z_c
    rbindlist(lapply(1:nrow(symmetry_operations), function(i) {
      new_x <- apply_operation(symmetry_operations[i, x], x_val, y_val, z_val)
      new_y <- apply_operation(symmetry_operations[i, y], x_val, y_val, z_val)
      new_z <- apply_operation(symmetry_operations[i, z], x_val, y_val, z_val)
      data.table(Label = paste(row$Label, i, sep = "_"), x_a = new_x %% 1, y_b = new_y %% 1, z_c = new_z %% 1)
    }))
  }
  transformed_coords <- rbindlist(lapply(1:nrow(atomic_coordinates), function(i) expand_coordinates(atomic_coordinates[i,])))

  # Round to handle floating point errors before finding unique rows
  precision <- 6
  transformed_coords[, `:=`(x_a = round(x_a, precision), y_b = round(y_b, precision), z_c = round(z_c, precision))]
  transformed_coords <- unique(transformed_coords, by = c("x_a", "y_b", "z_c"))
  return(transformed_coords)
}


#' @title Expand Coordinates into a Supercell
#' @description Takes a set of atomic coordinates within a single unit cell and
#'   replicates them into a 3x3x3 supercell grid.
#' @param transformed_coords A `data.table` of atom positions within one unit cell,
#'   typically from `apply_symmetry_operations`.
#' @return A `data.table` containing all atomic positions in the expanded supercell.
#' @family coordinate processors
#' @export
expand_transformed_coords <- function(transformed_coords) {
  n_cells <- 1
  cell_indices <- as.data.table(expand.grid(x = -n_cells:n_cells, y = -n_cells:n_cells, z = -n_cells:n_cells))

  expanded_coords <- rbindlist(lapply(1:nrow(cell_indices), function(i) {
    cell_shift <- cell_indices[i]
    transformed_coords[, .(Label = paste(Label, paste(cell_shift$x, cell_shift$y, cell_shift$z, sep = "_"), sep = "_"),
                           x_a = x_a + cell_shift$x, y_b = y_b + cell_shift$y, z_c = z_c + cell_shift$z)]
  }))

  # Round to handle floating point errors before finding unique rows
  precision <- 6
  expanded_coords[, `:=`(x_a = round(x_a, precision), y_b = round(y_b, precision), z_c = round(z_c, precision))]
  expanded_coords <- unique(expanded_coords, by = c("x_a", "y_b", "z_c"))
  return(expanded_coords)
}
