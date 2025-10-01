# Internal helper function to parse a single symmetry operation string component
# (e.g., "x-1/2", "y+z", "-z").
# It returns a list of two elements: 'W_row' (numeric vector of length 3 for x,y,z coefficients)
# and 'w_scalar' (numeric scalar for translation).
#' @noRd
.parse_sym_op_component <- function(op_str) {
  W_row <- c(0, 0, 0) # Coefficients for x, y, z
  w_scalar <- 0

  # Clean spaces and ensure consistent sign representation for fractions
  op_str <- gsub("\\s+", "", op_str)
  # Normalize signs for easier splitting, e.g., "x-1/2" -> "x+-1/2"
  op_str <- gsub("(?<=[^\\d/])-", "+-", op_str, perl = TRUE)
  # Split into terms like "x", "-y", "+1/2". Handles cases like "1/2+x" too.
  terms <- stringr::str_split(op_str, "(?<=[^/\\d])(?=[+-])", simplify = TRUE)
  terms <- terms[terms != ""]

  for (term in terms) {
    if (term == "") next

    # Handle fractional translation (e.g., +1/2, -1/3) or integer/decimal numbers
    if (grepl("^[+-]?(\\d+/\\d+|\\d*\\.\\d+|\\d+)$", term)) {
      w_scalar <- w_scalar + eval(parse(text = term)) # Safely evaluate numerical fractions/constants
      next
    }

    # Handle x, y, z terms
    if (grepl("x", term, fixed = TRUE)) {
      coeff_str <- gsub("x", "", term, fixed = TRUE)
      if (coeff_str == "" || coeff_str == "+") W_row[1] <- W_row[1] + 1
      else if (coeff_str == "-") W_row[1] <- W_row[1] - 1
      else W_row[1] <- W_row[1] + eval(parse(text = coeff_str)) # Handles "1/2x", "-1/3x"
    }
    if (grepl("y", term, fixed = TRUE)) {
      coeff_str <- gsub("y", "", term, fixed = TRUE)
      if (coeff_str == "" || coeff_str == "+") W_row[2] <- W_row[2] + 1
      else if (coeff_str == "-") W_row[2] <- W_row[2] - 1
      else W_row[2] <- W_row[2] + eval(parse(text = coeff_str))
    }
    if (grepl("z", term, fixed = TRUE)) {
      coeff_str <- gsub("z", "", term, fixed = TRUE)
      if (coeff_str == "" || coeff_str == "+") W_row[3] <- W_row[3] + 1
      else if (coeff_str == "-") W_row[3] <- W_row[3] - 1
      else W_row[3] <- W_row[3] + eval(parse(text = coeff_str))
    }
  }
  list(W_row = W_row, w_scalar = w_scalar)
}

# Main internal function to parse all symmetry operations into matrices
#' @noRd
.parse_symmetry_operations_to_matrices <- function(symmetry_operations_dt) {
  if (nrow(symmetry_operations_dt) == 0) return(NULL)

  num_sym_ops <- nrow(symmetry_operations_dt)
  W_matrices <- vector("list", num_sym_ops)
  w_vectors <- matrix(0, nrow = num_sym_ops, ncol = 3)

  for (i in 1:num_sym_ops) {
    op_row <- symmetry_operations_dt[i]
    W_current <- matrix(0, 3, 3)
    w_current <- c(0, 0, 0)

    # Parse x component
    parsed_x <- .parse_sym_op_component(op_row$x)
    W_current[1, ] <- parsed_x$W_row
    w_current[1] <- parsed_x$w_scalar

    # Parse y component
    parsed_y <- .parse_sym_op_component(op_row$y)
    W_current[2, ] <- parsed_y$W_row
    w_current[2] <- parsed_y$w_scalar

    # Parse z component
    parsed_z <- .parse_sym_op_component(op_row$z)
    W_current[3, ] <- parsed_z$W_row
    w_current[3] <- parsed_z$w_scalar

    W_matrices[[i]] <- W_current
    w_vectors[i, ] <- w_current
  }
  list(W_matrices = W_matrices, w_vectors = w_vectors)
}

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
#' cif_file <- system.file("extdata", "ICSD422.cif", package = "crystract")
#' if (file.exists(cif_file)) {
#'   cif_content <- read_cif_files(cif_file)[[1]]
#'   atoms <- extract_atomic_coordinates(cif_content)
#'   sym_ops <- extract_symmetry_operations(cif_content)
#'   full_cell_coords <- apply_symmetry_operations(atoms, sym_ops)
#'   print(full_cell_coords)
#' }
apply_symmetry_operations <- function(atomic_coordinates, symmetry_operations) {
  if (nrow(atomic_coordinates) == 0 || nrow(symmetry_operations) == 0) {
    return(data.table(Label = character(), x_a = numeric(), y_b = numeric(), z_c = numeric()))
  }

  # Parse all symmetry operations into matrix form once
  parsed_ops <- .parse_symmetry_operations_to_matrices(symmetry_operations)
  W_matrices <- parsed_ops$W_matrices
  w_vectors <- parsed_ops$w_vectors

  num_atoms <- nrow(atomic_coordinates)
  num_sym_ops <- length(W_matrices)
  total_generated_coords <- num_atoms * num_sym_ops

  # Pre-allocate output data.table for efficiency
  result_dt <- data.table(
    Label = character(total_generated_coords),
    x_a = numeric(total_generated_coords),
    y_b = numeric(total_generated_coords),
    z_c = numeric(total_generated_coords)
  )

  # Prepare original coordinates as a matrix (each row is a coordinate vector)
  coords_mat <- as.matrix(atomic_coordinates[, .(x_a, y_b, z_c)])

  # Fill labels in a vectorized way
  result_dt$Label <- paste(
    rep(atomic_coordinates$Label, each = num_sym_ops),
    rep(1:num_sym_ops, times = num_atoms), # Use operation index for unique labeling
    sep = "_"
  )

  # Apply transformations and fill coordinates using set() for speed
  idx_start <- 1
  for (i in 1:num_atoms) {
    original_coord <- coords_mat[i, ] # 1x3 numeric vector
    for (j in 1:num_sym_ops) {
      W_mat <- W_matrices[[j]] # 3x3 matrix
      w_vec <- w_vectors[j, ]  # 1x3 numeric vector

      # Apply affine transformation: x' = Wx + w
      # original_coord is 1x3, so W_mat %*% original_coord might error if original_coord is not a column vector.
      # Ensure original_coord is a column vector (3x1) for standard matrix multiplication.
      new_coord <- (W_mat %*% as.numeric(original_coord)) + w_vec

      row_idx <- idx_start + j - 1
      set(result_dt, row_idx, 2L, new_coord[1]) # x_a
      set(result_dt, row_idx, 3L, new_coord[2]) # y_b
      set(result_dt, row_idx, 4L, new_coord[3]) # z_c
    }
    idx_start <- idx_start + num_sym_ops
  }

  # Normalize coordinates to [0, 1) range
  result_dt[, `:=`(
    x_a = x_a %% 1,
    y_b = y_b %% 1,
    z_c = z_c %% 1
  )]
  # Handle negative results from modulo for consistency (e.g., -0.5 %% 1 is -0.5, needs to be 0.5)
  result_dt[x_a < 0, x_a := x_a + 1]
  result_dt[y_b < 0, y_b := y_b + 1]
  result_dt[z_c < 0, z_c := z_c + 1]

  # Round to handle floating point errors before finding unique rows
  precision <- 6
  result_dt[, `:=`(
    x_a = round(x_a, precision),
    y_b = round(y_b, precision),
    z_c = round(z_c, precision)
  )]
  transformed_coords <- unique(result_dt, by = c("x_a", "y_b", "z_c"))
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
  cell_indices <- as.data.table(expand.grid(
    x = -n_cells:n_cells,
    y = -n_cells:n_cells,
    z = -n_cells:n_cells
  ))

  expanded_coords <- rbindlist(lapply(1:nrow(cell_indices), function(i) {
    cell_shift <- cell_indices[i]
    transformed_coords[, .(
      Label = paste(
        Label,
        paste(cell_shift$x, cell_shift$y, cell_shift$z, sep = "_"),
        sep = "_"
      ),
      x_a = x_a + cell_shift$x,
      y_b = y_b + cell_shift$y,
      z_c = z_c + cell_shift$z
    )]
  }))

  # Round to handle floating point errors before finding unique rows
  precision <- 6
  expanded_coords[, `:=`(
    x_a = round(x_a, precision),
    y_b = round(y_b, precision),
    z_c = round(z_c, precision)
  )]
  expanded_coords <- unique(expanded_coords, by = c("x_a", "y_b", "z_c"))
  return(expanded_coords)
}
