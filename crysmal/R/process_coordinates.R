# Internal helper function to scale error values based on value's decimal places.
# Not exported
scale_error_internal <- function(value_str, error_str) {
  if (is.na(error_str) || error_str == "") return(NA_real_)

  value_str <- as.character(value_str)
  error_str <- as.character(error_str)

  parts <- strsplit(value_str, "\\.")[[1]]
  decimal_places <- 0
  if (length(parts) == 2 && !is.na(parts[2])) { # Check parts[2] is not NA
    decimal_places <- nchar(parts[2])
  } else if (length(parts) == 1 && grepl("^[0-9]+$", parts[1])) {
    # Value is an integer, no decimal part
    decimal_places <- 0
  } # Other cases (e.g. scientific notation) might not be handled by this simple parser.

  # CIF error notation (e.g., 1.23(5) -> error 0.05; 1.23(15) -> error 0.15)
  # implies error is on the last significant figures of the value.
  if (decimal_places > 0) {
    # Error applies to the precision shown by the value's decimal places.
    # e.g. value 1.234, error(5) -> 0.005; value 1.234, error(15) -> 0.015
    scaled_error <- as.numeric(error_str) / (10^decimal_places)
  } else { # Integer value, e.g. 123(5) -> error 5
    scaled_error <- as.numeric(error_str)
  }
  return(scaled_error)
}


#' Extract Atomic Coordinates from CIF Content
#'
#' @param cif_content A data.table representing the CIF file content.
#' @return A data.table with atomic coordinates (Label, x_a, y_b, z_c) and errors,
#'         or NULL if coordinates are not found or parsing fails.
#' @export
#' @examples
#' # cif_lines <- c("loop_", "_atom_site_label", "_atom_site_fract_x",
#' #                "_atom_site_fract_y", "_atom_site_fract_z",
#' #                "_atom_site_occupancy", # This can be a marker for data start
#' #                "Si1 X X X 0.0 0.0 0.0 1.0") # Example line structure
#' # cif_data <- data.table(V1 = cif_lines)
#' # extract_atomic_coordinates(cif_data)
extract_atomic_coordinates <- function(cif_content) {
  if (!is.data.table(cif_content) || !"V1" %in% names(cif_content)) return(NULL)

  labels <- vector()
  fractional_x <- vector()
  fractional_y <- vector()
  fractional_z <- vector()
  x_errors <- vector()
  y_errors <- vector()
  z_errors <- vector()

  # Locate the start of the atom site loop data.
  # A robust parser would interpret all _atom_site_ tags in a loop_ block.
  # This version uses common tags to identify the data block.
  atom_loop_header_indices <- grep("^_atom_site_", cif_content$V1)
  if (length(atom_loop_header_indices) == 0) return(NULL)

  # Determine start of data values by finding a common last tag in the header.
  atom_data_start_line_idx_tag <- grep("_atom_site_occupancy", cif_content$V1)
  if (length(atom_data_start_line_idx_tag) == 0) {
    atom_data_start_line_idx_tag <- grep("_atom_site_U_iso_or_equiv", cif_content$V1)
    if (length(atom_data_start_line_idx_tag) == 0) {
      if(length(atom_loop_header_indices) > 0) {
        atom_data_start_line_idx_tag <- max(atom_loop_header_indices)
      } else {
        return(NULL) # No atom site information identifiable
      }
    }
  }
  atom_data_start_line_idx_tag <- atom_data_start_line_idx_tag[1] # Use the first match

  atom_data_start_content_line <- atom_data_start_line_idx_tag + 1

  # Find end of the atomic data section.
  # Looks for next 'loop_', a blank line, or another section marker.
  section_end_candidates <- c(
    grep("^loop_", cif_content$V1[(atom_data_start_content_line):nrow(cif_content)]),
    grep("^\\s*$", cif_content$V1[(atom_data_start_content_line):nrow(cif_content)]),
    grep("^(#|data_)", cif_content$V1[(atom_data_start_content_line):nrow(cif_content)])
  )

  atom_data_end_content_line <- nrow(cif_content) # Default to end of file
  if (length(section_end_candidates) > 0) {
    atom_data_end_content_line <- min(section_end_candidates) + (atom_data_start_content_line - 1) -1
  }

  if (atom_data_start_content_line > atom_data_end_content_line) return(NULL)

  atom_lines <- cif_content$V1[atom_data_start_content_line:atom_data_end_content_line]
  atom_lines <- atom_lines[nzchar(trimws(atom_lines))] # Remove empty lines
  if (length(atom_lines) == 0) return(NULL)

  # Parsing assumes a specific column order after splitting the line,
  # typically: Label, Type, FractX, FractY, FractZ, Occupancy, etc.
  # This implementation expects fractional coordinates at specific indices (5,6,7)
  # after strsplit, which depends on the number of preceding items in the atom loop.
  # A full CIF parser would map _atom_site_ tags to columns dynamically.
  # This version relies on a consistent structure in the input CIF files.

  for (line in atom_lines) {
    properties <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(properties) < 7 ||startsWith(properties[1], "_")) next # Basic check & skip headers

    # Assuming Label at properties[1], and x,y,z values at properties[5], [6], [7]
    # This indexing is based on common CIF structures where several _atom_site_
    # tags (like type_symbol, disorder_group, etc.) might precede coordinates.
    labels <- c(labels, as.character(properties[1]))

    fractional_x_val_str <- gsub("\\(.*\\)", "", properties[5])
    fractional_y_val_str <- gsub("\\(.*\\)", "", properties[6])
    fractional_z_val_str <- gsub("\\(.*\\)", "", properties[7])

    x_error_str <- ifelse(grepl("\\(", properties[5]), gsub(".*\\((.*)\\).*", "\\1", properties[5]), NA_character_)
    y_error_str <- ifelse(grepl("\\(", properties[6]), gsub(".*\\((.*)\\).*", "\\1", properties[6]), NA_character_)
    z_error_str <- ifelse(grepl("\\(", properties[7]), gsub(".*\\((.*)\\).*", "\\1", properties[7]), NA_character_)

    fractional_x <- c(fractional_x, suppressWarnings(as.numeric(fractional_x_val_str)))
    fractional_y <- c(fractional_y, suppressWarnings(as.numeric(fractional_y_val_str)))
    fractional_z <- c(fractional_z, suppressWarnings(as.numeric(fractional_z_val_str)))

    x_errors <- c(x_errors, scale_error_internal(fractional_x_val_str, x_error_str))
    y_errors <- c(y_errors, scale_error_internal(fractional_y_val_str, y_error_str))
    z_errors <- c(z_errors, scale_error_internal(fractional_z_val_str, z_error_str))
  }

  if (length(labels) == 0) return(NULL) # No atoms found

  atomic_coordinates <- data.table(
    Label = labels,
    x_a = as.numeric(fractional_x),
    y_b = as.numeric(fractional_y),
    z_c = as.numeric(fractional_z),
    x_error = as.numeric(x_errors),
    y_error = as.numeric(y_errors),
    z_error = as.numeric(z_errors)
  )

  # Remove rows where coordinates are NA (can happen from parsing issues)
  atomic_coordinates <- atomic_coordinates[!is.na(x_a) & !is.na(y_b) & !is.na(z_c)]
  if(nrow(atomic_coordinates) == 0) return(NULL)

  return(atomic_coordinates)
}


#' Extract Symmetry Operations from CIF Content
#'
#' @param cif_content A data.table representing the CIF file content.
#' @return A data.table with symmetry operations (x, y, z components),
#'         or NULL if not found.
#' @export
#' @examples
#' # cif_lines <- c("_space_group_symop_operation_xyz", "1 'x, y, z'", "2 '-x, -y, -z'")
#' # cif_data <- data.table(V1 = cif_lines)
#' # extract_symmetry_operations(cif_data)
extract_symmetry_operations <- function(cif_content) {
  if (!is.data.table(cif_content) || !"V1" %in% names(cif_content)) return(NULL)

  symmetry_start_line_idx <- grep("^_space_group_symop_operation_xyz", cif_content$V1)
  if (length(symmetry_start_line_idx) == 0) {
    # Fallback for older CIFs or the alternative common tag
    symmetry_start_line_idx <- grep("^_symmetry_equiv_pos_as_xyz", cif_content$V1)
    if (length(symmetry_start_line_idx) == 0) return(NULL)
  }
  symmetry_start_line_idx <- symmetry_start_line_idx[1]

  symmetry_data_start_content_line <- symmetry_start_line_idx + 1

  # Find end of the symmetry operations section
  section_end_candidates <- c(
    grep("^loop_", cif_content$V1[symmetry_data_start_content_line:nrow(cif_content)]),
    grep("^\\s*$", cif_content$V1[symmetry_data_start_content_line:nrow(cif_content)]),
    grep("^(#|data_|_)", cif_content$V1[symmetry_data_start_content_line:nrow(cif_content)]) # Stop at next tag
  )

  symmetry_data_end_content_line <- nrow(cif_content) # Default to end of file
  if (length(section_end_candidates) > 0) {
    symmetry_data_end_content_line <- min(section_end_candidates) + (symmetry_data_start_content_line - 1) - 1
  }

  if (symmetry_data_start_content_line > symmetry_data_end_content_line) return(NULL)

  symmetry_lines_raw <- cif_content$V1[symmetry_data_start_content_line:symmetry_data_end_content_line]
  # Clean lines: remove leading numbers (e.g., "1 'x,y,z'") and quotes
  symmetry_lines_cleaned <- sapply(symmetry_lines_raw, function(line) {
    line <- trimws(line)
    line <- sub("^[0-9]+\\s+", "", line) # Remove leading integer and space if present
    line <- gsub("'", "", line) # Remove single quotes
    return(line)
  })

  symmetry_lines_cleaned <- symmetry_lines_cleaned[nzchar(symmetry_lines_cleaned)] # Remove empty lines
  if (length(symmetry_lines_cleaned) == 0) return(NULL)

  # Split each line into x, y, z components by comma
  symmetry_operations_matrix <- stringr::str_split_fixed(symmetry_lines_cleaned, ",\\s*", n = 3)
  if (ncol(symmetry_operations_matrix) != 3) return(NULL) # Expect x,y,z parts

  symmetry_operations <- data.table(
    x = trimws(symmetry_operations_matrix[, 1]),
    y = trimws(symmetry_operations_matrix[, 2]),
    z = trimws(symmetry_operations_matrix[, 3])
  )

  symmetry_operations <- symmetry_operations[x != "" & y != "" & z != ""] # Ensure no empty ops
  if(nrow(symmetry_operations) == 0) return(NULL)

  return(symmetry_operations)
}

# Internal helper function to apply a single operation component
# Not exported
apply_operation_component <- function(operation_str, x_coord, y_coord, z_coord) {
  # Replace x, y, z placeholders with their numeric values for evaluation.
  # Care is needed for order and partial matches (e.g., "x" in "+2*x").
  eval_env <- new.env()
  eval_env$x <- x_coord
  eval_env$y <- y_coord
  eval_env$z <- z_coord

  # Attempt to parse and evaluate the expression.
  # This can be unsafe if operation_str is from an untrusted source,
  # though basic arithmetic operations typical in symmetry ops should be fine.
  tryCatch({
    eval(parse(text = operation_str), envir = eval_env)
  }, error = function(e) {
    warning(paste("Error evaluating symmetry operation component:", operation_str, "-", e$message))
    NA_real_
  })
}


#' Apply Symmetry Operations to Atomic Coordinates
#'
#' @param atomic_coordinates A data.table with columns Label, x_a, y_b, z_c.
#' @param symmetry_operations A data.table with columns x, y, z representing symmetry operations.
#' @return A data.table with transformed coordinates, duplicates removed. Returns NULL on error or if inputs are invalid.
#' @export
#' @importFrom dplyr distinct
#' @examples
#' # atomic_coords <- data.table(Label="Si1", x_a=0.1, y_b=0.2, z_c=0.3)
#' # sym_ops_dt <- data.table(x = c("x", "-x"), y = c("y", "-y"), z = c("z", "-z"))
#' # apply_symmetry_operations(atomic_coords, sym_ops_dt)
apply_symmetry_operations <- function(atomic_coordinates, symmetry_operations) {
  if (is.null(atomic_coordinates) || nrow(atomic_coordinates) == 0 ||
      is.null(symmetry_operations) || nrow(symmetry_operations) == 0) {
    return(atomic_coordinates) # Return original if no ops or coords
  }

  # Helper to expand coordinates for a single atom row
  expand_coords_for_row <- function(atom_row, sym_ops) {
    original_x <- atom_row$x_a
    original_y <- atom_row$y_b
    original_z <- atom_row$z_c
    original_label <- atom_row$Label

    new_coords_list <- vector("list", nrow(sym_ops))

    for (i in 1:nrow(sym_ops)) {
      op_x_str <- sym_ops$x[i]
      op_y_str <- sym_ops$y[i]
      op_z_str <- sym_ops$z[i]

      new_x <- apply_operation_component(op_x_str, original_x, original_y, original_z)
      new_y <- apply_operation_component(op_y_str, original_x, original_y, original_z)
      new_z <- apply_operation_component(op_z_str, original_x, original_y, original_z)

      # Wrap coordinates to [0, 1) range.
      # x_frac = x - floor(x) maps to [0,1) for positive x, and (e.g.) [-0.9, 0) for negative x.
      # A robust way for all x: x_frac = x - floor(x).
      # E.g. -0.1 - floor(-0.1) = -0.1 - (-1) = 0.9.
      # E.g. -1.1 - floor(-1.1) = -1.1 - (-2) = 0.9.
      # This is generally effective for bringing coordinates into the unit cell.
      # However, values like 1.0 become 0.0, and -1.0 becomes 0.0.
      # Using a small tolerance for floor can help with floating point issues near integers.
      tol <- 1e-8
      new_x <- new_x - floor(new_x + tol)
      new_y <- new_y - floor(new_y + tol)
      new_z <- new_z - floor(new_z + tol)

      # An alternative for strict [0,1) range: new_coord %% 1, then adjust negatives if R's % gives them.
      # new_x = new_x %% 1; if (new_x < 0) new_x = new_x + 1

      new_coords_list[[i]] <- data.table(
        Label = paste(original_label, i, sep = "_sym"), # Unique label for symmetrically generated atom
        x_a = new_x,
        y_b = new_y,
        z_c = new_z
      )
    }
    return(rbindlist(new_coords_list))
  }

  transformed_coords_list <- lapply(1:nrow(atomic_coordinates), function(i) {
    expand_coords_for_row(atomic_coordinates[i, ], symmetry_operations)
  })
  transformed_coords <- rbindlist(transformed_coords_list)

  if (is.null(transformed_coords) || nrow(transformed_coords) == 0) return(NULL)

  # Remove duplicate coordinates. Rounding before checking for uniqueness is crucial for floating-point numbers.
  precision_digits <- 6
  transformed_coords_rounded <- copy(transformed_coords)
  transformed_coords_rounded[, `:=`(x_a = round(x_a, precision_digits),
                                    y_b = round(y_b, precision_digits),
                                    z_c = round(z_c, precision_digits))]

  unique_indices <- !duplicated(transformed_coords_rounded[, .(x_a, y_b, z_c)])
  transformed_coords <- transformed_coords[unique_indices]

  return(transformed_coords)
}


#' Expand Transformed Coordinates to Neighboring Unit Cells
#'
#' @param transformed_coords A data.table of atomic coordinates, typically from \code{apply_symmetry_operations}.
#'                           Must have Label, x_a, y_b, z_c.
#' @param n_cells Integer, number of unit cells to expand in each direction (e.g., 1 means -1 to +1).
#' @return A data.table with coordinates expanded into neighboring cells. Returns NULL if input is invalid.
#' @export
#' @importFrom dplyr distinct
#' @examples
#' # coords <- data.table(Label="Si1_sym1", x_a=0.1, y_b=0.2, z_c=0.3)
#' # expand_transformed_coords(coords)
expand_transformed_coords <- function(transformed_coords, n_cells = 1) {
  if (is.null(transformed_coords) || nrow(transformed_coords) == 0) {
    return(NULL)
  }

  # Generate a grid of unit cell translation vectors
  cell_shifts_xyz <- -n_cells:n_cells
  cell_indices <- as.data.table(expand.grid(
    dx = cell_shifts_xyz,
    dy = cell_shifts_xyz,
    dz = cell_shifts_xyz
  ))

  expanded_coords_list <- vector("list", nrow(cell_indices))

  for (i in 1:nrow(cell_indices)) {
    shift_vector <- cell_indices[i]

    temp_coords <- copy(transformed_coords)
    # Generate unique labels for atoms translated to different cells.
    temp_coords[, `:=`(
      x_a = x_a + shift_vector$dx,
      y_b = y_b + shift_vector$dy,
      z_c = z_c + shift_vector$dz,
      Label = paste(Label, shift_vector$dx, shift_vector$dy, shift_vector$dz, sep = "_cell")
    )]
    expanded_coords_list[[i]] <- temp_coords
  }
  expanded_coords <- rbindlist(expanded_coords_list)

  if (is.null(expanded_coords) || nrow(expanded_coords) == 0) return(NULL)

  # Remove duplicate coordinates if any might arise (e.g. from n_cells=0 or complex scenarios).
  # Rounding before distinct is good practice for floating point numbers, though simple additions here
  # are less prone to creating numerically identical but bitwise different coordinates.
  precision_digits <- 6
  expanded_coords_rounded <- copy(expanded_coords)
  expanded_coords_rounded[, `:=`(x_a = round(x_a, precision_digits),
                                 y_b = round(y_b, precision_digits),
                                 z_c = round(z_c, precision_digits))]
  unique_indices <- !duplicated(expanded_coords_rounded[, .(x_a, y_b, z_c)])
  expanded_coords <- expanded_coords[unique_indices]

  return(expanded_coords)
}
