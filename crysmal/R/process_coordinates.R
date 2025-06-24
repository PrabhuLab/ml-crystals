# Internal helper function to scale error values from the notebook.
# This version is less robust but matches the original logic.
# Not exported
scale_error_notebook <- function(value_str, error_str) {
  if (is.na(error_str) || error_str == "" || is.na(value_str)) return(NA_real_)
  value_str <- as.character(value_str); error_str <- as.character(error_str)
  split_val <- strsplit(value_str, "\\.")[[1]]
  if (length(split_val) < 2 || is.na(split_val[2])) return(NA_real_)
  decimal_places <- nchar(split_val[2])
  if (nchar(error_str) > decimal_places) return(NA_real_)
  scaled_error <- as.numeric(paste0("0.", strrep("0", decimal_places - nchar(error_str)), error_str))
  return(scaled_error)
}

# --- Documentation Group 1: Coordinate and Symmetry Extraction ---

#' Extract Coordinate and Symmetry Data from CIF Content
#'
#' Functions to extract atomic coordinates and symmetry operations from CIF file content.
#'
#' @param cif_content A data.table representing the CIF file content, typically
#'   from \code{fread(..., sep = "\\n", header = FALSE)}.
#' @return A data.table with the extracted information, or NULL if not found.
#' @name extract_coordinate_data
NULL

#' @rdname extract_coordinate_data
#' @export
#' @examples
#' # cif_lines <- c("loop_", "_atom_site_label", "_atom_site_fract_x",
#' #                "_atom_site_occupancy",
#' #                "Si1 X X X 0.0 0.0 0.0 1.0")
#' # cif_data <- data.table(V1 = cif_lines)
#' # extract_atomic_coordinates(cif_data)
extract_atomic_coordinates <- function(cif_content) {
  if (!is.data.table(cif_content) || !"V1" %in% names(cif_content)) return(NULL)
  atom_start <- grep("_atom_site_occupancy", cif_content$V1)
  if (length(atom_start) == 0) {
    atom_start <- grep("_atom_site_U_iso_or_equiv", cif_content$V1)
    if(length(atom_start) == 0) return(NULL)
  }
  atom_start <- atom_start[1]
  section_end_candidates <- c(
    grep("^loop_", cif_content$V1[(atom_start + 1):nrow(cif_content)]),
    grep("^#", cif_content$V1[(atom_start + 1):nrow(cif_content)]),
    grep("^\\s*$", cif_content$V1[(atom_start + 1):nrow(cif_content)])
  )
  atom_end <- nrow(cif_content) + 1
  if (length(section_end_candidates) > 0) atom_end <- min(section_end_candidates) + atom_start
  atom_lines <- cif_content$V1[(atom_start + 1):(atom_end - 1)]
  atom_lines <- atom_lines[nzchar(trimws(atom_lines))]
  if (length(atom_lines) == 0) return(NULL)
  parsed_atoms <- lapply(atom_lines, function(line) {
    properties <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(properties) < 7 || startsWith(properties[1], "_")) return(NULL)
    label <- as.character(properties[1])
    fx_str <- gsub("\\(.*\\)", "", properties[5]); fy_str <- gsub("\\(.*\\)", "", properties[6]); fz_str <- gsub("\\(.*\\)", "", properties[7])
    x_err_str <- ifelse(grepl("\\(", properties[5]), gsub(".*\\((.*)\\).*", "\\1", properties[5]), NA_character_)
    y_err_str <- ifelse(grepl("\\(", properties[6]), gsub(".*\\((.*)\\).*", "\\1", properties[6]), NA_character_)
    z_err_str <- ifelse(grepl("\\(", properties[7]), gsub(".*\\((.*)\\).*", "\\1", properties[7]), NA_character_)
    data.table(Label = label, x_a = suppressWarnings(as.numeric(fx_str)), y_b = suppressWarnings(as.numeric(fy_str)), z_c = suppressWarnings(as.numeric(fz_str)),
               x_error = scale_error_notebook(fx_str, x_err_str), y_error = scale_error_notebook(fy_str, y_err_str), z_error = scale_error_notebook(fz_str, z_err_str))
  })
  atomic_coordinates <- rbindlist(Filter(Negate(is.null), parsed_atoms))
  if (nrow(atomic_coordinates) == 0) return(NULL)
  return(atomic_coordinates)
}

#' @rdname extract_coordinate_data
#' @export
#' @examples
#' # cif_lines <- c("_space_group_symop_operation_xyz", "1 'x, y, z'", "2 '-x, -y, -z'")
#' # cif_data <- data.table(V1 = cif_lines)
#' # extract_symmetry_operations(cif_data)
extract_symmetry_operations <- function(cif_content) {
  if (!is.data.table(cif_content) || !"V1" %in% names(cif_content)) return(NULL)
  symmetry_start <- grep("_space_group_symop_operation_xyz", cif_content$V1)
  if (length(symmetry_start) == 0) return(NULL)
  symmetry_start <- symmetry_start[1]
  section_end_candidates <- c(
    grep("^loop_", cif_content$V1[(symmetry_start + 1):nrow(cif_content)]),
    grep("^#", cif_content$V1[(symmetry_start + 1):nrow(cif_content)]),
    grep("^\\s*$", cif_content$V1[(symmetry_start + 1):nrow(cif_content)])
  )
  symmetry_end <- nrow(cif_content) + 1
  if (length(section_end_candidates) > 0) symmetry_end <- min(section_end_candidates) + symmetry_start
  symmetry_lines <- cif_content$V1[(symmetry_start + 1):(symmetry_end - 1)]
  symmetry_lines <- symmetry_lines[nzchar(trimws(symmetry_lines))]
  if (length(symmetry_lines) == 0) return(NULL)
  symmetry_lines_cleaned <- gsub("^[0-9]+\\s+", "", symmetry_lines)
  symmetry_ops_matrix <- stringr::str_split_fixed(symmetry_lines_cleaned, ",\\s*", n = 3)
  symmetry_operations <- data.table(
    x = gsub("'", "", trimws(symmetry_ops_matrix[, 1])),
    y = gsub("'", "", trimws(symmetry_ops_matrix[, 2])),
    z = gsub("'", "", trimws(symmetry_ops_matrix[, 3]))
  )
  return(symmetry_operations)
}


# --- Documentation Group 2: Coordinate Transformation and Expansion ---

#' Transform and Expand Atomic Coordinates
#'
#' Functions to apply symmetry operations and expand the resulting unique coordinates
#' into neighboring unit cells.
#'
#' @param atomic_coordinates A data.table with columns Label, x_a, y_b, z_c.
#' @param symmetry_operations A data.table with columns x, y, z representing symmetry operations.
#' @param transformed_coords A data.table of atomic coordinates, typically from \code{apply_symmetry_operations}.
#' @param n_cells Integer, number of unit cells to expand in each direction (e.g., 1 means -1 to +1).
#' @return A data.table with the transformed or expanded coordinates.
#' @name transform_coords
NULL

#' @rdname transform_coords
#' @export
#' @examples
#' # atomic_coords <- data.table(Label="Si1", x_a=0.1, y_b=0.2, z_c=0.3)
#' # sym_ops_dt <- data.table(x = c("x", "-x"), y = c("y", "-y"), z = c("z", "-z"))
#' # apply_symmetry_operations(atomic_coords, sym_ops_dt)
apply_symmetry_operations <- function(atomic_coordinates, symmetry_operations) {
  if (is.null(atomic_coordinates) || nrow(atomic_coordinates) == 0 ||
      is.null(symmetry_operations) || nrow(symmetry_operations) == 0) {
    return(atomic_coordinates)
  }
  apply_operation <- function(operation, x, y, z) {
    operation <- gsub("x", sprintf("(%f)", x), operation, fixed = TRUE)
    operation <- gsub("y", sprintf("(%f)", y), operation, fixed = TRUE)
    operation <- gsub("z", sprintf("(%f)", z), operation, fixed = TRUE)
    eval(parse(text = operation))
  }
  expand_coordinates <- function(row) {
    x <- row$x_a; y <- row$y_b; z <- row$z_c
    rbindlist(lapply(1:nrow(symmetry_operations), function(i) {
      new_x <- apply_operation(symmetry_operations$x[i], x, y, z)
      new_y <- apply_operation(symmetry_operations$y[i], x, y, z)
      new_z <- apply_operation(symmetry_operations$z[i], x, y, z)
      data.table(
        Label = paste(row$Label, i, sep = "_"),
        x_a = ifelse(new_x > 1, new_x - 1, ifelse(new_x < 0, new_x + 1, new_x)),
        y_b = ifelse(new_y > 1, new_y - 1, ifelse(new_y < 0, new_y + 1, new_y)),
        z_c = ifelse(new_z > 1, new_z - 1, ifelse(new_z < 0, new_z + 1, new_z))
      )
    }))
  }
  transformed_coords <- rbindlist(lapply(1:nrow(atomic_coordinates), function(i) expand_coordinates(atomic_coordinates[i,])))

  # Replicating notebook's use of dplyr::distinct by using data.table::unique on rounded columns
  temp_rounded <- copy(transformed_coords)
  temp_rounded[, `:=` (x_a_round = round(x_a, 5), y_b_round = round(y_b, 5), z_c_round = round(z_c, 5))]
  unique_coords <- unique(temp_rounded, by = c("x_a_round", "y_b_round", "z_c_round"))

  return(unique_coords[, .(Label, x_a, y_b, z_c)])
}

#' @rdname transform_coords
#' @export
#' @examples
#' # coords <- data.table(Label="Si1_sym1", x_a=0.1, y_b=0.2, z_c=0.3)
#' # expand_transformed_coords(coords)
expand_transformed_coords <- function(transformed_coords, n_cells = 1) {
  if (is.null(transformed_coords) || nrow(transformed_coords) == 0) return(NULL)
  cell_indices <- as.data.table(expand.grid(dx = -n_cells:n_cells, dy = -n_cells:n_cells, dz = -n_cells:n_cells))
  expanded_coords <- rbindlist(lapply(1:nrow(cell_indices), function(i) {
    cell_shift <- cell_indices[i,]
    new_coords <- copy(transformed_coords)
    new_coords[, `:=`(Label = paste(Label, cell_shift$dx, cell_shift$dy, cell_shift$dz, sep = "_"),
                      x_a = x_a + cell_shift$dx, y_b = y_b + cell_shift$dy, z_c = z_c + cell_shift$dz)]
    return(new_coords)
  }))

  # Replicating notebook's use of dplyr::distinct
  temp_rounded <- copy(expanded_coords)
  temp_rounded[, `:=` (x_a_round = round(x_a, 5), y_b_round = round(y_b, 5), z_c_round = round(z_c, 5))]
  unique_coords <- unique(temp_rounded, by = c("x_a_round", "y_b_round", "z_c_round"))

  return(unique_coords[, .(Label, x_a, y_b, z_c)])
}
