#' Generic Value Extractor (Internal)
#'
#' @description Extracts a single value from CIF content based on a matching text pattern.
#' @param cif_content A data.table containing the lines of a CIF file.
#' @param pattern The text pattern (e.g., "_database_code_") to search for.
#' @param remove_pattern A boolean indicating whether to remove the search pattern from the result.
#' @return A character string of the cleaned value, or NA if not found.
#' @noRd
extract_value <- function(cif_content, pattern, remove_pattern = TRUE) {
  # Use base R's grepl for more robust pattern matching
  matching_lines <- cif_content[grepl(pattern, V1, fixed = TRUE)]

  if (nrow(matching_lines) > 0) {
    # Take the first matching line
    value <- matching_lines$V1[1]

    if (remove_pattern) {
      value <- gsub(pattern, "", value, fixed = TRUE)
    }
    value <- gsub("'", "", value)
    value <- trimws(value)
    return(value)
  } else {
    return(NA)
  }
}

#' @title Extract Database Code from CIF Content
#' @description Extracts the database code (e.g., ICSD number) from CIF data.
#' @param cif_content A `data.table` where each row is a line from a CIF file.
#' @return A character string of the database code, or `NA` if not found.
#' @family extractors
#' @export
#' @examples
#' cif_file <- system.file("extdata", "ICSD422.cif", package = "crysmal")
#' if (file.exists(cif_file)) {
#'   cif_content <- data.table::fread(cif_file, sep = "\n", header = FALSE)
#'   extract_database_code(cif_content)
#' }
extract_database_code <- function (cif_content) {
  extract_value(cif_content, "_database_code_")
}

#' @title Extract Chemical Formula from CIF Content
#' @param cif_content A `data.table` where each row is a line from a CIF file.
#' @return A character string of the chemical formula, or `NA` if not found.
#' @family extractors
#' @export
extract_chemical_formula <- function(cif_content) {
  extract_value(cif_content, "_chemical_formula_sum")
}

#' @title Extract Structure Type from CIF Content
#' @param cif_content A `data.table` where each row is a line from a CIF file.
#' @return A character string of the structure type, or `NA` if not found.
#' @family extractors
#' @export
extract_structure_type <- function(cif_content) {
  extract_value(cif_content, "_chemical_name_structure_type")
}

#' @title Extract Space Group Name from CIF Content
#' @param cif_content A `data.table` where each row is a line from a CIF file.
#' @return A character string of the space group name, or `NA` if not found.
#' @family extractors
#' @export
extract_space_group_name <- function(cif_content) {
  extract_value(cif_content, "_space_group_name_H-M_alt")
}

#' @title Extract Space Group Number from CIF Content
#' @param cif_content A `data.table` where each row is a line from a CIF file.
#' @return A character string of the space group number, or `NA` if not found.
#' @family extractors
#' @export
extract_space_group_number <- function(cif_content) {
  extract_value(cif_content, "_space_group_IT_number")
}

#' @title Extract Unit Cell Metrics
#' @description Parses the unit cell parameters (lengths a, b, c and angles
#'   alpha, beta, gamma) and their associated standard uncertainties from CIF
#'   content.
#' @param cif_content A `data.table` containing the lines of a CIF file.
#' @return A one-row `data.table` with columns for each cell parameter and its
#'   error.
#' @family extractors
#' @export
#' @examples
#' cif_file <- system.file("extdata", "ICSD422.cif", package = "crysmal")
#' if (file.exists(cif_file)) {
#'   cif_content <- data.table::fread(cif_file, sep = "\n", header = FALSE)
#'   metrics <- extract_unit_cell_metrics(cif_content)
#'   print(metrics)
#' }
extract_unit_cell_metrics <- function(cif_content) {
  cell_parameters <- c(
    "_cell_length_a",
    "_cell_length_b",
    "_cell_length_c",
    "_cell_angle_alpha",
    "_cell_angle_beta",
    "_cell_angle_gamma"
  )
  values <- list()
  errors <- list()

  scale_error <- function(value_str, error_str) {
    if (is.na(error_str) || error_str == "")
      return(NA_real_)
    decimal_pos <- regexpr("\\.", value_str)
    if (decimal_pos == -1) {
      as.numeric(error_str)
    }
    else {
      as.numeric(error_str) * 10^-(nchar(value_str) - decimal_pos)
    }
  }

  for (param in cell_parameters) {
    line <- cif_content[grepl(param, V1, fixed = TRUE)]$V1
    if (length(line) > 0) {
      match <- stringr::str_match(line[1], "\\s+([0-9\\.]+)(?:\\(([0-9]+)\\))?")
      if (!is.na(match[1, 1])) {
        value_str <- match[1, 2]
        error_str <- match[1, 3]
        values[[param]] <- as.numeric(value_str)
        errors[[paste0(param, "_error")]] <- scale_error(value_str, error_str)
      } else {
        values[[param]] <- NA_real_
        errors[[paste0(param, "_error")]] <- NA_real_
      }
    } else {
      values[[param]] <- NA_real_
      errors[[paste0(param, "_error")]] <- NA_real_
    }
  }
  return(as.data.table(c(values, errors)))
}

#' @title Extract Atomic Coordinates
#' @description Parses the atomic site information to extract the label, fractional
#'   coordinates (x, y, z), Wyckoff symbol, multiplicity, occupancy, and their standard
#'   uncertainties for each asymmetric atom in the unit cell.
#' @param cif_content A `data.table` containing the lines of a CIF file.
#' @return A `data.table` with atomic coordinate data. Returns `NULL` if not found.
#' @family extractors
#' @export
extract_atomic_coordinates <- function(cif_content) {
  # --- 1. Find the start of the atom site loop ---
  first_header_line_idx <- grep("^_atom_site_fract_x", cif_content$V1)
  if (is.na(first_header_line_idx)) {
    first_header_line_idx <- grep("^_atom_site_label", cif_content$V1)
    if (is.na(first_header_line_idx))
      return(NULL)
  }
  loop_start_line_idx <- max(grep("^loop_", cif_content$V1[1:first_header_line_idx]))
  if (is.infinite(loop_start_line_idx))
    return(NULL)

  # --- 2. Read the headers and find column indices ---
  line_indices <- (loop_start_line_idx + 1):nrow(cif_content)
  headers <- character()
  first_data_line_idx <- 0
  for (i in line_indices) {
    line <- cif_content$V1[i]
    if (startsWith(line, "_")) {
      headers <- c(headers, trimws(line))
    } else {
      first_data_line_idx <- i
      break
    }
  }
  tags_to_find <- c(
    label = "_atom_site_label",
    x = "_atom_site_fract_x",
    y = "_atom_site_fract_y",
    z = "_atom_site_fract_z",
    occupancy = "_atom_site_occupancy",
    wyckoff = "_atom_site_Wyckoff_symbol",
    multiplicity = "_atom_site_symmetry_multiplicity"
  )
  col_indices <- sapply(tags_to_find, function(tag) {
    idx <- which(headers == tag)
    if (length(idx) == 0)
      NA
    else
      idx
  })
  if (anyNA(col_indices[c("label", "x", "y", "z")])) {
    warning("CIF file is missing essential atom site tags (_label, _fract_x, _y, _z).")
    return(NULL)
  }

  # --- 3. Find the end of the data block ---
  end_candidates <- c(grep("^loop_|^_|^#", cif_content$V1[first_data_line_idx:nrow(cif_content)]),
                      grep("^\\s*$", cif_content$V1[first_data_line_idx:nrow(cif_content)]))
  last_data_line_idx <- if (length(end_candidates) > 0) {
    first_data_line_idx + min(end_candidates) - 2
  } else {
    nrow(cif_content)
  }
  if (first_data_line_idx > last_data_line_idx)
    return(NULL)

  # --- 4. Read data with fread ---
  data_lines <- cif_content$V1[first_data_line_idx:last_data_line_idx]
  atom_data <- fread(
    text = paste(data_lines, collapse = "\n"),
    header = FALSE,
    sep = "auto",
    quote = ""
  )

  # --- 5. VECTORIZED parsing of coordinates and errors ---
  parse_vector_with_error <- function(coord_vector) {
    matches <- stringr::str_match(coord_vector, "([0-9\\.\\-]+)(?:\\(([0-9]+)\\))?")
    value_str <- matches[, 2]
    error_str <- matches[, 3]
    decimal_pos <- regexpr("\\.", value_str)
    decimal_places <- ifelse(decimal_pos == -1, 0, nchar(value_str) - decimal_pos)
    scaled_error <- as.numeric(error_str) * 10^(-decimal_places)
    return(list(value = as.numeric(value_str), error = scaled_error))
  }

  x_data <- parse_vector_with_error(atom_data[[col_indices["x"]]])
  y_data <- parse_vector_with_error(atom_data[[col_indices["y"]]])
  z_data <- parse_vector_with_error(atom_data[[col_indices["z"]]])

  # Parse occupancy, assuming 1.0 if the column is missing
  occ_data <- if (!is.na(col_indices["occupancy"])) {
    parse_vector_with_error(atom_data[[col_indices["occupancy"]]])
  } else {
    list(value = rep(1.0, nrow(atom_data)),
         error = rep(NA_real_, nrow(atom_data)))
  }


  # --- 6. Assemble the final data table ---
  wyckoff_multiplicity <- if (!is.na(col_indices["multiplicity"])) {
    as.numeric(atom_data[[col_indices["multiplicity"]]])
  } else {
    rep(NA_real_, nrow(atom_data))
  }

  wyckoff_symbol <- if (!is.na(col_indices["wyckoff"])) {
    as.character(atom_data[[col_indices["wyckoff"]]])
  } else {
    rep(NA_character_, nrow(atom_data))
  }

  atomic_coordinates <- data.table(
    Label = atom_data[[col_indices["label"]]],
    WyckoffSymbol = wyckoff_symbol,
    WyckoffMultiplicity = wyckoff_multiplicity,
    Occupancy = occ_data$value,
    OccupancyError = occ_data$error,
    x_a = x_data$value,
    y_b = y_data$value,
    z_c = z_data$value,
    x_error = x_data$error,
    y_error = y_data$error,
    z_error = z_data$error
  )

  return(atomic_coordinates)
}

#' @title Extract Symmetry Operations
#' @description Parses the symmetry operation definitions from the CIF content.
#' @param cif_content A `data.table` containing the lines of a CIF file.
#' @return A `data.table` with symmetry operations. Returns `NULL` if not found.
#' @family extractors
#' @export
#' @examples
#' cif_file <- system.file("extdata", "ICSD422.cif", package = "crysmal")
#' if (file.exists(cif_file)) {
#'   cif_content <- data.table::fread(cif_file, sep = "\n", header = FALSE)
#'   sym_ops <- extract_symmetry_operations(cif_content)
#'   print(sym_ops)
#' }
extract_symmetry_operations <- function(cif_content) {
  # --- 1. Find the start of the symmetry loop ---
  symop_tag <- "_space_group_symop_operation_xyz"
  first_header_line_idx <- grep(symop_tag, cif_content$V1)[1]
  if (is.na(first_header_line_idx)) {
    symop_tag <- "_symmetry_equiv_pos_as_xyz"
    first_header_line_idx <- grep(symop_tag, cif_content$V1)[1]
    if (is.na(first_header_line_idx))
      return(NULL)
  }
  loop_start_line_idx <- max(grep("^loop_", cif_content$V1[1:first_header_line_idx]))
  if (is.infinite(loop_start_line_idx))
    return(NULL)

  # --- 2. Find the range of data lines ---
  line_indices <- (loop_start_line_idx + 1):nrow(cif_content)
  first_data_line_idx <- 0
  for (i in line_indices) {
    if (!startsWith(cif_content$V1[i], "_")) {
      first_data_line_idx <- i
      break
    }
  }
  if (first_data_line_idx == 0)
    return(NULL)

  end_candidates <- c(grep("^loop_|^_|^#", cif_content$V1[first_data_line_idx:nrow(cif_content)]),
                      grep("^\\s*$", cif_content$V1[first_data_line_idx:nrow(cif_content)]))
  last_data_line_idx <- if (length(end_candidates) > 0) {
    first_data_line_idx + min(end_candidates) - 2
  } else {
    nrow(cif_content)
  }
  if (first_data_line_idx > last_data_line_idx)
    return(NULL)

  # --- 3. Extract and parse raw data lines ---
  data_lines <- cif_content$V1[first_data_line_idx:last_data_line_idx]
  cleaned_lines <- trimws(gsub("^'|^\"|'$|\"$", "", trimws(sub(
    "^[0-9]+\\s+", "", data_lines
  ))))
  symmetry_matrix <- stringr::str_split_fixed(cleaned_lines, ",", n = 3)

  symmetry_dt <- data.table(
    x = trimws(symmetry_matrix[, 1]),
    y = trimws(symmetry_matrix[, 2]),
    z = trimws(symmetry_matrix[, 3])
  )

  return(symmetry_dt)
}
