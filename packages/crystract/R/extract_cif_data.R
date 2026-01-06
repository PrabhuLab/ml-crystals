#' Generic Value Extractor (Internal)
#'
#' @description Extracts a single value from CIF content based on a matching text pattern.
#' @param cif_content A data.table containing the lines of a CIF file.
#' @param pattern The text pattern (e.g., "_database_code_") to search for.
#' @param remove_pattern A boolean indicating whether to remove the search pattern from the result.
#' @return A character string of the cleaned value, or NA if not found.
#' @noRd
extract_value <- function(cif_content, pattern, remove_pattern = TRUE) {
  # fixed = TRUE automatically handles indented lines as it searches for the substring anywhere
  matching_lines <- cif_content[grepl(pattern, V1, fixed = TRUE)]
  if (nrow(matching_lines) > 0) {
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
#' @description Extracts the database code identifier (e.g., from `_database_code_`) from the CIF.
#' @param cif_content A `data.table` containing the lines of a CIF file.
#' @return A character string of the database code, or `NA` if not found.
#' @family extractors
#' @export
extract_database_code <- function (cif_content) {
  extract_value(cif_content, "_database_code_")
}

#' @title Extract Chemical Formula from CIF Content
#' @description Extracts the chemical sum formula (e.g., from `_chemical_formula_sum`).
#' @param cif_content A `data.table` containing the lines of a CIF file.
#' @return A character string of the chemical formula, or `NA` if not found.
#' @family extractors
#' @export
extract_chemical_formula <- function(cif_content) {
  extract_value(cif_content, "_chemical_formula_sum")
}

#' @title Extract Structure Type from CIF Content
#' @description Extracts the structure type name (e.g., from `_chemical_name_structure_type`).
#' @param cif_content A `data.table` containing the lines of a CIF file.
#' @return A character string of the structure type, or `NA` if not found.
#' @family extractors
#' @export
extract_structure_type <- function(cif_content) {
  extract_value(cif_content, "_chemical_name_structure_type")
}

#' @title Extract Space Group Name from CIF Content
#' @description Extracts the Hermann-Mauguin space group name (e.g., from `_space_group_name_H-M_alt`).
#' @param cif_content A `data.table` containing the lines of a CIF file.
#' @return A character string of the space group name, or `NA` if not found.
#' @family extractors
#' @export
extract_space_group_name <- function(cif_content) {
  extract_value(cif_content, "_space_group_name_H-M_alt")
}

#' @title Extract Space Group Number from CIF Content
#' @description Extracts the International Tables space group number (e.g., from `_space_group_IT_number`).
#' @param cif_content A `data.table` containing the lines of a CIF file.
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
    } else {
      as.numeric(error_str) * 10^-(nchar(value_str) - decimal_pos)
    }
  }

  for (param in cell_parameters) {
    # fixed = TRUE handles indented lines automatically
    line <- cif_content[grepl(param, V1, fixed = TRUE)]$V1
    if (length(line) > 0) {
      match <- stringr::str_match(line[1], "\\s+([0-9\\.]+)(?:\\(([0-9]+)\\))?")
      if (!is.na(match[1, 1])) {
        values[[param]] <- as.numeric(match[1, 2])
        errors[[paste0(param, "_error")]] <- scale_error(match[1, 2], match[1, 3])
      } else {
        values[[param]] <- NA_real_
        errors[[paste0(param, "_error")]] <- NA_real_
      }
    } else {
      values[[param]] <- NA_real_
      errors[[paste0(param, "_error")]] <- NA_real_
    }
  }
  return(data.table::as.data.table(c(values, errors)))
}

#' @title Extract Atomic Coordinates
#' @description Parses atomic site info, including labels, fractional coordinates,
#'   and other properties. It efficiently normalizes atom labels to a
#'   'SymbolNumber' format and validates them against the chemical formula.
#' @param cif_content A `data.table` containing the lines of a CIF file.
#' @param chemical_formula The chemical formula string for validation.
#' @return A `data.table` with atomic coordinate data, or `NULL` if not found.
#' @family extractors
#' @export
extract_atomic_coordinates <- function(cif_content, chemical_formula = NA) {
  # --- 1. Find Header Block ---
  # Regex: Matches start of line (^), optional whitespace (\\s*), then the tag
  first_header_line_idx <- grep("^\\s*_atom_site_fract_x", cif_content$V1)
  if (length(first_header_line_idx) == 0) {
    first_header_line_idx <- grep("^\\s*_atom_site_label", cif_content$V1)
    if (length(first_header_line_idx) == 0)
      return(NULL)
  }
  first_header_line_idx <- first_header_line_idx[1]

  # --- 2. Find Loop Start ---
  # Regex: Matches start of line (^), optional whitespace (\\s*), then "loop_"
  loop_start_line_idx <- max(grep("^\\s*loop_", cif_content$V1[1:first_header_line_idx]))
  if (is.infinite(loop_start_line_idx))
    return(NULL)

  # --- 3. Parse Headers ---
  line_indices <- (loop_start_line_idx + 1):nrow(cif_content)
  headers <- character()
  first_data_line_idx <- 0

  for (i in line_indices) {
    line <- trimws(cif_content$V1[i])
    if (line == "" || startsWith(line, "#"))
      next
    # trimws ensures this check works even if indented
    if (startsWith(line, "_")) {
      headers <- c(headers, line)
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

  if (anyNA(col_indices[c("label", "x", "y", "z")]))
    return(NULL)

  # --- 4. Determine Data Block Range ---
  # Regex: Find the next loop, new tag (_), or comment (#), allowing indentation
  end_candidates <- c(
    grep("^\\s*loop_|^\\s*_|^\\s*#", cif_content$V1[first_data_line_idx:nrow(cif_content)]),
    grep("^\\s*$", cif_content$V1[first_data_line_idx:nrow(cif_content)])
  )

  last_data_line_idx <- if (length(end_candidates) > 0)
    first_data_line_idx + min(end_candidates) - 2
  else
    nrow(cif_content)

  if (first_data_line_idx > last_data_line_idx)
    return(NULL)

  # --- 5. Read Data ---
  data_lines <- cif_content$V1[first_data_line_idx:last_data_line_idx]

  # quote = "" prevents crashes on unquoted atoms like Ag or N
  atom_data <- data.table::fread(
    text = paste(data_lines, collapse = "\n"),
    header = FALSE,
    sep = "auto",
    quote = ""
  )

  # Clean strings (remove manual quotes if present)
  clean_str <- function(col_data) {
    gsub("^'|'$|^\"|\"$", "", as.character(col_data))
  }

  parse_vector_with_error <- function(coord_vector) {
    coord_vector <- clean_str(coord_vector)
    # Group 1: Number (allows e/E notation)
    # Group 2: Error (optional digits in parens)
    matches <- stringr::str_match(coord_vector,
                                  "([0-9\\.\\-]+(?:[eE][+\\-]?[0-9]+)?)(?:\\(([0-9]+)\\))?")
    value_str <- matches[, 2]
    error_str <- matches[, 3]

    # Handle numeric conversion
    values <- as.numeric(value_str)

    # Handle error scaling
    errors <- mapply(function(val_s, err_s) {
      if (is.na(err_s) || err_s == "")
        return(0) # Default to 0 if NA
      if (is.na(val_s))
        return(NA_real_)

      # Check if scientific notation is present in the VALUE string
      if (grepl("[eE]", val_s)) {
        # Scientific notation logic:
        # If val is 4e-4 and err is 2 (representing 0.0002), we need to adjust.
        # Standard CIF practice for error in parens usually implies the last significant digit.
        # However, purely numerical parsing is safer if not strictly defined.
        # Fallback: Treat error as raw numeric if scientific notation is used,
        # or assume CIF standard doesn't use sci-notation with parenthesis errors often.
        # Simplest valid approach:
        return(as.numeric(err_s)) # Returns raw number, often effectively 0 if not scaled
      }

      decimal_pos <- regexpr("\\.", val_s)
      decimal_places <- ifelse(decimal_pos == -1, 0, nchar(val_s) - decimal_pos)
      as.numeric(err_s) * 10^(-decimal_places)
    }, value_str, error_str)

    return(list(value = values, error = errors))
  }

  # Safely extract columns
  safe_extract <- function(idx) {
    if (idx > ncol(atom_data))
      return(rep(NA, nrow(atom_data)))
    return(atom_data[[idx]])
  }

  x_data <- parse_vector_with_error(safe_extract(col_indices["x"]))
  y_data <- parse_vector_with_error(safe_extract(col_indices["y"]))
  z_data <- parse_vector_with_error(safe_extract(col_indices["z"]))

  occ_data <- if (!is.na(col_indices["occupancy"]))
    parse_vector_with_error(safe_extract(col_indices["occupancy"]))
  else
    list(value = rep(1.0, nrow(atom_data)),
         error = rep(NA_real_, nrow(atom_data)))

  atomic_coordinates <- data.table::data.table(
    Label = clean_str(safe_extract(col_indices["label"])),
    WyckoffSymbol = if (!is.na(col_indices["wyckoff"]))
      clean_str(safe_extract(col_indices["wyckoff"]))
    else
      NA_character_,
    WyckoffMultiplicity = if (!is.na(col_indices["multiplicity"]))
      as.numeric(clean_str(safe_extract(col_indices["multiplicity"])))
    else
      NA_real_,
    Occupancy = occ_data$value,
    OccupancyError = occ_data$error,
    x_a = x_data$value,
    y_b = y_data$value,
    z_c = z_data$value,
    x_error = x_data$error,
    y_error = y_data$error,
    z_error = z_data$error
  )

  # Remove rows that failed to parse (NAs in coordinates)
  atomic_coordinates <- atomic_coordinates[!is.na(x_a) &
                                             !is.na(y_b) &
                                             !is.na(z_c)]
  if (nrow(atomic_coordinates) == 0)
    return(NULL)

  # --- 6. Heuristic Label Correction ---
  valid_elements <- c(
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og"
  )

  base_symbols <- stringr::str_extract(atomic_coordinates$Label, "^[A-Za-z]+")
  unique_base_symbols <- unique(base_symbols)
  non_standard_symbols <- setdiff(unique_base_symbols, valid_elements)
  corrections_made <- list()

  if (length(non_standard_symbols) > 0) {
    elements_in_formula <- if (!is.na(chemical_formula))
      unique(stringr::str_extract_all(chemical_formula, "[A-Z][a-z]?")[[1]])
    else
      character(0)
    for (sym in non_standard_symbols) {
      corrected_sym <- NA_character_
      if (toupper(sym) %in% c("OW", "WAT", "OH")) {
        if ("O" %in% elements_in_formula ||
            length(elements_in_formula) == 0)
          corrected_sym <- "O"
      } else {
        prefix2 <- substr(sym, 1, 2)
        prefix1 <- substr(sym, 1, 1)
        if (prefix2 %in% valid_elements &&
            (prefix2 %in% elements_in_formula ||
             length(elements_in_formula) == 0)) {
          corrected_sym <- prefix2
        } else if (prefix1 %in% valid_elements &&
                   (prefix1 %in% elements_in_formula ||
                    length(elements_in_formula) == 0)) {
          corrected_sym <- prefix1
        }
      }
      if (!is.na(corrected_sym)) {
        indices_to_fix <- which(base_symbols == sym)
        atomic_coordinates[indices_to_fix, Label := sub("^[A-Za-z]+", corrected_sym, Label)]
        corrections_made[[sym]] <- corrected_sym
      }
    }
  }

  if (length(corrections_made) > 0) {
    correction_messages <- sapply(names(corrections_made), function(key)
      paste0("'", key, "' -> '", corrections_made[[key]], "'"))
    warning(paste0(
      "Corrected non-standard atom labels: ",
      paste(correction_messages, collapse = ", ")
    ))
  }

  # Normalize Labels
  symbols <- sub("^([A-Z][a-z]?).*", "\\1", atomic_coordinates$Label, perl = TRUE)
  numbers <- stringr::str_remove_all(atomic_coordinates$Label, "[^0-9]")
  atomic_coordinates[, Label := paste0(symbols, numbers)]

  return(atomic_coordinates)
}

#' @title Extract Symmetry Operations
#' @description Parses the symmetry operation definitions from the CIF content.
#'   This function uses manual token parsing to avoid issues with inconsistent quoting and spacing in CIF files.
#' @param cif_content A `data.table` containing the lines of a CIF file.
#' @return A `data.table` with symmetry operations. Returns `NULL` if not found.
#' @family extractors
#' @export
extract_symmetry_operations <- function(cif_content) {
  # --- 1. Define valid tags ---
  target_tags <- c("_space_group_symop_operation_xyz",
                   "_symmetry_equiv_pos_as_xyz")
  header_matches <- grep(paste(target_tags, collapse = "|"), cif_content$V1)
  if (length(header_matches) == 0)
    return(NULL)
  first_header_line_idx <- header_matches[1]

  # --- 2. Find Loop Start ---
  # Regex: Matches start of line (^), optional whitespace (\\s*), then "loop_"
  loop_start_line_idx <- max(grep("^\\s*loop_", cif_content$V1[1:first_header_line_idx]))
  if (is.infinite(loop_start_line_idx))
    return(NULL)

  # --- 3. Parse Headers to find Index ---
  line_indices <- (loop_start_line_idx + 1):nrow(cif_content)
  headers <- character()
  first_data_line_idx <- 0
  for (i in line_indices) {
    line <- trimws(cif_content$V1[i])
    if (line == "" || startsWith(line, "#"))
      next
    # trimws ensures this check works even if indented
    if (startsWith(line, "_")) {
      headers <- c(headers, line)
    } else {
      first_data_line_idx <- i
      break
    }
  }

  xyz_col_idx <- which(headers %in% target_tags)
  if (length(xyz_col_idx) == 0)
    return(NULL)
  xyz_col_idx <- xyz_col_idx[1]

  # --- 4. Determine Data Range ---
  # Regex: Find the next loop, new tag (_), or comment (#), allowing indentation
  end_candidates <- c(
    grep("^\\s*loop_|^\\s*_|^\\s*#", cif_content$V1[first_data_line_idx:nrow(cif_content)]),
    grep("^\\s*$", cif_content$V1[first_data_line_idx:nrow(cif_content)])
  )
  last_data_line_idx <- if (length(end_candidates) > 0)
    first_data_line_idx + min(end_candidates) - 2
  else
    nrow(cif_content)
  if (first_data_line_idx > last_data_line_idx)
    return(NULL)

  # --- 5. Manual Line Parsing ---
  data_lines <- cif_content$V1[first_data_line_idx:last_data_line_idx]

  extract_nth_token <- function(line, n) {
    line <- trimws(line)
    pattern <- "(?:'[^']*'|\"[^\"]*\"|\\S+)"
    matches <- gregexpr(pattern, line)
    tokens <- regmatches(line, matches)[[1]]
    if (length(tokens) < n)
      return(NA)
    return(tokens[n])
  }

  raw_strings <- sapply(data_lines, extract_nth_token, n = xyz_col_idx)

  # Filter out failures
  raw_strings <- raw_strings[!is.na(raw_strings)]
  if (length(raw_strings) == 0)
    return(NULL)

  # --- 6. Aggressive Cleaning ---
  cleaned_strings <- gsub("\\s+", "", raw_strings)
  cleaned_strings <- gsub("['\"]", "", cleaned_strings)
  cleaned_strings <- tolower(cleaned_strings)

  symmetry_matrix <- stringr::str_split_fixed(cleaned_strings, ",", n = 3)
  symmetry_dt <- data.table::data.table(x = symmetry_matrix[, 1], y = symmetry_matrix[, 2], z = symmetry_matrix[, 3])

  symmetry_dt <- symmetry_dt[x != "" & y != "" & z != ""]

  if (nrow(symmetry_dt) == 0)
    return(NULL)

  return(symmetry_dt)
}
