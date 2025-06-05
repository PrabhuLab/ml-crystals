#' Read CIF Files from a List of Paths
#'
#' Reads multiple CIF files into a list of data.tables.
#' Each data.table represents one CIF file with lines as rows.
#'
#' @param cif_file_paths A character vector of paths to CIF files.
#' @return A list of data.tables. Each data.table has one column 'V1' containing lines from a CIF file.
#'         Returns an empty list if no files are found or readable.
#' @export
#' @examples
#' # temp_cif_file <- tempfile(fileext = ".cif")
#' # writeLines(c("data_test", "_cell_length_a 10"), temp_cif_file)
#' # cif_contents <- read_cif_files(temp_cif_file)
#' # unlink(temp_cif_file)
#' # print(cif_contents)
read_cif_files <- function(cif_file_paths) {
  if (length(cif_file_paths) == 0) {
    warning("No CIF file paths provided.")
    return(list())
  }

  cif_contents_list <- lapply(cif_file_paths, function(file_path) {
    tryCatch({
      # fread with sep="\n" reads each line into column V1.
      fread(file_path, sep = "\n", header = FALSE, strip.white = FALSE, quote="", data.table = TRUE)
    }, error = function(e) {
      warning(paste("Failed to read or process CIF file:", file_path, "-", e$message))
      NULL # Return NULL for files that failed to read
    })
  })

  successful_reads <- !sapply(cif_contents_list, is.null)
  cif_contents_list <- cif_contents_list[successful_reads]
  names(cif_contents_list) <- basename(cif_file_paths[successful_reads])

  return(cif_contents_list)
}


#' Process a Single CIF File Content
#'
#' Extracts crystallographic information, calculates geometric properties,
#' and determines bonding and coordination for a single CIF file's content.
#'
#' @param cif_content A data.table representing the content of one CIF file,
#'                    as returned by `fread(filepath, sep="\\n", header=FALSE)`.
#' @param bonding_method A character string specifying the bonding algorithm to use.
#'                       Options: "min_dist", "brunner", "hoppe". Default is "min_dist".
#' @param min_dist_delta Numeric, delta parameter for `minimum_distance` method.
#' @param brunner_delta Numeric, delta parameter for `brunner` method.
#' @param hoppe_bs_threshold Numeric, bond strength threshold for `hoppe` method.
#' @param hoppe_tolerance Numeric, tolerance for `hoppe` method convergence.
#' @param expand_n_cells Integer, number of unit cells to expand for distance calculations.
#' @return A data.table row (or list that can be row-bound) containing all
#'         extracted and calculated data for the CIF file. Returns NULL if essential
#'         data (e.g. coordinates or cell parameters) is missing.
#' @export
#' @examples
#' # dummy_cif_lines <- c(
#' #   "data_test", "_database_code_ '00000'", "_chemical_formula_sum 'H2 O1'",
#' #   "_cell_length_a   5.0", "_cell_length_b   5.0", "_cell_length_c   5.0",
#' #   "_cell_angle_alpha 90.0", "_cell_angle_beta  90.0", "_cell_angle_gamma 90.0",
#' #   "_space_group_symop_operation_xyz", "  'x, y, z'",
#' #   "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_foo", "_atom_site_bar",
#' #   "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy",
#' #   "O1 O X X 0.0 0.0 0.0 1.0" # Example line for coord parsing
#' # )
#' # cif_dt <- data.table(V1 = dummy_cif_lines)
#' # result <- process_single_cif_data(cif_dt)
#' # print(result)
process_single_cif_data <- function(cif_content,
                                    bonding_method = "min_dist",
                                    min_dist_delta = 0.1,
                                    brunner_delta = 0.0001,
                                    hoppe_bs_threshold = 0.5,
                                    hoppe_tolerance = 0.001,
                                    expand_n_cells = 1) {
  if (is.null(cif_content) || nrow(cif_content) == 0) {
    warning("Empty CIF content provided.")
    return(NULL)
  }

  # --- Extraction ---
  db_code <- extract_database_code(cif_content)
  chem_formula <- extract_chemical_formula(cif_content)
  struct_type <- extract_structure_type(cif_content)
  sg_name <- extract_space_group_name(cif_content)
  sg_number <- extract_space_group_number(cif_content)

  unit_cell_metrics <- extract_unit_cell_metrics(cif_content)
  if (is.null(unit_cell_metrics) || nrow(unit_cell_metrics) == 0 || any(is.na(unlist(unit_cell_metrics)))) {
    warning(paste("Essential unit cell metrics missing or NA for CIF. DB Code:", db_code %||% "N/A", ". Skipping detailed processing."))
    return(data.table(
        database_code = db_code, chemical_formula = chem_formula, structure_type = struct_type,
        space_group_name = sg_name, space_group_number = sg_number,
        error_message = "Missing unit cell metrics"
    ))
  }

  atomic_coords <- extract_atomic_coordinates(cif_content)
  if (is.null(atomic_coords) || nrow(atomic_coords) == 0) {
    warning(paste("Atomic coordinates missing or unparseable for CIF. DB Code:", db_code %||% "N/A", ". Skipping detailed processing."))
    return(data.table(
        database_code = db_code, chemical_formula = chem_formula, structure_type = struct_type,
        space_group_name = sg_name, space_group_number = sg_number,
        unit_cell_metrics = list(unit_cell_metrics),
        error_message = "Missing atomic coordinates"
    ))
  }

  sym_ops <- extract_symmetry_operations(cif_content)
  if (is.null(sym_ops) || nrow(sym_ops) == 0) {
      warning(paste("No symmetry operations found/parsed for CIF. Assuming P1 (identity only). DB Code:", db_code %||% "N/A"))
      sym_ops <- data.table(x="x", y="y", z="z") # Default to identity operation
  }

  # --- Processing Coordinates ---
  transformed_coords <- apply_symmetry_operations(atomic_coords, sym_ops)
  if (is.null(transformed_coords) || nrow(transformed_coords) == 0) {
      warning(paste("Failed to apply symmetry operations. DB Code:", db_code %||% "N/A", ". Using input coordinates."))
      transformed_coords <- atomic_coords # Fallback to original coordinates
  }

  # `coords_in_cell_and_neighbors` will be the search space for `calculate_distances`.
  # It includes atoms generated by symmetry within the primary cell (`transformed_coords`)
  # and their translations to neighboring cells.
  coords_in_cell_and_neighbors <- expand_transformed_coords(transformed_coords, n_cells = expand_n_cells)

  if (is.null(coords_in_cell_and_neighbors) || nrow(coords_in_cell_and_neighbors) == 0) {
      warning(paste("Failed to expand coordinates to neighboring cells. DB Code:", db_code %||% "N/A", ". Using only in-cell transformed coordinates."))
      coords_in_cell_and_neighbors <- transformed_coords # Fallback, might miss inter-cell bonds
  }

  # --- Calculate Properties ---
  # `atomic_coords` (asymmetric unit) are the reference atoms (Atom1)
  # `coords_in_cell_and_neighbors` are the atoms to search against (Atom2)
  distances <- calculate_distances(atomic_coords, coords_in_cell_and_neighbors, unit_cell_metrics)

  bonded_pairs <- NULL
  if (!is.null(distances) && nrow(distances) > 0) {
    if (bonding_method == "min_dist") {
      bonded_pairs <- minimum_distance(distances, delta = min_dist_delta)
    } else if (bonding_method == "brunner") {
      bonded_pairs <- brunner(distances, delta = brunner_delta)
    } else if (bonding_method == "hoppe") {
      bonded_pairs <- hoppe(distances, bond_strength_threshold = hoppe_bs_threshold, tolerance = hoppe_tolerance)
    } else {
      warning(paste("Unknown bonding method:", bonding_method, ". Defaulting to 'min_dist'."))
      bonded_pairs <- minimum_distance(distances, delta = min_dist_delta)
    }
  } else {
     warning(paste("No distances calculated or distances table is empty. Cannot determine bonds. DB Code:", db_code %||% "N/A"))
  }

  neighbor_counts <- calculate_neighbor_counts(bonded_pairs)

  bond_angles <- NULL
  if (!is.null(bonded_pairs) && nrow(bonded_pairs) > 0) {
    # For `calculate_angles`, `atomic_coords` defines the central atoms for angles.
    # `coords_in_cell_and_neighbors` is the comprehensive set to find all atom positions.
    bond_angles <- calculate_angles(bonded_pairs, atomic_coords, coords_in_cell_and_neighbors, unit_cell_metrics)
  }

  # --- Compile Results ---
  result_dt <- data.table(
    database_code = db_code,
    chemical_formula = chem_formula,
    structure_type = struct_type,
    space_group_name = sg_name,
    space_group_number = sg_number,
    unit_cell_metrics = list(unit_cell_metrics),
    atomic_coordinates_input = list(atomic_coords),
    symmetry_operations_used = list(sym_ops),
    expanded_coordinates_searched = list(coords_in_cell_and_neighbors),
    distances_calculated = list(distances),
    bonded_pairs_identified = list(bonded_pairs),
    neighbor_counts_calculated = list(neighbor_counts),
    bond_angles_calculated = list(bond_angles),
    error_message = NA_character_
  )
  return(result_dt)
}

# Helper for warning messages when db_code might be NULL
"%||%" <- function(a, b) if (!is.null(a)) a else b


#' Analyze Multiple CIF Files from a Folder or List of Paths
#'
#' Reads CIF files, processes each to extract crystallographic data,
#' calculate geometric properties, and determine bonding, then aggregates
#' the results into a single data.table.
#'
#' @param cif_input Either a character string path to a folder containing CIF files,
#'                  or a character vector of direct paths to CIF files.
#' @param pattern A glob pattern to match CIF files within the folder (e.g., "*.cif").
#'                Used only if `cif_input` is a folder path.
#' @param bonding_method Character, bonding algorithm: "min_dist", "brunner", "hoppe".
#' @param ... Additional parameters passed to `process_single_cif_data` (e.g., deltas for bonding methods).
#' @return A data.table where each row summarizes one CIF file.
#' @export
#' @examples
#' # temp_dir <- tempdir()
#' # cif1_path <- file.path(temp_dir, "test1.cif")
#' # cif_content_minimal <- c(
#' #   "data_test", "_database_code_ '00001'", "_chemical_formula_sum 'A1'",
#' #   "_cell_length_a 1", "_cell_length_b 1", "_cell_length_c 1",
#' #   "_cell_angle_alpha 90", "_cell_angle_beta 90", "_cell_angle_gamma 90",
#' #   "_space_group_symop_operation_xyz 'x,y,z'",
#' #   "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_foo", "_atom_site_bar",
#' #   "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy",
#' #   "A1 A X X 0.1 0.1 0.1 1.0"
#' # )
#' # writeLines(cif_content_minimal, cif1_path)
#' # results <- analyze_cif_files(temp_dir, pattern = "*.cif")
#' # print(results)
#' # unlink(cif1_path); unlink(temp_dir, recursive=TRUE)
analyze_cif_files <- function(cif_input, pattern = "*.cif", bonding_method = "min_dist", ...) {
  cif_file_paths <- character(0)
  if (length(cif_input) == 1 && dir.exists(cif_input)) {
    cif_file_paths <- list.files(path = cif_input, pattern = pattern, full.names = TRUE)
    if (length(cif_file_paths) == 0) {
      warning(paste("No files matching pattern", pattern, "found in folder", cif_input))
      return(data.table())
    }
  } else if (is.character(cif_input) && all(file.exists(cif_input))) {
    cif_file_paths <- cif_input
  } else {
    stop("cif_input must be a valid folder path or a vector of valid file paths.")
  }

  message(paste("Found", length(cif_file_paths), "CIF files to process."))

  cif_contents_list <- read_cif_files(cif_file_paths)
  if (length(cif_contents_list) == 0) {
    warning("No CIF files could be read successfully.")
    return(data.table())
  }

  all_results_list <- lapply(names(cif_contents_list), function(file_name) {
    cif_content <- cif_contents_list[[file_name]]
    message(paste("Processing CIF file:", file_name))

    single_result <- process_single_cif_data(cif_content,
                                             bonding_method = bonding_method,
                                             ...)
    if (!is.null(single_result)) {
      single_result[, source_file := file_name] # Add source file name
      return(single_result)
    }
    return(NULL) # Return NULL if processing failed significantly
  })

  all_results_list <- all_results_list[!sapply(all_results_list, is.null)]

  if (length(all_results_list) == 0) {
    warning("No CIF files were processed successfully.")
    return(data.table())
  }

  # Combine all results; fill=TRUE handles cases where some files might have
  # error messages and not all data columns.
  main_table <- rbindlist(all_results_list, fill = TRUE)

  return(main_table)
}
