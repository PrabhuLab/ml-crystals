#' Read CIF Files from a List of Paths
#' @param cif_file_paths A character vector of paths to CIF files.
#' @return A list of data.tables.
#' @export
read_cif_files <- function(cif_file_paths) {
  if (length(cif_file_paths) == 0) return(list())
  cif_contents_list <- lapply(cif_file_paths, function(file_path) {
    tryCatch({
      fread(file_path, sep = "\n", header = FALSE, strip.white = FALSE, quote="", data.table = TRUE)
    }, error = function(e) {
      warning(paste("Failed to read or process CIF file:", file_path, "-", e$message))
      NULL
    })
  })
  successful_reads <- !sapply(cif_contents_list, is.null)
  cif_contents_list <- cif_contents_list[successful_reads]
  names(cif_contents_list) <- basename(cif_file_paths[successful_reads])
  return(cif_contents_list)
}

"%||%" <- function(a, b) if (!is.null(a)) a else b

#' Process a Single CIF File Content
#' @param cif_content A data.table representing the content of one CIF file.
#' @param bonding_method A character string specifying the bonding algorithm ("min_dist", "brunner", "hoppe").
#' @param ... Additional parameters passed to the bonding function.
#' @return A data.table row containing all extracted and calculated data for the CIF file.
#' @export
process_single_cif_data <- function(cif_content, bonding_method = "min_dist", ...) {
  if (is.null(cif_content) || nrow(cif_content) == 0) return(NULL)
  db_code <- extract_database_code(cif_content)
  # --- Extraction from Notebook Logic ---
  unit_cell_metrics <- extract_unit_cell_metrics(cif_content)
  if (is.null(unit_cell_metrics) || nrow(unit_cell_metrics) == 0 || any(is.na(unlist(unit_cell_metrics)))) {
    return(data.table(database_code = db_code, error_message = "Missing unit cell metrics"))
  }
  atomic_coords <- extract_atomic_coordinates(cif_content)
  if (is.null(atomic_coords) || nrow(atomic_coords) == 0) {
    return(data.table(database_code = db_code, error_message = "Missing atomic coordinates"))
  }
  sym_ops <- extract_symmetry_operations(cif_content)

  # --- Processing ---
  transformed_coords <- apply_symmetry_operations(atomic_coords, sym_ops)
  expanded_coords <- expand_transformed_coords(transformed_coords)
  distances <- calculate_distances(atomic_coords, expanded_coords, unit_cell_metrics)

  bonding_args <- list(distances = distances, ...)
  bonded_pairs <- switch(bonding_method,
                         "min_dist" = do.call(minimum_distance, bonding_args),
                         "brunner" = do.call(brunner, bonding_args),
                         "hoppe" = do.call(hoppe, bonding_args),
                         { warning("Unknown bonding method. Defaulting to 'min_dist'.");
                           do.call(minimum_distance, bonding_args) })

  neighbor_counts <- calculate_neighbor_counts(bonded_pairs)
  bond_angles <- calculate_angles(bonded_pairs, atomic_coords, expanded_coords, unit_cell_metrics)

  # --- Compile Results with clear and consistent names ---
  data.table(
    database_code = db_code,
    chemical_formula = extract_chemical_formula(cif_content),
    structure_type = extract_structure_type(cif_content),
    space_group_name = extract_space_group_name(cif_content),
    space_group_number = extract_space_group_number(cif_content),
    unit_cell_metrics = list(unit_cell_metrics),
    atomic_coordinates_input = list(atomic_coords),      # Original coords
    symmetry_operations = list(sym_ops),                 # The symops found
    atomic_coordinates_transformed = list(transformed_coords), # Coords after symm
    atomic_coordinates_expanded = list(expanded_coords), # Coords after expanding cells
    distances_calculated = list(distances),
    bonded_pairs_identified = list(bonded_pairs),
    neighbor_counts_calculated = list(neighbor_counts),
    bond_angles_calculated = list(bond_angles),
    error_message = NA_character_
  )
}

#' Analyze Multiple CIF Files
#' @param cif_input Path to a folder or a vector of file paths.
#' @param pattern Glob pattern for files if `cif_input` is a folder.
#' @param bonding_method Bonding algorithm to use.
#' @param ... Additional parameters passed to `process_single_cif_data`.
#' @return A data.table where each row summarizes one CIF file.
#' @export
analyze_cif_files <- function(cif_input, pattern = "*.cif", bonding_method = "min_dist", ...) {
  cif_file_paths <- character(0)
  if (length(cif_input) == 1 && dir.exists(cif_input)) {
    cif_file_paths <- list.files(path = cif_input, pattern = pattern, full.names = TRUE)
  } else if (is.character(cif_input) && all(file.exists(cif_input))) {
    cif_file_paths <- cif_input
  } else {
    stop("cif_input must be a valid folder path or a vector of valid file paths.")
  }
  if (length(cif_file_paths) == 0) {
    warning("No CIF files found or specified.")
    return(data.table())
  }
  message(paste("Found", length(cif_file_paths), "CIF files to process."))
  cif_contents_list <- read_cif_files(cif_file_paths)

  all_results_list <- lapply(names(cif_contents_list), function(file_name) {
    message(paste("Processing CIF file:", file_name))
    single_result <- process_single_cif_data(cif_contents_list[[file_name]], bonding_method = bonding_method, ...)
    if (!is.null(single_result)) single_result[, source_file := file_name]
    return(single_result)
  })

  main_table <- rbindlist(Filter(Negate(is.null), all_results_list), fill = TRUE)
  return(main_table)
}
