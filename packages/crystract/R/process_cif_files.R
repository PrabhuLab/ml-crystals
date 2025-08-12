#' @title Read CIF Files into Memory
#' @description Reads one or more CIF files from disk and loads each into a
#'   `data.table`.
#' @param file_paths A character vector of paths to the CIF files.
#' @return A list of `data.table` objects.
#' @export
#' @examples
#' cif_file <- system.file("extdata", "ICSD422.cif", package = "crysmal")
#' if (file.exists(cif_file)) {
#'   cif_data_list <- read_cif_files(cif_file)
#'   print(cif_data_list[[1]][1:5, ])
#' }
read_cif_files <- function(file_paths) {
  # Set the names of the list to the base filenames
  cif_list <- lapply(
    file_paths,
    fread,
    sep = "\n",
    header = FALSE,
    strip.white = FALSE
  )
  names(cif_list) <- basename(file_paths)
  return(cif_list)
}

#' @title Process a Single CIF Data Object
#' @description This function orchestrates the entire analysis pipeline for a
#'   single crystal structure.
#' @param cif_content A `data.table` containing the lines of a single CIF file.
#' @param file_name The name of the original CIF file.
#' @return A one-row `data.table` with results. Returns `NULL` on failure.
#' @export
#' @examples
#' cif_file <- system.file("extdata", "ICSD422.cif", package = "crysmal")
#' if (file.exists(cif_file)) {
#'   cif_content <- read_cif_files(cif_file)[[1]]
#'   processed_data <- process_single_cif_data(cif_content, basename(cif_file))
#'   str(processed_data, max.level = 1)
#' }
process_single_cif_data <- function(cif_content, file_name = "unknown") {
  database_code <- extract_database_code(cif_content)
  chemical_formula <- extract_chemical_formula(cif_content)
  structure_type <- extract_structure_type(cif_content)
  space_group_name <- extract_space_group_name(cif_content)
  space_group_number <- extract_space_group_number(cif_content)
  unit_cell_metrics <- extract_unit_cell_metrics(cif_content)
  atomic_coordinates <- extract_atomic_coordinates(cif_content)
  symmetry_operations <- extract_symmetry_operations(cif_content)

  if (is.null(atomic_coordinates) ||
      is.null(symmetry_operations) || is.null(unit_cell_metrics)) {
    warning(paste(
      "Could not process",
      file_name,
      "due to missing essential data."
    ))
    return(NULL)
  }

  transformed_coords <- apply_symmetry_operations(atomic_coordinates, symmetry_operations)
  expanded_coords <- expand_transformed_coords(transformed_coords)
  distances <- calculate_distances(atomic_coordinates, expanded_coords, unit_cell_metrics)
  bonded_pairs <- minimum_distance(distances)
  brunner_pairs <- brunner(distances)
  hoppe_pairs <- hoppe(distances)
  neighbor_counts <- calculate_neighbor_counts(bonded_pairs)
  bond_angles <- calculate_angles(bonded_pairs,
                                  atomic_coordinates,
                                  expanded_coords,
                                  unit_cell_metrics)
  bonded_pairs <- propagate_distance_error(bonded_pairs, atomic_coordinates, unit_cell_metrics)
  bond_angles <- propagate_angle_error(bond_angles,
                                       atomic_coordinates,
                                       expanded_coords,
                                       unit_cell_metrics)

  return(
    data.table(
      file_name = file_name,
      database_code = database_code,
      chemical_formula = chemical_formula,
      structure_type = structure_type,
      space_group_name = space_group_name,
      space_group_number = space_group_number,
      unit_cell_metrics = list(unit_cell_metrics),
      atomic_coordinates = list(atomic_coordinates),
      symmetry_operations = list(symmetry_operations),
      transformed_coords = list(transformed_coords),
      expanded_coords = list(expanded_coords),
      distances = list(distances),
      bonded_pairs = list(bonded_pairs),
      brunner_pairs = list(brunner_pairs),
      hoppe_pairs = list(hoppe_pairs),
      neighbor_counts = list(neighbor_counts),
      bond_angles = list(bond_angles)
    )
  )
}

#' @title Analyze a Batch of CIF Files
#' @description A high-level wrapper function that reads, processes, and analyzes
#'   one or more CIF files.
#' @param file_paths A character vector of paths to the CIF files to be analyzed.
#' @return A `data.table` where each row summarizes the analysis of one CIF file.
#' @export
#' @examples
#' cif_file <- system.file("extdata", "ICSD422.cif", package = "crysmal")
#' if (file.exists(cif_file)) {
#'   analysis_results <- analyze_cif_files(cif_file)
#'   print(analysis_results[, .(file_name, database_code, chemical_formula)])
#' }
analyze_cif_files <- function(file_paths) {
  cif_contents_list <- read_cif_files(file_paths)

  # Use mapply to iterate over both the content and the names (filenames)
  results_list <- mapply(function(cif_content, file_name) {
    tryCatch({
      process_single_cif_data(cif_content, file_name)
    }, error = function(e) {
      warning(paste("Failed to process file '", file_name, "': ", e$message, sep = ""))
      return(NULL)
    })
  },
  cif_contents_list,
  names(cif_contents_list),
  SIMPLIFY = FALSE)

  successful_results <- results_list[!sapply(results_list, is.null)]

  if (length(successful_results) > 0) {
    return(rbindlist(successful_results, fill = TRUE))
  } else {
    warning("No CIF files were processed successfully.")
    return(data.table())
  }
}
