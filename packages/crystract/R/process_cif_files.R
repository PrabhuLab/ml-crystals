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
#' @description This function orchestrates the analysis pipeline for a
#'   single crystal structure with user-configurable steps.
#' @param cif_content A `data.table` containing the lines of a single CIF file.
#' @param file_name The name of the original CIF file.
#' @param perform_extraction Logical. If `TRUE`, extracts all metadata and basic
#'   structural data. This is required for all subsequent steps.
#' @param perform_calcs_and_transforms Logical. If `TRUE`, generates the full
#'   unit cell, expands it to a supercell, and calculates interatomic distances.
#' @param bonding_algorithms A character vector specifying one or more bonding
#'   algorithms to use. Results are stored in separate columns. The first
#'   algorithm in the vector is used for subsequent neighbor and angle calculations.
#'   Options: `"minimum_distance"`, `"brunner"`, `"hoppe"`. Use `"none"` to skip.
#' @param calculate_bond_angles Logical. If `TRUE`, computes bond angles based on
#'   the primary bonding algorithm.
#' @param perform_error_propagation Logical. If `TRUE`, calculates uncertainties
#'   for all computed distances and angles.
#' @return A one-row `data.table` with results. Returns `NULL` on failure.
#' @export
process_single_cif_data <- function(cif_content,
                                    file_name = "unknown",
                                    perform_extraction = TRUE,
                                    perform_calcs_and_transforms = TRUE,
                                    bonding_algorithms = c("minimum_distance"),
                                    calculate_bond_angles = TRUE,
                                    perform_error_propagation = TRUE) {
  # --- Step 1: Data Extraction ---
  if (!perform_extraction) {
    return(data.table(file_name = file_name))
  }
  database_code <- extract_database_code(cif_content)
  chemical_formula <- extract_chemical_formula(cif_content)
  structure_type <- extract_structure_type(cif_content)
  space_group_name <- extract_space_group_name(cif_content)
  space_group_number <- extract_space_group_number(cif_content)
  unit_cell_metrics <- extract_unit_cell_metrics(cif_content)
  atomic_coordinates <- extract_atomic_coordinates(cif_content)
  symmetry_operations <- extract_symmetry_operations(cif_content)

  # Initialize remaining result variables to NULL
  transformed_coords <- NULL
  expanded_coords <- NULL
  distances <- NULL
  bonded_pairs_md <- NULL
  bonded_pairs_brunner <- NULL
  bonded_pairs_hoppe <- NULL
  primary_bonded_pairs <- NULL
  neighbor_counts <- NULL
  bond_angles <- NULL

  # Essential data check
  if (is.null(atomic_coordinates) || is.null(symmetry_operations) || is.null(unit_cell_metrics)) {
    warning(paste("Could not process", file_name, "due to missing essential data."))
    return(data.table(
      file_name = file_name,
      database_code = database_code,
      chemical_formula = chemical_formula,
      structure_type = structure_type,
      space_group_name = space_group_name,
      space_group_number = space_group_number,
      unit_cell_metrics = list(unit_cell_metrics),
      atomic_coordinates = list(atomic_coordinates),
      symmetry_operations = list(symmetry_operations)
    ))
  }

  # --- Step 2: Calculations, Transformations, and Expansions ---
  if (perform_calcs_and_transforms) {
    transformed_coords <- apply_symmetry_operations(atomic_coordinates, symmetry_operations)
    expanded_coords <- expand_transformed_coords(transformed_coords)
    distances <- calculate_distances(atomic_coordinates, expanded_coords, unit_cell_metrics)
  }

  # --- Step 3: Choose and Apply Bonding Algorithm(s) ---
  if (!is.null(distances) && length(bonding_algorithms) > 0 && !"none" %in% bonding_algorithms) {
    primary_algo <- bonding_algorithms[1]

    for (algo in unique(bonding_algorithms)) {
      current_bonds <- switch(
        algo,
        "minimum_distance" = minimum_distance(distances),
        "brunner" = brunner(distances),
        "hoppe" = hoppe(distances),
        {
          warning(paste("Invalid bonding algorithm '", algo, "' ignored.", sep = ""))
          NULL
        }
      )
      if (is.null(current_bonds)) next
      if (algo == "minimum_distance") bonded_pairs_md <- current_bonds
      if (algo == "brunner") bonded_pairs_brunner <- current_bonds
      if (algo == "hoppe") bonded_pairs_hoppe <- current_bonds
      if (algo == primary_algo) primary_bonded_pairs <- current_bonds
    }
    if (!is.null(primary_bonded_pairs)) {
      neighbor_counts <- calculate_neighbor_counts(primary_bonded_pairs)
    }
  }

  # --- Step 4: Calculate Bond Angles (based on primary algorithm) ---
  if (calculate_bond_angles && !is.null(primary_bonded_pairs)) {
    bond_angles <- calculate_angles(primary_bonded_pairs, atomic_coordinates, expanded_coords, unit_cell_metrics)
  }

  # --- Step 5: Error Propagation ---
  if (perform_error_propagation) {
    if (!is.null(bonded_pairs_md)) {
      bonded_pairs_md <- propagate_distance_error(bonded_pairs_md, atomic_coordinates, unit_cell_metrics)
    }
    if (!is.null(bonded_pairs_brunner)) {
      bonded_pairs_brunner <- propagate_distance_error(bonded_pairs_brunner, atomic_coordinates, unit_cell_metrics)
    }
    if (!is.null(bonded_pairs_hoppe)) {
      bonded_pairs_hoppe <- propagate_distance_error(bonded_pairs_hoppe, atomic_coordinates, unit_cell_metrics)
    }
    if (!is.null(bond_angles)) {
      bond_angles <- propagate_angle_error(bond_angles, atomic_coordinates, expanded_coords, unit_cell_metrics)
    }
  }

  # --- Assemble Final Results ---
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
      bonded_pairs_minimum_distance = list(bonded_pairs_md),
      bonded_pairs_brunner = list(bonded_pairs_brunner),
      bonded_pairs_hoppe = list(bonded_pairs_hoppe),
      neighbor_counts = list(neighbor_counts),
      bond_angles = list(bond_angles)
    )
  )
}


#' @title Analyze a Batch of CIF Files
#' @description A high-level wrapper function that reads, processes, and analyzes
#'   one or more CIF files, with parameters to control the workflow.
#' @param file_paths A character vector of paths to the CIF files to be analyzed.
#' @param perform_extraction Logical. If `TRUE` (default), extracts all metadata
#'   and base crystallographic data. Required for all other steps.
#' @param perform_calcs_and_transforms Logical. If `TRUE` (default), generates
#'   the full unit cell, expands it to a supercell, and calculates distances.
#' @param bonding_algorithms A character vector specifying which bonding algorithms
#'   to apply. The results for each are stored in separate columns (e.g.,
#'   `bonded_pairs_minimum_distance`). The first algorithm listed is considered
#'   primary and is used for subsequent `neighbor_counts` and `bond_angles`
#'   calculations.
#'   Default: `c("minimum_distance")`.
#'   Options: `"minimum_distance"`, `"brunner"`, `"hoppe"`.
#'   Use `"none"` to skip all bonding analysis.
#' @param calculate_bond_angles Logical. If `TRUE` (default), calculates bond
#'   angles based on the primary bonding algorithm.
#' @param perform_error_propagation Logical. If `TRUE` (default), propagates
#'   experimental uncertainties to all calculated distances and angles.
#' @return A `data.table` where each row summarizes the analysis of one CIF file.
#' @export
#' @examples
#' cif_file <- system.file("extdata", "ICSD422.cif", package = "crysmal")
#' if (file.exists(cif_file)) {
#'   # Example 1: Run two bonding algorithms, but no error propagation.
#'   analysis_results <- analyze_cif_files(
#'     cif_file,
#'     bonding_algorithms = c("minimum_distance", "hoppe"),
#'     perform_error_propagation = FALSE
#'   )
#'   # The output will have 'bonded_pairs_minimum_distance' and 'bonded_pairs_hoppe' columns.
#'   # 'neighbor_counts' and 'bond_angles' will be based on the minimum_distance method.
#'   # str(analysis_results, max.level = 1)
#'
#'   # Example 2: Run only the initial extraction and transformations.
#'   analysis_subset <- analyze_cif_files(
#'     cif_file,
#'     bonding_algorithms = "none",
#'     calculate_bond_angles = FALSE,
#'     perform_error_propagation = FALSE
#'   )
#'   # str(analysis_subset, max.level = 1)
#' }
analyze_cif_files <- function(file_paths,
                              perform_extraction = TRUE,
                              perform_calcs_and_transforms = TRUE,
                              bonding_algorithms = c("minimum_distance"),
                              calculate_bond_angles = TRUE,
                              perform_error_propagation = TRUE) {
  cif_contents_list <- read_cif_files(file_paths)

  # Use mapply to iterate over both the content and the names (filenames)
  results_list <- mapply(function(cif_content, file_name) {
    tryCatch({
      process_single_cif_data(
        cif_content,
        file_name,
        perform_extraction = perform_extraction,
        perform_calcs_and_transforms = perform_calcs_and_transforms,
        bonding_algorithms = bonding_algorithms,
        calculate_bond_angles = calculate_bond_angles,
        perform_error_propagation = perform_error_propagation
      )
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
