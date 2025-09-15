#' @title Read CIF Files into Memory
#' @description Reads one or more CIF files from disk and loads each into a
#'   `data.table`.
#' @param file_paths A character vector of paths to the CIF files.
#' @return A list of `data.table` objects.
#' @export
#' @examples
#' cif_file <- system.file("extdata", "ICSD422.cif", package = "crystract")
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

#' @title Analyze the Content of a Single CIF File
#' @description The core worker function that orchestrates the analysis pipeline
#'   for a single crystal structure's data. It is called by `analyze_cif_files`
#'   for batch processing.
#' @param cif_content A `data.table` containing the lines of a single CIF file.
#' @param file_name The name of the original CIF file, used for labeling output.
#' @param perform_extraction Logical. If `TRUE`, extracts all metadata and basic
#'   structural data.
#' @param perform_calcs_and_transforms Logical. If `TRUE`, generates the full
#'   unit cell, expands it, and calculates interatomic distances.
#' @param bonding_algorithms A character vector of bonding algorithms to use.
#'   The first algorithm is used for subsequent angle/neighbor calculations.
#'   Options: `"minimum_distance"`, `"brunner"`, `"hoppe"`. Use `"none"` to skip.
#' @param calculate_bond_angles Logical. If `TRUE`, computes bond angles.
#' @param perform_error_propagation Logical. If `TRUE`, calculates uncertainties.
#'   Any missing error values in the CIF file are treated as zero.
#' @return A one-row `data.table` with results. If essential data for calculations
#'   is missing, a warning is issued and only partial results are returned.
#' @export
#' @examples
#' cif_file <- system.file("extdata", "ICSD422.cif", package = "crystract")
#' if (file.exists(cif_file)) {
#'   cif_content <- read_cif_files(cif_file)[[1]]
#'   # Using the single-file analyzer directly:
#'   single_result <- analyze_single_cif(cif_content, basename(cif_file))
#'   str(single_result, max.level = 1)
#' }
analyze_single_cif <- function(cif_content,
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

  # --- Step 2: Validate Essential Data for Calculations ---
  failure_reasons <- c()
  if (is.null(atomic_coordinates))
    failure_reasons <- c(failure_reasons, "missing atomic coordinates")
  if (is.null(symmetry_operations))
    failure_reasons <- c(failure_reasons, "missing symmetry operations")
  if (is.null(unit_cell_metrics)) {
    failure_reasons <- c(failure_reasons, "missing unit cell metrics block")
  } else {
    required_metrics <- c(
      "_cell_length_a",
      "_cell_length_b",
      "_cell_length_c",
      "_cell_angle_alpha",
      "_cell_angle_beta",
      "_cell_angle_gamma"
    )
    if (any(is.na(unit_cell_metrics[, ..required_metrics]))) {
      failure_reasons <- c(failure_reasons, "incomplete unit cell parameters")
    }
  }

  if (length(failure_reasons) > 0) {
    warning(
      paste(
        "Skipping calculations for file '",
        file_name,
        "' due to: ",
        paste(failure_reasons, collapse = "; "),
        ".",
        sep = ""
      )
    )
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
        symmetry_operations = list(symmetry_operations)
      )
    )
  }

  # --- Step 3: Initialize Result Variables ---
  transformed_coords <- NULL
  expanded_coords <- NULL
  distances <- NULL
  bonded_pairs_md <- NULL
  bonded_pairs_brunner <- NULL
  bonded_pairs_hoppe <- NULL
  primary_bonded_pairs <- NULL
  neighbor_counts <- NULL
  bond_angles <- NULL

  # --- Step 4: Calculations, Transformations, and Expansions ---
  if (perform_calcs_and_transforms) {
    transformed_coords <- apply_symmetry_operations(atomic_coordinates, symmetry_operations)
    expanded_coords <- expand_transformed_coords(transformed_coords)
    distances <- calculate_distances(atomic_coordinates, expanded_coords, unit_cell_metrics)
  }

  # --- Step 5: Bonding Algorithms ---
  if (!is.null(distances) &&
      length(bonding_algorithms) > 0 &&
      !"none" %in% bonding_algorithms) {
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
      if (is.null(current_bonds))
        next
      if (algo == "minimum_distance")
        bonded_pairs_md <- current_bonds
      if (algo == "brunner")
        bonded_pairs_brunner <- current_bonds
      if (algo == "hoppe")
        bonded_pairs_hoppe <- current_bonds
      if (algo == primary_algo)
        primary_bonded_pairs <- current_bonds
    }
    if (!is.null(primary_bonded_pairs)) {
      neighbor_counts <- calculate_neighbor_counts(primary_bonded_pairs)
    }
  }

  # --- Step 6: Bond Angles ---
  if (calculate_bond_angles && !is.null(primary_bonded_pairs)) {
    bond_angles <- calculate_angles(primary_bonded_pairs,
                                    atomic_coordinates,
                                    expanded_coords,
                                    unit_cell_metrics)
  }

  # --- Step 7: Error Propagation ---
  if (perform_error_propagation) {
    # These functions are designed to treat NA errors as 0.
    if (!is.null(bonded_pairs_md))
      bonded_pairs_md <- propagate_distance_error(bonded_pairs_md, atomic_coordinates, unit_cell_metrics)
    if (!is.null(bonded_pairs_brunner))
      bonded_pairs_brunner <- propagate_distance_error(bonded_pairs_brunner,
                                                       atomic_coordinates,
                                                       unit_cell_metrics)
    if (!is.null(bonded_pairs_hoppe))
      bonded_pairs_hoppe <- propagate_distance_error(bonded_pairs_hoppe,
                                                     atomic_coordinates,
                                                     unit_cell_metrics)
    if (!is.null(bond_angles))
      bond_angles <- propagate_angle_error(bond_angles,
                                           atomic_coordinates,
                                           expanded_coords,
                                           unit_cell_metrics)
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
#' @description A high-level wrapper function that reads and analyzes one or more
#'   CIF files in batch. It calls `analyze_single_cif` for each file.
#' @param file_paths A character vector of paths to the CIF files to be analyzed.
#' @param perform_extraction Logical. If `TRUE` (default), extracts all metadata
#'   and base crystallographic data.
#' @param perform_calcs_and_transforms Logical. If `TRUE` (default), generates
#'   the full unit cell, expands it to a supercell, and calculates distances.
#' @param bonding_algorithms A character vector specifying which bonding algorithms
#'   to apply. The first algorithm listed is used for subsequent calculations.
#'   Default: `c("minimum_distance")`. Options: `"minimum_distance"`, `"brunner"`,
#'   `"hoppe"`. Use `"none"` to skip bonding analysis.
#' @param calculate_bond_angles Logical. If `TRUE` (default), calculates bond
#'   angles based on the primary bonding algorithm.
#' @param perform_error_propagation Logical. If `TRUE` (default), propagates
#'   experimental uncertainties. Any missing uncertainty values in a CIF file
#'   will be treated as zero during this step.
#' @return A `data.table` where each row summarizes the analysis of one CIF file.
#'   Files with missing essential data will generate warnings and have `NULL`
#'   values in the calculated columns.
#' @export
#' @examples
#' cif_file <- system.file("extdata", "ICSD422.cif", package = "crystract")
#' if (file.exists(cif_file)) {
#'   # This will run the full analysis on all specified files.
#'   analysis_results <- analyze_cif_files(cif_file)
#'   print(analysis_results[, .(file_name, database_code, chemical_formula)])
#' }
analyze_cif_files <- function(file_paths,
                              perform_extraction = TRUE,
                              perform_calcs_and_transforms = TRUE,
                              bonding_algorithms = c("minimum_distance"),
                              calculate_bond_angles = TRUE,
                              perform_error_propagation = TRUE) {
  cif_contents_list <- read_cif_files(file_paths)

  results_list <- mapply(function(cif_content, file_name) {
    tryCatch({
      # Call the single-file worker function
      analyze_single_cif(
        cif_content,
        file_name,
        perform_extraction = perform_extraction,
        perform_calcs_and_transforms = perform_calcs_and_transforms,
        bonding_algorithms = bonding_algorithms,
        calculate_bond_angles = calculate_bond_angles,
        perform_error_propagation = perform_error_propagation
      )
    }, error = function(e) {
      warning(
        paste(
          "A critical error occurred while processing file '",
          file_name,
          "': ",
          e$message,
          sep = ""
        )
      )
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

#' @title Export Analysis Results to a Directory of CSVs
#' @description Takes the output `data.table` from `analyze_cif_files` or
#'   `analyze_single_cif` and exports its contents into a structured directory
#'   of CSV files. A 'meta' folder is created for top-level data, and
#'   separate folders are created for each nested table (e.g., atomic_coordinates,
#'   bond_angles).
#'
#' @details
#' The function operates as follows:
#' 1.  It creates the main `output_dir`. If the directory already exists,
#'     the function will stop unless `overwrite = TRUE`.
#' 2.  It identifies which columns in the input `analysis_results` are standard
#'     data (metadata) and which are list-columns containing nested `data.table`s.
#' 3.  The metadata is saved as a single `meta_summary.csv` file inside a `meta` sub-directory.
#' 4.  For each list-column (e.g., `unit_cell_metrics`), it creates a sub-directory
#'     with that name (e.g., `output_dir/unit_cell_metrics/`).
#' 5.  Inside each sub-directory, it iterates through every row of the original
#'     `analysis_results` table. For each row, it saves the corresponding nested
#'     `data.table` as a CSV file. The CSV is named after the CIF file it
#'     originated from (e.g., `ICSD422.csv`).
#'
#' This structure makes it easy to access all data of a specific type (e.g., all
#' bond angle tables) or all data related to a single original CIF file.
#'
#' @param analysis_results A `data.table` object, typically the output from
#'   `analyze_cif_files`.
#' @param output_dir A character string specifying the path to the main output
#'   directory where the folders and files will be created.
#' @param overwrite A logical value. If `TRUE`, any existing directory at
#'   `output_dir` will be removed and recreated. If `FALSE` (the default), the
#'   function will stop with an error if the directory already exists.
#' @return Invisibly returns the path to the `output_dir`.
#' @family post-processing
#' @export
#' @examples
#' # This is a full workflow example.
#'
#' # 1. Define path to an example CIF file
#' cif_file <- system.file("extdata", "ICSD422.cif", package = "crystract")
#' if (file.exists(cif_file)) {
#'   # 2. Run the analysis
#'   analysis_results <- analyze_cif_files(cif_file)
#'
#'   # 3. Define a temporary output directory for this example
#'   export_path <- file.path(tempdir(), "crystract_export")
#'
#'   # 4. Export the results, overwriting if the directory exists
#'   export_analysis_to_csv(analysis_results, export_path, overwrite = TRUE)
#'
#'   # 5. List the created files and folders to verify
#'   cat("Exported directory structure:\n")
#'   print(list.files(export_path, recursive = TRUE))
#'
#'   # 6. Clean up the temporary directory
#'   unlink(export_path, recursive = TRUE)
#' }
export_analysis_to_csv <- function(analysis_results, output_dir, overwrite = FALSE) {
  # --- Input Validation ---
  if (!inherits(analysis_results, "data.table") ||
      nrow(analysis_results) == 0) {
    stop("`analysis_results` must be a non-empty data.table from analyze_cif_files().")
  }
  if (!"file_name" %in% names(analysis_results)) {
    stop("`analysis_results` must contain a 'file_name' column.")
  }

  # --- Directory Management ---
  if (dir.exists(output_dir)) {
    if (overwrite) {
      unlink(output_dir, recursive = TRUE, force = TRUE)
      dir.create(output_dir, recursive = TRUE)
    } else {
      stop(
        paste(
          "Output directory '",
          output_dir,
          "' already exists. Use overwrite = TRUE to replace it."
        )
      )
    }
  } else {
    dir.create(output_dir, recursive = TRUE)
  }

  # --- Identify Column Types ---
  is_nested_col <- sapply(analysis_results, is.list)
  nested_col_names <- names(is_nested_col)[is_nested_col]
  meta_col_names <- names(is_nested_col)[!is_nested_col]

  # --- Export Metadata ---
  meta_dir <- file.path(output_dir, "meta")
  dir.create(meta_dir)
  meta_data <- analysis_results[, ..meta_col_names]
  data.table::fwrite(meta_data, file.path(meta_dir, "meta_summary.csv"))

  # --- Export Nested Tables ---
  for (col_name in nested_col_names) {
    # Skip creating empty folders for columns that are all NULL
    if (all(sapply(analysis_results[[col_name]], is.null))) {
      next
    }

    # Create a subdirectory for the nested table type
    nested_dir <- file.path(output_dir, col_name)
    dir.create(nested_dir)

    # Iterate through each file's results (each row)
    for (i in 1:nrow(analysis_results)) {
      nested_table <- analysis_results[[col_name]][[i]]

      if (!is.null(nested_table) &&
          inherits(nested_table, "data.table") && nrow(nested_table) > 0) {
        cif_name <- analysis_results$file_name[i]
        csv_name <- sub("(?i)\\.cif$", ".csv", cif_name, perl = TRUE)
        output_path <- file.path(nested_dir, csv_name)
        data.table::fwrite(nested_table, output_path)
      }
    }
  }

  message(paste(
    "Analysis successfully exported to:",
    normalizePath(output_dir)
  ))
  invisible(output_dir)
}
