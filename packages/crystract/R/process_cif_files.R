#' Internal, Safe Worker for Parallel CIF Analysis
#'
#' @description A robust wrapper for `analyze_single_cif` that includes error
#'   handling. Designed to be called by `future.apply` without capturing a
#'   large parent environment.
#'
#' @param cif_content A data.table of CIF content.
#' @param file_name The name of the file.
#' @param ... Additional arguments passed directly to `analyze_single_cif`.
#' @return A one-row data.table on success, or NULL on failure.
#' @noRd
.analyze_single_cif_safe <- function(cif_content, file_name, ...) {
  tryCatch({
    analyze_single_cif(cif_content = cif_content, file_name = file_name, ...)
  }, error = function(e) {
    warning(
      sprintf(
        "A critical error occurred while processing '%s': %s",
        file_name,
        e$message
      )
    )
    return(NULL)
  })
}

#' @title Read CIF Files into Memory
#' @description Reads one or more CIF files from disk and loads each into a
#'   `data.table`. It includes encoding sanitization to handle legacy characters.
#' @param file_paths A character vector of paths to the CIF files.
#' @return A list of `data.table` objects.
#' @export
read_cif_files <- function(file_paths) {
  # Set the names of the list to the base filenames
  cif_list <- lapply(file_paths, function(fp) {
    # Read the file
    dt <- data.table::fread(
      fp,
      sep = "\n",
      header = FALSE,
      strip.white = FALSE,
      encoding = "UTF-8" # Attempt UTF-8 first
    )

    # Sanitize: Convert strings to valid UTF-8, replacing invalid bytes with empty strings
    if (nrow(dt) > 0) {
      dt[, V1 := iconv(V1, to = "UTF-8", sub = "")]
    }

    return(dt)
  })

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
  atomic_coordinates <- extract_atomic_coordinates(cif_content, chemical_formula)
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
    return(NULL)
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
#' @description A high-level wrapper that analyzes CIF files in batch, supporting
#'   parallel processing and batch-to-disk operations for large datasets.
#'
#' @details
#' This function is optimized for large datasets. It processes files in chunks
#' (`batch_size`) and uses a memory-safe parallel backend to avoid data transfer
#' bottlenecks.
#'
#' It operates in two modes:
#' 1.  **In-Memory (default):** If `output_dir` is `NULL`, results are aggregated
#'     in memory and returned as one `data.table`. Ideal for smaller jobs.
#' 2.  **Batch-to-Disk:** If `output_dir` is a path, each processed batch is
#'     saved as an RDS file to the specified directory. This is the recommended
#'     mode for very large datasets as it prevents memory overload.
#'
#' @param file_paths A character vector of CIF file paths, or a named list of
#'   `data.table`s each containing CIF content.
#' @param workers An integer for the number of parallel workers. Defaults to `1` (sequential).
#' @param output_dir Optional path to a directory for saving batched results.
#'   If `NULL` (default), results are returned in memory.
#' @param batch_size Number of files per batch. Defaults to 1000.
#' @param ... Additional arguments passed directly to `analyze_single_cif`, such
#'   as `bonding_algorithms`, `calculate_bond_angles`, etc.
#'
#' @return If `output_dir` is `NULL`, a single `data.table` with all results.
#'   If `output_dir` is provided, invisibly returns the path to the output directory.
#' @export
analyze_cif_files <- function(file_paths,
                              workers = 1,
                              output_dir = NULL,
                              batch_size = 1000,
                              ...) {
  # --- 1. Handle Input: File Paths or Pre-parsed List ---
  if (is.character(file_paths)) {
    # Check if the input is a single directory path
    if (length(file_paths) == 1 && dir.exists(file_paths)) {
      message(sprintf(
        "Input is a directory. Searching for .cif files in '%s'...",
        file_paths
      ))
      found_files <- list.files(
        file_paths,
        pattern = "\\.cif$",
        full.names = TRUE,
        ignore.case = TRUE
      )

      if (length(found_files) == 0) {
        warning(sprintf("No .cif files found in directory '%s'.", file_paths))
        return(data.table())
      }
      file_paths <- found_files
    }

    cif_contents_list <- read_cif_files(file_paths)
  } else if (is.list(file_paths) &&
             (length(file_paths) == 0 ||
              inherits(file_paths[[1]], "data.table"))) {
    cif_contents_list <- file_paths
    if (is.null(names(cif_contents_list))) {
      names(cif_contents_list) <- paste0("unnamed_", seq_along(cif_contents_list))
    }
  } else {
    stop(
      "`file_paths` must be a character vector of paths, a directory path, or a list of data.tables."
    )
  }

  # --- 2. Set Up Parallel Plan if specified ---
  if (workers > 1) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop(
        "Packages 'future' and 'future.apply' are required for parallel processing.",
        call. = FALSE
      )
    }
    # Set up the plan and ensure it's reverted on exit
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    message(paste0("Setting parallel plan to use ", workers, " workers."))
    future::plan(future::multisession, workers = workers)
  }

  # --- 3. Main Logic: Batch Processing ---
  total_files <- length(cif_contents_list)
  if (total_files == 0) {
    message("No files to process.")
    return(data.table())
  }

  batch_starts <- seq(1, total_files, by = batch_size)
  num_batches <- length(batch_starts)

  # Store all results here if not writing to disk
  all_results_list <- if (is.null(output_dir))
    list()
  else
    NULL

  message(sprintf(
    "Starting analysis of %d files in %d batches.",
    total_files,
    num_batches
  ))

  # List of extra arguments for the worker
  more_args <- list(...)

  for (i in seq_along(batch_starts)) {
    start_index <- batch_starts[i]
    end_index <- min(start_index + batch_size - 1, total_files)

    # Create the small batch to be processed
    cif_batch <- cif_contents_list[start_index:end_index]

    message(
      sprintf(
        "\n--- Processing Batch %d of %d (Files %d to %d) ---",
        i,
        num_batches,
        start_index,
        end_index
      )
    )

    # This is the key change: we call our safe, top-level worker.
    # The 'future.apply' function sends ONLY the small 'cif_batch' and its names
    # to the workers, avoiding the memory overload.
    results_list <- if (workers > 1) {
      future.apply::future_mapply(
        FUN = .analyze_single_cif_safe,
        cif_content = cif_batch,
        file_name = names(cif_batch),
        MoreArgs = more_args,
        SIMPLIFY = FALSE,
        future.seed = TRUE
      )
    } else {
      # Fallback to sequential mapply if workers = 1
      mapply(
        FUN = .analyze_single_cif_safe,
        cif_content = cif_batch,
        file_name = names(cif_batch),
        MoreArgs = more_args,
        SIMPLIFY = FALSE
      )
    }

    successful_results <- results_list[!sapply(results_list, is.null)]

    if (length(successful_results) > 0) {
      batch_dt <- rbindlist(successful_results, fill = TRUE)

      if (is.null(output_dir)) {
        # Mode 1: Append to list for in-memory aggregation
        all_results_list[[i]] <- batch_dt
      } else {
        # Mode 2: Save batch to its own file on disk
        if (!dir.exists(output_dir))
          dir.create(output_dir, recursive = TRUE)
        batch_filename <- file.path(output_dir, paste0("batch_", i, ".rds"))
        saveRDS(batch_dt, batch_filename)
        message(
          sprintf(
            "Batch %d complete. Saved %d results to '%s'.",
            i,
            nrow(batch_dt),
            batch_filename
          )
        )
      }
    } else {
      message(sprintf(
        "Batch %d complete. No structures processed successfully.\n",
        i
      ))
    }

    # Clean up memory
    rm(cif_batch, results_list, successful_results)
    gc()
  }

  # --- 4. Finalize and Return ---
  message("----------------------------------\nAnalysis Complete!")

  if (is.null(output_dir)) {
    # In-memory mode: combine all batch results and return
    if (length(all_results_list) > 0) {
      final_dt <- rbindlist(all_results_list, fill = TRUE)
      message(sprintf(
        "Successfully processed a total of %d structures.",
        nrow(final_dt)
      ))
      return(final_dt)
    } else {
      warning("No structures were processed successfully.")
      return(data.table())
    }
  } else {
    # Batch-to-disk mode: return the path to the results directory
    message(sprintf(
      "Batch results have been saved in the '%s/' directory.",
      normalizePath(output_dir)
    ))
    return(invisible(normalizePath(output_dir)))
  }
}

#' @title Aggregate Batched Analysis Results
#' @description Reads all `batch_*.rds` files from a directory and combines them.
#'   For very large datasets, it supports selecting specific columns to reduce
#'   memory usage during aggregation.
#'
#' @param input_dir The directory containing RDS files from `analyze_cif_files`.
#' @param cols_to_keep An optional character vector of column names to keep.
#'   If provided, only these columns will be loaded from each batch file,
#'   dramatically reducing memory consumption. If `NULL` (default), all columns
#'   are loaded.
#' @return A single `data.table` containing the aggregated results.
#' @family post-processing
#' @export
aggregate_batch_results <- function(input_dir, cols_to_keep = NULL) {
  batch_files <- list.files(path = input_dir,
                            pattern = "^batch_\\d+\\.rds$",
                            full.names = TRUE)
  if (length(batch_files) == 0) {
    stop("No batch result files (e.g., 'batch_1.rds') found in the directory.")
  }

  message(sprintf("Found %d batch files to aggregate.", length(batch_files)))

  # Use lapply to read and optionally subset each file
  all_results_list <- lapply(batch_files, function(file) {
    batch_dt <- readRDS(file)
    if (!is.null(cols_to_keep)) {
      # Check which of the requested columns actually exist in the data.table
      valid_cols <- intersect(cols_to_keep, names(batch_dt))
      if (length(valid_cols) > 0) {
        return(batch_dt[, ..valid_cols])
      }
    }
    return(batch_dt)
  })

  aggregated_results <- rbindlist(all_results_list, fill = TRUE)

  message(sprintf(
    "Aggregation complete. Total rows: %d.",
    nrow(aggregated_results)
  ))

  return(aggregated_results)
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
          inherits(nested_table, "data.table") &&
          nrow(nested_table) > 0) {
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
