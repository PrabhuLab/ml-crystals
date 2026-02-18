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
#'   for a single crystal structure's data. It dynamically adjusts supercell
#'   size based on the requested bonding algorithms to ensure Voronoi accuracy.
#' @param cif_content Either a `data.table` containing the lines of a CIF file,
#'   OR a character string specifying the file path to a CIF file.
#' @param file_name The name of the original CIF file, used for labeling output.
#' @param perform_extraction Logical. If `TRUE`, extracts all metadata.
#' @param perform_calcs_and_transforms Logical. If `TRUE`, generates unit cell
#'   and supercell.
#' @param bonding_algorithms Character vector. Options: `"minimum_distance"`,
#'   `"brunner"`, `"econ"`, `"crystal_nn"`, `"voronoi"`.
#' @param calculate_bond_angles Logical.
#' @param perform_error_propagation Logical.
#' @param tolerance Numeric. Cutoff for floating-point noise and atom merging (default 1e-4).
#' @param minimum_distance_delta Numeric. Tolerance for min-dist algorithm.
#' @return A one-row `data.table` with results.
#' @export
analyze_single_cif <- function(cif_content,
                               file_name = NULL,
                               perform_extraction = TRUE,
                               perform_calcs_and_transforms = TRUE,
                               bonding_algorithms = c("minimum_distance"),
                               calculate_bond_angles = TRUE,
                               perform_error_propagation = TRUE,
                               tolerance = 1e-4,
                               minimum_distance_delta = 0.1) {
  # --- Step 0: Handle Input Type ---
  if (is.character(cif_content)) {
    if (length(cif_content) != 1)
      stop("`analyze_single_cif` accepts exactly one string.")
    if (!file.exists(cif_content))
      stop(paste("File not found:", cif_content))
    if (is.null(file_name))
      file_name <- basename(cif_content)
    cif_content <- read_cif_files(cif_content)[[1]]
  } else if (inherits(cif_content, "data.table")) {
    if (is.null(file_name))
      file_name <- "unknown"
  } else {
    stop("`cif_content` must be a file path or loaded CIF data.")
  }

  # --- Step 1: Data Extraction ---
  if (!perform_extraction)
    return(data.table::data.table(file_name = file_name))

  database_code <- extract_database_code(cif_content)
  chemical_formula <- extract_chemical_formula(cif_content)
  structure_type <- extract_structure_type(cif_content)
  space_group_name <- extract_space_group_name(cif_content)
  space_group_number <- extract_space_group_number(cif_content)
  unit_cell_metrics <- extract_unit_cell_metrics(cif_content)
  atomic_coordinates <- extract_atomic_coordinates(cif_content, chemical_formula)
  symmetry_operations <- extract_symmetry_operations(cif_content)

  # --- Step 2: Validate Essential Data ---
  if (is.null(atomic_coordinates) ||
      is.null(symmetry_operations) ||
      is.null(unit_cell_metrics)) {
    warning(paste(
      "Skipping '",
      file_name,
      "': missing essential structural data.",
      sep = ""
    ))
    return(NULL)
  }

  # --- Step 3: Initialize Result Variables ---
  transformed_coords <- NULL
  expanded_coords <- NULL
  distances <- NULL
  algo_results <- list()

  # --- Step 4: Calculations, Transformations, and Expansions ---
  if (perform_calcs_and_transforms) {
    transformed_coords <- apply_symmetry_operations(atomic_coordinates,
                                                    symmetry_operations,
                                                    unit_cell_metrics,
                                                    tolerance = tolerance)

    # DYNAMIC EXPANSION LOGIC
    # Constrain logic: only expand > 3x3x3 if Voronoi/CrystalNN are requested
    expansion_factors <- c(1, 1, 1) # Default 3x3x3 (-1:1)

    needs_voronoi <- any(c("voronoi", "crystal_nn") %in% bonding_algorithms)

    if (needs_voronoi) {
      # Voronoi needs typically 13.0 Angstroms coverage
      voronoi_cutoff <- 13.0
      expansion_factors <- calculate_expansion_factors(unit_cell_metrics, voronoi_cutoff)
    }

    expanded_coords <- expand_transformed_coords(transformed_coords, expansion_factors)

    distances <- calculate_distances(atomic_coordinates,
                                     expanded_coords,
                                     unit_cell_metrics,
                                     tolerance = tolerance)
  }

  # --- Step 5: Bonding Algorithms & Step 6: Angles & Errors ---
  if (!is.null(distances) &&
      length(bonding_algorithms) > 0 &&
      !"none" %in% bonding_algorithms) {
    for (algo in unique(bonding_algorithms)) {
      current_bonds <- switch(
        algo,
        "minimum_distance" = minimum_distance(distances, delta = minimum_distance_delta),
        "brunner" = brunner_nn_reciprocal(distances),
        "econ" = econ_nn(distances, atomic_coordinates),
        "voronoi" = voronoi_nn(atomic_coordinates, expanded_coords, unit_cell_metrics),
        "crystal_nn" = crystal_nn(
          distances,
          atomic_coordinates,
          expanded_coords,
          unit_cell_metrics
        ),
        {
          warning(paste("Invalid algorithm '", algo, "' ignored.", sep = ""))
          NULL
        }
      )

      if (!is.null(current_bonds) && nrow(current_bonds) > 0) {
        out_algo_name <- algo

        if (perform_error_propagation) {
          current_bonds <- propagate_distance_error(current_bonds,
                                                    atomic_coordinates,
                                                    unit_cell_metrics)
        }

        algo_results[[paste0("bonds_", out_algo_name)]] <- list(current_bonds)
        cn_table <- calculate_neighbor_counts(current_bonds)
        algo_results[[paste0("cn_", out_algo_name)]] <- list(cn_table)

        if (calculate_bond_angles) {
          current_angles <- calculate_angles(current_bonds,
                                             atomic_coordinates,
                                             expanded_coords,
                                             unit_cell_metrics)

          if (perform_error_propagation &&
              !is.null(current_angles) &&
              nrow(current_angles) > 0) {
            current_angles <- propagate_angle_error(
              current_angles,
              atomic_coordinates,
              expanded_coords,
              unit_cell_metrics
            )
          }
          algo_results[[paste0("angles_", out_algo_name)]] <- list(current_angles)
        }
      }
    }
  }

  # --- Step 7: Construct Return Data Table ---
  final_dt <- data.table::data.table(
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
    distances = list(distances)
  )

  if (length(algo_results) > 0) {
    for (res_name in names(algo_results)) {
      set(final_dt, j = res_name, value = algo_results[[res_name]])
    }
  }

  return(final_dt)
}

#' @title Analyze a Batch of CIF Files
#' @description A high-level wrapper that analyzes CIF files in batch, supporting
#'   parallel processing and batch-to-disk operations.
#' @param file_paths Character vector of paths or list of data.tables.
#' @param workers Integer. Number of parallel workers.
#' @param output_dir Path to output directory (optional).
#' @param batch_size Integer.
#' @param ... Args passed to `analyze_single_cif`.
#' @return Data.table or output directory path.
#' @export
analyze_cif_files <- function(file_paths,
                              workers = 1,
                              output_dir = NULL,
                              batch_size = 1000,
                              ...) {
  if (is.character(file_paths)) {
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
        return(data.table::data.table())
      }
      file_paths <- found_files
    }
    cif_contents_list <- read_cif_files(file_paths)
  } else if (is.list(file_paths) &&
             (length(file_paths) == 0 ||
              inherits(file_paths[[1]], "data.table"))) {
    cif_contents_list <- file_paths
    if (is.null(names(cif_contents_list)))
      names(cif_contents_list) <- paste0("unnamed_", seq_along(cif_contents_list))
  } else {
    stop(
      "`file_paths` must be a character vector of paths, a directory path, or a list of data.tables."
    )
  }

  if (workers > 1) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop(
        "Packages 'future' and 'future.apply' are required for parallel processing.",
        call. = FALSE
      )
    }
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    message(paste0("Setting parallel plan to use ", workers, " workers."))
    future::plan(future::multisession, workers = workers)
  }

  total_files <- length(cif_contents_list)
  if (total_files == 0) {
    message("No files to process.")
    return(data.table::data.table())
  }

  batch_starts <- seq(1, total_files, by = batch_size)
  num_batches <- length(batch_starts)
  all_results_list <- if (is.null(output_dir))
    list()
  else
    NULL
  message(sprintf(
    "Starting analysis of %d files in %d batches.",
    total_files,
    num_batches
  ))
  more_args <- list(...)

  for (i in seq_along(batch_starts)) {
    start_index <- batch_starts[i]
    end_index <- min(start_index + batch_size - 1, total_files)
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
      batch_dt <- data.table::rbindlist(successful_results, fill = TRUE)
      if (is.null(output_dir)) {
        all_results_list[[i]] <- batch_dt
      } else {
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
    rm(cif_batch, results_list, successful_results)
    gc()
  }

  message("----------------------------------\nAnalysis Complete!")
  if (is.null(output_dir)) {
    if (length(all_results_list) > 0) {
      final_dt <- data.table::rbindlist(all_results_list, fill = TRUE)
      message(sprintf(
        "Successfully processed a total of %d structures.",
        nrow(final_dt)
      ))
      return(final_dt)
    } else {
      warning("No structures were processed successfully.")
      return(data.table::data.table())
    }
  } else {
    message(sprintf(
      "Batch results have been saved in the '%s/' directory.",
      normalizePath(output_dir)
    ))
    return(invisible(normalizePath(output_dir)))
  }
}

#' @title Aggregate Batched Analysis Results
#' @description Reads all `batch_*.rds` files from a directory and combines them.
#' @param input_dir The directory containing RDS files from `analyze_cif_files`.
#' @param cols_to_keep Optional character vector of column names to keep.
#' @return A single `data.table` containing the aggregated results.
#' @export
aggregate_batch_results <- function(input_dir, cols_to_keep = NULL) {
  batch_files <- list.files(path = input_dir,
                            pattern = "^batch_\\d+\\.rds$",
                            full.names = TRUE)
  if (length(batch_files) == 0) {
    stop("No batch result files (e.g., 'batch_1.rds') found in the directory.")
  }

  message(sprintf("Found %d batch files to aggregate.", length(batch_files)))

  all_results_list <- lapply(batch_files, function(file) {
    batch_dt <- readRDS(file)
    if (!is.null(cols_to_keep)) {
      valid_cols <- intersect(cols_to_keep, names(batch_dt))
      if (length(valid_cols) > 0) {
        return(batch_dt[, ..valid_cols])
      }
    }
    return(batch_dt)
  })

  aggregated_results <- data.table::rbindlist(all_results_list, fill = TRUE)
  message(sprintf(
    "Aggregation complete. Total rows: %d.",
    nrow(aggregated_results)
  ))
  return(aggregated_results)
}

#' @title Export Analysis Results to a Directory of CSVs
#' @description Exports analysis results to CSV structure.
#' @param analysis_results A `data.table` object.
#' @param output_dir Path to main output directory.
#' @param overwrite Logical.
#' @return Invisibly returns output_dir.
#' @export
export_analysis_to_csv <- function(analysis_results, output_dir, overwrite = FALSE) {
  if (!inherits(analysis_results, "data.table") ||
      nrow(analysis_results) == 0) {
    stop("`analysis_results` must be a non-empty data.table.")
  }
  if (!"file_name" %in% names(analysis_results)) {
    stop("`analysis_results` must contain a 'file_name' column.")
  }

  if (dir.exists(output_dir)) {
    if (overwrite) {
      unlink(output_dir, recursive = TRUE, force = TRUE)
      dir.create(output_dir, recursive = TRUE)
    } else {
      stop(paste(
        "Output directory '",
        output_dir,
        "' already exists. Use overwrite = TRUE."
      ))
    }
  } else {
    dir.create(output_dir, recursive = TRUE)
  }

  is_nested_col <- sapply(analysis_results, is.list)
  nested_col_names <- names(is_nested_col)[is_nested_col]
  meta_col_names <- names(is_nested_col)[!is_nested_col]

  meta_dir <- file.path(output_dir, "meta")
  dir.create(meta_dir)
  meta_data <- analysis_results[, ..meta_col_names]
  data.table::fwrite(meta_data, file.path(meta_dir, "meta_summary.csv"))

  for (col_name in nested_col_names) {
    if (all(sapply(analysis_results[[col_name]], is.null)))
      next
    nested_dir <- file.path(output_dir, col_name)
    dir.create(nested_dir)
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
