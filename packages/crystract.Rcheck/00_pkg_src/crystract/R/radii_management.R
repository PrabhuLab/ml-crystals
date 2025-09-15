#' @title Set or Reset a Custom Atomic Radii Table
#' @description Allows the user to provide their own table of atomic radii for
#' the current R session. This custom table will be used by functions like
#' `filter_ghost_distances`.
#'
#' @details
#' The provided `data.table` or `data.frame` must contain at least three columns:
#' \itemize{
#'   \item `Symbol` (character): The chemical symbol of the element (e.g., "Si").
#'   \item `Radius` (numeric): The atomic radius in Angstroms (Å).
#'   \item `Type` (character): A descriptor for the radius type (e.g., "covalent", "ionic").
#' }
#' This allows for storing multiple types of radii in the same table, which can
#' be selected using the `radii_type` argument in relevant functions.
#'
#' To revert to using the package's default radii table, call the function
#' with `NULL` or without any arguments.
#'
#' @param radii_data A `data.table` or `data.frame` containing the custom radii
#' data, or `NULL` to reset to the package default.
#' @return Invisibly returns `NULL`. A message is printed to the console
#' confirming the action.
#' @export
#' @family post-processing
#' @examples
#' # 1. Create a custom radii table with both covalent and ionic radii
#' my_radii <- data.table::data.table(
#'   Symbol = c("Si", "O", "O"),
#'   Radius = c(1.11, 0.66, 1.40),
#'   Type   = c("covalent", "covalent", "ionic")
#' )
#'
#' # 2. Set this table for the current session
#' set_radii_data(my_radii)
#'
#' # Now, any subsequent calls to filter_ghost_distances() will use this table.
#' # You could, for example, filter based on ionic radii by setting radii_type = "ionic".
#'
#' # 3. Reset to the package's default covalent radii table
#' set_radii_data(NULL)
#'
set_radii_data <- function(radii_data = NULL) {
  if (is.null(radii_data)) {
    # If NULL, remove the custom data from the environment
    if (exists("user_radii_data", envir = .crystract_env)) {
      rm("user_radii_data", envir = .crystract_env)
    }
    message("Session radii data reset to package default.")
  } else {
    # --- Input Validation ---
    if (!is.data.frame(radii_data)) {
      stop("`radii_data` must be a data.frame or data.table.")
    }

    required_cols <- c("Symbol", "Radius", "Type")
    if (!all(required_cols %in% names(radii_data))) {
      stop(paste(
        "`radii_data` must contain the columns:",
        paste(required_cols, collapse = ", ")
      ))
    }

    # Coerce to data.table and store in the package environment
    custom_dt <- data.table::as.data.table(radii_data)

    if (!is.character(custom_dt$Symbol) ||
        !is.numeric(custom_dt$Radius) || !is.character(custom_dt$Type)) {
      stop(
        "Please check column types: Symbol (character), Radius (numeric), Type (character)."
      )
    }

    assign("user_radii_data", custom_dt, envir = .crystract_env)
    message("Custom radii table has been set for the current R session.")
  }
  invisible(NULL)
}

#' @title Get Radii Data (Internal)
#' @description An internal helper function that retrieves the active radii table.
#' It prioritizes a user-defined table set via `set_radii_data`, falling
#' back to the package's default `covalent_radii` data if none is set.
#' @return A `data.table` of atomic radii.
#' @noRd
get_radii_data <- function() {
  if (exists("user_radii_data", envir = .crystract_env)) {
    return(get("user_radii_data", envir = .crystract_env))
  } else {
    # Returns the package's default data
    return(crystract::covalent_radii)
  }
}
