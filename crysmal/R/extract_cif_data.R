#' @import data.table
#' @importFrom stringr str_split_fixed str_extract

# Internal helper function to extract a specific value
# Based on the original notebook logic. Not exported.
extract_value <- function(cif_content, pattern, remove_pattern = TRUE) {
  if (!is.data.table(cif_content) || !"V1" %in% names(cif_content)) {
    return(NA_character_)
  }
  lines_dt <- cif_content[V1 %like% pattern]
  if (nrow(lines_dt) > 0) {
    value <- lines_dt$V1[1] # Take the first match
    if (remove_pattern) {
      value <- gsub(pattern, "", value)
    }
    value <- gsub("'", "", value)
    value <- trimws(value)
    return(value)
  } else {
    return(NA_character_)
  }
}

#' Extract Crystallographic Information from CIF Content
#'
#' A suite of functions to extract specific metadata tags from CIF file content.
#'
#' @param cif_content A data.table representing the CIF file content, typically from \code{fread(..., sep = "\\n", header = FALSE)}.
#' @return A character string with the requested value, or NA if not found.
#' @name extract_cif
NULL

#' @rdname extract_cif
#' @export
#' @examples
#' # cif_data <- data.table(V1 = c("data_example", "_database_code_ '12345'"))
#' # extract_database_code(cif_data)
extract_database_code <- function (cif_content) {
  extract_value(cif_content, "_database_code_")
}

#' @rdname extract_cif
#' @export
#' @examples
#' # cif_data <- data.table(V1 = c("_chemical_formula_sum 'H2 O'"))
#' # extract_chemical_formula(cif_data)
extract_chemical_formula <- function(cif_content) {
  extract_value(cif_content, "_chemical_formula_sum")
}

#' @rdname extract_cif
#' @export
#' @examples
#' # cif_data <- data.table(V1 = c("_chemical_name_structure_type 'NaCl'"))
#' # extract_structure_type(cif_data)
extract_structure_type <- function(cif_content) {
  extract_value(cif_content, "_chemical_name_structure_type")
}

#' @rdname extract_cif
#' @export
#' @examples
#' # cif_data <- data.table(V1 = c("_space_group_name_H-M_alt 'P 1'"))
#' # extract_space_group_name(cif_data)
extract_space_group_name <- function(cif_content) {
  extract_value(cif_content, "_space_group_name_H-M_alt")
}

#' @rdname extract_cif
#' @export
#' @examples
#' # cif_data <- data.table(V1 = c("_space_group_IT_number 1"))
#' # extract_space_group_number(cif_data)
extract_space_group_number <- function(cif_content) {
  extract_value(cif_content, "_space_group_IT_number")
}

#' Extract Unit Cell Metrics from CIF Content
#'
#' @param cif_content A data.table representing the CIF file content.
#' @return A data.table with unit cell parameters (a, b, c, alpha, beta, gamma),
#'         or a data.table with NAs if parameters are not found.
#' @export
#' @examples
#' # cif_lines <- c("_cell_length_a   10.0", "_cell_length_b   10.0(1)",
#' #                "_cell_length_c   10.0", "_cell_angle_alpha 90.0",
#' #                "_cell_angle_beta  90.0", "_cell_angle_gamma 90.0")
#' # cif_data <- data.table(V1 = cif_lines)
#' # extract_unit_cell_metrics(cif_data)
extract_unit_cell_metrics <- function(cif_content) {
  cell_parameters <- c(
    "_cell_length_a", "_cell_length_b", "_cell_length_c",
    "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"
  )

  # Logic from the original notebook
  values <- sapply(cell_parameters, function(param) {
    line_dt <- cif_content[V1 %like% paste0("^", param)]
    if (nrow(line_dt) > 0) {
      # This simple regex from notebook doesn't handle errors in parentheses like 8.11(6)
      value <- gsub(".*\\s+([0-9\\.]+).*", "\\1", line_dt$V1[1])
      return(suppressWarnings(as.numeric(value)))
    } else {
      return(NA_real_)
    }
  })

  unit_cell_metrics <- as.data.table(t(values))
  setnames(unit_cell_metrics, cell_parameters)

  return(unit_cell_metrics)
}
