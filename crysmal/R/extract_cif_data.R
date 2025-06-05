# Internal helper function to extract a specific value
# Not exported
extract_value <- function(cif_content, pattern, remove_pattern = TRUE) {
  # Ensure cif_content is a data.table and has V1 column
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

#' Extract Database Code from CIF Content
#'
#' @param cif_content A data.table representing the CIF file content, typically from \code{fread(..., sep = "\\n", header = FALSE)}.
#' @return A character string with the database code, or NA if not found.
#' @export
#' @examples
#' # cif_data <- data.table(V1 = c("data_example", "_database_code_ '12345'"))
#' # extract_database_code(cif_data)
extract_database_code <- function (cif_content) {
  extract_value(cif_content, "_database_code_")
}

#' Extract Chemical Formula Sum from CIF Content
#'
#' @param cif_content A data.table representing the CIF file content.
#' @return A character string with the chemical formula sum, or NA if not found.
#' @export
#' @examples
#' # cif_data <- data.table(V1 = c("_chemical_formula_sum 'H2 O'"))
#' # extract_chemical_formula(cif_data)
extract_chemical_formula <- function(cif_content) {
  extract_value(cif_content, "_chemical_formula_sum")
}

#' Extract Structure Type from CIF Content
#'
#' @param cif_content A data.table representing the CIF file content.
#' @return A character string with the structure type, or NA if not found.
#' @export
#' @examples
#' # cif_data <- data.table(V1 = c("_chemical_name_structure_type 'NaCl'"))
#' # extract_structure_type(cif_data)
extract_structure_type <- function(cif_content) {
  extract_value(cif_content, "_chemical_name_structure_type")
}

#' Extract Space Group Name (H-M) from CIF Content
#'
#' @param cif_content A data.table representing the CIF file content.
#' @return A character string with the space group name, or NA if not found.
#' @export
#' @examples
#' # cif_data <- data.table(V1 = c("_space_group_name_H-M_alt 'P 1'"))
#' # extract_space_group_name(cif_data)
extract_space_group_name <- function(cif_content) {
  extract_value(cif_content, "_space_group_name_H-M_alt")
}

#' Extract Space Group IT Number from CIF Content
#'
#' @param cif_content A data.table representing the CIF file content.
#' @return A character string with the space group IT number, or NA if not found.
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
#' # cif_lines <- c("_cell_length_a   10.0", "_cell_length_b   10.0",
#' #                "_cell_length_c   10.0", "_cell_angle_alpha 90.0",
#' #                "_cell_angle_beta  90.0", "_cell_angle_gamma 90.0")
#' # cif_data <- data.table(V1 = cif_lines)
#' # extract_unit_cell_metrics(cif_data)
extract_unit_cell_metrics <- function(cif_content) {
  cell_parameters <- c(
    "_cell_length_a", "_cell_length_b", "_cell_length_c",
    "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"
  )

  values <- sapply(cell_parameters, function(param) {
    line_dt <- cif_content[V1 %like% paste0("^", param)] # Match from start of line
    if (nrow(line_dt) > 0) {
      line <- line_dt$V1[1]
      # Extract numeric value, removing uncertainty in parentheses if present
      value_str <- stringr::str_extract(line, "[0-9\\.]+\\(?[0-9]*\\)?$")
      value_str <- gsub("\\(.*\\)", "", value_str) # Remove parentheses and content
      value <- suppressWarnings(as.numeric(value_str))
      return(value)
    } else {
      return(NA_real_)
    }
  })

  unit_cell_metrics <- as.data.table(t(values))
  if (nrow(unit_cell_metrics) > 0) {
    setnames(unit_cell_metrics, cell_parameters)
  } else { # Handle case where no parameters are found (empty input cif_content)
    empty_dt <- data.table(matrix(ncol = length(cell_parameters), nrow = 0))
    setnames(empty_dt, cell_parameters)
    return(empty_dt)
  }
  # If unit_cell_metrics was created from sapply, it will have 1 row.
  # If all values were NA, it's a row of NAs, which is the desired outcome.
  return(unit_cell_metrics)
}
