#' Function to Extract a Value From a CIF File
#' @import data.table
#' @param cif_content a string containing the content of a CIF file
#' @param pattern a regular expression pattern to match the desired value
#' @param remove_pattern a logical value indicating whether to remove the pattern from the value
#'
#' @return
#' @export
#' @examples
#' \dontrun{
#'  extract_value(cif_content, "data_ICSD_code")
#'  }
extract_value <- function(cif_content, pattern, remove_pattern = TRUE) {
  # Bind V1 locally
  V1 <-  cif_content$V1
  # Find lines that match the given pattern
  lines <- cif_content[V1 %like% pattern]
  if (nrow(lines) > 0) {
    value <- lines$V1
    if (remove_pattern) {
      value <- gsub(pattern, "", value)
    }
    # Clean the value by removing quotes and trimming whitespace
    value <- gsub("'", "", value)
    value <- trimws(value)
    return(value)
  } else {
    return(NA)
  }
}

#' Uses the extract_value function to extract the ICSD code from a CIF file
#'
#' @param cif_content a string containing the content of a CIF file
#'
#' @return The code from the ICSD database
#' @export
#'
#' @examples
#' \dontrun{
#' extract_ICSD_code(cif_content)
#' }
extract_ICSD_code <- function (cif_content) {
  extract_value(cif_content, "_database_code_ICSD")
}

#' Uses the extract_value function to extract the chemical formula from a CIF file
#'
#' @param cif_content a string containing the content of a CIF file
#'
#' @return The chemical formula of the compound
#' @export
#'
#' @examples
#' \dontrun{
#' extract_chemical_formula(cif_content)
#' }
ex_cf <- function(cif_content) {
  extract_value(cif_content, "_chemical_formula_sum")
}

#' Uses the extract_value function to extract the structure type from a CIF file
#'
#' @param cif_content a string containing the content of a CIF file
#'
#' @return The structure type of the compound
#' @export
#'
#' @examples
#' \dontrun{
#' extract_structure_type(cif_content)
#' }
ex_st <- function(cif_content) {
  extract_value(cif_content, "_chemical_name_structure_type")
}

#' Uses the extract_value function to extract the space group name from a CIF file
#'
#' @param cif_content a string containing the content of a CIF file
#'
#' @return The space group name
#' @export
#'
#' @examples
#' \dontrun{
#' extract_space_group_name(cif_content)
#' }
ex_sg_name <- function(cif_content) {
  extract_value(cif_content, "_space_group_name_H-M_alt")
}

#' Uses the extract_value function to extract the space group number from a CIF file
#'
#' @param cif_content a string containing the content of a CIF file
#'
#' @return The space group number
#' @export
#'
#' @examples
#' \dontrun{
#' extract_space_group_number(cif_content)
#' }
ex_sg_num <- function(cif_content) {
  extract_value(cif_content, "_space_group_IT_number")
}
