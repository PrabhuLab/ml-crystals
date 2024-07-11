#' Extract symmetry operations from CIF file
#'
#' @param cif_content a string containing the content of a CIF file
#' @import data.table
#' @return A data table with the extracted symmetry operations
#' @export
#'
#' @examples
#' \dontrun{
#' extract_symmetry_operations(cif_content)
#' }
extract_symmetry_operations <- function(cif_content) {
  symmetry_start <- grep("_space_group_symop_operation_xyz", cif_content$V1)
  symmetry_end <- grep("^loop_", cif_content$V1[symmetry_start:length(cif_content$V1)]) + symmetry_start - 1

  if (length(symmetry_end) > 1) {
    symmetry_end <- symmetry_end[1]
  }

  symmetry_lines <- cif_content$V1[(symmetry_start + 1):(symmetry_end - 1)]
  symmetry_operations <- str_split_fixed(symmetry_lines, ",\\s*", n = 3)

  symmetry_operations <- data.table(
    x = gsub("^[0-9]+\\s+'|'", "", symmetry_operations[, 1]),
    y = gsub("^[0-9]+\\s+'|'", "", symmetry_operations[, 2]),
    z = gsub("^[0-9]+\\s+'|'", "", symmetry_operations[, 3])
  )

  return(symmetry_operations)
}
