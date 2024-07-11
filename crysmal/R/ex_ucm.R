#' Function to extract unit cell metrics
#'
#' @param cif_content a string containing the content of a CIF file
#' @import data.table
#' @return A data table with the extracted unit cell metrics
#' @export
#'
#' @examples
#' \dontrun{
#' extract_unit_cell_metrics(cif_content)
#' }
ex_ucm <- function(cif_content) {
  # Locally bind V1
  V1 <- cif_content$V1
  # Define the cell parameters to extract
  cell_parameters <- c(
    "_cell_length_a",
    "_cell_length_b",
    "_cell_length_c",
    "_cell_angle_alpha",
    "_cell_angle_beta",
    "_cell_angle_gamma"
  )

  # Extract and clean values for each cell parameter
  values <- sapply(cell_parameters, function(param) {
    line <- cif_content[V1 %like% param]$V1
    if (length(line) > 0) {
      value <- gsub(".*\\s+([0-9\\.]+).*", "\\1", line)
      return(as.numeric(value))
    } else {
      return(NA)
    }
  })
  # Create a data table with the extracted values
  unit_cell_metrics <- as.data.table(t(values))
  setnames(unit_cell_metrics, cell_parameters)

  return(unit_cell_metrics)
}
