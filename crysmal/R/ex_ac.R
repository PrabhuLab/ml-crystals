#' Extracting atomic coordinates
#'
#' @param cif_content a string containing the content of a CIF file
#' @import data.table
#' @return A data table with the extracted atomic coordinates
#' @export
#'
#' @examples
#' \dontrun{
#' extract_atomic_coordinates(cif_content)
#' }
extract_atomic_coordinates <- function(cif_content) {
  # Locally bind V1
  V1 <- cif_content$V1
  # Initialize lists to store atom properties
  labels <- vector()
  fractional_x <- vector()
  fractional_y <- vector()
  fractional_z <- vector()
  x_errors <- vector()
  y_errors <- vector()
  z_errors <- vector()

  # Find the start and end of the atom section
  atom_start <- grep("_atom_site_occupancy", cif_content$V1)
  atom_end <- grep("^loop_|^#End", cif_content$V1[atom_start:length(cif_content$V1)]) + atom_start - 1

  # If atom_end has multiple matches, take the first match after atom_start
  if(length(atom_end) > 1) {
    atom_end <- atom_end[1]
  }

  # Iterate over each line containing atomic coordinates
  for (line in cif_content$V1[(atom_start + 1):(atom_end - 1)]) {
    # Extract atom properties from the line
    properties <- strsplit(line, "\\s+")[[1]]

    # Extract and store relevant information
    labels <- c(labels, properties[1])

    # Extract fractional coordinates and errors separately
    fractional_x_val <- gsub("\\(.*\\)", "", properties[5])
    fractional_y_val <- gsub("\\(.*\\)", "", properties[6])
    fractional_z_val <- gsub("\\(.*\\)", "", properties[7])

    x_error <- ifelse(grepl("\\(", properties[5]), gsub(".*\\((.*)\\).*", "\\1", properties[5]), NA)
    y_error <- ifelse(grepl("\\(", properties[6]), gsub(".*\\((.*)\\).*", "\\1", properties[6]), NA)
    z_error <- ifelse(grepl("\\(", properties[7]), gsub(".*\\((.*)\\).*", "\\1", properties[7]), NA)

    fractional_x <- c(fractional_x, fractional_x_val)
    fractional_y <- c(fractional_y, fractional_y_val)
    fractional_z <- c(fractional_z, fractional_z_val)
    x_errors <- c(x_errors, x_error)
    y_errors <- c(y_errors, y_error)
    z_errors <- c(z_errors, z_error)
  }

  # Create a data table to store the atomic coordinates
  atomic_coordinates <- data.table(
    Label = labels,
    x_a = as.numeric(fractional_x),
    y_b = as.numeric(fractional_y),
    z_c = as.numeric(fractional_z),
    x_error = as.numeric(x_errors),
    y_error = as.numeric(y_errors),
    z_error = as.numeric(z_errors)
  )
  return(atomic_coordinates)
}
