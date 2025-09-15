#' Atomic Radii Data for Bond-Length Estimation
#'
#' A data.table containing atomic radii for elements. This data is used for
#' estimating plausible bond lengths to filter out non-physical "ghost" distances
#' that can occur in disordered structures. Radii are in Angstroms (Å).
#'
#' @format A data table with three columns:
#' \describe{
#' \item{Symbol}{The chemical symbol of the element.}
#' \item{Radius}{The atomic radius in Angstroms.}
#' \item{Type}{The type of radius, e.g., "covalent". The default table only contains covalent radii.}
#' }
#' @source J. Emsley. The Elements. Third edition 1998, Oxford University Press.
#' As provided by Julia-Maria Huebner.
#' @export
"covalent_radii"

#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # This makes the data available for use within the package's internal functions.
  utils::data("covalent_radii", package = pkgname, envir = parent.env(environment()))
}
