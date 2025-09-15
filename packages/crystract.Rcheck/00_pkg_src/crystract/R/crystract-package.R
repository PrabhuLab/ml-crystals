#' crystract: Tools to Extract Data from CIF Files for Crystallography
#'
#' @description
#' Provides a suite of functions to parse Crystallographic Information
#' Files (.cif), extracting essential data such as chemical formulas, unit cell
#' parameters, atomic coordinates, and symmetry operations. It also includes
#' tools to calculate interatomic distances, identify bonded pairs using various
#' algorithms (Minimum Distance, Brunner's, Hoppe's), determine nearest
#' neighbor counts, and calculate bond angles. The package is designed to
#' facilitate the preparation of crystallographic data for further analysis,
#' including machine learning applications in materials science.
#'
#' @keywords internal
"_PACKAGE"
#'
#' @importFrom data.table .N := .SD as.data.table copy data.table fread rbindlist set setkey setnames
#' @importFrom dplyr distinct
#' @importFrom stats na.omit
#' @importFrom stringr str_match str_split_fixed
#' @importFrom utils combn globalVariables
NULL
