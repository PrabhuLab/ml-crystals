#' crysmal: Tools to Extract Data From CIF Files for Crystallography
#'
#' Provides a suite of functions to parse Crystallographic Information
#' Files (.cif), extracting essential data such as chemical formulas, unit cell
#' parameters, atomic coordinates, and symmetry operations. It also includes
#' tools to calculate interatomic distances, identify bonded pairs using various
#' algorithms (Minimum Distance, Brunner's, Hoppe's), determine nearest
#' neighbor counts, and calculate bond angles. The package is designed to
#' facilitate the preparation of crystallographic data for further analysis,
#' including machine learning applications in materials science.
#'
#' @docType package
#' @name crysmal-package
#' @aliases crysmal
#' @import data.table
#' @importFrom stringr str_split_fixed str_extract
#' @importFrom dplyr distinct
#' @importFrom utils combn
NULL
