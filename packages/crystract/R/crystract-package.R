#' crystract: Tools to Extract Data from CIF Files for Crystallography
#'
#' @description
#' Parse Crystallographic Information Files (.cif) and extract essential
#' data such as chemical formulas, unit cell parameters, atomic coordinates,
#' and symmetry operations. Calculate interatomic distances, identify bonded
#' pairs using various algorithms (Minimum Distance, Brunner's, EconNN,
#' CrystalNN, Voronoi), determine nearest neighbor counts, and calculate bond
#' angles. Facilitates the preparation of crystallographic data for further
#' analysis, including machine learning applications in materials science.
#'
#' @keywords internal
"_PACKAGE"
#'
#' @importFrom data.table .N := .SD as.data.table copy data.table fread rbindlist set setkey setnames
#' @importFrom dplyr distinct
#' @importFrom stats as.dist cutree hclust na.omit
#' @importFrom geometry delaunayn
#' @importFrom stringr str_match str_split_fixed
#' @importFrom utils combn globalVariables
NULL
