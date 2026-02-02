#' @title Load Internal Data (Internal)
#' @description Loads data from inst/extdata or package data.
#' @param filename String. Name of file in inst/extdata.
#' @noRd
load_ext_data <- function(filename) {
  path <- system.file("extdata", filename, package = "crystract")
  if (path == "") {
    # Fallback for dev mode
    path <- file.path("inst", "extdata", filename)
  }
  if (!file.exists(path))
    return(NULL)
  data.table::fread(path)
}

#' @title Get Electronegativity
#' @description Looks up Pauling electronegativity. Matches pymatgen `Element.X`.
#' @param symbols Character vector of element symbols.
#' @return Numeric vector of electronegativities.
#' @noRd
get_electronegativity <- function(symbols) {
  # Load data from the provided CSV
  en_data <- load_ext_data("pauling_electronegativity_stable_elements.csv")
  if (is.null(en_data)) {
    warning("pauling_electronegativity_stable_elements.csv not found.")
    return(rep(NA_real_, length(symbols)))
  }

  vals <- en_data$PaulingElectronegativity[match(symbols, en_data$Symbol)]
  return(vals)
}

#' @title Get Atomic/Ionic Radius
#' @description Mimics pymatgen `_get_radius` logic.
#' 1. If oxidation state is provided and non-zero: Look up Shannon Radius.
#' 2. If fractional oxidation state: Linear interpolation between floor/ceil states.
#' 3. If no specific match: Use average cationic/anionic radius for that element.
#' 4. Fallback (or if OS is 0/NA): Use Covalent/Atomic Radius.
#'
#' @param symbol Character. Element symbol.
#' @param oxidation_state Numeric. Formal charge.
#' @param cn Numeric. Coordination number (optional, default 6).
#' @return Numeric radius in Angstroms.
#' @noRd
get_radius <- function(symbol,
                       oxidation_state = NA,
                       cn = 6) {
  # 1. Load Covalent/Atomic Radii (Ultimate Fallback)
  atomic_rads <- get_radii_data()

  # Helper for default radius
  get_default_radius <- function(sym) {
    r_val <- atomic_rads[Symbol == sym & Type == "covalent"]$Radius
    if (length(r_val) > 0)
      return(r_val[1])
    return(0)
  }

  # 2. Handle 0 or NA Oxidation State -> Return Default immediately
  if (is.na(oxidation_state) || oxidation_state == 0) {
    return(get_default_radius(symbol))
  }

  # 3. Load Shannon Radii
  shannon_rads <- load_ext_data("Shannon_Radii.csv")
  if (is.null(shannon_rads)) {
    return(get_default_radius(symbol))
  }

  # Normalize columns
  data.table::setnames(
    shannon_rads,
    old = c("Element", "Charge", "Ionic Radius"),
    new = c("Symbol", "OxidationState", "Radius"),
    skip_absent = TRUE
  )

  # Helper to fetch radius for a specific integer/float OS
  fetch_radius_for_os <- function(sym, os) {
    # Try exact match for Symbol + OS + CN
    match <- shannon_rads[Symbol == sym &
                            OxidationState == os &
                            Coordination == cn]
    if (nrow(match) > 0)
      return(match$Radius[1])

    # Relax CN: Find matches for Symbol + OS, pick closest CN
    matches <- shannon_rads[Symbol == sym &
                              OxidationState == os]
    if (nrow(matches) > 0) {
      matches[, Diff := abs(as.numeric(Coordination) - cn)]
      return(matches[order(Diff)]$Radius[1])
    }
    return(NA_real_)
  }

  # 4. Check Exact Match
  r_exact <- fetch_radius_for_os(symbol, oxidation_state)
  if (!is.na(r_exact))
    return(r_exact)

  # 5. Handle Fractional Oxidation State (Linear Interpolation)
  # Matches pymatgen logic: if math.floor(oxi) and math.ceil(oxi) in radii...
  os_floor <- floor(oxidation_state)
  os_ceil <- ceiling(oxidation_state)

  if (os_floor != os_ceil) {
    r_floor <- fetch_radius_for_os(symbol, os_floor)
    r_ceil <- fetch_radius_for_os(symbol, os_ceil)

    if (!is.na(r_floor) && !is.na(r_ceil)) {
      x <- oxidation_state - os_floor
      return((1 - x) * r_floor + x * r_ceil)
    }
  }

  # 6. Fallback: Average Cationic / Anionic Radius
  # Matches pymatgen logic for `average_cationic_radius` / `average_anionic_radius`
  if (oxidation_state > 0) {
    avg_cat <- mean(shannon_rads[Symbol == symbol &
                                   OxidationState > 0]$Radius, na.rm = TRUE)
    if (!is.nan(avg_cat))
      return(avg_cat)
  } else if (oxidation_state < 0) {
    avg_an <- mean(shannon_rads[Symbol == symbol &
                                  OxidationState < 0]$Radius, na.rm = TRUE)
    if (!is.nan(avg_an))
      return(avg_an)
  }

  # 7. Final Fallback
  return(get_default_radius(symbol))
}
