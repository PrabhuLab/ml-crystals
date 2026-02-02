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
#' @description Looks up Pauling electronegativity. matches pymatgen `Element.X`.
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
#' @description mimic pymatgen `_get_radius` logic.
#' 1. If oxidation state is provided and non-zero: Look up Shannon Radius.
#' 2. If Shannon radius not found for that CN, try to average or find closest.
#' 3. If oxidation state is 0 or NA: Look up Covalent/Atomic Radius.
#'
#' @param symbol Character. Element symbol.
#' @param oxidation_state Numeric. Formal charge.
#' @param cn Numeric. Coordination number (optional).
#' @return Numeric radius in Angstroms.
#' @noRd
get_radius <- function(symbol,
                       oxidation_state = NA,
                       cn = 6) {
  # Atomic/Covalent (Fallback)
  atomic_rads <- get_radii_data()

  # Shannon (Ionic) - Load from the provided CSV
  shannon_rads <- load_ext_data("Shannon_Radii.csv")
  if (!is.null(shannon_rads)) {
    # Rename columns to match internal function expectations
    data.table::setnames(
      shannon_rads,
      old = c("Element", "Charge", "Ionic Radius"),
      new = c("Symbol", "OxidationState", "Radius"),
      skip_absent = TRUE
    )
  }


  radius <- 1.0 # Ultimate fallback

  if (!is.na(oxidation_state) &&
      oxidation_state != 0 && !is.null(shannon_rads)) {
    # Try exact match for Symbol + OS + CN
    match <- shannon_rads[Symbol == symbol &
                            OxidationState == oxidation_state &
                            Coordination == cn]

    if (nrow(match) > 0) {
      radius <- match$Radius[1]
    } else {
      # Relax CN: Find matches for Symbol + OS, pick closest CN
      matches <- shannon_rads[Symbol == symbol &
                                OxidationState == oxidation_state]
      if (nrow(matches) > 0) {
        matches[, Diff := abs(as.numeric(Coordination) - cn)]
        radius <- matches[order(Diff)]$Radius[1]
      } else {
        # No match for this OS? Fallback to atomic.
        r_atom <- atomic_rads[Symbol == symbol &
                                Type == "covalent"]$Radius
        if (length(r_atom) > 0)
          radius <- r_atom[1]
      }
    }
  } else {
    # Neutral or no OS provided
    r_atom <- atomic_rads[Symbol == symbol &
                            Type == "covalent"]$Radius
    if (length(r_atom) > 0)
      radius <- r_atom[1]
  }

  return(radius)
}
