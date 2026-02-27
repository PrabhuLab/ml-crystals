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

  # fill=TRUE handles the footer warning
  data.table::fread(path, fill = TRUE)
}

#' @title Get Cached Shannon Radii (Internal)
#' @description Loads, cleans, and caches the Shannon Radii table.
#' @noRd
get_shannon_radii <- function() {
  # 1. return cached version if exists
  if (exists("shannon_radii_cache", envir = .crystract_env)) {
    return(get("shannon_radii_cache", envir = .crystract_env))
  }

  # 2. Load from disk
  dt <- load_ext_data("Shannon_Radii.csv")

  if (is.null(dt))
    return(NULL)

  # 3. Robust Column Renaming
  # Remove spaces/punctuation to handle "Ionic Radius" vs "IonicRadius"
  current_names <- names(dt)
  clean_names <- gsub("[^A-Za-z0-9]", "", current_names)

  # Define mapping based on Clean Names -> Standard Names
  target_map <- c("Element" = "Symbol",
                  "Charge" = "OxidationState",
                  "IonicRadius" = "Radius")

  for (clean_n in names(target_map)) {
    idx <- which(clean_names == clean_n)
    if (length(idx) > 0) {
      setnames(dt, current_names[idx[1]], target_map[[clean_n]])
    }
  }

  # 4. Clean Coordination Column (Roman -> Integer)
  # Many Shannon tables use VI, IV, etc.
  if ("Coordination" %in% names(dt)) {
    # Create a mapping for common Roman numerals found in Shannon data
    roman_map <- c(
      "I" = 1,
      "II" = 2,
      "III" = 3,
      "IV" = 4,
      "V" = 5,
      "VI" = 6,
      "VII" = 7,
      "VIII" = 8,
      "IX" = 9,
      "X" = 10,
      "XI" = 11,
      "XII" = 12,
      "XIII" = 13,
      "XIV" = 14
    )

    # Function to convert single value
    clean_coord <- function(x) {
      # Remove distinct markers like 'sq', 'py' sometimes found in Shannon datasets
      x_clean <- gsub("[^A-Za-z0-9]", "", as.character(x))

      # Check if it's already a number
      if (grepl("^[0-9]+$", x_clean))
        return(as.numeric(x_clean))

      # Check Roman map
      if (x_clean %in% names(roman_map))
        return(roman_map[[x_clean]])

      # Fallback: Try removing letters and parsing (e.g. "4sq" -> 4)
      num_part <- gsub("[^0-9]", "", as.character(x))
      if (nchar(num_part) > 0)
        return(as.numeric(num_part))

      return(NA_real_)
    }

    # Apply conversion
    dt[, Coordination := sapply(Coordination, clean_coord)]
  }

  # 4. Cache and return
  assign("shannon_radii_cache", dt, envir = .crystract_env)
  return(dt)
}

#' @title Calculate Cross Product (Internal)
#' @description Calculates the cross product of two 3D vectors.
#' @param v1 Numeric vector of length 3.
#' @param v2 Numeric vector of length 3.
#' @return Numeric vector of length 3.
#' @noRd
cross_product <- function(v1, v2) {
  c(v1[2] * v2[3] - v1[3] * v2[2], v1[3] * v2[1] - v1[1] * v2[3], v1[1] * v2[2] - v1[2] * v2[1])
}

#' @title Calculate Circumcenter (Internal)
#' @description Calculates the circumcenter of a tetrahedron defined by 4 points.
#'   Includes robust checks for coplanar/degenerate vertices which cause
#'   numerical instability in Voronoi tessellation.
#' @param pts A 4x3 matrix of coordinates.
#' @return Numeric vector of length 3. Returns `rep(NA_real_, 3)` if the
#'   tetrahedron is degenerate (volume ~ 0).
#' @noRd
calculate_circumcenter <- function(pts) {
  # Relative vectors from point 1
  r2 <- pts[2, ] - pts[1, ]
  r3 <- pts[3, ] - pts[1, ]
  r4 <- pts[4, ] - pts[1, ]

  # Construct matrix A (rows are difference vectors)
  A <- rbind(r2, r3, r4)

  # Check for degeneracy using determinant
  # The determinant is 6 * Volume of tetrahedron.
  det_A <- det(A)

  # Tolerance matches pymatgen/scipy qhull sensitivity defaults
  if (abs(det_A) < 1e-10) {
    return(rep(NA_real_, 3))
  }

  # Squared magnitudes
  sq2 <- sum(r2^2)
  sq3 <- sum(r3^2)
  sq4 <- sum(r4^2)

  b <- 0.5 * c(sq2, sq3, sq4)

  # Solve linear system Ax = b with fallback
  x <- tryCatch(
    solve(A, b),
    error = function(e)
      rep(NA_real_, 3)
  )

  if (anyNA(x))
    return(rep(NA_real_, 3))

  return(pts[1, ] + x)
}

#' @title Calculate Solid Angle (Internal)
#' @description Calculates the solid angle subtended by a polygon face from a center point,
#' following the method used in pymatgen (Oosterom & Strackee).
#' @param center Numeric vector of length 3 (center atom coords).
#' @param face_coords Matrix of face vertices (Nx3).
#' @return Numeric solid angle in steradians.
#' @noRd
calculate_solid_angle <- function(center, face_coords) {
  # Calculate displacement vectors from center to vertices
  r <- t(t(face_coords) - center)
  r_norms <- sqrt(rowSums(r^2))

  angle <- 0
  n <- nrow(face_coords)

  if (n < 3)
    return(0)

  # Iterate through the triangle fan (0, i, i+1 logic)
  for (i in 2:(n - 1)) {
    r0 <- r[1, ]
    ri <- r[i, ]
    rj <- r[i + 1, ]

    n0 <- r_norms[1]
    ni <- r_norms[i]
    nj <- r_norms[i + 1]

    # Scalar triple product: r0 . (ri x rj)
    cp <- cross_product(ri, rj)
    tp <- abs(sum(r0 * cp))

    # Denominator (Oosterom & Strackee)
    denom <- n0 * ni * nj +
      nj * sum(r0 * ri) +
      ni * sum(r0 * rj) +
      n0 * sum(ri * rj)

    if (denom == 0) {
      term <- if (tp > 0)
        0.5 * pi
      else-0.5 * pi
    } else {
      term <- atan(tp / denom)
    }

    contrib <- if (term > 0)
      term
    else
      term + pi
    angle <- angle + (contrib * 2)
  }
  return(angle)
}

#' @title Check if Elements are Metals
#' @description Checks a vector of symbols against a standard list of metals.
#' Based on Science Notes list.
#' @param symbols Character vector of element symbols.
#' @return Logical vector.
#' @noRd
is_metal <- function(symbols) {
  metals <- c(
    "Li",
    "Be",
    "Na",
    "Mg",
    "Al",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og"
  )
  return(symbols %in% metals)
}

#' @title Get Electronegativity
#' @description Looks up Pauling electronegativity. Matches pymatgen `Element.X`.
#' @param symbols Character vector of element symbols.
#' @return Numeric vector of electronegativities.
#' @noRd
get_electronegativity <- function(symbols) {
  vals <- pauling_en$PaulingElectronegativity[match(symbols, pauling_en$Symbol)]
  return(vals)
}

#' @title Get Atomic/Ionic Radius (Pymatgen Style - For EconNN)
#' @description Implements Pymatgen's logic:
#' 1. Ionic (Shannon) if Oxi != 0.
#' 2. Covalent (Cordero) if Oxi == 0 or Ionic not found.
#' 3. Atomic if Covalent not found.
#' Uses cached data to prevent disk I/O on every call.
#' @param symbol Character. Element symbol.
#' @param oxidation_state Numeric. Formal charge.
#' @param cn Numeric. Coordination number.
#' @return Numeric radius in Angstroms.
#' @noRd
get_radius_pymatgen_style <- function(symbol,
                                      oxidation_state = NA,
                                      cn = 6) {
  atomic_rads <- get_radii_data()
  # USE CACHED LOADER
  shannon_rads <- get_shannon_radii()

  # --- CHECK 1: Oxidation State ---
  if (!is.na(oxidation_state) &&
      oxidation_state != 0 && !is.null(shannon_rads)) {
    # Exact match for Symbol + OS + CN
    match <- shannon_rads[Symbol == symbol &
                            OxidationState == oxidation_state &
                            Coordination == cn]

    if (nrow(match) > 0)
      return(match$Radius[1])

    # Relax CN: Find matches for Symbol + OS, pick closest CN
    matches <- shannon_rads[Symbol == symbol &
                              OxidationState == oxidation_state]
    if (nrow(matches) > 0) {
      # The coordination column is now numeric, so subtraction is safe
      matches[, Diff := abs(Coordination - cn)]
      return(matches[order(Diff)]$Radius[1])
    }
  }

  # --- CHECK 2: Covalent (Cordero) ---
  r_atom <- atomic_rads[Symbol == symbol &
                          Type == "covalent"]$Radius
  if (length(r_atom) > 0)
    return(r_atom[1])

  # --- CHECK 3: Atomic ---
  radius_fallback <- atomic_rads[Symbol == symbol &
                                   Type == "atomic"]$Radius
  if (length(radius_fallback) > 0)
    return(radius_fallback[1])

  # Pymatgen _get_radius returns 0 if no match found
  return(0)
}

#' @title Get Atomic/Ionic Radius (Decision Tree - For CrystalNN)
#' @description Implements Julia-Maria Huebner's specific Decision Tree logic,
#' aligned with Pymatgen's _get_radius fallback.
#' Uses cached data to prevent disk I/O on every call.
#' @param symbol Character. Element symbol.
#' @param oxidation_state Numeric. Formal charge.
#' @param cn Numeric. Coordination number.
#' @return Numeric radius in Angstroms.
#' @noRd
get_radius_decision_tree <- function(symbol,
                                     oxidation_state = NA,
                                     cn = 6) {
  atomic_rads <- get_radii_data()
  # USE CACHED LOADER
  shannon_rads <- get_shannon_radii()

  # --- CHECK 1: Oxidation State ---
  if (!is.na(oxidation_state) &&
      oxidation_state != 0 && !is.null(shannon_rads)) {
    match <- shannon_rads[Symbol == symbol &
                            OxidationState == oxidation_state &
                            Coordination == cn]
    if (nrow(match) > 0)
      return(match$Radius[1])

    matches <- shannon_rads[Symbol == symbol &
                              OxidationState == oxidation_state]
    if (nrow(matches) > 0) {
      # The coordination column is now numeric, so subtraction is safe
      matches[, Diff := abs(Coordination - cn)]
      return(matches[order(Diff)]$Radius[1])
    }
  }

  # --- CHECK 2: Covalent (Cordero) ---
  # Pymatgen logic falls back to covalent/atomic if ionic not found or oxi=0
  radius_match <- atomic_rads[Symbol == symbol &
                                Type == "covalent"]$Radius
  if (length(radius_match) > 0)
    return(radius_match[1])

  # --- CHECK 3: Atomic ---
  radius_fallback <- atomic_rads[Symbol == symbol &
                                   Type == "atomic"]$Radius
  if (length(radius_fallback) > 0)
    return(radius_fallback[1])

  return(1.0)
}

#' @title Snap to Common Crystallographic Fractions (Internal)
#' @description Vectorized function that snaps decimal values (like 0.3333) to exact
#'   fractions (like 1/3) to avoid floating point precision issues. Common denominators:
#'   1, 2, 3, 4, 6, 8, 12, 24.
#' @param vals Numeric vector of values to snap.
#' @param tolerance Numeric tolerance.
#' @param wrap Logical. If TRUE, wraps values to[0, 1) and safely binds exactly 1.0 to 0.0.
#' @return Numeric vector of snapped values.
#' @noRd
snap_to_fraction <- function(vals,
                             tolerance = 1e-4,
                             wrap = TRUE) {
  if (length(vals) == 0)
    return(vals)

  vals <- as.numeric(vals)

  if (wrap) {
    vals <- vals %% 1.0
  }

  snapped_vals <- vals
  valid_mask <- !is.na(vals)
  snapped_mask <- rep(FALSE, length(vals))

  denominators <- c(1, 2, 3, 4, 6, 8, 12, 24)

  for (d in denominators) {
    to_check <- valid_mask & !snapped_mask
    if (!any(to_check))
      break

    diffs <- abs(vals[to_check] * d - round(vals[to_check] * d))
    matches <- diffs < tolerance

    if (any(matches)) {
      match_idx <- which(to_check)[matches]
      snapped_vals[match_idx] <- round(vals[match_idx] * d) / d
      snapped_mask[match_idx] <- TRUE
    }
  }

  if (wrap) {
    wrap_mask <- valid_mask &
      (abs(snapped_vals - 1.0) < 1e-7 | abs(snapped_vals) < 1e-7)
    snapped_vals[wrap_mask] <- 0.0
  }

  return(snapped_vals)
}

#' @title Merge Partial Occupancy Sites (Internal)
#' @description Identifies atoms in the asymmetric unit that share the same
#'   fractional coordinates and merges them into a single representative site.
#'   This is required for Voronoi and CrystalNN to function with disordered structures.
#' @param atomic_coordinates A `data.table` of atomic coordinates.
#' @param tolerance Numeric. Snap tolerance for identifying overlapping sites.
#' @return A `data.table` with overlapping atoms combined into single labels.
#' @noRd
merge_partial_occupancy_sites <- function(atomic_coordinates, tolerance = 1e-4) {
  if (is.null(atomic_coordinates) ||
      nrow(atomic_coordinates) == 0)
    return(atomic_coordinates)

  # 1. Create grouping keys based on snapped coordinates
  dt <- copy(atomic_coordinates)
  dt[, `:=`(
    gx = snap_to_fraction(x_a, tolerance = tolerance, wrap = TRUE),
    gy = snap_to_fraction(y_b, tolerance = tolerance, wrap = TRUE),
    gz = snap_to_fraction(z_c, tolerance = tolerance, wrap = TRUE)
  )]

  # 2. Check if there are actually any overlapping sites
  overlap_check <- dt[, .N, by = .(gx, gy, gz)]
  if (all(overlap_check$N == 1))
    return(atomic_coordinates)

  # 3. Merge overlapping sites
  collapsed <- dt[, .(
    Label = paste(Label, collapse = "_"),
    TypeSymbol = if (all(is.na(TypeSymbol)))
      NA_character_
    else
      paste(na.omit(TypeSymbol), collapse = "/"),
    WyckoffSymbol = WyckoffSymbol[1],
    WyckoffMultiplicity = WyckoffMultiplicity[1],
    OxidationState = mean(OxidationState, na.rm = TRUE),
    ThermalParam = mean(ThermalParam, na.rm = TRUE),
    Occupancy = sum(Occupancy),
    OccupancyError = if (all(is.na(OccupancyError)))
      NA_real_
    else
      sqrt(sum(OccupancyError^2, na.rm = TRUE)),
    x_a = x_a[1],
    y_b = y_b[1],
    z_c = z_c[1],
    x_error = x_error[1],
    y_error = y_error[1],
    z_error = z_error[1]
  ), by = .(gx, gy, gz)]

  collapsed[, c("gx", "gy", "gz") := NULL]

  return(collapsed)
}
