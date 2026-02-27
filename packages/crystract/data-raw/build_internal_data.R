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
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
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
#' Uses internal sysdata to prevent disk I/O on every call.
#' @param symbol Character. Element symbol.
#' @param oxidation_state Numeric. Formal charge.
#' @param cn Numeric. Coordination number.
#' @return Numeric radius in Angstroms.
#' @noRd
get_radius_pymatgen_style <- function(symbol,
                                      oxidation_state = NA,
                                      cn = 6) {
  atomic_rads <- get_radii_data()
  shannon_rads <- shannon_radii

  # --- CHECK 1: Oxidation State ---
  if (!is.na(oxidation_state) && oxidation_state != 0) {
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
#' Uses internal sysdata to prevent disk I/O on every call.
#' @param symbol Character. Element symbol.
#' @param oxidation_state Numeric. Formal charge.
#' @param cn Numeric. Coordination number.
#' @return Numeric radius in Angstroms.
#' @noRd
get_radius_decision_tree <- function(symbol,
                                     oxidation_state = NA,
                                     cn = 6) {
  atomic_rads <- get_radii_data()
  shannon_rads <- shannon_radii

  # --- CHECK 1: Oxidation State ---
  if (!is.na(oxidation_state) && oxidation_state != 0) {
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
