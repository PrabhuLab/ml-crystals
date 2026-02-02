#' @title Calculate Cross Product (Internal)
#' @description Calculates the cross product of two 3D vectors.
#' @param v1 Numeric vector of length 3.
#' @param v2 Numeric vector of length 3.
#' @return Numeric vector of length 3.
#' @noRd
cross_product <- function(v1, v2) {
  c(v1[2] * v2[3] - v1[3] * v2[2], v1[3] * v2[1] - v1[1] * v2[3], v1[1] * v2[2] - v1[2] * v2[1])
}

#' @title Calculate Tetrahedron Volume (Internal)
#' @description Calculates the volume of a tetrahedron given 4 vertices.
#' Used for reconstruction of Voronoi face areas.
#' Matches pymatgen `vol_tetra`.
#' @param vt1,vt2,vt3,vt4 Numeric vectors of length 3 (vertices).
#' @return Numeric volume.
#' @noRd
vol_tetra <- function(vt1, vt2, vt3, vt4) {
  # Formula: |dot((v1 - v4), cross((v2 - v4), (v3 - v4)))| / 6
  diff1 <- vt1 - vt4
  diff2 <- vt2 - vt4
  diff3 <- vt3 - vt4
  cp <- cross_product(diff2, diff3)
  abs(sum(diff1 * cp)) / 6
}

#' @title Calculate Circumcenter (Internal)
#' @description Calculates the circumcenter of a tetrahedron defined by 4 points.
#' @param pts A 4x3 matrix of coordinates.
#' @return Numeric vector of length 3.
#' @noRd
calculate_circumcenter <- function(pts) {
  # Relative vectors from point 1
  r2 <- pts[2, ] - pts[1, ]
  r3 <- pts[3, ] - pts[1, ]
  r4 <- pts[4, ] - pts[1, ]

  # Squared magnitudes
  sq2 <- sum(r2^2)
  sq3 <- sum(r3^2)
  sq4 <- sum(r4^2)

  # Linear system Ax = b
  A <- rbind(r2, r3, r4)
  b <- 0.5 * c(sq2, sq3, sq4)

  # Solve with error handling for degenerate tetrahedra
  x <- tryCatch(
    solve(A, b),
    error = function(e)
      rep(NaN, 3)
  )

  return(pts[1, ] + x)
}

#' @title Calculate Solid Angle (Internal)
#' @description Calculates the solid angle subtended by a polygon face from a center point,
#' following the Oosterom & Strackee method used in pymatgen `solid_angle`.
#' @param center Numeric vector of length 3 (center atom coords).
#' @param face_coords Matrix of face vertices (Nx3).
#' @return Numeric solid angle in steradians.
#' @noRd
calculate_solid_angle <- function(center, face_coords) {
  # Compute displacement from center
  disp <- t(t(face_coords) - center)
  r_norm <- sqrt(rowSums(disp^2))

  angle <- 0
  n <- nrow(face_coords)

  if (n < 3)
    return(0)

  # Iterate through the triangle fan (1, i, i+1 logic)
  # Pymatgen logic: "for ii in range(1, len(disp) - 1): jj = ii + 1"
  # This corresponds to indices 0, ii, jj in 0-indexed Python arrays.
  # In R (1-based), this is indices 1, ii, jj where ii runs 2 to n-1.

  for (ii in 2:(n - 1)) {
    jj <- ii + 1

    r0 <- disp[1, ] # Python index 0
    ri <- disp[ii, ]
    rj <- disp[jj, ]

    n0 <- r_norm[1]
    ni <- r_norm[ii]
    nj <- r_norm[jj]

    # Scalar triple product: r0 . (ri x rj)
    cp <- cross_product(ri, rj)
    tp <- abs(sum(r0 * cp))

    # Denominator (Oosterom & Strackee)
    de <- n0 * ni * nj +
      nj * sum(r0 * ri) +
      ni * sum(r0 * rj) +
      n0 * sum(ri * rj)

    if (de == 0) {
      term <- if (tp > 0)
        0.5 * pi
      else
        - 0.5 * pi
    } else {
      term <- atan(tp / de)
    }

    contrib <- if (term > 0)
      term
    else
      term + pi
    angle <- angle + (contrib * 2)
  }
  return(angle)
}

#' @title Semicircle Integral (Internal)
#' @description Integral between two bounds of a unit semicircle.
#' Used in CrystalNN for fingerprint/CN probability calculation.
#' @param x1 Lower bound.
#' @param x2 Upper bound.
#' @param r Radius (default 1).
#' @return Numeric integral value.
#' @noRd
semicircle_integral <- function(x1, x2, r = 1) {
  # Helper for indefinite integral area
  calc_area <- function(x, rad) {
    if (abs(x) >= rad)
      return(0.25 * pi * rad^2 * sign(x))
    0.5 * (x * sqrt(rad^2 - x^2) + rad^2 * atan(x / sqrt(rad^2 - x^2)))
  }

  # Handle full quadrant cases
  area1 <- if (x1 >= r)
    0.25 * pi * r^2
  else
    calc_area(x1, r)
  area2 <- if (x2 >= r)
    0.25 * pi * r^2
  else
    calc_area(x2, r)

  return((area1 - area2) / (0.25 * pi * r^2))
}
