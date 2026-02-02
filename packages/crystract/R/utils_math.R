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
  # Pymatgen implementation iterates 0, 1, 2... we use 1-based indexing
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
