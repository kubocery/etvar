# ==============================================================================
# Generalized Euler / Givens rotation matrix
# ==============================================================================

#' Construct an orthogonal rotation matrix from Givens angles
#'
#' Builds an `n`-dimensional orthogonal rotation matrix by composing
#' `n*(n-1)/2` elementary Givens rotations, one for each pair of coordinate
#' planes. This parameterization maps a vector of unconstrained angles to the
#' space of orthogonal matrices, which is used by [etvar()] to search for the
#' structural identification.
#'
#' @param theta Numeric vector of length `n*(n-1)/2`. Rotation angles (radians).
#' @param n     Integer. Dimension of the rotation matrix.
#' @return Matrix (n x n). Orthogonal rotation matrix with determinant +1.
#' @export
#' @examples
#' Q <- euler_rotation(c(0.5, -0.3, 1.2), n = 3)
#' round(t(Q) %*% Q, 10)  # Should be identity
euler_rotation <- function(theta, n) {
  n <- .etvar_integer_scalar(n, "n", min = 1L)
  if (!is.numeric(theta) || anyNA(theta) || any(!is.finite(theta))) {
    stop("`theta` must be a finite numeric vector.", call. = FALSE)
  }

  expected <- n * (n - 1) / 2
  if (length(theta) != expected) {
    stop("Length of `theta` must be n * (n - 1) / 2 = ", expected,
         " for n = ", n, "; got ", length(theta), ".", call. = FALSE)
  }

  R_mat     <- diag(nrow = n)
  angle_idx <- 1L

  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      Rij       <- diag(nrow = n)
      ang       <- theta[angle_idx]
      Rij[i, i] <-  cos(ang)
      Rij[i, j] <- -sin(ang)
      Rij[j, i] <-  sin(ang)
      Rij[j, j] <-  cos(ang)
      R_mat     <- R_mat %*% Rij
      angle_idx <- angle_idx + 1L
    }
  }

  R_mat
}
