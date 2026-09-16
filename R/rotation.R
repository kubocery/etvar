# ==============================================================================
# Generalized Euler / Givens rotation matrix
# ==============================================================================

#' Construct an orthogonal rotation matrix from Givens angles
#'
#' Builds an `n`-dimensional orthogonal rotation matrix by composing
#' `n*(n-1)/2` elementary Givens rotations, one for each pair of coordinate
#' planes. This parameterization maps a vector of unconstrained angles to the
#' space of orthogonal matrices. It remains available for constructing legacy
#' full-rotation starting values; [etvar()] itself optimises only the selected
#' target-shock direction.
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


# Target-direction parameterisation ------------------------------------------

# Convert a non-zero target direction to the canonical angles used by the
# target-specific Givens sequence below.
.etvar_angles_from_direction <- function(direction, shock_index) {
  if (!is.numeric(direction) || length(direction) < 1L || anyNA(direction) ||
      any(!is.finite(direction))) {
    stop("`direction` must be a finite numeric vector.", call. = FALSE)
  }

  direction_norm <- sqrt(sum(direction^2))
  if (!is.finite(direction_norm) ||
      direction_norm <= sqrt(.Machine$double.eps)) {
    stop("`direction` must have non-zero length.", call. = FALSE)
  }
  direction <- direction / direction_norm
  n <- length(direction)
  shock_index <- .etvar_integer_scalar(shock_index, "shock_index",
                                        min = 1L, max = n)
  if (n == 1L) return(numeric())

  other_indices <- setdiff(seq_len(n), shock_index)
  theta <- numeric(n - 1L)
  theta[1L] <- atan2(
    direction[other_indices[1L]], direction[shock_index]
  )
  if (n > 2L) {
    for (i in seq.int(2L, n - 1L)) {
      previous_indices <- c(shock_index, other_indices[seq_len(i - 1L)])
      previous_norm <- sqrt(sum(direction[previous_indices]^2))
      theta[i] <- atan2(direction[other_indices[i]], previous_norm)
    }
  }
  theta
}


.etvar_target_rotation <- function(theta, n, shock_index) {
  n <- .etvar_integer_scalar(n, "n", min = 1L)
  shock_index <- .etvar_integer_scalar(shock_index, "shock_index",
                                        min = 1L, max = n)
  expected <- n - 1L
  if (!is.numeric(theta) || anyNA(theta) || any(!is.finite(theta)) ||
      length(theta) != expected) {
    stop("Target-direction angles must be a finite numeric vector of length ",
         expected, ".", call. = FALSE)
  }
  if (n == 1L) return(matrix(1, 1L, 1L))

  rotation <- diag(n)
  other_indices <- setdiff(seq_len(n), shock_index)
  for (i in seq_along(other_indices)) {
    other <- other_indices[i]
    givens <- diag(n)
    givens[shock_index, shock_index] <- cos(theta[i])
    givens[shock_index, other] <- -sin(theta[i])
    givens[other, shock_index] <- sin(theta[i])
    givens[other, other] <- cos(theta[i])
    rotation <- rotation %*% givens
  }
  rotation
}


.etvar_direction_from_angles <- function(theta, n, shock_index) {
  rotation <- .etvar_target_rotation(theta, n, shock_index)
  drop(rotation[, shock_index])
}


# Construct the deterministic special-orthogonal completion associated with a
# target direction. The non-target columns are not identified by the objective.
.etvar_complete_target_rotation <- function(direction, shock_index) {
  theta <- .etvar_angles_from_direction(direction, shock_index)
  .etvar_target_rotation(theta, length(direction), shock_index)
}


# Accept the new n - 1 target-direction angles and legacy full Euler angles.
# In the only ambiguous non-trivial case (n = 2), "auto" prefers the new
# target-direction representation; callers can request "full" explicitly.
.etvar_prepare_target_start <- function(start, n, shock_index,
                                        parameterization = c(
                                          "auto", "target", "full"
                                        )) {
  parameterization <- match.arg(parameterization)
  n <- .etvar_integer_scalar(n, "n", min = 1L)
  shock_index <- .etvar_integer_scalar(shock_index, "shock_index",
                                        min = 1L, max = n)
  if (!is.numeric(start) || anyNA(start) || any(!is.finite(start))) {
    stop("`start` must be a finite numeric vector.", call. = FALSE)
  }

  n_target <- n - 1L
  n_full <- n * (n - 1L) / 2L
  selected <- parameterization
  if (selected == "auto") {
    if (length(start) == n_target) {
      selected <- "target"
    } else if (length(start) == n_full) {
      selected <- "full"
    } else {
      stop("`start` must have length ", n_target,
           " (target direction) or ", n_full,
           " (legacy full rotation).", call. = FALSE)
    }
  }

  expected <- if (selected == "target") n_target else n_full
  if (length(start) != expected) {
    stop("`start` must have length ", expected, " for `start_parameterization = \"",
         selected, "\"`.", call. = FALSE)
  }

  if (selected == "target") return(as.numeric(start))

  legacy_rotation <- euler_rotation(start, n)
  .etvar_angles_from_direction(legacy_rotation[, shock_index], shock_index)
}
