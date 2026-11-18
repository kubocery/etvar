#' Orthogonal rotation from Givens angles
#'
#' Builds an n x n rotation matrix as the product of the `n * (n - 1) / 2`
#' elementary Givens rotations, one for each pair of coordinates, in the order
#' (1, 2), (1, 3), ..., (n - 1, n). This full parameterisation is what the
#' original MATLAB implementation of ETVAR optimised. [etvar()] itself only
#' estimates the target column of the rotation, but rotations built here can
#' still be used as starting values (`start_parameterization = "full"`).
#'
#' @param theta Numeric vector of `n * (n - 1) / 2` angles, in radians.
#' @param n Integer. Dimension of the rotation.
#' @return An n x n orthogonal matrix with determinant 1.
#' @examples
#' Q <- euler_rotation(c(0.5, -0.3, 1.2), n = 3)
#' round(crossprod(Q), 10)
#' @export
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

  rotation <- diag(nrow = n)
  k <- 0L
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      k <- k + 1L
      rotation <- rotation %*% givens(n, i, j, theta[k])
    }
  }
  rotation
}


# Target-direction parameterisation --------------------------------------------
#
# ETVAR identifies only the target column q of the rotation. q is described by
# n - 1 angles: starting from the unit vector of the target shock, it is
# rotated towards each of the other coordinates in turn. The remaining columns
# come from the same product of rotations; they complete q to a rotation but
# are not identified by the event criterion.

givens <- function(n, i, j, angle) {
  g <- diag(nrow = n)
  g[i, i] <- cos(angle)
  g[i, j] <- -sin(angle)
  g[j, i] <- sin(angle)
  g[j, j] <- cos(angle)
  g
}

.etvar_target_rotation <- function(theta, n, shock_index) {
  n <- .etvar_integer_scalar(n, "n", min = 1L)
  shock_index <- .etvar_integer_scalar(shock_index, "shock_index",
                                       min = 1L, max = n)
  if (!is.numeric(theta) || anyNA(theta) || any(!is.finite(theta)) ||
      length(theta) != n - 1L) {
    stop("Target-direction angles must be a finite numeric vector of length ",
         n - 1L, ".", call. = FALSE)
  }

  rotation <- diag(nrow = n)
  others <- setdiff(seq_len(n), shock_index)
  for (i in seq_along(others)) {
    rotation <- rotation %*% givens(n, shock_index, others[i], theta[i])
  }
  rotation
}

.etvar_direction_from_angles <- function(theta, n, shock_index) {
  drop(.etvar_target_rotation(theta, n, shock_index)[, shock_index])
}

# Inverse of .etvar_direction_from_angles(): the angles of a (non-zero)
# direction, which is normalised to unit length first.
.etvar_angles_from_direction <- function(direction, shock_index) {
  if (!is.numeric(direction) || length(direction) < 1L || anyNA(direction) ||
      any(!is.finite(direction))) {
    stop("`direction` must be a finite numeric vector.", call. = FALSE)
  }
  size <- sqrt(sum(direction^2))
  if (size <= sqrt(.Machine$double.eps)) {
    stop("`direction` must have non-zero length.", call. = FALSE)
  }
  direction <- direction / size
  n <- length(direction)
  shock_index <- .etvar_integer_scalar(shock_index, "shock_index",
                                       min = 1L, max = n)
  if (n == 1L) {
    return(numeric())
  }

  others <- setdiff(seq_len(n), shock_index)
  theta <- numeric(n - 1L)
  theta[1L] <- atan2(direction[others[1L]], direction[shock_index])
  for (i in seq_len(n - 1L)[-1L]) {
    done <- c(shock_index, others[seq_len(i - 1L)])
    theta[i] <- atan2(direction[others[i]], sqrt(sum(direction[done]^2)))
  }
  theta
}

.etvar_complete_target_rotation <- function(direction, shock_index) {
  theta <- .etvar_angles_from_direction(direction, shock_index)
  .etvar_target_rotation(theta, length(direction), shock_index)
}

# Turns user-supplied angles into target-direction angles. Angles of a full
# rotation are converted by taking the target column. "auto" decides from the
# length; the only ambiguous case, n = 2, is read as target angles.
.etvar_prepare_target_start <- function(start, n, shock_index,
                                        parameterization = c("auto", "target",
                                                             "full")) {
  parameterization <- match.arg(parameterization)
  n <- .etvar_integer_scalar(n, "n", min = 1L)
  shock_index <- .etvar_integer_scalar(shock_index, "shock_index",
                                       min = 1L, max = n)
  if (!is.numeric(start) || anyNA(start) || any(!is.finite(start))) {
    stop("`start` must be a finite numeric vector.", call. = FALSE)
  }

  n_target <- n - 1L
  n_full <- n * (n - 1L) / 2L
  if (parameterization == "auto") {
    parameterization <- if (length(start) == n_target) {
      "target"
    } else if (length(start) == n_full) {
      "full"
    } else {
      stop("`start` must have length ", n_target, " (target direction) or ",
           n_full, " (full rotation).", call. = FALSE)
    }
  }

  expected <- if (parameterization == "target") n_target else n_full
  if (length(start) != expected) {
    stop("`start` must have length ", expected,
         " for `start_parameterization = \"", parameterization, "\"`.",
         call. = FALSE)
  }

  if (parameterization == "target") {
    return(as.numeric(start))
  }
  full_rotation <- euler_rotation(start, n)
  .etvar_angles_from_direction(full_rotation[, shock_index], shock_index)
}
