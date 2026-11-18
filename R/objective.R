#' ETVAR objective function
#'
#' Measures how far the movements implied by a candidate target shock are
#' from the observed movements on the event date(s). For each event date and
#' variable, the change implied by the deterministic component plus the target
#' shock is compared with the actual change, both divided by the standard
#' deviation of that variable's first differences; the objective is the sum
#' of the absolute deviations over variables and event dates.
#'
#' The candidate shock enters only through its unit direction q in the
#' Cholesky-whitened space, and its contribution is a quadratic form in q. The
#' quadratic forms are computed once, so evaluating the objective does not
#' require a historical decomposition. [etvar_optimize()] minimises this
#' function.
#'
#' @param theta Numeric vector of angles. By default the `n - 1` angles of the
#'   target direction; with `parameterization = "full"`, the `n * (n - 1) / 2`
#'   angles of a full Givens rotation (see [euler_rotation()]).
#' @param var_est An `etvar_var` object from [estimate_var()].
#' @param target_obs Integer vector. Row(s) of the event date(s) in the
#'   effective sample, i.e. in `var_est$Y`. With several rows the deviations
#'   are summed over all of them.
#' @param shock_index Integer. Which structural shock to target. Default 1.
#' @param parameterization Character. How to read `theta`: `"target"` (n - 1
#'   angles), `"full"` (a full rotation), or `"auto"` to decide from its
#'   length. With two variables both have length one and `"auto"` reads
#'   `theta` as `"target"`.
#' @return The objective value, a single non-negative number (lower is
#'   better).
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3)
#' var_est <- estimate_var(y, lags = 2, nimp_wold = 10)
#' # Event on raw row 80, i.e. effective-sample row 80 - lags
#' etvar_objective(c(0, 0), var_est, target_obs = 78)
#' etvar_objective(c(0.5, -0.3), var_est, target_obs = 78)
#' @export
etvar_objective <- function(theta, var_est, target_obs, shock_index = 1L,
                            parameterization = c("auto", "target", "full")) {
  parameterization <- match.arg(parameterization)
  context <- .etvar_prepare_objective(var_est, target_obs, shock_index)
  .etvar_objective_prepared(theta, context, parameterization)
}


# Precomputes everything that does not depend on the target direction q.
#
# With A0 = -L Q, the target shock adds L q q' z_t to the reduced-form
# residual at time t (z_t = L^{-1} u_t). Running these innovations through the
# VAR dynamics gives, at every date and for every variable, a linear function
# of the n^2 products q_a q_b. For each event date we keep the coefficients
# of the implied first difference, scaled like the data.
.etvar_prepare_objective <- function(var_est, target_obs, shock_index = 1L) {
  if (!inherits(var_est, "etvar_var")) {
    stop("`var_est` must be an `etvar_var` object from `estimate_var()`.",
         call. = FALSE)
  }
  n <- var_est$n
  lags <- var_est$lags
  shock_index <- .etvar_integer_scalar(shock_index, "shock_index",
                                       min = 1L, max = n)
  target_obs <- .etvar_index_vector(target_obs, "target_obs", nrow(var_est$Y))

  real_diff <- diff(var_est$data)
  scales <- apply(real_diff, 2L, stats::sd)
  if (any(!is.finite(scales)) || any(scales <= sqrt(.Machine$double.eps))) {
    stop("Every variable must have non-zero finite first-difference ",
         "variability.", call. = FALSE)
  }

  deterministic <- deterministic_path(var_est$ALPHA, var_est$X, n, lags,
                                      var_est$trend, var_est$constant)
  # Row i of diff() is the change from raw row i to raw row i + 1, so the
  # change into effective row t sits in row t + lags - 1.
  diff_rows <- target_obs + lags - 1L
  baseline_scaled <- sweep(
    diff(deterministic)[diff_rows, , drop = FALSE] -
      real_diff[diff_rows, , drop = FALSE],
    2L, scales, "/"
  )

  chol_lower <- t(chol(var_est$SIGMA))
  whitened <- var_est$eps

  # Column a + (b - 1) n of the innovation is L[, a] * z_t[b], matching
  # as.vector(q %o% q).
  innovation_base <- chol_lower[, rep(seq_len(n), times = n), drop = FALSE]
  transitions <- lapply(seq_len(lags), function(k) {
    t(var_est$ALPHA[var_est$constant + var_est$trend + (k - 1L) * n +
                      seq_len(n), , drop = FALSE])
  })

  # Only the last p coefficient matrices are needed, so keep them in a ring
  # buffer instead of storing the whole history.
  n_terms <- n * n
  buffer <- array(0, dim = c(lags + 1L, n, n_terms))
  slot <- function(raw_row) (raw_row - 1L) %% (lags + 1L) + 1L
  event_rows <- target_obs + lags
  event_coefficients <- array(0, dim = c(length(target_obs), n, n_terms))

  for (raw_row in seq.int(lags + 1L, max(event_rows))) {
    current <- sweep(innovation_base, 2L,
                     rep(whitened[raw_row - lags, ], each = n), "*")
    for (k in seq_len(lags)) {
      previous <- matrix(buffer[slot(raw_row - k), , ], nrow = n)
      current <- current + transitions[[k]] %*% previous
    }

    hits <- which(event_rows == raw_row)
    if (length(hits) > 0L) {
      previous <- matrix(buffer[slot(raw_row - 1L), , ], nrow = n)
      change <- sweep(current - previous, 1L, scales, "/")
      for (i in hits) event_coefficients[i, , ] <- change
    }
    buffer[slot(raw_row), , ] <- current
  }

  structure(
    list(
      n = n,
      shock_index = shock_index,
      target_obs = target_obs,
      baseline_scaled = baseline_scaled,
      event_coefficients_scaled = event_coefficients,
      scales = scales
    ),
    class = "etvar_objective_context"
  )
}


# Evaluates the objective from a prepared context. This is what optim() calls,
# so it only builds the direction and a few small matrix products.
.etvar_objective_prepared <- function(theta, context,
                                      parameterization = c("auto", "target",
                                                           "full")) {
  if (!inherits(context, "etvar_objective_context")) {
    stop("`context` must be prepared by `.etvar_prepare_objective()`.",
         call. = FALSE)
  }
  parameterization <- match.arg(parameterization)
  n <- context$n

  angles <- .etvar_prepare_target_start(
    theta, n, context$shock_index, parameterization
  )
  q <- .etvar_direction_from_angles(angles, n, context$shock_index)
  q_products <- as.vector(tcrossprod(q))

  deviations <- context$baseline_scaled
  for (i in seq_len(nrow(deviations))) {
    coefficients <- matrix(context$event_coefficients_scaled[i, , ], nrow = n)
    deviations[i, ] <- deviations[i, ] + drop(coefficients %*% q_products)
  }
  sum(abs(deviations))
}
