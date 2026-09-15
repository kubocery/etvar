# ==============================================================================
# ETVAR objective function
# ==============================================================================

#' ETVAR objective function
#'
#' Given a vector of rotation angles π, this function:
#' 1. Constructs the structural impact matrix A0 = chol(Σ)' Q(π).
#' 2. Runs a historical decomposition.
#' 3. Isolates the counterfactual path driven by the target structural shock.
#' 4. Computes the scaled absolute deviation between the counterfactual and
#'    realised growth rates on the event date.
#'
#' The optimiser minimises this function to find the rotation that makes the
#' target shock best explain the event-date movements.
#'
#' @param theta       Numeric vector. Rotation angles (length n*(n-1)/2).
#' @param var_est     An `etvar_var` object from [estimate_var()].
#' @param target_obs  Integer or integer vector. Row index (indices) of the
#'                    event date(s) in the effective sample (i.e., in
#'                    `var_est$Y`). Passing multiple indices defines a target
#'                    window: the objective is the sum across all event dates.
#' @param shock_index Integer. Which structural shock to target. Default 1.
#' @return Scalar. Sum of absolute scaled deviations across all target dates
#'   (lower is better).
#' @export
etvar_objective <- function(theta, var_est, target_obs, shock_index = 1L) {
  if (!inherits(var_est, "etvar_var")) {
    stop("`var_est` must be an `etvar_var` object from `estimate_var()`.",
         call. = FALSE)
  }

  SIGMA <- var_est$SIGMA
  data  <- var_est$data
  n     <- var_est$n
  shock_index <- .etvar_integer_scalar(shock_index, "shock_index",
                                        min = 1L, max = n)
  target_obs <- .etvar_index_vector(target_obs, "target_obs", nrow(var_est$Y))

  R_mat <- euler_rotation(theta, n)
  A0    <- -t(chol(SIGMA)) %*% R_mat

  decomp <- hist_decomposition(
    A0       = A0,
    ALPHA    = var_est$ALPHA,
    Y        = var_est$Y,
    X        = var_est$X,
    n        = n,
    lags     = var_est$lags,
    trend    = var_est$trend,
    constant = var_est$constant
  )

  # Permute (Time x Variable x Shock) -> (Time x Shock x Variable)
  hist_perm   <- aperm(decomp$hist_matrix, c(1, 3, 2))
  shock_contr <- hist_perm[, shock_index, ]  # target shock contribution

  CF_data <- decomp$deterministic + shock_contr

  # Growth-rate deviations at the event date

  CF_diff   <- diff(CF_data)
  real_diff <- diff(data)

  sc <- apply(real_diff, 2, stats::sd)
  if (any(!is.finite(sc)) || any(sc <= sqrt(.Machine$double.eps))) {
    stop("Every variable must have non-zero finite first-difference ",
         "variability.", call. = FALSE)
  }

  # `target_obs` indexes the effective sample, whereas `diff(data)` indexes
  # changes in the raw sample. Add the lag offset, then account for the row
  # removed by `diff()`.
  idx <- target_obs + var_est$lags - 1L

  lib_CF   <- CF_diff[idx, , drop = FALSE]
  lib_real <- real_diff[idx, , drop = FALSE]

  sum(abs(sweep(lib_CF, 2, sc, "/") - sweep(lib_real, 2, sc, "/")))
}
