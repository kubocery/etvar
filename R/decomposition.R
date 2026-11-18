#' Historical decomposition of a structural VAR
#'
#' Splits every variable into a deterministic component (pre-sample values
#' carried forward by the constant and trend) and the contribution of each
#' structural shock. Because the VAR is linear, the contributions and the
#' deterministic component add up to the data.
#'
#' @param A0 Matrix (n x n). Structural impact matrix; the structural errors
#'   are `solve(A0, u_t)` for reduced-form residuals `u_t`.
#' @param ALPHA Matrix (k x n). VAR coefficients, e.g. from [estimate_var()].
#' @param Y Matrix (T_eff x n). Dependent variables.
#' @param X Matrix (T_eff x k). Regressors.
#' @param n Integer. Number of variables.
#' @param lags Integer. Number of lags.
#' @param trend Integer. 0 = none, 1 = linear, 2 = quadratic.
#' @param constant Integer. 1 = intercept, 0 = no intercept.
#' @return A list with
#'   \describe{
#'     \item{hist_matrix}{Array (T x n x n). Element `[t, v, s]` is the
#'       contribution of shock `s` to variable `v` at time `t`, where
#'       T = T_eff + lags covers the raw sample. The first `lags` rows are
#'       zero.}
#'     \item{deterministic}{Matrix (T x n). Deterministic component.}
#'     \item{errors_str}{Matrix (T_eff x n). Structural errors.}
#'   }
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3)
#' var_est <- estimate_var(y, lags = 2, nimp_wold = 10)
#' A0 <- t(chol(var_est$SIGMA))  # recursive (Cholesky) identification
#' hd <- hist_decomposition(A0, var_est$ALPHA, var_est$Y, var_est$X,
#'                          n = 3, lags = 2, trend = 0, constant = 1)
#'
#' # contributions plus the deterministic part reproduce the data
#' total <- hd$deterministic + rowSums(hd$hist_matrix, dims = 2)
#' all.equal(total, y, check.attributes = FALSE)
#' @export
hist_decomposition <- function(A0, ALPHA, Y, X, n, lags, trend, constant) {
  inputs <- .etvar_validate_decomposition(
    A0, ALPHA, Y, X, n, lags, trend, constant
  )
  hd <- with(inputs, decompose_history(A0, ALPHA, Y, X, n, lags, trend,
                                       constant, datestar = NULL))
  hd[c("hist_matrix", "deterministic", "errors_str")]
}


#' Historical decomposition with an intervention path
#'
#' Like [hist_decomposition()], but also decomposes an intervention path on
#' which all structural errors after a given date are set to zero. ETVAR uses
#' it to separate the target shocks up to and including the event date from
#' the target shocks that arrive afterwards.
#'
#' @inheritParams hist_decomposition
#' @param datestar Integer. Number of periods at the end of the sample whose
#'   structural errors are set to zero on the intervention path. In [etvar()]
#'   this is the number of observations after the first event date.
#' @return A list with the elements of [hist_decomposition()] and
#'   \describe{
#'     \item{hist_matrix_stop}{Array (T x n x n). Shock contributions on the
#'       intervention path.}
#'     \item{all_data}{Array (T x n x n). `all_data[, , s]` is the path the
#'       model implies when only shock `s` operates: the deterministic
#'       component plus the contribution of shock `s`.}
#'   }
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3)
#' var_est <- estimate_var(y, lags = 2, nimp_wold = 10)
#' A0 <- t(chol(var_est$SIGMA))  # recursive (Cholesky) identification
#' hd <- hist_decomposition_3d(A0, var_est$ALPHA, var_est$Y, var_est$X,
#'                             n = 3, lags = 2, trend = 0, constant = 1,
#'                             datestar = 20)
#' dim(hd$hist_matrix_stop)
#' @export
hist_decomposition_3d <- function(A0, ALPHA, Y, X, n, lags, trend, constant,
                                  datestar) {
  inputs <- .etvar_validate_decomposition(
    A0, ALPHA, Y, X, n, lags, trend, constant
  )
  datestar <- .etvar_integer_scalar(datestar, "datestar", min = 0L,
                                    max = nrow(inputs$Y))
  with(inputs, decompose_history(A0, ALPHA, Y, X, n, lags, trend, constant,
                                 datestar = datestar))
}


#' Rearrange shock contributions from a historical decomposition
#'
#' Pulls the pieces that the ETVAR plots need out of the (time x variable x
#' shock) arrays returned by [hist_decomposition_3d()]: the contribution of
#' the target shock, the sum over all other shocks and, when an intervention
#' path is supplied, the split of the target contribution into the part up to
#' the event and the part after it.
#'
#' @param hist_matrix Array (T x n x n). Shock contributions from
#'   [hist_decomposition()] or [hist_decomposition_3d()].
#' @param hist_matrix_stop Array (T x n x n) of intervention-path
#'   contributions from [hist_decomposition_3d()], or `NULL`.
#' @param T_obs Integer. Number of observations T.
#' @param n Integer. Number of variables (and shocks).
#' @param target_shock Integer. The shock to isolate. Default 1.
#' @return A list of T x n matrices:
#'   \describe{
#'     \item{total}{Sum of all shock contributions.}
#'     \item{target}{Contribution of the target shock.}
#'     \item{target_stop}{Target-shock contribution on the intervention path,
#'       or `NULL`.}
#'     \item{target_rest}{`target - target_stop`: the effect of target shocks
#'       after the intervention date, or `NULL`.}
#'     \item{others}{Contribution of all other shocks combined.}
#'   }
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3)
#' var_est <- estimate_var(y, lags = 2, nimp_wold = 10)
#' A0 <- t(chol(var_est$SIGMA))  # recursive (Cholesky) identification
#' hd <- hist_decomposition_3d(A0, var_est$ALPHA, var_est$Y, var_est$X,
#'                             n = 3, lags = 2, trend = 0, constant = 1,
#'                             datestar = 20)
#' contributions <- rearrange_shocks(hd$hist_matrix, hd$hist_matrix_stop,
#'                                   T_obs = nrow(y), n = 3)
#' head(contributions$target)
#' @export
rearrange_shocks <- function(hist_matrix, hist_matrix_stop = NULL,
                             T_obs, n, target_shock = 1L) {
  T_obs <- .etvar_integer_scalar(T_obs, "T_obs", min = 1L)
  n <- .etvar_integer_scalar(n, "n", min = 1L)
  target_shock <- .etvar_integer_scalar(target_shock, "target_shock",
                                        min = 1L, max = n)
  hist_matrix <- .etvar_numeric_array(hist_matrix, "hist_matrix",
                                      c(T_obs, n, n))
  if (!is.null(hist_matrix_stop)) {
    hist_matrix_stop <- .etvar_numeric_array(
      hist_matrix_stop, "hist_matrix_stop", c(T_obs, n, n)
    )
  }

  total <- rowSums(hist_matrix, dims = 2L)
  target <- matrix(hist_matrix[, , target_shock], nrow = T_obs, ncol = n)
  target_stop <- NULL
  target_rest <- NULL
  if (!is.null(hist_matrix_stop)) {
    target_stop <- matrix(hist_matrix_stop[, , target_shock],
                          nrow = T_obs, ncol = n)
    target_rest <- target - target_stop
  }

  list(
    total = total,
    target = target,
    target_stop = target_stop,
    target_rest = target_rest,
    others = total - target
  )
}


# Shared engine for the two decomposition functions (inputs already
# validated). Each shock's contribution follows the VAR recursion from a zero
# start, driven only by A0[, s] times that shock's structural errors.
decompose_history <- function(A0, ALPHA, Y, X, n, lags, trend, constant,
                              datestar = NULL) {
  n_obs <- nrow(Y)
  errors_str <- t(solve(A0, t(Y - X %*% ALPHA)))
  deterministic <- deterministic_path(ALPHA, X, n, lags, trend, constant)
  slopes <- lag_coefficients(ALPHA, n, lags, constant, trend)
  zeros <- matrix(0, nrow = lags, ncol = n)

  contributions <- function(errors) {
    out <- array(0, dim = c(n_obs + lags, n, n))
    for (s in seq_len(n)) {
      out[-seq_len(lags), , s] <- var_recursion(
        slopes, zeros, outer(errors[, s], A0[, s])
      )
    }
    out
  }

  hist_matrix <- contributions(errors_str)
  out <- list(
    hist_matrix = hist_matrix,
    deterministic = deterministic,
    errors_str = errors_str
  )
  if (is.null(datestar)) {
    return(out)
  }

  errors_stop <- errors_str
  if (datestar > 0L) {
    errors_stop[seq.int(n_obs - datestar + 1L, n_obs), ] <- 0
  }
  out$hist_matrix_stop <- contributions(errors_stop)
  out$all_data <- sweep(hist_matrix, c(1L, 2L), deterministic, "+")
  out
}

# Contribution of a single structural error -- shock `shock` at effective
# row `row` -- to all variables over the raw sample: zero before that date,
# then the impulse response scaled by the size of the error.
single_error_contribution <- function(var_est, A0, errors_str, shock, row) {
  n <- var_est$n
  lags <- var_est$lags
  offset <- matrix(0, nrow = nrow(errors_str), ncol = n)
  offset[row, ] <- A0[, shock] * errors_str[row, shock]
  slopes <- lag_coefficients(var_est$ALPHA, n, lags, var_est$constant,
                             var_est$trend)
  path <- var_recursion(slopes, matrix(0, nrow = lags, ncol = n), offset)
  rbind(matrix(0, nrow = lags, ncol = n), path)
}
