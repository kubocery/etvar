# ==============================================================================
# Historical decomposition
# ==============================================================================

# ---------- internal helper: build deterministic regressor vector ----------
.build_reg_vec <- function(constant, trend, tt, tt2, auxi) {
  if (constant == 1L) {
    if (trend == 1L)      return(c(1, tt, auxi))
    if (trend == 2L)      return(c(1, tt, tt2, auxi))
    c(1, auxi)
  } else {
    if (trend == 1L)      return(c(tt, auxi))
    if (trend == 2L)      return(c(tt, tt2, auxi))
    auxi
  }
}


#' Historical decomposition of a structural VAR
#'
#' Decomposes each variable into the contribution of each structural shock and
#' a deterministic (trend + constant) component. For every shock, a
#' counterfactual path is built recursively using only that shock's
#' reduced-form innovations; the deterministic component is then subtracted.
#'
#' @param A0       Matrix (n x n). Structural impact matrix.
#' @param ALPHA    Matrix (k x n). VAR coefficient matrix.
#' @param Y        Matrix (T_eff x n). Dependent variable matrix.
#' @param X        Matrix (T_eff x k). Regressor matrix.
#' @param n        Integer. Number of variables.
#' @param lags     Integer. Number of lags.
#' @param trend    Integer. 0 = none, 1 = linear, 2 = quadratic.
#' @param constant Integer. 1 = include intercept, 0 = omit.
#' @return Named list:
#'   \describe{
#'     \item{hist_matrix}{Array (Traw x n x n). Detrended shock contributions
#'       where `[t, v, s]` is the contribution of shock `s` to variable `v`
#'       at time `t`.}
#'     \item{deterministic}{Matrix (Traw x n). Deterministic component path.}
#'     \item{errors_str}{Matrix (T_eff x n). Structural errors.}
#'   }
#' @export
hist_decomposition <- function(A0, ALPHA, Y, X, n, lags, trend, constant) {
  inputs <- .etvar_validate_decomposition(
    A0, ALPHA, Y, X, n, lags, trend, constant
  )
  A0 <- inputs$A0
  ALPHA <- inputs$ALPHA
  Y <- inputs$Y
  X <- inputs$X
  n <- inputs$n
  lags <- inputs$lags
  trend <- inputs$trend
  constant <- inputs$constant
  Traw <- nrow(Y) + lags

  errors_red <- Y - X %*% ALPHA
  errors_str <- t(solve(A0, t(errors_red)))

  Y1_shock_detrend <- array(0, dim = c(Traw, n, n))
  errors_str10     <- array(0, dim = c(Traw - lags, n, n))
  Y1_shock         <- array(0, dim = c(Traw, n, n))
  Y1_shock_det     <- matrix(0, nrow = Traw, ncol = n)

  nc <- as.integer(constant == 1L)
  nt <- trend

  # --- Initialise from pre-sample lags ---
  for (i_lag in seq_len(lags)) {
    idx_end   <- (lags - (i_lag - 1)) * n + (nc + nt)
    idx_start <- idx_end - n + 1
    val       <- X[1, idx_start:idx_end]
    Y1_shock_det[i_lag, ] <- val
    for (i in seq_len(n)) Y1_shock[i_lag, , i] <- val
  }

  auxi_det <- numeric(lags * n)
  for (i in seq_len(lags)) {
    idx_s <- 1 + (i - 1) * n
    idx_e <- n + (i - 1) * n
    auxi_det[idx_s:idx_e] <- Y1_shock_det[(lags + 1) - i, ]
  }

  # --- Map structural errors to reduced-form ---
  for (i in seq_len(n)) {
    lala      <- matrix(0, nrow = Traw - lags, ncol = n)
    lala[, i] <- errors_str[, i]
    errors_str10[, , i] <- t(A0 %*% t(lala))
  }

  # --- Main recursive loop ---
  tt  <- 1
  tt2 <- 1

  for (i_time in seq.int(lags + 1L, Traw)) {
    for (iauxi in seq_len(lags)) {
      idx_s <- 1 + (iauxi - 1) * n
      idx_e <- n + (iauxi - 1) * n
      auxi_det[idx_s:idx_e] <- Y1_shock_det[i_time - iauxi, ]
    }

    reg_vec <- .build_reg_vec(constant, trend, tt, tt2, auxi_det)
    Y1_shock_det[i_time, ] <- matrix(reg_vec, nrow = 1) %*% ALPHA

    for (i_shock in seq_len(n)) {
      auxi_shock <- numeric(lags * n)
      for (iauxi2 in seq_len(lags)) {
        idx_s <- 1 + (iauxi2 - 1) * n
        idx_e <- n + (iauxi2 - 1) * n
        auxi_shock[idx_s:idx_e] <- Y1_shock[i_time - iauxi2, , i_shock]
      }

      reg_shock <- .build_reg_vec(constant, trend, tt, tt2, auxi_shock)

      term1 <- matrix(reg_shock, nrow = 1) %*% ALPHA
      term2 <- errors_str10[i_time - lags, , i_shock]

      Y1_shock[i_time, , i_shock]         <- term1 + term2
      Y1_shock_detrend[i_time, , i_shock] <-
        Y1_shock[i_time, , i_shock] - Y1_shock_det[i_time, ]
    }

    tt  <- tt + 1
    tt2 <- tt^2
  }

  list(
    hist_matrix   = Y1_shock_detrend,
    deterministic = Y1_shock_det,
    errors_str    = errors_str
  )
}


#' Historical decomposition with intervention (stop) path
#'
#' Extends [hist_decomposition()] by simultaneously computing a counterfactual
#' path where structural errors are zeroed out in the last `datestar` periods.
#' This separates the contribution of shocks occurring on/after the event from
#' those that occurred before.
#'
#' @inheritParams hist_decomposition
#' @param datestar Integer. Number of periods from the end of the sample in
#'   which structural errors are set to zero (intervention window).
#' @return Named list:
#'   \describe{
#'     \item{hist_matrix}{Array (Traw x n x n). Detrended shock contributions.}
#'     \item{deterministic}{Matrix (Traw x n). Deterministic component.}
#'     \item{errors_str}{Matrix (T_eff x n). Structural errors.}
#'     \item{hist_matrix_stop}{Array (Traw x n x n). Detrended contributions
#'       on the intervention path.}
#'     \item{all_data}{Array (Traw x n x n). Full (non-detrended) shock paths.}
#'   }
#' @export
hist_decomposition_3d <- function(A0, ALPHA, Y, X, n, lags, trend, constant,
                                  datestar) {
  inputs <- .etvar_validate_decomposition(
    A0, ALPHA, Y, X, n, lags, trend, constant
  )
  A0 <- inputs$A0
  ALPHA <- inputs$ALPHA
  Y <- inputs$Y
  X <- inputs$X
  n <- inputs$n
  lags <- inputs$lags
  trend <- inputs$trend
  constant <- inputs$constant
  datestar <- .etvar_integer_scalar(datestar, "datestar", min = 0L,
                                    max = nrow(Y))
  Traw <- nrow(Y) + lags

  errors_red <- Y - X %*% ALPHA
  errors_str <- t(solve(A0, t(errors_red)))

  errors_str_stop <- errors_str
  if (datestar > 0) {
    idx_s <- nrow(errors_str_stop) - datestar + 1
    idx_e <- nrow(errors_str_stop)
    errors_str_stop[idx_s:idx_e, ] <- 0
  }

  Y1_shock_detrend      <- array(0, dim = c(Traw, n, n))
  Y1_shock_detrend_stop <- array(0, dim = c(Traw, n, n))
  errors_str10          <- array(0, dim = c(Traw - lags, n, n))
  errors_str10_stop     <- array(0, dim = c(Traw - lags, n, n))
  Y1_shock              <- array(0, dim = c(Traw, n, n))
  Y1_shock_stop         <- array(0, dim = c(Traw, n, n))
  Y1_shock_det          <- matrix(0, nrow = Traw, ncol = n)

  nc <- as.integer(constant == 1L)
  nt <- trend

  for (i_lag in seq_len(lags)) {
    idx_end   <- (lags - (i_lag - 1)) * n + (nc + nt)
    idx_start <- idx_end - n + 1
    val <- X[1, idx_start:idx_end]
    Y1_shock_det[i_lag, ] <- val
    for (i in seq_len(n)) {
      Y1_shock[i_lag, , i]      <- val
      Y1_shock_stop[i_lag, , i] <- val
    }
  }

  auxi_det <- numeric(lags * n)
  for (i in seq_len(lags)) {
    idx_s <- 1 + (i - 1) * n
    idx_e <- n + (i - 1) * n
    auxi_det[idx_s:idx_e] <- Y1_shock_det[(lags + 1) - i, ]
  }

  for (i in seq_len(n)) {
    lala      <- matrix(0, nrow = Traw - lags, ncol = n)
    lala[, i] <- errors_str[, i]
    errors_str10[, , i] <- t(A0 %*% t(lala))
    lala[, i] <- errors_str_stop[, i]
    errors_str10_stop[, , i] <- t(A0 %*% t(lala))
  }

  auxi_shock      <- numeric(lags * n)
  auxi_shock_stop <- numeric(lags * n)
  tt  <- 1
  tt2 <- 1

  for (i_time in seq.int(lags + 1L, Traw)) {
    for (iauxi in seq_len(lags)) {
      idx_s <- 1 + (iauxi - 1) * n
      idx_e <- n + (iauxi - 1) * n
      auxi_det[idx_s:idx_e] <- Y1_shock_det[i_time - iauxi, ]
    }

    reg_det <- .build_reg_vec(constant, trend, tt, tt2, auxi_det)
    Y1_shock_det[i_time, ] <- reg_det %*% ALPHA

    for (i_shock in seq_len(n)) {
      for (iauxi2 in seq_len(lags)) {
        idx_s <- 1 + (iauxi2 - 1) * n
        idx_e <- n + (iauxi2 - 1) * n
        auxi_shock[idx_s:idx_e]      <- Y1_shock[i_time - iauxi2, , i_shock]
        auxi_shock_stop[idx_s:idx_e] <-
          Y1_shock_stop[i_time - iauxi2, , i_shock]
      }

      err_idx <- i_time - lags

      reg_shock <- .build_reg_vec(constant, trend, tt, tt2, auxi_shock)
      Y1_shock[i_time, , i_shock] <-
        (reg_shock %*% ALPHA) + errors_str10[err_idx, , i_shock]
      Y1_shock_detrend[i_time, , i_shock] <-
        Y1_shock[i_time, , i_shock] - Y1_shock_det[i_time, ]

      reg_stop <- .build_reg_vec(constant, trend, tt, tt2, auxi_shock_stop)
      Y1_shock_stop[i_time, , i_shock] <-
        (reg_stop %*% ALPHA) + errors_str10_stop[err_idx, , i_shock]
      Y1_shock_detrend_stop[i_time, , i_shock] <-
        Y1_shock_stop[i_time, , i_shock] - Y1_shock_det[i_time, ]
    }

    tt  <- tt + 1
    tt2 <- tt^2
  }

  list(
    hist_matrix      = Y1_shock_detrend,
    deterministic    = Y1_shock_det,
    errors_str       = errors_str,
    hist_matrix_stop = Y1_shock_detrend_stop,
    all_data         = Y1_shock
  )
}


#' Rearrange shock contributions from a historical decomposition
#'
#' Permutes the 3-D decomposition array from (Time × Variable × Shock) to
#' (Time × Shock × Variable) and extracts convenient per-shock matrices.
#' Generalises beyond a single "tariff shock": specify `target_shock` to
#' choose which shock to isolate.
#'
#' @param hist_matrix      Array (T x n x n). Detrended shock contributions
#'   from [hist_decomposition()] or [hist_decomposition_3d()].
#' @param hist_matrix_stop Array (T x n x n). Intervention-path contributions
#'   (from [hist_decomposition_3d()]). Pass `NULL` if no intervention path
#'   was computed.
#' @param T_obs            Integer. Number of observations.
#' @param n                Integer. Number of variables / shocks.
#' @param target_shock     Integer. Which shock to isolate. Default 1.
#' @return Named list:
#'   \describe{
#'     \item{total}{Matrix (T x n). Sum of all shock contributions per
#'       variable.}
#'     \item{target}{Matrix (T x n). Contribution of the target shock.}
#'     \item{target_stop}{Matrix (T x n). Target shock on the intervention path
#'       (or `NULL`).}
#'     \item{target_rest}{Matrix (T x n). Continued target shocks
#'       (target - target_stop), or `NULL`.}
#'     \item{others}{Matrix (T x n). All shocks except the target.}
#'   }
#' @export
rearrange_shocks <- function(hist_matrix, hist_matrix_stop = NULL,
                             T_obs, n, target_shock = 1L) {
  T_obs <- .etvar_integer_scalar(T_obs, "T_obs", min = 1L)
  n <- .etvar_integer_scalar(n, "n", min = 1L)
  target_shock <- .etvar_integer_scalar(target_shock, "target_shock",
                                         min = 1L, max = n)
  hist_matrix <- .etvar_numeric_array(
    hist_matrix, "hist_matrix", c(T_obs, n, n)
  )
  if (!is.null(hist_matrix_stop)) {
    hist_matrix_stop <- .etvar_numeric_array(
      hist_matrix_stop, "hist_matrix_stop", c(T_obs, n, n)
    )
  }

  hist_all      <- array(0, dim = c(T_obs, n, n))
  has_stop      <- !is.null(hist_matrix_stop)
  if (has_stop) hist_all_stop <- array(0, dim = c(T_obs, n, n))

  for (i in seq_len(n)) {
    hist_all[, , i] <- hist_matrix[, i, ]
    if (has_stop) hist_all_stop[, , i] <- hist_matrix_stop[, i, ]
  }

  total       <- matrix(0, nrow = T_obs, ncol = n)
  target_mat  <- matrix(0, nrow = T_obs, ncol = n)
  target_stop <- if (has_stop) matrix(0, nrow = T_obs, ncol = n) else NULL

  for (iv in seq_len(n)) {
    target_mat[, iv] <- hist_all[, target_shock, iv]
    total[, iv]      <- rowSums(hist_all[, , iv])
    if (has_stop) target_stop[, iv] <- hist_all_stop[, target_shock, iv]
  }

  others      <- total - target_mat
  target_rest <- if (has_stop) target_mat - target_stop else NULL

  list(
    total       = total,
    target      = target_mat,
    target_stop = target_stop,
    target_rest = target_rest,
    others      = others
  )
}
