# ==============================================================================
# VAR estimation internals and main estimation function
# ==============================================================================

#' Trim data to the effective sample
#'
#' Removes the first `lags` rows from the data matrix so that it aligns with
#' the regressor matrix built by [build_var_regressors()].
#'
#' @param data Matrix (T x n). Raw data.
#' @param lags Integer. Number of lags to drop.
#' @return Matrix ((T - lags) x n).
#' @keywords internal
trim_to_effective <- function(data, lags) {
  data[(lags + 1):nrow(data), , drop = FALSE]
}

#' Build the VAR regressor matrix
#'
#' Constructs a matrix `X = [deterministic | Y_{t-1} | ... | Y_{t-p}]` where
#' the deterministic block may include a constant, linear trend, or quadratic
#' trend depending on the flags.
#'
#' @param data    Matrix (T x n). Raw data.
#' @param Tr      Integer. Total rows in raw data.
#' @param n       Integer. Number of variables.
#' @param lags    Integer. Number of lags.
#' @param trend   Integer. 0 = none, 1 = linear, 2 = quadratic.
#' @param constant Integer. 1 = include intercept, 0 = omit.
#' @return Matrix ((Tr - lags) x k). Regressor matrix.
#' @keywords internal
build_var_regressors <- function(data, Tr, n, lags, trend, constant) {
  Ylag <- matrix(0, nrow = Tr, ncol = n * lags)

  for (i in 1:lags) {
    col_start <- n * (i - 1) + 1
    col_end   <- n * i
    src_start <- lags + 1 - i
    src_end   <- Tr - i
    Ylag[(lags + 1):Tr, col_start:col_end] <- data[src_start:src_end, 1:n]
  }

  Ylag_trimmed <- Ylag[(lags + 1):Tr, , drop = FALSE]
  T_eff        <- Tr - lags
  ones_vec     <- matrix(1, nrow = T_eff, ncol = 1)
  tt           <- matrix(1:T_eff, ncol = 1)
  tt2          <- tt^2

  if (constant == 1) {
    if (trend == 1)      Xlag <- cbind(ones_vec, tt, Ylag_trimmed)
    else if (trend == 2) Xlag <- cbind(ones_vec, tt, tt2, Ylag_trimmed)
    else                 Xlag <- cbind(ones_vec, Ylag_trimmed)
  } else {
    if (trend == 1)      Xlag <- cbind(tt, Ylag_trimmed)
    else if (trend == 2) Xlag <- cbind(tt, tt2, Ylag_trimmed)
    else                 Xlag <- Ylag_trimmed
  }

  return(Xlag)
}

#' Build the VAR companion matrix
#'
#' Stacks a VAR(p) as a VAR(1) system in companion form. Strips deterministic
#' coefficient rows (constant + trend) from `ALPHA` and builds:
#' \deqn{F = \begin{pmatrix} A_1 & A_2 & \cdots & A_p \\
#' I & 0 & \cdots & 0 \\ & \ddots & & \vdots \\ 0 & \cdots & I & 0
#' \end{pmatrix}}
#'
#' @param ALPHA    Matrix (k x n). OLS coefficient matrix.
#' @param n        Integer. Number of endogenous variables.
#' @param lags     Integer. Number of lags.
#' @param constant Integer. 1 = intercept present, 0 = none.
#' @param trend    Integer. Trend order (0, 1, or 2).
#' @return Matrix (n*lags x n*lags). Companion matrix.
#' @keywords internal
build_companion <- function(ALPHA, n, lags, constant, trend) {
  n_det     <- constant + trend
  lag_coefs <- t(ALPHA[(n_det + 1):nrow(ALPHA), , drop = FALSE])
  bottom    <- cbind(diag(n * (lags - 1)),
                     matrix(0, nrow = n * (lags - 1), ncol = n))
  rbind(lag_coefs, bottom)
}

#' Compute Wold impulse responses
#'
#' Computes the moving-average (Wold) representation by powering up the
#' companion matrix, and the Cholesky-orthogonalized IRFs.
#'
#' @param BB    Matrix (nL x nL). Companion matrix.
#' @param SIGMA Matrix (n x n). Reduced-form error covariance.
#' @param n     Integer. Number of endogenous variables.
#' @param nimp  Integer. Number of horizons.
#' @return Named list with `wimpu` (companion powers) and `cimpu` (Cholesky IRFs).
#' @keywords internal
compute_wold <- function(BB, SIGMA, n, nimp) {
  rows_bb    <- nrow(BB)
  cols_bb    <- ncol(BB)
  wimpu      <- array(0, dim = c(rows_bb, cols_bb, nimp))
  cimpu      <- array(0, dim = c(n, n, nimp))
  chol_lower <- t(chol(SIGMA))

  for (j in 1:nimp) {
    wimpu[, , j] <- if (j == 1) diag(rows_bb) else wimpu[, , j - 1] %*% BB
    cimpu[, , j] <- wimpu[1:n, 1:n, j] %*% chol_lower
  }

  list(cimpu = cimpu, wimpu = wimpu)
}


#' Estimate a reduced-form VAR model
#'
#' Estimates a VAR(p) by OLS, computes the error covariance matrix (with
#' degrees-of-freedom correction), Cholesky-orthogonalized shocks, the
#' companion matrix, and the Wold impulse response representation.
#'
#' @param data     Matrix (T x n). Observed time-series data.
#' @param lags     Integer. Number of lags p.
#' @param constant Logical or integer. Include an intercept? Default `TRUE`.
#' @param trend    Integer. Trend specification: 0 = none, 1 = linear,
#'                 2 = quadratic. Default 0.
#' @param nimp_wold Integer. Number of horizons for the Wold decomposition.
#'                  Default 100.
#' @param estimation_rows Integer row indices or a logical vector selecting
#'   which rows of the effective sample (`Y` and `X`) are used to estimate the
#'   VAR coefficients and innovation covariance. Residuals are still evaluated
#'   for the full effective sample. `NULL` uses every effective-sample row.
#' @return A named list (class `"etvar_var"`) with elements:
#'   \describe{
#'     \item{Y}{Dependent variable matrix (T_eff x n).}
#'     \item{X}{Regressor matrix (T_eff x k).}
#'     \item{ALPHA}{OLS coefficients (k x n).}
#'     \item{SIGMA}{Error covariance matrix (n x n).}
#'     \item{eps}{Cholesky-orthogonalized residuals (T_eff x n).}
#'     \item{residuals}{OLS residuals evaluated over the full effective sample
#'       (T_eff x n).}
#'     \item{residuals_estimation}{OLS residuals for the estimation rows.}
#'     \item{estimation_rows}{Indices of the effective-sample rows used for
#'       estimation.}
#'     \item{companion}{Companion matrix (n*lags x n*lags).}
#'     \item{wimpu}{Wold multiplier array (nL x nL x nimp_wold).}
#'     \item{cimpu}{Cholesky IRFs (n x n x nimp_wold).}
#'     \item{n}{Number of variables.}
#'     \item{lags}{Number of lags.}
#'     \item{constant}{Constant flag (integer).}
#'     \item{trend}{Trend order (integer).}
#'     \item{data}{Original data matrix.}
#'   }
#' @export
#' @examples
#' # Estimate a VAR(4) on simulated data
#' set.seed(1)
#' Y <- matrix(rnorm(500), ncol = 2)
#' var_est <- estimate_var(Y, lags = 4)
#' var_est$n
estimate_var <- function(data, lags, constant = TRUE, trend = 0,
                         nimp_wold = 100, estimation_rows = NULL) {
  const_int <- as.integer(constant)
  data      <- as.matrix(data)
  Tr        <- nrow(data)
  n         <- ncol(data)

  Y <- trim_to_effective(data, lags)
  X <- build_var_regressors(data, Tr, n, lags, trend, const_int)

  Tr_eff <- nrow(Y)

  if (is.null(estimation_rows)) {
    estimation_rows <- seq_len(Tr_eff)
  } else if (is.logical(estimation_rows)) {
    if (length(estimation_rows) != Tr_eff || anyNA(estimation_rows)) {
      stop("Logical `estimation_rows` must have one non-missing value ",
           "per effective-sample row.")
    }
    estimation_rows <- which(estimation_rows)
  } else {
    if (!is.numeric(estimation_rows) || anyNA(estimation_rows) ||
        any(!is.finite(estimation_rows)) ||
        any(estimation_rows != as.integer(estimation_rows))) {
      stop("`estimation_rows` must contain finite, non-missing row indices.")
    }
    estimation_rows <- sort(unique(as.integer(estimation_rows)))
    if (any(estimation_rows < 1L | estimation_rows > Tr_eff)) {
      stop("`estimation_rows` contains indices outside the effective sample.")
    }
  }

  X_est <- X[estimation_rows, , drop = FALSE]
  Y_est <- Y[estimation_rows, , drop = FALSE]
  n_parameters <- ncol(X_est)
  dof <- nrow(Y_est) - n_parameters

  if (dof <= 0L) {
    stop("The VAR estimation window must contain more usable observations ",
         "than coefficients per equation (", n_parameters, ").")
  }
  if (qr(X_est)$rank < n_parameters) {
    stop("The VAR regressor matrix is rank deficient in the estimation window.")
  }

  ALPHA    <- qr.solve(X_est, Y_est)
  resid_est <- Y_est - X_est %*% ALPHA
  resid     <- Y - X %*% ALPHA
  sse       <- t(resid_est) %*% resid_est
  SIGMA <- sse / as.numeric(dof)

  eps <- resid %*% solve(chol(SIGMA))

  BB   <- build_companion(ALPHA, n, lags, const_int, trend)
  wold <- compute_wold(BB, SIGMA, n, nimp_wold)

  structure(list(
    Y         = Y,
    X         = X,
    ALPHA     = ALPHA,
    SIGMA     = SIGMA,
    eps       = eps,
    residuals = resid,
    residuals_estimation = resid_est,
    estimation_rows = estimation_rows,
    nobs_estimation = length(estimation_rows),
    companion = BB,
    wimpu     = wold$wimpu,
    cimpu     = wold$cimpu,
    n         = n,
    lags      = lags,
    constant  = const_int,
    trend     = trend,
    data      = data
  ), class = "etvar_var")
}
