#' Estimate a reduced-form VAR
#'
#' Fits a VAR(p) by OLS, equation by equation. Besides the coefficients and the
#' residual covariance matrix (with a degrees-of-freedom correction), the
#' result holds the companion matrix and the reduced-form moving-average
#' (Wold) multipliers that the rest of the package builds on.
#'
#' @param data Matrix or data frame (T x n) of observed series.
#' @param lags Integer. Number of lags p.
#' @param constant Logical. Include an intercept? Default `TRUE`.
#' @param trend Integer. Deterministic trend: 0 = none, 1 = linear,
#'   2 = quadratic. Default 0.
#' @param nimp_wold Integer. Number of horizons (including horizon zero) for
#'   the Wold multipliers. Default 100.
#' @param estimation_rows Integer row indices or a logical vector selecting
#'   the rows of the effective sample (`Y` and `X`) used to estimate the
#'   coefficients and the covariance matrix. Residuals are still computed for
#'   the whole effective sample. `NULL` (the default) uses every row.
#' @return A list of class `"etvar_var"` with elements:
#'   \describe{
#'     \item{Y}{Dependent variables over the effective sample (T_eff x n),
#'       i.e. `data` without its first `lags` rows.}
#'     \item{X}{Regressors (T_eff x k): deterministic terms followed by the
#'       lags of all variables.}
#'     \item{ALPHA}{OLS coefficients (k x n).}
#'     \item{SIGMA}{Residual covariance matrix (n x n).}
#'     \item{eps}{Cholesky-orthogonalised residuals (T_eff x n).}
#'     \item{residuals}{Residuals over the whole effective sample (T_eff x n).}
#'     \item{residuals_estimation}{Residuals of the estimation rows.}
#'     \item{estimation_rows}{Effective-sample rows used for estimation.}
#'     \item{nobs_estimation}{Number of estimation rows.}
#'     \item{companion}{Companion matrix (n p x n p).}
#'     \item{wimpu}{Powers of the companion matrix (n p x n p x
#'       `nimp_wold`); the top-left n x n block of slice h + 1 is the
#'       horizon-h Wold multiplier.}
#'     \item{cimpu}{Cholesky impulse responses (n x n x `nimp_wold`).}
#'     \item{n, lags, constant, trend}{The model specification.}
#'     \item{data}{The data as a numeric matrix.}
#'   }
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(500), ncol = 2)
#' var_est <- estimate_var(y, lags = 4)
#' var_est$ALPHA
#' @export
estimate_var <- function(data, lags, constant = TRUE, trend = 0,
                         nimp_wold = 100, estimation_rows = NULL) {
  data <- .etvar_numeric_matrix(data)
  lags <- .etvar_integer_scalar(lags, "lags", min = 1L)
  constant <- as.integer(.etvar_flag(constant, "constant"))
  trend <- .etvar_integer_scalar(trend, "trend", min = 0L, max = 2L)
  nimp_wold <- .etvar_integer_scalar(nimp_wold, "nimp_wold", min = 1L)
  n <- ncol(data)

  if (nrow(data) <= lags) {
    stop("`data` must contain more rows than `lags`.", call. = FALSE)
  }

  Y <- trim_to_effective(data, lags)
  X <- build_var_regressors(data, lags, trend, constant)
  estimation_rows <- check_estimation_rows(estimation_rows, nrow(Y))

  X_est <- X[estimation_rows, , drop = FALSE]
  Y_est <- Y[estimation_rows, , drop = FALSE]
  dof <- nrow(X_est) - ncol(X_est)
  if (dof <= 0L) {
    stop("The VAR estimation window must contain more usable observations ",
         "than coefficients per equation (", ncol(X_est), ").", call. = FALSE)
  }
  qr_x <- qr(X_est)
  if (qr_x$rank < ncol(X_est)) {
    stop("The VAR regressor matrix is rank deficient in the estimation window.",
         call. = FALSE)
  }

  ALPHA <- qr.coef(qr_x, Y_est)
  resid_est <- Y_est - X_est %*% ALPHA
  resid <- Y - X %*% ALPHA
  SIGMA <- crossprod(resid_est) / dof

  chol_sigma <- tryCatch(chol(SIGMA), error = function(e) {
    stop("The estimated innovation covariance matrix is not positive ",
         "definite. Use a longer estimation window or remove collinear ",
         "variables.", call. = FALSE)
  })

  companion <- build_companion(ALPHA, n, lags, constant, trend)
  wold <- compute_wold(companion, SIGMA, n, nimp_wold)

  structure(list(
    Y = Y,
    X = X,
    ALPHA = ALPHA,
    SIGMA = SIGMA,
    eps = resid %*% solve(chol_sigma),
    residuals = resid,
    residuals_estimation = resid_est,
    estimation_rows = estimation_rows,
    nobs_estimation = length(estimation_rows),
    companion = companion,
    wimpu = wold$wimpu,
    cimpu = wold$cimpu,
    n = n,
    lags = lags,
    constant = constant,
    trend = trend,
    data = data
  ), class = "etvar_var")
}


# Internal helpers -------------------------------------------------------------

check_estimation_rows <- function(rows, n_obs) {
  if (is.null(rows)) {
    return(seq_len(n_obs))
  }
  if (is.logical(rows)) {
    if (length(rows) != n_obs || anyNA(rows)) {
      stop("Logical `estimation_rows` must have one non-missing value ",
           "per effective-sample row.", call. = FALSE)
    }
    return(which(rows))
  }
  if (!is.numeric(rows) || anyNA(rows) || any(!is.finite(rows)) ||
      any(rows != round(rows))) {
    stop("`estimation_rows` must contain finite, non-missing row indices.",
         call. = FALSE)
  }
  rows <- sort(unique(as.integer(rows)))
  if (any(rows < 1L | rows > n_obs)) {
    stop("`estimation_rows` contains indices outside the effective sample.",
         call. = FALSE)
  }
  rows
}

# The rows of `data` that have a full set of lags.
trim_to_effective <- function(data, lags) {
  data[seq.int(lags + 1L, nrow(data)), , drop = FALSE]
}

# Constant, linear trend and squared trend, in that order. The trend counts
# effective-sample periods, so it is 1 in the first period with a full set of
# lags. Returns NULL when the model has no deterministic terms.
deterministic_regressors <- function(n_obs, constant, trend) {
  t <- seq_len(n_obs)
  cbind(
    if (constant == 1L) rep(1, n_obs),
    if (trend >= 1L) t,
    if (trend >= 2L) t^2
  )
}

# X = [deterministic terms | y_{t-1} | ... | y_{t-p}] for the effective sample.
build_var_regressors <- function(data, lags, trend, constant) {
  n_raw <- nrow(data)
  lagged <- lapply(seq_len(lags), function(k) {
    data[seq.int(lags + 1L - k, n_raw - k), , drop = FALSE]
  })
  unname(cbind(
    deterministic_regressors(n_raw - lags, constant, trend),
    do.call(cbind, lagged)
  ))
}

# The rows of ALPHA that multiply the lagged variables, (n p) x n.
lag_coefficients <- function(ALPHA, n, lags, constant, trend) {
  ALPHA[constant + trend + seq_len(n * lags), , drop = FALSE]
}

# Companion form of the VAR: the lag matrices A_1, ..., A_p in the first n
# rows and an identity block that shifts the older lags down.
build_companion <- function(ALPHA, n, lags, constant, trend) {
  top <- t(lag_coefficients(ALPHA, n, lags, constant, trend))
  shift <- cbind(diag(n * (lags - 1L)), matrix(0, n * (lags - 1L), n))
  rbind(top, shift)
}

# Powers of the companion matrix (wimpu) and the Cholesky impulse responses
# (cimpu) for horizons 0, ..., nimp - 1.
compute_wold <- function(companion, SIGMA, n, nimp) {
  m <- nrow(companion)
  wimpu <- array(0, dim = c(m, m, nimp))
  cimpu <- array(0, dim = c(n, n, nimp))
  chol_lower <- t(chol(SIGMA))

  power <- diag(m)
  for (h in seq_len(nimp)) {
    if (h > 1L) power <- power %*% companion
    wimpu[, , h] <- power
    cimpu[, , h] <- power[seq_len(n), seq_len(n)] %*% chol_lower
  }
  list(wimpu = wimpu, cimpu = cimpu)
}

# Runs the VAR forward, y_t = offset_t + A_1 y_{t-1} + ... + A_p y_{t-p},
# for every row of `offset`. `initial` holds the p pre-sample values (oldest
# first) and `offset` everything that does not depend on lagged y:
# deterministic terms and shocks. Only the simulated periods are returned.
#
# Shock contributions in the historical decomposition, the deterministic path
# and the bootstrap samples are all produced by this one loop.
var_recursion <- function(lag_coefficients, initial, offset) {
  lags <- nrow(initial)
  n <- ncol(offset)
  transition <- t(lag_coefficients)
  offset <- t(offset)
  path <- matrix(0, nrow = n, ncol = ncol(offset))

  # state = (y_{t-1}, ..., y_{t-p}) stacked, matching the order of the lag
  # columns in X.
  state <- as.vector(t(initial[lags:1, , drop = FALSE]))
  keep <- seq_len(n * (lags - 1L))
  for (i in seq_len(ncol(offset))) {
    y <- offset[, i] + transition %*% state
    path[, i] <- y
    state <- c(y, state[keep])
  }
  t(path)
}

# The path the model implies with every shock switched off: the pre-sample
# values (read from the first row of X) carried forward by the lag dynamics
# and the deterministic terms.
deterministic_path <- function(ALPHA, X, n, lags, trend, constant) {
  n_det <- constant + trend
  first_lags <- matrix(X[1L, n_det + seq_len(n * lags)], nrow = lags,
                       byrow = TRUE)
  initial <- first_lags[lags:1, , drop = FALSE]

  offset <- matrix(0, nrow = nrow(X), ncol = n)
  if (n_det > 0L) {
    offset <- deterministic_regressors(nrow(X), constant, trend) %*%
      ALPHA[seq_len(n_det), , drop = FALSE]
  }

  path <- var_recursion(lag_coefficients(ALPHA, n, lags, constant, trend),
                        initial, offset)
  rbind(initial, path)
}
