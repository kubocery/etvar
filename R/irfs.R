#' Structural impulse responses
#'
#' Multiplies the reduced-form moving-average (Wold) multipliers by the impact
#' matrix. The response at horizon h is `Psi_h %*% A0`, so column j traces
#' the effect of a positive one-unit structural shock j (one standard
#' deviation, since the structural shocks have unit variance). The sign
#' matches the structural errors used in the historical decomposition.
#'
#' @param var_est An `etvar_var` object from [estimate_var()].
#' @param A0 Matrix (n x n). Structural impact matrix.
#' @param nimp Integer. Number of horizons, including horizon zero. Default
#'   30; at most `var_est`'s `nimp_wold`.
#' @return Array (n x n x nimp). Element `[i, j, h]` is the response of
#'   variable i to shock j at horizon h - 1.
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3)
#' var_est <- estimate_var(y, lags = 2, nimp_wold = 10)
#' A0 <- t(chol(var_est$SIGMA))  # recursive (Cholesky) identification
#' irfs <- compute_structural_irfs(var_est, A0, nimp = 10)
#' irfs[, 1, 1:3]  # responses to shock 1 at horizons 0-2
#' @export
compute_structural_irfs <- function(var_est, A0, nimp = 30L) {
  if (!inherits(var_est, "etvar_var")) {
    stop("`var_est` must be an `etvar_var` object from `estimate_var()`.",
         call. = FALSE)
  }
  n <- var_est$n
  nimp <- .etvar_integer_scalar(nimp, "nimp", min = 1L,
                                max = dim(var_est$wimpu)[3L])
  A0 <- .etvar_numeric_matrix(A0, "A0")
  if (!identical(dim(A0), c(n, n))) {
    stop("`A0` must be an `n` by `n` matrix.", call. = FALSE)
  }

  irfs <- array(0, dim = c(n, n, nimp))
  for (h in seq_len(nimp)) {
    irfs[, , h] <- var_est$wimpu[seq_len(n), seq_len(n), h] %*% A0
  }
  irfs
}


#' Bootstrap confidence bands for structural impulse responses
#'
#' Residual bootstrap: each draw runs the fitted VAR over its estimation
#' sample, starting from the observed values just before it, with residuals
#' resampled (whole rows, so the correlation between variables within a
#' period is kept), re-estimates the VAR on the simulated data, and
#' recomputes the impulse responses with the rotation of the fitted model
#' held fixed. Holding the rotation fixed keeps
#' the target shock the same shock in every draw. The bands are percentile
#' intervals across the draws.
#'
#' @param object An `etvar` object.
#' @param level Numeric. Confidence level. Default 0.95.
#' @param n_boot Integer. Number of successful draws. Default 500.
#' @param seed Integer or `NULL`. Random seed for reproducible bands. The
#'   caller's random-number state is restored afterwards.
#' @param nimp Integer or `NULL`. Number of horizons. Defaults to the number
#'   stored in `object$irfs`.
#' @param max_attempts Integer. Maximum number of draws to try. A draw is
#'   discarded and redrawn when its covariance matrix is not positive definite
#'   or the simulated data explode. Default `3 * n_boot`.
#' @return A list with `lower` and `upper` bounds (n x n x `nimp` arrays,
#'   shaped like `object$irfs` with the default `nimp`), the individual
#'   bootstrap responses in `draws` (n x n x `nimp` x `n_boot`), and the
#'   settings (`level`, `n_boot` successful draws, `n_attempts`, `method`).
#'   [plot_etvar_irfs()] uses `draws` to normalise each draw on its own when
#'   the responses are rescaled with `scaling_index`.
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3,
#'             dimnames = list(NULL, c("x1", "x2", "x3")))
#' dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 100)
#' fit <- etvar(y, dates, target_date = dates[80], lags = 2, nimp = 10)
#' ci <- compute_irf_confint(fit, level = 0.9, n_boot = 50, seed = 1)
#' ci$lower[, 1, 1:3]
#' plot(fit, type = "irfs", irf_ci = ci)
#' @export
compute_irf_confint <- function(object, level = 0.95, n_boot = 500L,
                                seed = NULL, nimp = NULL,
                                max_attempts = n_boot * 3L) {
  if (!inherits(object, "etvar")) {
    stop("`object` must be an `etvar` object.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be a number between 0 and 1.", call. = FALSE)
  }
  n_boot <- .etvar_integer_scalar(n_boot, "n_boot", min = 1L)
  if (is.null(nimp)) nimp <- dim(object$irfs)[3L]
  nimp <- .etvar_integer_scalar(nimp, "nimp", min = 1L)
  max_attempts <- .etvar_integer_scalar(max_attempts, "max_attempts",
                                        min = n_boot)

  n <- object$var$n
  Q <- object$rotation
  if (!is.numeric(Q) || !identical(dim(Q), c(n, n))) {
    stop("`object$rotation` must be an `n` by `n` matrix.", call. = FALSE)
  }

  if (!is.null(seed)) {
    seed <- .etvar_integer_scalar(seed, "seed", min = 0L,
                                  max = .Machine$integer.max)
    old_rng <- .etvar_get_rng()
    on.exit(.etvar_restore_rng(old_rng), add = TRUE)
    set.seed(seed)
  }

  var_est <- object$var
  draws <- array(NA_real_, dim = c(n, n, nimp, n_boot))
  n_good <- 0L
  n_try <- 0L
  while (n_good < n_boot && n_try < max_attempts) {
    n_try <- n_try + 1L
    boot_irfs <- tryCatch({
      boot_data <- simulate_var_residual_bootstrap(var_est)
      if (any(!is.finite(boot_data))) stop("explosive bootstrap sample")
      bootstrap_irfs(boot_data, var_est, Q, nimp)
    }, error = function(e) NULL)

    if (is.null(boot_irfs) || any(!is.finite(boot_irfs))) next
    n_good <- n_good + 1L
    draws[, , , n_good] <- boot_irfs
  }

  if (n_good == 0L) {
    stop("No valid bootstrap IRF draws were produced.", call. = FALSE)
  }
  if (n_good < n_boot) {
    warning("Only ", n_good, " valid bootstrap draws were produced after ",
            n_try, " attempts.", call. = FALSE)
    draws <- draws[, , , seq_len(n_good), drop = FALSE]
  }

  alpha <- 1 - level
  flat <- matrix(draws, ncol = n_good)
  list(
    lower = array(row_quantiles(flat, alpha / 2), dim = c(n, n, nimp)),
    upper = array(row_quantiles(flat, 1 - alpha / 2), dim = c(n, n, nimp)),
    draws = draws,
    level = level,
    n_boot = n_good,
    n_attempts = n_try,
    method = "residual bootstrap",
    fixed_theta = TRUE
  )
}


# Bootstrap helpers ------------------------------------------------------------

# One bootstrap sample. The fitted VAR is run over the estimation window,
# starting from the observed values just before it, with centred residuals
# from the window drawn with replacement. Rows outside the window keep the
# observed data; the re-estimation does not use them.
simulate_var_residual_bootstrap <- function(var_est) {
  data <- var_est$data
  lags <- var_est$lags
  n_obs <- nrow(data) - lags
  rows <- var_est$estimation_rows
  if (is.null(rows)) rows <- seq_len(n_obs)
  simulated <- seq.int(min(rows), max(rows))

  pool <- var_est$residuals_estimation
  if (is.null(pool)) pool <- var_est$residuals
  pool <- sweep(pool, 2L, colMeans(pool))
  offset <- pool[sample.int(nrow(pool), size = length(simulated),
                            replace = TRUE), , drop = FALSE]

  n_det <- var_est$constant + var_est$trend
  if (n_det > 0L) {
    regressors <- deterministic_regressors(n_obs, var_est$constant,
                                           var_est$trend)
    offset <- offset + regressors[simulated, , drop = FALSE] %*%
      var_est$ALPHA[seq_len(n_det), , drop = FALSE]
  }

  # Effective row r is raw row r + lags, so its lags are raw rows
  # r, ..., r + lags - 1.
  initial <- data[min(rows) - 1L + seq_len(lags), , drop = FALSE]
  slopes <- lag_coefficients(var_est$ALPHA, var_est$n, lags,
                             var_est$constant, var_est$trend)
  out <- data
  out[simulated + lags, ] <- var_recursion(slopes, initial, offset)
  out
}

# Impulse responses for one bootstrap sample with the rotation Q held fixed.
# A slimmed-down estimate_var(): one QR decomposition, and only the first n
# rows of the companion powers.
bootstrap_irfs <- function(boot_data, var_est, Q, nimp) {
  n <- var_est$n
  lags <- var_est$lags
  constant <- var_est$constant
  trend <- var_est$trend

  Y <- trim_to_effective(boot_data, lags)
  X <- build_var_regressors(boot_data, lags, trend, constant)
  rows <- var_est$estimation_rows
  if (!is.null(rows)) {
    X <- X[rows, , drop = FALSE]
    Y <- Y[rows, , drop = FALSE]
  }

  qr_x <- qr(X)
  if (qr_x$rank < ncol(X)) stop("rank-deficient bootstrap regressors")
  ALPHA <- qr.coef(qr_x, Y)
  SIGMA <- crossprod(Y - X %*% ALPHA) / (nrow(X) - ncol(X))
  A0 <- -t(chol(SIGMA)) %*% Q

  companion <- build_companion(ALPHA, n, lags, constant, trend)
  power_rows <- diag(1, nrow = n, ncol = n * lags)
  irfs <- array(0, dim = c(n, n, nimp))
  for (h in seq_len(nimp)) {
    if (h > 1L) power_rows <- power_rows %*% companion
    irfs[, , h] <- power_rows[, seq_len(n), drop = FALSE] %*% A0
  }
  irfs
}

# Type-7 quantiles (the stats::quantile() default) of every row of `x`, from
# a single sort instead of one quantile() call per row.
row_quantiles <- function(x, prob) {
  n_draws <- ncol(x)
  sorted <- matrix(x[order(row(x), x)], nrow = nrow(x), ncol = n_draws,
                   byrow = TRUE)
  index <- 1 + (n_draws - 1) * prob
  lo <- floor(index)
  hi <- ceiling(index)
  h <- index - lo
  out <- sorted[, lo]
  upper <- sorted[, hi]
  interpolate <- which(index > lo & upper != out)
  out[interpolate] <- (1 - h) * out[interpolate] + h * upper[interpolate]
  out
}
