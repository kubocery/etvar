# ==============================================================================
# Structural impulse response functions
# ==============================================================================

#' Compute structural impulse response functions
#'
#' Multiplies each Wold impulse-response matrix by the structural impact
#' matrix A0 to obtain the response of every variable to every structural
#' shock at each horizon.
#'
#' @param var_est An `etvar_var` object from [estimate_var()].
#' @param A0      Matrix (n x n). Structural impact matrix.
#' @param nimp    Integer. Number of impulse-response horizons, including
#'   horizon zero. Default 30.
#' @return Array (n x n x nimp) where element `[i, j, h]` is the response
#'         of variable `i` to structural shock `j` at array position `h`;
#'         position 1 is horizon zero.
#' @export
compute_structural_irfs <- function(var_est, A0, nimp = 30L) {
  if (!inherits(var_est, "etvar_var")) {
    stop("`var_est` must be an `etvar_var` object from `estimate_var()`.",
         call. = FALSE)
  }
  n    <- var_est$n
  nimp <- .etvar_integer_scalar(nimp, "nimp", min = 1L,
                                max = dim(var_est$wimpu)[3L])
  A0 <- .etvar_numeric_matrix(A0, "A0")
  if (!identical(dim(A0), c(n, n))) {
    stop("`A0` must be an `n` by `n` matrix.", call. = FALSE)
  }
  irfs <- array(0, dim = c(n, n, nimp))

  for (j in seq_len(nimp)) {
    irfs[, , j] <- -var_est$wimpu[1:n, 1:n, j] %*% A0
  }

  irfs
}


#' Simulate one residual-bootstrap VAR sample
#'
#' @param var_est An `etvar_var` object from [estimate_var()].
#' @return Matrix with the same dimensions as `var_est$data`.
#' @keywords internal
simulate_var_residual_bootstrap <- function(var_est) {
  data     <- var_est$data
  lags     <- var_est$lags
  T_raw    <- nrow(data)
  T_eff    <- T_raw - lags
  n        <- var_est$n
  residual_pool <- if (!is.null(var_est$residuals_estimation)) {
    var_est$residuals_estimation
  } else {
    var_est$residuals
  }
  residuals <- sweep(residual_pool, 2, colMeans(residual_pool), "-")

  if (lags < 1L) {
    stop("Residual-bootstrap IRF intervals require `lags` >= 1.")
  }

  boot_data <- matrix(NA_real_, nrow = T_raw, ncol = n)
  boot_data[seq_len(lags), ] <- data[seq_len(lags), , drop = FALSE]

  draw_idx <- sample.int(nrow(residuals), size = T_eff, replace = TRUE)
  boot_residuals <- residuals[draw_idx, , drop = FALSE]

  # Split ALPHA into deterministic rows (constant, trend, squared trend) and
  # lag rows. The deterministic part and the resampled residuals do not
  # depend on the simulated path, so they are added for all periods at once.
  ALPHA <- var_est$ALPHA
  n_det <- var_est$constant + var_est$trend
  lag_coefficients <- ALPHA[n_det + seq_len(n * lags), , drop = FALSE]
  offset <- boot_residuals
  if (n_det > 0L) {
    tt <- seq_len(T_eff)
    deterministic <- cbind(
      if (var_est$constant == 1L) rep(1, T_eff),
      if (var_est$trend >= 1L) tt,
      if (var_est$trend >= 2L) tt^2
    )
    offset <- offset +
      deterministic %*% ALPHA[seq_len(n_det), , drop = FALSE]
  }

  # Companion-form recursion: `state` stacks (y_{t-1}, ..., y_{t-p}), matching
  # the column order of the lag block in the regressor matrix. The path is
  # built column-wise (n x T_eff), which is much cheaper than row assignment.
  transition <- t(lag_coefficients)
  offset <- t(offset)
  path <- matrix(0, nrow = n, ncol = T_eff)
  state <- as.vector(t(data[lags:1L, , drop = FALSE]))
  keep <- seq_len(n * (lags - 1L))
  for (tt in seq_len(T_eff)) {
    y <- offset[, tt] + transition %*% state
    path[, tt] <- y
    state <- c(y, state[keep])
  }
  boot_data[lags + seq_len(T_eff), ] <- t(path)

  colnames(boot_data) <- colnames(data)
  boot_data
}


# Re-estimate the VAR on one bootstrap sample and return its structural IRFs
# for a fixed rotation `Q`. This is a lean version of estimate_var() plus
# compute_structural_irfs(): it performs a single QR decomposition and builds
# only the first n rows of the companion-matrix powers needed for the IRFs.
.etvar_bootstrap_irfs <- function(boot_data, lags, constant, trend,
                                  estimation_rows, Q, nimp) {
  n <- ncol(boot_data)
  Y <- trim_to_effective(boot_data, lags)
  X <- build_var_regressors(boot_data, nrow(boot_data), n, lags, trend,
                            constant)
  if (!is.null(estimation_rows)) {
    X <- X[estimation_rows, , drop = FALSE]
    Y <- Y[estimation_rows, , drop = FALSE]
  }
  n_parameters <- ncol(X)
  dof <- nrow(Y) - n_parameters
  if (dof <= 0L) stop("Too few observations in the bootstrap sample.")

  qr_x <- qr(X)
  if (qr_x$rank < n_parameters) stop("Rank-deficient bootstrap regressors.")
  ALPHA <- qr.coef(qr_x, Y)
  residuals <- Y - X %*% ALPHA
  SIGMA <- crossprod(residuals) / dof
  A0 <- -t(chol(SIGMA)) %*% Q

  companion <- build_companion(ALPHA, n, lags, constant, trend)
  irfs <- array(0, dim = c(n, n, nimp))
  power_rows <- diag(1, nrow = n, ncol = n * lags)
  for (h in seq_len(nimp)) {
    if (h > 1L) power_rows <- power_rows %*% companion
    irfs[, , h] <- -power_rows[, seq_len(n), drop = FALSE] %*% A0
  }
  irfs
}


# Type-7 quantiles (the stats::quantile() default) of every row of `x`,
# computed with one global sort instead of one quantile() call per row.
.etvar_row_quantiles <- function(x, prob) {
  n_draws <- ncol(x)
  sorted <- matrix(x[order(row(x), x)], nrow = nrow(x), ncol = n_draws,
                   byrow = TRUE)
  index <- 1 + (n_draws - 1) * prob
  lo <- floor(index)
  hi <- ceiling(index)
  out <- sorted[, lo]
  h <- index - lo
  upper <- sorted[, hi]
  interpolate <- which(index > lo & upper != out)
  out[interpolate] <- (1 - h) * out[interpolate] + h * upper[interpolate]
  out
}


#' Compute bootstrap confidence intervals for structural IRFs
#'
#' Builds residual-bootstrap samples from the fitted reduced-form VAR,
#' re-estimates the VAR on each bootstrap sample, and recomputes structural
#' IRFs using the fitted ETVAR target direction and its stored orthogonal
#' completion. The rotation is held fixed so the structural shock labels stay
#' aligned with the fitted model.
#'
#' @param object An `etvar` object.
#' @param level Numeric. Confidence level. Default `0.95`.
#' @param n_boot Integer. Number of successful bootstrap draws. Default `500`.
#' @param seed Integer or `NULL`. Optional random seed. The caller's
#'   random-number state is restored on exit.
#' @param nimp Integer or `NULL`. Number of IRF horizons. Default: the number
#'   of horizons stored in `object$irfs`.
#' @param max_attempts Integer. Maximum bootstrap attempts. Extra attempts are
#'   used when a draw fails because the bootstrap covariance matrix is not
#'   positive definite.
#' @return A list with `lower` and `upper` arrays, each with the same
#'   dimensions as `object$irfs`, plus bootstrap metadata.
#' @export
compute_irf_confint <- function(object, level = 0.95, n_boot = 500L,
                                seed = NULL, nimp = NULL,
                                max_attempts = n_boot * 3L) {
  if (!inherits(object, "etvar")) {
    stop("`object` must be an `etvar` object.")
  }
  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be a numeric scalar between 0 and 1.")
  }

  n_boot <- .etvar_integer_scalar(n_boot, "n_boot", min = 1L)

  if (is.null(nimp)) nimp <- dim(object$irfs)[3]
  nimp <- .etvar_integer_scalar(nimp, "nimp", min = 1L)

  max_attempts <- .etvar_integer_scalar(max_attempts, "max_attempts",
                                         min = n_boot)
  if (!is.null(seed)) {
    seed <- .etvar_integer_scalar(seed, "seed", min = 0L,
                                  max = .Machine$integer.max)
  }

  if (!is.null(seed)) {
    old_rng <- .etvar_get_rng()
    on.exit(.etvar_restore_rng(old_rng), add = TRUE)
    set.seed(seed)
  }

  n <- object$var$n
  if (!is.null(object$rotation) &&
      is.numeric(object$rotation) &&
      identical(dim(object$rotation), c(n, n))) {
    Q <- object$rotation
  } else if (identical(object$parameterization, "target_direction")) {
    Q <- .etvar_target_rotation(object$theta, n, object$shock_index)
  } else {
    # Compatibility with etvar objects created before target-direction
    # parameterisation was introduced.
    Q <- euler_rotation(object$theta, n)
  }
  draws  <- array(NA_real_, dim = c(n, n, nimp, n_boot))
  n_good <- 0L
  n_try  <- 0L

  while (n_good < n_boot && n_try < max_attempts) {
    n_try <- n_try + 1L

    boot_irfs <- tryCatch({
      boot_data <- simulate_var_residual_bootstrap(object$var)
      if (any(!is.finite(boot_data))) stop("Non-finite bootstrap sample.")
      .etvar_bootstrap_irfs(
        boot_data       = boot_data,
        lags            = object$var$lags,
        constant        = object$var$constant,
        trend           = object$var$trend,
        estimation_rows = object$var$estimation_rows,
        Q               = Q,
        nimp            = nimp
      )
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
  draw_matrix <- matrix(draws, nrow = n * n * nimp, ncol = n_good)
  lower <- array(.etvar_row_quantiles(draw_matrix, alpha / 2),
                 dim = c(n, n, nimp))
  upper <- array(.etvar_row_quantiles(draw_matrix, 1 - alpha / 2),
                 dim = c(n, n, nimp))

  list(
    lower        = lower,
    upper        = upper,
    level        = level,
    n_boot       = n_good,
    n_attempts   = n_try,
    method       = "residual bootstrap",
    fixed_theta  = TRUE
  )
}
