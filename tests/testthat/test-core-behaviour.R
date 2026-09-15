make_core_test_data <- function(n_obs = 90L) {
  set.seed(202)
  innovations <- matrix(stats::rnorm(n_obs * 2L), ncol = 2L)
  data <- matrix(0, nrow = n_obs, ncol = 2L)
  for (i in 3:n_obs) {
    data[i, 1L] <- 0.4 * data[i - 1L, 1L] -
      0.1 * data[i - 2L, 1L] + innovations[i, 1L]
    data[i, 2L] <- 0.2 * data[i - 1L, 1L] +
      0.3 * data[i - 1L, 2L] + innovations[i, 2L]
  }
  colnames(data) <- c("y1", "y2")
  data
}

test_that("the objective targets the requested raw-sample event", {
  data <- make_core_test_data()
  var_est <- estimate_var(data, lags = 2L, nimp_wold = 5L)
  target_eff <- 50L
  theta <- 0

  A0 <- -t(chol(var_est$SIGMA)) %*% euler_rotation(theta, var_est$n)
  decomp <- hist_decomposition(
    A0, var_est$ALPHA, var_est$Y, var_est$X, var_est$n,
    var_est$lags, var_est$trend, var_est$constant
  )
  shock_contribution <- decomp$hist_matrix[, , 1L]
  counterfactual <- decomp$deterministic + shock_contribution
  scales <- apply(diff(data), 2L, stats::sd)
  raw_diff_row <- target_eff + var_est$lags - 1L
  expected <- sum(abs(
    diff(counterfactual)[raw_diff_row, ] / scales -
      diff(data)[raw_diff_row, ] / scales
  ))

  expect_equal(etvar_objective(theta, var_est, target_eff), expected)
})

test_that("historical decomposition reconstructs the observed data", {
  data <- make_core_test_data()
  var_est <- estimate_var(data, lags = 2L, nimp_wold = 5L)
  A0 <- -t(chol(var_est$SIGMA))
  decomp <- hist_decomposition(
    A0, var_est$ALPHA, var_est$Y, var_est$X, var_est$n,
    var_est$lags, var_est$trend, var_est$constant
  )

  shock_sum <- apply(decomp$hist_matrix, c(1L, 2L), sum)
  expect_equal(unname(decomp$deterministic + shock_sum), unname(data),
               tolerance = 1e-8)
})

test_that("public estimators reject invalid scalar and data inputs", {
  data <- make_core_test_data()

  expect_error(estimate_var(data, lags = 0), "greater than or equal to 1")
  expect_error(estimate_var(data, lags = 1.5), "must be one integer")
  expect_error(estimate_var(cbind(data, NA_real_), lags = 2), "finite")
  expect_error(estimate_var(data.frame(x = 1:10, y = letters[1:10]), 2),
               "Every column")
  expect_error(euler_rotation(0, 3), "Length of `theta`")
  expect_identical(euler_rotation(numeric(), 1L), matrix(1, 1L, 1L))
})

test_that("etvar computes as many Wold horizons as requested", {
  data <- make_core_test_data()
  dates <- seq.Date(as.Date("2021-01-01"), by = "day", length.out = nrow(data))

  fit <- etvar(
    data, dates, target_date = dates[70L], lags = 2L,
    optimize = FALSE, nimp = 105L, verbose = FALSE
  )

  expect_equal(dim(fit$irfs)[3L], 105L)
  expect_equal(dim(fit$var$wimpu)[3L], 105L)
})

test_that("bootstrap intervals and plotters handle fitted output", {
  skip_if_not_installed("ggplot2")
  data <- make_core_test_data()
  dates <- seq.Date(as.Date("2021-01-01"), by = "day", length.out = nrow(data))
  fit <- etvar(
    data, dates, target_date = dates[70L], lags = 2L,
    optimize = FALSE, nimp = 4L, verbose = FALSE
  )

  intervals <- compute_irf_confint(fit, n_boot = 2L, seed = 1L,
                                   max_attempts = 6L)
  expect_identical(dim(intervals$lower), dim(fit$irfs))
  expect_s3_class(
    plot_etvar_counterfactual(dates, data[, 1L],
                              fit$hd$all_data[, 1L, 1L]),
    "ggplot"
  )
  expect_s3_class(
    plot_etvar_errors(dates[-seq_len(fit$var$lags)],
                      fit$hd$errors_str[, 1L]),
    "ggplot"
  )

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)
  old_mar <- graphics::par("mar")
  plot_etvar_irfs(fit$irfs, fit$variable_names)
  expect_equal(graphics::par("mar"), old_mar)
  plot_etvar_hd(
    dates, data, fit$hd$deterministic, fit$shocks$target_stop,
    fit$shocks$target_rest, fit$shocks$others, fit$datestar,
    fit$variable_names, ask = FALSE
  )
  expect_equal(graphics::par("mar"), old_mar)
})
