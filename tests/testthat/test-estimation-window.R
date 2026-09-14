make_window_test_data <- function(n_obs = 120L) {
  set.seed(101)
  innovations <- matrix(stats::rnorm(n_obs * 2L), ncol = 2L)
  data <- matrix(0, nrow = n_obs, ncol = 2L)

  for (i in 3:n_obs) {
    data[i, 1] <- 0.35 * data[i - 1L, 1] -
      0.10 * data[i - 2L, 1] + innovations[i, 1]
    data[i, 2] <- 0.20 * data[i - 1L, 2] +
      0.15 * data[i - 1L, 1] + innovations[i, 2]
  }

  colnames(data) <- c("y1", "y2")
  data
}

test_that("VAR window estimates parameters but retains full-sample errors", {
  data <- make_window_test_data()
  dates <- seq.Date(as.Date("2020-01-01"), by = "day",
                    length.out = nrow(data))
  lags <- 2L
  start_date <- dates[25L]
  end_date <- dates[70L]

  fit <- etvar(
    data = data,
    dates = dates,
    target_date = dates[90L],
    lags = lags,
    optimize = FALSE,
    start_angles = 0,
    verbose = FALSE,
    start_date = start_date,
    end_date = end_date
  )

  expected_rows <- which(
    dates[(lags + 1L):length(dates)] >= start_date &
      dates[(lags + 1L):length(dates)] <= end_date
  )
  X_est <- fit$var$X[expected_rows, , drop = FALSE]
  Y_est <- fit$var$Y[expected_rows, , drop = FALSE]
  expected_alpha <- qr.solve(X_est, Y_est)
  expected_residuals <- Y_est - X_est %*% expected_alpha
  expected_sigma <- crossprod(expected_residuals) /
    (nrow(Y_est) - ncol(X_est))

  expect_equal(fit$var$estimation_rows, expected_rows)
  expect_equal(fit$estimation_dates, dates[25L:70L])
  expect_equal(fit$var$ALPHA, expected_alpha)
  expect_equal(fit$var$SIGMA, expected_sigma)
  expect_equal(nrow(fit$var$residuals), nrow(data) - lags)
  expect_equal(nrow(fit$hd$errors_str), nrow(data) - lags)
  expect_equal(nrow(fit$shocks$total), nrow(data))
  expect_true(all(is.finite(fit$hd$errors_str)))
})

test_that("NULL date bounds preserve full-sample estimation", {
  data <- make_window_test_data()
  dates <- seq.Date(as.Date("2020-01-01"), by = "day",
                    length.out = nrow(data))

  fit <- etvar(
    data = data,
    dates = dates,
    target_date = dates[90L],
    lags = 2L,
    optimize = FALSE,
    start_angles = 0,
    verbose = FALSE
  )

  expect_equal(fit$var$estimation_rows, seq_len(nrow(data) - 2L))
  expect_equal(fit$var$residuals_estimation, fit$var$residuals)
})

test_that("too-short estimation windows fail clearly", {
  data <- make_window_test_data()
  dates <- seq.Date(as.Date("2020-01-01"), by = "day",
                    length.out = nrow(data))

  expect_error(
    etvar(
      data = data,
      dates = dates,
      target_date = dates[90L],
      lags = 2L,
      optimize = FALSE,
      start_angles = 0,
      verbose = FALSE,
      start_date = dates[20L],
      end_date = dates[22L]
    ),
    "more usable observations"
  )
})
