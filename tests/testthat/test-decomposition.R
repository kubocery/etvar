make_decomposition_data <- function() {
  set.seed(611)
  data <- matrix(stats::rnorm(240L), ncol = 3L)
  data[, 1L] <- cumsum(data[, 1L]) / 5
  data
}

test_that("decompositions add up to the data for every deterministic term", {
  data <- make_decomposition_data()
  specifications <- expand.grid(constant = c(FALSE, TRUE), trend = 0:2)

  for (i in seq_len(nrow(specifications))) {
    var_est <- estimate_var(data, lags = 3L,
                            constant = specifications$constant[i],
                            trend = specifications$trend[i], nimp_wold = 2L)
    hd <- hist_decomposition(
      -t(chol(var_est$SIGMA)), var_est$ALPHA, var_est$Y, var_est$X, 3L, 3L,
      var_est$trend, var_est$constant
    )
    expect_equal(hd$deterministic + rowSums(hd$hist_matrix, dims = 2L),
                 data, ignore_attr = TRUE)
  }
})

test_that("the intervention path only removes errors after the cut-off", {
  data <- make_decomposition_data()
  var_est <- estimate_var(data, lags = 2L, nimp_wold = 2L)
  A0 <- -t(chol(var_est$SIGMA))
  args <- list(A0, var_est$ALPHA, var_est$Y, var_est$X, 3L, 2L, 0L, 1L)
  n_obs <- nrow(data)

  plain <- do.call(hist_decomposition, args)
  cut <- do.call(hist_decomposition_3d, c(args, datestar = 10L))
  uncut <- do.call(hist_decomposition_3d, c(args, datestar = 0L))

  expect_equal(cut[names(plain)], plain)
  expect_equal(uncut$hist_matrix_stop, uncut$hist_matrix)
  before <- seq_len(n_obs - 10L)
  expect_equal(cut$hist_matrix_stop[before, , ], cut$hist_matrix[before, , ])
  expect_false(isTRUE(all.equal(cut$hist_matrix_stop[n_obs, , ],
                                cut$hist_matrix[n_obs, , ])))
  expect_equal(cut$all_data[, , 2L],
               cut$deterministic + cut$hist_matrix[, , 2L])
})

test_that("rearrange_shocks() splits the contributions consistently", {
  data <- make_decomposition_data()
  var_est <- estimate_var(data, lags = 2L, nimp_wold = 2L)
  hd <- hist_decomposition_3d(
    -t(chol(var_est$SIGMA)), var_est$ALPHA, var_est$Y, var_est$X, 3L, 2L,
    0L, 1L, datestar = 10L
  )

  parts <- rearrange_shocks(hd$hist_matrix, hd$hist_matrix_stop,
                            T_obs = nrow(data), n = 3L, target_shock = 2L)
  expect_equal(parts$target, hd$hist_matrix[, , 2L])
  expect_equal(parts$total, rowSums(hd$hist_matrix, dims = 2L))
  expect_equal(parts$target + parts$others, parts$total)
  expect_equal(parts$target_stop + parts$target_rest, parts$target)

  no_stop <- rearrange_shocks(hd$hist_matrix, NULL, nrow(data), 3L, 2L)
  expect_null(no_stop$target_stop)
  expect_null(no_stop$target_rest)
})

test_that("the event-day shock is separated from earlier target shocks", {
  data <- make_decomposition_data()
  dates <- seq(as.Date("2022-01-01"), by = "day", length.out = nrow(data))
  fit <- etvar(data, dates, target_date = dates[60L], lags = 2L, nimp = 3L)
  parts <- fit$shock_contributions
  event_row <- fit$target_obs
  error <- fit$hd$errors_str[event_row - fit$var$lags, fit$shock_index]

  expect_equal(sum(abs(parts$target_event[seq_len(event_row - 1L), ])), 0)
  expect_equal(parts$target_event[event_row, ],
               unname(fit$A0[, fit$shock_index]) * error)
  expect_equal(parts$target_before + parts$target_event, parts$target_stop)
  expect_equal(
    parts$target_event + parts$target_rest + parts$others +
      parts$target_before,
    unname(fit$data - fit$hd$deterministic)
  )
})
