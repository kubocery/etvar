make_irf_fit <- function() {
  set.seed(510)
  data <- matrix(stats::rnorm(240L), ncol = 3L,
                 dimnames = list(NULL, c("a", "b", "c")))
  dates <- seq(as.Date("2022-01-01"), by = "day", length.out = nrow(data))
  etvar(data, dates, target_date = dates[60L], lags = 2L, nimp = 6L)
}

test_that("impulse responses are responses to a positive structural shock", {
  fit <- make_irf_fit()
  event_row <- fit$target_obs - fit$var$lags

  # A0 maps the structural errors into the residuals ...
  expect_equal(
    drop(fit$A0 %*% fit$hd$errors_str[event_row, ]),
    fit$var$residuals[event_row, ],
    ignore_attr = TRUE
  )
  # ... and the impact response to shock j is A0[, j]: a positive error moves
  # the variables the way the impulse responses say.
  expect_equal(fit$irfs[, , 1L], unname(fit$A0))
})

test_that("the horizon-one response applies the first lag matrix to A0", {
  fit <- make_irf_fit()
  first_lag <- t(fit$var$ALPHA[1L + seq_len(3L), ])

  expect_equal(fit$irfs[, , 2L], unname(first_lag %*% fit$A0))
})

test_that("bootstrap bands are reproducible and leave the caller's RNG alone", {
  fit <- make_irf_fit()
  set.seed(1)
  expected_next <- stats::runif(1L)

  set.seed(1)
  first <- compute_irf_confint(fit, n_boot = 20L, seed = 42L)
  expect_identical(stats::runif(1L), expected_next)

  second <- compute_irf_confint(fit, n_boot = 20L, seed = 42L)
  expect_identical(first$lower, second$lower)
  expect_identical(first$upper, second$upper)
  expect_identical(dim(first$lower), dim(fit$irfs))
  expect_equal(sum(first$lower > first$upper), 0)
  expect_identical(first$n_boot, 20L)
})

test_that("bootstrap quantiles match stats::quantile()", {
  set.seed(3)
  draws <- matrix(stats::rnorm(200L), nrow = 10L)
  draws[1L, ] <- 1
  draws[2L, 1:10] <- 0

  for (p in c(0, 0.05, 0.5, 0.975, 1)) {
    expect_equal(
      etvar:::row_quantiles(draws, p),
      apply(draws, 1L, stats::quantile, probs = p, names = FALSE)
    )
  }
})

test_that("normalised bands are computed from normalised draws", {
  skip_if_not_installed("ggplot2")
  fit <- make_irf_fit()
  ci <- compute_irf_confint(fit, n_boot = 30L, seed = 4L)

  pl <- plot_etvar_irfs(fit$irfs, fit$variable_names,
                        scaling_index = c(2, 1, 1),
                        scaling_multiplier = 0.5, irf_ci = ci)
  anchor <- pl$data[pl$data$Variable == "b" & pl$data$Horizon == 0L, ]

  expect_equal(anchor$Response, 0.5)
  expect_equal(c(anchor$Lower, anchor$Upper), c(0.5, 0.5))
  expect_equal(sum(pl$data$Lower > pl$data$Upper), 0)
})

test_that("bootstrap samples start from the data at the estimation window", {
  set.seed(515)
  data <- matrix(stats::rnorm(240L), ncol = 2L)
  var_est <- estimate_var(data, lags = 2L, nimp_wold = 2L,
                          estimation_rows = 41:100)
  sample <- etvar:::simulate_var_residual_bootstrap(var_est)

  untouched <- c(seq_len(42L), 103:120)
  expect_identical(sample[untouched, ], data[untouched, ])
  expect_false(isTRUE(all.equal(sample[43:102, ], data[43:102, ])))
})
