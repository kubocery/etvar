make_fevd_test_object <- function(irfs, target_shock = 1L) {
  structure(
    list(
      irfs = irfs,
      shock_index = target_shock,
      variable_names = paste0("y", seq_len(dim(irfs)[1L])),
      var = list(n = dim(irfs)[1L])
    ),
    class = "etvar"
  )
}

make_fevd_var_data <- function(n_obs = 90L) {
  set.seed(303)
  innovations <- matrix(stats::rnorm(n_obs * 2L), ncol = 2L)
  data <- matrix(0, nrow = n_obs, ncol = 2L)
  for (i in 3:n_obs) {
    data[i, 1L] <- 0.4 * data[i - 1L, 1L] + innovations[i, 1L]
    data[i, 2L] <- 0.2 * data[i - 1L, 1L] +
      0.3 * data[i - 1L, 2L] + innovations[i, 2L]
  }
  colnames(data) <- c("y1", "y2")
  data
}

test_that("structural FEVD matches cumulative squared IRF shares", {
  irfs <- array(
    c(
      1, 2, 2, 1,
      2, 0, 1, 3,
      1, 1, 2, 2
    ),
    dim = c(2L, 2L, 3L)
  )
  object <- make_fevd_test_object(irfs)

  result <- compute_fevd(object)

  expected_h1 <- sweep(
    irfs[, , 1L]^2,
    1L,
    rowSums(irfs[, , 1L]^2),
    "/"
  )
  cumulative_h2 <- irfs[, , 1L]^2 + irfs[, , 2L]^2
  expected_h2 <- sweep(cumulative_h2, 1L, rowSums(cumulative_h2), "/")

  expect_s3_class(result, "etvar_fevd")
  expect_identical(dim(result$shares), c(2L, 2L, 3L))
  expect_equal(unname(result$shares[, , 1L]), expected_h1)
  expect_equal(unname(result$shares[, , 2L]), expected_h2)
  expect_equal(unname(result$target), unname(result$shares[, 1L, ]))
  expect_equal(unname(result$target + result$others), matrix(1, 2L, 3L))
  expect_equal(
    unname(apply(result$shares, c(1L, 3L), sum)),
    matrix(1, 2L, 3L)
  )
})

test_that("FEVD can isolate a non-default target and is sign invariant", {
  irfs <- array(seq_len(18L), dim = c(3L, 3L, 2L))
  object <- make_fevd_test_object(irfs, target_shock = 2L)

  result <- compute_fevd(object, nimp = 1L)
  sign_reversed <- object
  sign_reversed$irfs <- -object$irfs

  expect_identical(result$target_shock, 2L)
  expect_identical(result$horizons, 1L)
  expect_equal(
    unname(result$target),
    matrix(unname(result$shares[, 2L, 1L]), ncol = 1L)
  )
  expect_equal(result$shares, compute_fevd(sign_reversed, nimp = 1L)$shares)
})

test_that("one-shock FEVD assigns zero to other shocks", {
  irfs <- array(c(2, 1, 0.5), dim = c(1L, 1L, 3L))
  result <- compute_fevd(make_fevd_test_object(irfs))

  expect_equal(unname(result$target), matrix(1, 1L, 3L))
  expect_equal(unname(result$others), matrix(0, 1L, 3L))
})

test_that("FEVD validates horizons and plot input", {
  skip_if_not_installed("ggplot2")
  irfs <- array(seq_len(12L), dim = c(2L, 2L, 3L))
  object <- make_fevd_test_object(irfs)
  result <- compute_fevd(object)

  expect_error(compute_fevd(object, nimp = 4L), "between 1 and 3")
  expect_error(compute_fevd(object, target_shock = 3L), "between 1 and 2")
  expect_s3_class(plot_etvar_fevd(result), "ggplot")
  expect_s3_class(
    plot_etvar_fevd(result, variables = 2L, title = NA),
    "ggplot"
  )
  expect_error(plot_etvar_fevd(result, variables = 3L), "between 1 and 2")
})

test_that("plot.etvar dispatches to the FEVD plot", {
  skip_if_not_installed("ggplot2")
  data <- make_fevd_var_data()
  dates <- seq.Date(as.Date("2021-01-01"), by = "day", length.out = nrow(data))
  fit <- etvar(
    data, dates, target_date = dates[70L], lags = 2L,
    optimize = FALSE, nimp = 4L, verbose = FALSE
  )

  result <- compute_fevd(fit)
  expect_equal(unname(result$target + result$others), matrix(1, 2L, 4L))
  expect_s3_class(plot(fit, type = "fevd"), "ggplot")
  expect_s3_class(plot(fit, type = "fevd", variables = 1L), "ggplot")
})
