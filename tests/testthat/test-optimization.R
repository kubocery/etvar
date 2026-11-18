test_that("two-variable models are optimised globally and without warnings", {
  # On this series the earlier Nelder-Mead/BFGS search stopped in a local
  # minimum (0.29); the global minimum is about 0.077.
  set.seed(11)
  data <- apply(matrix(stats::rnorm(400L), ncol = 2L), 2L, cumsum) / 10
  var_est <- estimate_var(data, lags = 2L, nimp_wold = 3L)

  expect_no_warning(
    result <- etvar_optimize(var_est = var_est, target_obs = 148L,
                             shock_index = 2L)
  )

  context <- etvar:::.etvar_prepare_objective(var_est, 148L, shock_index = 2L)
  grid <- seq(-pi / 2, pi / 2, length.out = 2001L)
  grid_values <- vapply(grid, etvar:::.etvar_objective_prepared, numeric(1),
                        context = context, parameterization = "target")
  expect_lte(result$value, min(grid_values) + 1e-8)
  expect_identical(result$optimization$objective_by_start_stage$method,
                   c("grid search", "Brent"))
})

test_that("a direction and its negative have the same objective", {
  set.seed(512)
  var_est <- estimate_var(matrix(stats::rnorm(300L), ncol = 3L), lags = 1L,
                          nimp_wold = 2L)
  theta <- c(0.4, -1.1)
  direction <- etvar:::.etvar_direction_from_angles(theta, 3L, 1L)
  flipped <- etvar:::.etvar_angles_from_direction(-direction, 1L)

  expect_equal(etvar_objective(flipped, var_est, 50L),
               etvar_objective(theta, var_est, 50L))
})

test_that("the optimised target shock is positive on the first event date", {
  set.seed(512)
  var_est <- estimate_var(matrix(stats::rnorm(300L), ncol = 3L), lags = 1L,
                          nimp_wold = 2L)

  result <- suppressWarnings(etvar_optimize(
    var_est = var_est, target_obs = c(70L, 40L), shock_index = 2L
  ))
  A0 <- -t(chol(var_est$SIGMA)) %*% result$rotation
  errors <- t(solve(A0, t(var_est$residuals)))

  expect_gt(errors[40L, 2L], 0)
  expect_equal(result$rotation[, 2L], result$target_direction)
})

test_that("starts are reproducible and the caller's RNG is restored", {
  set.seed(513)
  var_est <- estimate_var(matrix(stats::rnorm(300L), ncol = 3L), lags = 1L,
                          nimp_wold = 2L)
  set.seed(2)
  expected_next <- stats::runif(1L)

  set.seed(2)
  first <- suppressWarnings(etvar_optimize(
    var_est = var_est, target_obs = 50L, n_starts = 3L, seed = 9L,
    control_nm = list(maxit = 50L), control_bfgs = list(maxit = 5L)
  ))
  expect_identical(stats::runif(1L), expected_next)

  second <- suppressWarnings(etvar_optimize(
    var_est = var_est, target_obs = 50L, n_starts = 3L, seed = 9L,
    control_nm = list(maxit = 50L), control_bfgs = list(maxit = 5L)
  ))
  expect_identical(first$optimization$start_values,
                   second$optimization$start_values)
  expect_identical(first$value, second$value)
  expect_identical(nrow(first$optimization$start_values), 3L)
})
