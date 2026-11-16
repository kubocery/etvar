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

test_that("verbose sequential optimisation enables live optim tracing", {
  data <- make_core_test_data(n_obs = 45L)
  var_est <- estimate_var(data, lags = 2L, nimp_wold = 3L)

  output <- capture.output(
    suppressMessages(suppressWarnings(
      etvar_optimize(
        start = 0,
        var_est = var_est,
        target_obs = 30L,
        control_nm = list(reltol = 1e-6, maxit = 1L, trace = 0),
        control_bfgs = list(
          reltol = 1e-6, maxit = 1L, trace = 0, REPORT = 1L
        ),
        verbose = TRUE
      )
    ))
  )

  expect_true(any(grepl("Nelder-Mead direct search", output, fixed = TRUE)))
  expect_true(any(grepl("initial", output, fixed = TRUE)))
  expect_true(any(grepl("value", output, fixed = TRUE)))
})

test_that("target-direction parameterisation replaces the full rotation", {
  set.seed(404)
  data <- matrix(stats::rnorm(240L), ncol = 3L)
  colnames(data) <- c("y1", "y2", "y3")
  dates <- seq.Date(as.Date("2022-01-01"), by = "day", length.out = nrow(data))
  var_est <- estimate_var(data, lags = 2L, nimp_wold = 4L)

  legacy_angles <- c(0.3, -0.7, 1.1)
  legacy_rotation <- euler_rotation(legacy_angles, 3L)
  target_angles <- etvar:::.etvar_prepare_target_start(
    legacy_angles, n = 3L, shock_index = 2L, parameterization = "full"
  )
  target_rotation <- etvar:::.etvar_target_rotation(
    target_angles, n = 3L, shock_index = 2L
  )

  expect_length(target_angles, 2L)
  expect_equal(target_rotation[, 2L], legacy_rotation[, 2L], tolerance = 1e-12)
  expect_equal(crossprod(target_rotation), diag(3L), tolerance = 1e-12)
  expect_equal(det(target_rotation), 1, tolerance = 1e-12)
  expect_equal(
    etvar_objective(
      legacy_angles, var_est, target_obs = 50L, shock_index = 2L,
      parameterization = "full"
    ),
    etvar_objective(
      target_angles, var_est, target_obs = 50L, shock_index = 2L,
      parameterization = "target"
    ),
    tolerance = 1e-10
  )

  fit <- etvar(
    data, dates, target_date = dates[60L], lags = 2L,
    shock_index = 2L, start_angles = legacy_angles,
    start_parameterization = "full", optimize = FALSE, nimp = 4L,
    verbose = FALSE
  )
  expect_length(fit$theta, 2L)
  expect_identical(fit$parameterization, "target_direction")
  expect_equal(fit$target_direction, legacy_rotation[, 2L], tolerance = 1e-12)
  expect_equal(fit$rotation[, 2L], fit$target_direction, tolerance = 1e-12)
})

test_that("prepared objective matches the full historical decomposition", {
  set.seed(406)
  data <- matrix(stats::rnorm(270L), ncol = 3L)
  target_obs <- c(40L, 70L)
  target_angles <- c(0.35, -0.6)

  slow_objective <- function(theta, var_est, target_obs, shock_index) {
    Q <- etvar:::.etvar_target_rotation(theta, var_est$n, shock_index)
    A0 <- -t(chol(var_est$SIGMA)) %*% Q
    decomp <- hist_decomposition(
      A0, var_est$ALPHA, var_est$Y, var_est$X, var_est$n,
      var_est$lags, var_est$trend, var_est$constant
    )
    counterfactual <- decomp$deterministic +
      decomp$hist_matrix[, , shock_index]
    counterfactual_diff <- diff(counterfactual)
    real_diff <- diff(var_est$data)
    scales <- apply(real_diff, 2L, stats::sd)
    indices <- target_obs + var_est$lags - 1L
    sum(abs(
      sweep(counterfactual_diff[indices, , drop = FALSE], 2L, scales, "/") -
        sweep(real_diff[indices, , drop = FALSE], 2L, scales, "/")
    ))
  }

  specifications <- expand.grid(
    constant = c(FALSE, TRUE),
    trend = 0:2,
    KEEP.OUT.ATTRS = FALSE
  )

  for (i in seq_len(nrow(specifications))) {
    var_est <- estimate_var(
      data,
      lags = 2L,
      constant = specifications$constant[i],
      trend = specifications$trend[i],
      nimp_wold = 3L
    )
    context <- etvar:::.etvar_prepare_objective(
      var_est, target_obs, shock_index = 2L
    )
    prepared <- etvar:::.etvar_objective_prepared(
      target_angles, context, parameterization = "target"
    )
    public <- etvar_objective(
      target_angles, var_est, target_obs, shock_index = 2L,
      parameterization = "target"
    )
    reference <- slow_objective(
      target_angles, var_est, target_obs, shock_index = 2L
    )

    expect_s3_class(context, "etvar_objective_context")
    expect_equal(prepared, reference, tolerance = 1e-9)
    expect_equal(public, reference, tolerance = 1e-9)
  }
})

test_that("optimiser returns n - 1 parameters and a valid completion", {
  set.seed(405)
  data <- matrix(stats::rnorm(150L), ncol = 3L)
  var_est <- estimate_var(data, lags = 1L, nimp_wold = 3L)

  result <- suppressWarnings(etvar_optimize(
    start = c(0.2, -0.3, 0.5),
    start_parameterization = "full",
    var_est = var_est,
    target_obs = 35L,
    shock_index = 2L,
    control_nm = list(reltol = 1e-4, maxit = 1L, trace = 0),
    control_bfgs = list(reltol = 1e-4, maxit = 1L, trace = 0, REPORT = 1L),
    verbose = FALSE
  ))

  expect_length(result$par, 2L)
  expect_length(result$target_direction, 3L)
  expect_equal(result$rotation[, 2L], result$target_direction,
               tolerance = 1e-12)
  expect_equal(crossprod(result$rotation), diag(3L), tolerance = 1e-12)
  expect_equal(det(result$rotation), 1, tolerance = 1e-12)
  expect_identical(result$parameterization, "target_direction")
  expect_named(result$optimization, c(
    "convergence", "message", "selected_stage", "selected_start",
    "objective_by_start_stage", "n_evaluations", "n_gradient_evaluations",
    "elapsed_seconds", "controls", "seed", "n_starts", "start_values"
  ))
  expect_s3_class(result$optimization$objective_by_start_stage, "data.frame")
  expect_setequal(
    names(result$optimization$objective_by_start_stage),
    c("start", "stage", "method", "objective", "convergence", "message",
      "fn_evaluations", "gr_evaluations", "elapsed_seconds")
  )
  expect_equal(nrow(result$optimization$objective_by_start_stage), 3L)
  expect_identical(
    result$optimization$selected_stage,
    result$stage
  )
  expect_identical(
    result$optimization$selected_start,
    result$start
  )
  expect_true(result$optimization$n_evaluations >= 1L)
  expect_true(is.numeric(result$optimization$elapsed_seconds))
  expect_equal(dim(result$optimization$start_values), c(1L, 2L))
})

test_that("etvar retains optimisation diagnostics and warns on nonzero convergence", {
  set.seed(407)
  data <- matrix(stats::rnorm(165L), ncol = 3L)
  dates <- seq.Date(as.Date("2021-01-01"), by = "day", length.out = nrow(data))

  fit <- NULL
  expect_warning(
    fit <- etvar(
      data, dates, target_date = dates[45L], lags = 2L,
      control_nm = list(reltol = 1e-4, maxit = 1L, trace = 0),
      control_bfgs = list(reltol = 1e-4, maxit = 1L, trace = 0,
                          REPORT = 1L),
      nimp = 4L, verbose = FALSE
    ),
    "nonzero convergence code"
  )

  expect_s3_class(fit, "etvar")
  selected_row <- fit$optimization$objective_by_start_stage$start ==
    fit$optimization$selected_start &
    fit$optimization$objective_by_start_stage$stage ==
      fit$optimization$selected_stage
  expect_true(any(selected_row))
  expect_equal(
    fit$optimization$objective_by_start_stage$objective[selected_row][1L],
    fit$objective
  )
  expect_identical(
    fit$optimization$convergence,
    fit$optimization$objective_by_start_stage$convergence[selected_row][1L]
  )
  expect_false(identical(fit$optimization$convergence, 0L))
  expect_true(nchar(fit$optimization$message) > 0L)
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
  irf_plot <- plot_etvar_irfs(
    fit$irfs, fit$variable_names, variables = 2L,
    irf_ci = intervals, facet_ncol = 1L
  )
  counterfactual_plot <- plot_etvar_counterfactual(
    dates, data, fit$hd$all_data[, , 1L],
    variable_names = fit$variable_names
  )
  hd_plot <- plot_etvar_hd(
    dates, data, fit$hd$deterministic,
    fit$shock_contributions$target_stop,
    fit$shock_contributions$target_rest,
    fit$shock_contributions$others, fit$datestar,
    fit$variable_names, ask = FALSE, variables = 2L
  )
  error_plot <- plot_etvar_errors(
    dates[-seq_len(fit$var$lags)], fit$hd$errors_str[, 1L]
  )

  plots <- list(irf_plot, counterfactual_plot, hd_plot, error_plot)
  for (pl in plots) {
    expect_s3_class(pl, "ggplot")
    expect_no_error(ggplot2::ggplot_build(pl))
  }
  expect_identical(levels(irf_plot$data$Variable), "y2")
  expect_identical(levels(counterfactual_plot$data$Variable), c("y1", "y2"))
  expect_identical(levels(hd_plot$data$Variable), "y2")
  expect_identical(levels(error_plot$data$Shock), "Target shock")
})

test_that("historical decomposition stacks the event shock from zero", {
  skip_if_not_installed("ggplot2")

  dates <- as.Date(c("2025-04-03", "2025-04-04"))
  zero <- matrix(0, nrow = 2L, ncol = 1L)
  event <- matrix(c(1, -1), ncol = 1L)
  continued <- matrix(c(2, -2), ncol = 1L)
  others <- matrix(c(3, -3), ncol = 1L)

  plot <- plot_etvar_hd(
    date_vec = dates,
    data = event + continued + others,
    deterministic = zero,
    target_stop = event,
    target_rest = continued,
    others = others,
    datestar = 1L,
    variable_names = "y1",
    n_tail = 2L
  )
  bars <- ggplot2::ggplot_build(plot)$data[[1L]]

  event_bars <- bars[bars$fill == "#D95319", , drop = FALSE]
  continued_bars <- bars[bars$fill == "#EDB120", , drop = FALSE]
  other_bars <- bars[bars$fill == "#4DBEEE", , drop = FALSE]

  expect_equal(event_bars$ymin, c(0, -1))
  expect_equal(event_bars$ymax, c(1, 0))
  expect_equal(continued_bars$ymin, c(1, -3))
  expect_equal(continued_bars$ymax, c(3, -1))
  expect_equal(other_bars$ymin, c(3, -6))
  expect_equal(other_bars$ymax, c(6, -3))
})

test_that("plot.etvar uses consistent faceted selection interfaces", {
  skip_if_not_installed("ggplot2")
  data <- make_core_test_data()
  dates <- seq.Date(as.Date("2021-01-01"), by = "day", length.out = nrow(data))
  fit <- etvar(
    data, dates, target_date = dates[70L], lags = 2L,
    optimize = FALSE, nimp = 4L, verbose = FALSE
  )

  for (type in c("irfs", "fevd", "counterfactual", "hd")) {
    pl <- plot(fit, type = type, variables = 2L, facet_ncol = 1L)
    expect_s3_class(pl, "ggplot")
    expect_no_error(ggplot2::ggplot_build(pl))
    expect_identical(levels(pl$data$Variable), "y2")
  }

  legacy_plot <- plot(fit, type = "counterfactual", variable = 1L)
  expect_s3_class(legacy_plot, "ggplot")
  expect_identical(levels(legacy_plot$data$Variable), "y1")

  target_error <- plot(fit, type = "errors")
  expect_s3_class(target_error, "ggplot")
  expect_no_error(ggplot2::ggplot_build(target_error))
  expect_identical(
    levels(target_error$data$Shock),
    "Target shock (1)"
  )
  expect_identical(target_error$labels$title, "Structural error")
  expect_identical(target_error$layers[[2L]]$aes_params$colour, "black")
  legacy_error <- plot(fit, type = "errors", shock = 1L)
  expect_identical(levels(legacy_error$data$Shock), "Target shock (1)")

  expect_error(
    plot(fit, type = "hd", variable = 1L, variables = 2L),
    "only one"
  )
})
