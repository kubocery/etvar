make_etvar_inputs <- function() {
  set.seed(711)
  data <- matrix(stats::rnorm(180L), ncol = 2L,
                 dimnames = list(NULL, c("y1", "y2")))
  dates <- seq(as.Date("2021-01-01"), by = "day", length.out = nrow(data))
  list(data = data, dates = dates)
}

test_that("etvar() rejects inconsistent event dates", {
  inputs <- make_etvar_inputs()
  fit_on <- function(...) {
    etvar(inputs$data, inputs$dates, lags = 2L, optimize = FALSE, ...)
  }

  expect_error(fit_on(target_date = "1999-01-01"), "not found in `dates`")
  expect_error(fit_on(target_date = inputs$dates[c(50L, 50L)]), "duplicate")
  expect_error(fit_on(target_date = inputs$dates[2L]), "after the first")
  expect_error(
    etvar(inputs$data, rev(inputs$dates), inputs$dates[50L], lags = 2L),
    "strictly increasing"
  )
  expect_error(fit_on(target_date = inputs$dates[50L],
                      variable_names = "only one"),
               "one non-empty name per variable")
})

test_that("event windows are sorted and anchored at the first date", {
  inputs <- make_etvar_inputs()
  fit <- etvar(inputs$data, inputs$dates,
               target_date = inputs$dates[c(70L, 55L)], lags = 2L, nimp = 3L)

  expect_identical(fit$target_obs, c(55L, 70L))
  expect_identical(fit$target_date, inputs$dates[55L])
  expect_identical(fit$datestar, nrow(inputs$data) - 55L)
})

test_that("print() and summary() describe the fit", {
  inputs <- make_etvar_inputs()
  fit <- etvar(inputs$data, inputs$dates, target_date = inputs$dates[60L],
               lags = 2L, nimp = 3L)

  expect_output(print(fit), "Event date: +2021-03-01")
  expect_output(summary(fit), "Largest companion eigenvalue")
  expect_output(summary(fit), "2\\. y2")
  expect_output(printed <- withVisible(print(fit)))
  expect_false(printed$visible)
  expect_identical(printed$value, fit)
})

test_that("a single series is handled as a one-shock model", {
  inputs <- make_etvar_inputs()
  fit <- etvar(inputs$data[, 1L, drop = FALSE], inputs$dates,
               target_date = inputs$dates[60L], lags = 2L, nimp = 4L)

  expect_length(fit$theta, 0L)
  expect_identical(dim(fit$irfs), c(1L, 1L, 4L))
  expect_equal(unname(compute_fevd(fit)$target), matrix(1, 1L, 4L))
})

test_that("FEVD plots label a non-target shock as such", {
  skip_if_not_installed("ggplot2")
  inputs <- make_etvar_inputs()
  fit <- etvar(inputs$data, inputs$dates, target_date = inputs$dates[60L],
               lags = 2L, nimp = 3L)

  expect_identical(levels(plot(fit, type = "fevd")$data$Component),
                   c("Target shock (1)", "Other shocks"))
  expect_identical(levels(plot(fit, type = "fevd", shock = 2L)$data$Component),
                   c("Shock 2", "Other shocks"))
})
