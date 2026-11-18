#' Estimate an event-targeted VAR
#'
#' Fits an ETVAR model in one call: estimates the reduced-form VAR, finds the
#' structural shock that best explains how all variables moved on the event
#' date(s), and computes the impulse responses, the historical decomposition
#' and everything the plotting methods need.
#'
#' The structural shocks are a rotation of the Cholesky shocks, `A0 = -L Q`
#' with `L` the lower Cholesky factor of the residual covariance and `Q` an
#' orthogonal matrix. The event criterion pins down only the target column of
#' `Q`, a unit vector described by `n - 1` angles, so only that column is
#' estimated. The other columns complete it to a rotation in a fixed way; the
#' individual non-target shocks are therefore not identified, although their
#' combined contribution is.
#'
#' A direction and its negative fit the event equally well. When the
#' direction is optimised, its sign is chosen so that the target shock is
#' positive on the first event date. The impulse responses describe a
#' positive one-standard-deviation shock, so on impact they move the same way
#' as the data did on the event date. Angles passed with `optimize = FALSE`
#' are used as given.
#'
#' @param data Matrix or data frame (T x n) of numeric series, one column per
#'   variable.
#' @param dates Vector of T dates (anything [as.Date()] accepts), strictly
#'   increasing.
#' @param target_date One or more event dates, each of which must appear in
#'   `dates`. With several dates the deviations are summed over all of them
#'   (a target window).
#' @param lags Integer. Number of VAR lags. Default 8.
#' @param constant Logical. Include an intercept? Default `TRUE`.
#' @param trend Integer. Deterministic trend: 0 = none, 1 = linear,
#'   2 = quadratic. Default 0.
#' @param shock_index Integer. Which structural shock is the target shock.
#'   Default 1.
#' @param start_angles Numeric vector or `NULL`. Starting angles for the
#'   optimiser, or the final angles when `optimize = FALSE`: either the `n - 1`
#'   angles of the target direction or the `n * (n - 1) / 2` angles of a full
#'   rotation as in [euler_rotation()]. `NULL` starts from the Cholesky
#'   direction.
#' @param optimize Logical. Optimise the direction? Default `TRUE`. With
#'   `FALSE`, the model is built from `start_angles` (or the Cholesky
#'   direction) as they are.
#' @param n_starts Integer. Number of starting points when
#'   `start_angles = NULL`; the Cholesky direction is always one of them.
#'   More starts make it less likely that the optimiser stops in a local
#'   minimum. Default 1.
#' @param seed Integer or `NULL`. Seed for the random starting points.
#' @param nimp Integer. Number of impulse-response horizons, including horizon
#'   zero. Default 30.
#' @param variable_names Character vector of n names. Defaults to the column
#'   names of `data`, or `V1`, `V2`, ... .
#' @param control_nm,control_bfgs Lists of [optim()] controls for the
#'   Nelder-Mead and BFGS stages of [etvar_optimize()].
#' @param n_cores Integer. Number of cores for running several starting
#'   points in parallel. Default 1.
#' @param verbose Logical. Print progress messages? Default `FALSE`.
#' @param start_date,end_date Dates or `NULL`. Restrict the sample used to
#'   estimate the VAR coefficients (inclusive bounds). The fitted VAR is still
#'   applied to the whole sample, so residuals, the decomposition and the
#'   event dates may lie outside the window. `NULL` means no bound.
#' @param start_parameterization Character. How to read `start_angles`:
#'   `"target"`, `"full"`, or `"auto"` to decide from its length (with two
#'   variables `"auto"` reads it as `"target"`).
#' @return An object of class `"etvar"`, a list with
#'   \describe{
#'     \item{var}{The reduced-form VAR, see [estimate_var()].}
#'     \item{theta}{The `n - 1` angles of the target direction.}
#'     \item{target_direction}{The target column of `rotation`, a unit
#'       vector.}
#'     \item{rotation}{The rotation `Q` (n x n).}
#'     \item{A0}{Structural impact matrix `-L Q` (n x n).}
#'     \item{irfs}{Structural impulse responses (n x n x `nimp`); see
#'       [compute_structural_irfs()].}
#'     \item{hd}{Historical decomposition from [hist_decomposition_3d()], with
#'       the intervention path starting after the first event date.}
#'     \item{shock_contributions}{Contributions to the data (T x n matrices):
#'       those of [rearrange_shocks()] (`total`, `target`, `target_stop` =
#'       target shocks up to and including the first event date,
#'       `target_rest` = target shocks after it, `others`), plus
#'       `target_event`, the contribution of the target shock on the first
#'       event date alone, and `target_before`, that of the target shocks
#'       before it.}
#'     \item{dates, data, variable_names}{The input data.}
#'     \item{target_dates}{All event dates, sorted.}
#'     \item{target_date}{The first event date.}
#'     \item{target_obs}{Rows of the event dates in `data`.}
#'     \item{datestar}{Number of observations after the first event date.}
#'     \item{start_date, end_date}{Bounds of the VAR estimation sample.}
#'     \item{estimation_dates}{Dates actually used to estimate the VAR.}
#'     \item{shock_index}{The target shock.}
#'     \item{objective}{The minimised objective, see [etvar_objective()].}
#'     \item{optimization}{Optimiser diagnostics, see [etvar_optimize()].}
#'   }
#' @references
#' Pinter, G., Smets, F. and Üslü, S. (2025). Market whiplash after the
#' 2025 tariff shock: an event-targeted VAR approach. *BIS Working Papers*,
#' No. 1282. Bank for International Settlements.
#' <https://www.bis.org/publications/working-paper-1282-market-whiplash-after-2025-tariff-shock-event-targeted-var-approach>
#' @seealso [plot.etvar()], [compute_fevd()], [compute_irf_confint()]
#' @examples
#' set.seed(42)
#' y <- matrix(rnorm(160), ncol = 2)
#' dates <- seq(as.Date("2020-01-01"), by = "day", length.out = 80)
#' fit <- etvar(y, dates, target_date = dates[60], lags = 2)
#' fit
#' @export
etvar <- function(data, dates, target_date, lags = 8L,
                  constant = TRUE, trend = 0L, shock_index = 1L,
                  start_angles = NULL, optimize = TRUE,
                  n_starts = 1L, seed = NULL,
                  nimp = 30L, variable_names = NULL,
                  control_nm = list(reltol = 1e-9, maxit = 20000, trace = 0),
                  control_bfgs = list(reltol = 1e-9, maxit = 4000,
                                      trace = 0, REPORT = 1),
                  n_cores = 1L,
                  verbose = FALSE,
                  start_date = NULL, end_date = NULL,
                  start_parameterization = c("auto", "target", "full")) {
  data <- .etvar_numeric_matrix(data)
  dates <- .etvar_as_date(dates, "dates")
  target_dates <- .etvar_as_date(target_date, "target_date")
  lags <- .etvar_integer_scalar(lags, "lags", min = 1L)
  constant <- .etvar_flag(constant, "constant")
  trend <- .etvar_integer_scalar(trend, "trend", min = 0L, max = 2L)
  shock_index <- .etvar_integer_scalar(shock_index, "shock_index", min = 1L,
                                       max = ncol(data))
  start_parameterization <- match.arg(start_parameterization)
  optimize <- .etvar_flag(optimize, "optimize")
  n_starts <- .etvar_integer_scalar(n_starts, "n_starts", min = 1L)
  nimp <- .etvar_integer_scalar(nimp, "nimp", min = 1L)
  n_cores <- .etvar_integer_scalar(n_cores, "n_cores", min = 1L)
  verbose <- .etvar_flag(verbose, "verbose")
  if (!is.list(control_nm) || !is.list(control_bfgs)) {
    stop("`control_nm` and `control_bfgs` must be lists.", call. = FALSE)
  }
  if (!is.null(seed)) {
    seed <- .etvar_integer_scalar(seed, "seed", min = 0L,
                                  max = .Machine$integer.max)
  }
  n <- ncol(data)

  start <- NULL
  if (!is.null(start_angles)) {
    start <- .etvar_prepare_target_start(start_angles, n, shock_index,
                                         start_parameterization)
  }

  if (is.null(variable_names)) {
    variable_names <- colnames(data)
    if (is.null(variable_names)) variable_names <- paste0("V", seq_len(n))
  }
  .etvar_check_labels(variable_names, n, "variable_names")

  if (length(dates) != nrow(data)) {
    stop("`dates` must have the same length as the number of rows in `data`.",
         call. = FALSE)
  }
  if (is.unsorted(dates, strictly = TRUE)) {
    stop("`dates` must be strictly increasing.", call. = FALSE)
  }
  if (nrow(data) <= lags) {
    stop("`data` must contain more rows than `lags`.", call. = FALSE)
  }
  if (anyDuplicated(target_dates)) {
    stop("`target_date` must not contain duplicate dates.", call. = FALSE)
  }

  target_obs <- match(target_dates, dates)
  if (anyNA(target_obs)) {
    stop("date(s) not found in `dates`: ",
         paste(format(target_dates[is.na(target_obs)], "%Y-%m-%d"),
               collapse = ", "),
         call. = FALSE)
  }
  target_obs <- sort(target_obs)
  target_dates <- dates[target_obs]
  if (any(target_obs <= lags)) {
    stop("Target dates must occur after the first `lags` observations.",
         call. = FALSE)
  }

  # Dates of the effective sample (rows with a full set of lags).
  effective_dates <- dates[-seq_len(lags)]
  start_date <- date_bound(start_date, "start_date", effective_dates[1L])
  end_date <- date_bound(end_date, "end_date",
                         effective_dates[length(effective_dates)])
  if (start_date > end_date) {
    stop("`start_date` must be on or before `end_date`.", call. = FALSE)
  }
  estimation_rows <- which(effective_dates >= start_date &
                             effective_dates <= end_date)
  if (length(estimation_rows) == 0L) {
    stop("No effective-sample observations fall between `start_date` and ",
         "`end_date`.", call. = FALSE)
  }
  estimation_dates <- effective_dates[estimation_rows]

  if (verbose) {
    message("Estimating VAR(", lags, ") from ",
            format(estimation_dates[1L], "%Y-%m-%d"), " to ",
            format(estimation_dates[length(estimation_dates)], "%Y-%m-%d"),
            " (", length(estimation_dates), " observations) ...")
  }
  var_est <- estimate_var(data, lags = lags, constant = constant,
                          trend = trend, nimp_wold = nimp,
                          estimation_rows = estimation_rows)
  var_est$estimation_dates <- estimation_dates
  var_est$estimation_start_date <- estimation_dates[1L]
  var_est$estimation_end_date <- estimation_dates[length(estimation_dates)]

  target_eff <- target_obs - lags
  if (optimize) {
    if (verbose) message("Optimising target-shock direction ...")
    opt <- etvar_optimize(
      start = start,
      var_est = var_est,
      target_obs = target_eff,
      shock_index = shock_index,
      start_parameterization = "target",
      n_starts = n_starts,
      seed = seed,
      control_nm = control_nm,
      control_bfgs = control_bfgs,
      n_cores = n_cores,
      verbose = verbose
    )
  } else {
    if (is.null(start)) start <- cholesky_angles(n, shock_index)
    opt <- fixed_direction(start, var_est, target_eff, shock_index,
                           control_nm, control_bfgs, n_cores, seed)
    if (verbose) message("Using the supplied angles. Objective: ", opt$value)
  }

  # The minus sign follows the original implementation. It only fixes the
  # sign convention of the structural errors; the impulse responses below
  # are always responses to a positive shock.
  A0 <- -t(chol(var_est$SIGMA)) %*% opt$rotation

  if (verbose) message("Computing impulse responses ...")
  irfs <- compute_structural_irfs(var_est, A0, nimp)

  if (verbose) message("Running the historical decomposition ...")
  datestar <- length(dates) - target_obs[1L]
  hd <- hist_decomposition_3d(
    A0 = A0, ALPHA = var_est$ALPHA, Y = var_est$Y, X = var_est$X,
    n = n, lags = lags, trend = trend, constant = as.integer(constant),
    datestar = datestar
  )
  shock_contributions <- rearrange_shocks(
    hd$hist_matrix, hd$hist_matrix_stop, nrow(data), n, shock_index
  )
  # target_stop holds every target shock up to the event, including those
  # from years before it. Separate the event-day shock from that history.
  shock_contributions$target_event <- single_error_contribution(
    var_est, A0, hd$errors_str, shock_index, target_eff[1L]
  )
  shock_contributions$target_before <- shock_contributions$target_stop -
    shock_contributions$target_event

  structure(list(
    var = var_est,
    theta = opt$par,
    target_direction = drop(opt$rotation[, shock_index]),
    rotation = opt$rotation,
    A0 = A0,
    irfs = irfs,
    hd = hd,
    shock_contributions = shock_contributions,
    dates = dates,
    data = data,
    target_dates = target_dates,
    target_date = target_dates[1L],
    target_obs = target_obs,
    start_date = start_date,
    end_date = end_date,
    estimation_dates = estimation_dates,
    datestar = datestar,
    shock_index = shock_index,
    variable_names = variable_names,
    objective = opt$value,
    optimization = opt$optimization
  ), class = "etvar")
}


date_bound <- function(value, arg, default) {
  if (is.null(value)) {
    return(default)
  }
  parsed <- tryCatch(as.Date(value), error = function(e) as.Date(NA))
  if (length(parsed) != 1L || is.na(parsed)) {
    stop("`", arg, "` must be one non-missing date or `NULL`.", call. = FALSE)
  }
  parsed
}

# The optimize = FALSE counterpart of etvar_optimize(): evaluates the
# objective at the given angles and records that nothing was optimised.
fixed_direction <- function(theta, var_est, target_obs, shock_index,
                            control_nm, control_bfgs, n_cores, seed) {
  started <- proc.time()[["elapsed"]]
  direction <- .etvar_direction_from_angles(theta, var_est$n, shock_index)
  theta <- .etvar_angles_from_direction(direction, shock_index)
  value <- etvar_objective(theta, var_est, target_obs, shock_index,
                           parameterization = "target")
  record <- stage_record(
    NA_integer_, "not_optimized", NA_character_,
    list(value = value, convergence = 0L, message = "Optimization skipped.",
         counts = c(`function` = 1L, gradient = 0L)),
    proc.time()[["elapsed"]] - started
  )
  controls <- list(nelder_mead = control_nm, bfgs = control_bfgs,
                   n_cores = n_cores, parallel = FALSE)
  summary <- optimization_summary(record, record, controls, seed,
                                  list(theta), started)
  summary$n_starts <- 0L

  list(
    par = theta,
    rotation = .etvar_target_rotation(theta, var_est$n, shock_index),
    value = value,
    optimization = summary
  )
}


#' Print an etvar object
#'
#' Prints a short description of a fitted ETVAR model.
#'
#' @param x An `etvar` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3,
#'             dimnames = list(NULL, c("x1", "x2", "x3")))
#' dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 100)
#' fit <- etvar(y, dates, target_date = dates[80], lags = 2, nimp = 10)
#' print(fit)
#' @export
print.etvar <- function(x, ...) {
  cat("Event-Targeted VAR (ETVAR)\n")
  print_field("Variables", x$var$n)
  print_field("Lags", x$var$lags)
  if (length(x$target_dates) == 1L) {
    print_field("Event date", format(x$target_date))
  } else {
    print_field("Event window", paste(
      format(min(x$target_dates)), "to", format(max(x$target_dates)),
      paste0("(", length(x$target_dates), " dates)")
    ))
  }
  print_field("VAR sample", estimation_sample(x$var))
  print_field("Target shock", x$shock_index)
  print_field("Objective", format(round(x$objective, 6)))
  print_field("IRF horizons", dim(x$irfs)[3L])
  invisible(x)
}

print_field <- function(label, value) {
  if (length(value) == 0L) {
    return(invisible())
  }
  cat(sprintf("  %-15s %s\n", paste0(label, ":"), value))
}

estimation_sample <- function(var) {
  if (is.null(var$estimation_start_date)) {
    return(NULL)
  }
  paste(format(var$estimation_start_date), "to",
        format(var$estimation_end_date),
        paste0("(", var$nobs_estimation, " observations)"))
}


#' Summarise an etvar object
#'
#' Prints the model specification, the identification settings and the
#' largest companion-matrix eigenvalue (a stability check: it should be below
#' one) of a fitted ETVAR model.
#'
#' @param object An `etvar` object.
#' @param ... Ignored.
#' @return `object`, invisibly.
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3,
#'             dimnames = list(NULL, c("x1", "x2", "x3")))
#' dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 100)
#' fit <- etvar(y, dates, target_date = dates[80], lags = 2, nimp = 10)
#' summary(fit)
#' @export
summary.etvar <- function(object, ...) {
  var <- object$var
  cat("Event-Targeted VAR (ETVAR)\n\nModel\n")
  print_field("Variables", var$n)
  print_field("Lags", var$lags)
  print_field("Constant", if (var$constant == 1L) "yes" else "no")
  print_field("Trend", c("none", "linear", "quadratic")[var$trend + 1L])
  print_field("Observations", paste0(nrow(object$data), " (", nrow(var$Y),
                                     " with a full set of lags)"))
  print_field("VAR sample", estimation_sample(var))

  cat("\nIdentification\n")
  plural <- if (length(object$target_dates) > 1L) "s" else ""
  print_field(paste0("Event date", plural),
              paste(format(object$target_dates), collapse = ", "))
  print_field("Target shock", object$shock_index)
  print_field("Angles", length(object$theta))
  print_field("Objective", format(round(object$objective, 6)))

  cat("\nVariables\n")
  cat(paste0("  ", seq_along(object$variable_names), ". ",
             object$variable_names, "\n"), sep = "")

  max_eigenvalue <- max(Mod(eigen(var$companion, only.values = TRUE)$values))
  cat("\nLargest companion eigenvalue:", format(round(max_eigenvalue, 4)),
      "(below 1 for a stable VAR)\n")
  invisible(object)
}


#' Plot an etvar object
#'
#' One entry point for the ETVAR plots; `type` chooses the plot and the other
#' arguments are passed on to the matching `plot_etvar_*()` function.
#'
#' @param x An `etvar` object.
#' @param type Character. One of
#'   \describe{
#'     \item{`"irfs"`}{Impulse responses to the target shock (or `shock`),
#'       see [plot_etvar_irfs()].}
#'     \item{`"fevd"`}{Forecast error variance shares of the target shock
#'       and of all other shocks, see [plot_etvar_fevd()].}
#'     \item{`"hd"`}{Historical decomposition around the event, see
#'       [plot_etvar_hd()]: the target shock on the first event date, target
#'       shocks after it, and everything else (other shocks and target shocks
#'       before the event).}
#'     \item{`"counterfactual"`}{Observed data against the path implied by
#'       the target shock alone, see [plot_etvar_counterfactual()].}
#'     \item{`"errors"`}{The structural errors of the target shock, see
#'       [plot_etvar_errors()].}
#'   }
#' @param variables Integer vector or `NULL`. Variables to show (all plots
#'   except `"errors"`). `NULL` shows every variable.
#' @param shock Integer or `NULL`. Shock to show in the `"irfs"`, `"fevd"`,
#'   `"counterfactual"` and `"errors"` plots. Defaults to the target shock.
#' @param scaling_index,scaling_multiplier,common_ylim Passed to
#'   [plot_etvar_irfs()]: normalise the responses by one IRF element, rescale
#'   them, and share the y axis across panels.
#' @param event_dates,event_labels Passed to [plot_etvar_errors()]: vertical
#'   markers and their labels.
#' @param n_tail Integer. Number of trailing observations shown in the
#'   `"hd"`, `"counterfactual"` and `"errors"` plots. Default 40.
#' @param x_date_breaks Character. Spacing of the date axis labels in the
#'   `"counterfactual"` and `"errors"` plots, e.g. `"2 weeks"`. Default
#'   `"1 week"`.
#' @param x_break_every Integer. In the `"hd"` plot, label every
#'   `x_break_every`-th date. Default 2.
#' @param merge_target_components Logical. In the `"hd"` plot, show the target
#'   shock as one component instead of splitting it at the event date.
#'   Default `FALSE`.
#' @param facet_ncol Integer. Number of panel columns. Default 2.
#' @param title Character, `NULL` or `NA`. Plot title; `NULL` uses the default
#'   title and `NA` removes it.
#' @param conf_int Logical. For `"irfs"`, compute bootstrap confidence bands
#'   with [compute_irf_confint()] and draw them. Default `FALSE`.
#' @param conf_level,n_boot,ci_seed Confidence level, number of bootstrap
#'   draws and seed used when `conf_int = TRUE`.
#' @param irf_ci Optional result of [compute_irf_confint()]. When supplied,
#'   its bands are drawn (whatever `conf_int` says) and no bootstrap is run.
#' @param ci_col,ci_border Fill and border colour of the confidence bands.
#' @param ... Not used.
#' @return A [ggplot2::ggplot] object.
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3,
#'             dimnames = list(NULL, c("x1", "x2", "x3")))
#' dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 100)
#' fit <- etvar(y, dates, target_date = dates[80], lags = 2, nimp = 10)
#' plot(fit, type = "irfs")
#' plot(fit, type = "fevd")
#' plot(fit, type = "hd", n_tail = 30)
#' plot(fit, type = "counterfactual", variables = 1, n_tail = 30)
#' plot(fit, type = "errors", n_tail = 40)
#' @export
plot.etvar <- function(
    x,
    type = c("irfs", "fevd", "hd", "counterfactual", "errors"),
    variables = NULL,
    shock = NULL,
    scaling_index = NULL,
    scaling_multiplier = 1,
    common_ylim = FALSE,
    event_dates = NULL,
    event_labels = NULL,
    n_tail = 40L,
    x_date_breaks = "1 week",
    x_break_every = 2L,
    merge_target_components = FALSE,
    facet_ncol = 2L,
    title = NULL,
    conf_int = FALSE,
    conf_level = 0.95,
    n_boot = 500L,
    ci_seed = NULL,
    irf_ci = NULL,
    ci_col = grDevices::adjustcolor("gray70", alpha.f = 0.45),
    ci_border = NA,
    ...) {
  type <- match.arg(type)
  n <- x$var$n
  if (type == "errors" && !is.null(variables)) {
    stop("`variables` does not apply to structural-error plots; use `shock`.",
         call. = FALSE)
  }
  if (!is.null(variables)) {
    variables <- .etvar_index_vector(variables, "variables", n)
  }
  if (is.null(shock)) shock <- x$shock_index
  shock <- .etvar_integer_scalar(shock, "shock", min = 1L, max = n)
  conf_int <- .etvar_flag(conf_int, "conf_int")

  switch(type,
    irfs = {
      if (conf_int && is.null(irf_ci)) {
        irf_ci <- compute_irf_confint(x, level = conf_level, n_boot = n_boot,
                                      seed = ci_seed)
      }
      plot_etvar_irfs(
        irfs = x$irfs,
        variable_names = x$variable_names,
        shock_index = shock,
        scaling_index = scaling_index,
        scaling_multiplier = scaling_multiplier,
        common_ylim = common_ylim,
        variables = variables,
        title = title,
        irf_ci = irf_ci,
        ci_col = ci_col,
        ci_border = ci_border,
        facet_ncol = facet_ncol
      )
    },
    fevd = plot_etvar_fevd(
      compute_fevd(x, target_shock = shock),
      variables = variables,
      title = title,
      facet_ncol = facet_ncol
    ),
    hd = plot_etvar_hd(
      date_vec = x$dates,
      data = x$data,
      deterministic = x$hd$deterministic,
      target_event = x$shock_contributions$target_event,
      target_rest = x$shock_contributions$target_rest,
      others = x$shock_contributions$others +
        x$shock_contributions$target_before,
      datestar = x$datestar,
      variable_names = x$variable_names,
      variables = variables,
      n_tail = n_tail,
      x_break_every = x_break_every,
      merge_target_components = merge_target_components,
      title = title,
      facet_ncol = facet_ncol
    ),
    counterfactual = plot_etvar_counterfactual(
      date_vec = x$dates,
      realized = x$data,
      counterfactual = x$hd$all_data[, , shock],
      variable_names = x$variable_names,
      variables = variables,
      n_tail = n_tail,
      x_date_breaks = x_date_breaks,
      title = title,
      facet_ncol = facet_ncol
    ),
    errors = plot_etvar_errors(
      date_vec = x$dates[-seq_len(x$var$lags)],
      errors = x$hd$errors_str[, shock],
      n_tail = n_tail,
      event_dates = event_dates,
      event_labels = event_labels,
      x_date_breaks = x_date_breaks,
      title = title,
      panel_title = if (shock == x$shock_index) {
        paste0("Target shock (", shock, ")")
      } else {
        paste("Shock", shock)
      }
    )
  )
}
