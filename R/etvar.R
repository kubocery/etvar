# ==============================================================================
# Main user-facing function and S3 methods
# ==============================================================================

#' Estimate an Event-Targeted VAR model
#'
#' This is the main entry point for the package. It runs the full ETVAR
#' pipeline in a single call:
#'
#' 1. Estimate a reduced-form VAR.
#' 2. Optimise the rotation angles to target the event date.
#' 3. Build the structural impact matrix A0.
#' 4. Compute impulse response functions.
#' 5. Perform historical decomposition (with intervention path).
#' 6. Rearrange shock contributions for easy plotting.
#'
#' @param data           Matrix or data.frame (T x n). Time-series data.
#'   Each column is a variable. Must be numeric.
#' @param dates          Date vector of length T. Observation dates.
#' @param target_date    Date, character, or a **vector** of dates coercible
#'   with `as.Date()`. A single date targets one event; a vector of dates
#'   defines a target window where the objective is summed across all event
#'   dates simultaneously. All supplied dates must exist in `dates`.
#' @param start_date     Date, character, or `NULL`. Inclusive start of the
#'   sample used to estimate the initial reduced-form VAR. The fitted VAR is
#'   still evaluated over the full `data` sample. `NULL` uses the first
#'   available effective-sample date.
#' @param end_date       Date, character, or `NULL`. Inclusive end of the
#'   sample used to estimate the initial reduced-form VAR. The fitted VAR is
#'   still evaluated over the full `data` sample. `NULL` uses the last
#'   available date.
#' @param lags           Integer. Number of VAR lags. Default 8.
#' @param constant       Logical. Include an intercept? Default `TRUE`.
#' @param trend          Integer. Trend order (0, 1, or 2). Default 0.
#' @param shock_index    Integer. Which structural shock to target. Default 1.
#' @param start_angles   Numeric vector of length n*(n-1)/2, or `NULL`.
#'   Pre-computed angles can be supplied to warm-start the optimisation.
#'   When `NULL` and `n_starts = 1`, zeros are used (Cholesky starting point).
#'   When `NULL` and `n_starts > 1`, a random multi-start search is performed
#'   (recommended when no prior solution is available).
#' @param optimize       Logical. Run the three-stage optimisation? Default
#'   `TRUE`. Set to `FALSE` when supplying final `start_angles` directly.
#' @param n_starts       Integer. Total number of optimiser starting points.
#'   Only used when `start_angles = NULL`; the zero vector is always included.
#'   A higher number gives a better chance of finding the global minimum at the
#'   cost of runtime. Default 1 (single start from zeros).
#' @param seed           Integer or `NULL`. Random seed for reproducible
#'   multi-start draws. Default `NULL`.
#' @param nimp           Integer. Number of IRF horizons, including horizon
#'   zero. Default 30.
#' @param variable_names Character vector of length n. Names for each variable.
#'   If `NULL`, uses column names of `data` or `V1, V2, ...`.
#' @param control_nm     List. Nelder-Mead control parameters. See [optim()].
#' @param control_bfgs   List. BFGS control parameters. See [optim()].
#' @param n_cores        Integer. Number of CPU cores for parallel multi-start
#'   optimisation. Only used when `start_angles = NULL` and `n_starts > 1`.
#'   Default 1 (sequential); the maximum is 2.
#' @param verbose        Logical. Print progress messages? Default `TRUE`.
#' @return An S3 object of class `"etvar"` containing:
#'   \describe{
#'     \item{var}{`etvar_var` object. Reduced-form VAR estimates.}
#'     \item{theta}{Numeric vector. Optimal rotation angles.}
#'     \item{A0}{Matrix (n x n). Structural impact matrix.}
#'     \item{irfs}{Array (n x n x nimp). Structural IRFs.}
#'     \item{hd}{Named list from [hist_decomposition_3d()].}
#'     \item{shocks}{Named list from [rearrange_shocks()].}
#'     \item{dates}{Date vector.}
#'     \item{data}{Original data matrix.}
#'     \item{target_dates}{Date vector. All targeted event dates.}
#'     \item{target_date}{Date. First (earliest) targeted event date.}
#'     \item{target_obs}{Integer vector. Row indices of the event dates in
#'       the raw data.}
#'     \item{start_date}{Date. Requested/resolved lower bound of the initial
#'       VAR estimation sample.}
#'     \item{end_date}{Date. Requested/resolved upper bound of the initial
#'       VAR estimation sample.}
#'     \item{estimation_dates}{Date vector. Effective observations actually
#'       used to estimate the initial VAR.}
#'     \item{datestar}{Offset from end for the earliest event date.}
#'     \item{shock_index}{Which shock was targeted.}
#'     \item{variable_names}{Character vector of variable names.}
#'     \item{objective}{Final value of the objective function.}
#'   }
#' @export
#' @examples
#' set.seed(42)
#' example_data <- matrix(rnorm(160), ncol = 2)
#' example_dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 80)
#' fit <- etvar(example_data, example_dates, target_date = example_dates[60],
#'              lags = 2, optimize = FALSE, verbose = FALSE)
#' print(fit)
etvar <- function(data, dates, target_date, lags = 8L,
                  constant = TRUE, trend = 0L, shock_index = 1L,
                  start_angles = NULL, optimize = TRUE,
                  n_starts = 1L, seed = NULL,
                  nimp = 30L, variable_names = NULL,
                  control_nm   = list(reltol = 1e-9, maxit = 20000, trace = 0),
                  control_bfgs = list(reltol = 1e-9, maxit = 4000,
                                      trace = 0, REPORT = 1),
                  n_cores = 1L,
                  verbose = TRUE,
                  start_date = NULL, end_date = NULL) {

  # --- Input validation ---
  data         <- .etvar_numeric_matrix(data)
  dates        <- .etvar_as_date(dates, "dates")
  target_dates <- .etvar_as_date(target_date, "target_date")
  lags <- .etvar_integer_scalar(lags, "lags", min = 1L)
  constant <- .etvar_flag(constant, "constant")
  trend <- .etvar_integer_scalar(trend, "trend", min = 0L, max = 2L)
  shock_index <- .etvar_integer_scalar(shock_index, "shock_index", min = 1L,
                                        max = ncol(data))
  optimize <- .etvar_flag(optimize, "optimize")
  n_starts <- .etvar_integer_scalar(n_starts, "n_starts", min = 1L)
  nimp <- .etvar_integer_scalar(nimp, "nimp", min = 1L)
  n_cores <- .etvar_integer_scalar(n_cores, "n_cores", min = 1L, max = 2L)
  verbose <- .etvar_flag(verbose, "verbose")
  if (!is.list(control_nm) || !is.list(control_bfgs)) {
    stop("`control_nm` and `control_bfgs` must be lists.", call. = FALSE)
  }
  if (!is.null(seed)) {
    seed <- .etvar_integer_scalar(seed, "seed", min = 0L,
                                  max = .Machine$integer.max)
  }
  n <- ncol(data)

  n_angles <- n * (n - 1L) / 2L
  if (!is.null(start_angles) &&
      (!is.numeric(start_angles) || anyNA(start_angles) ||
       any(!is.finite(start_angles)) || length(start_angles) != n_angles)) {
    stop("`start_angles` must be a finite numeric vector of length ",
         n_angles, ".", call. = FALSE)
  }

  if (is.null(variable_names)) {
    variable_names <- if (!is.null(colnames(data))) {
      colnames(data)
    } else {
      paste0("V", seq_len(n))
    }
  }
  if (!is.character(variable_names) || length(variable_names) != n ||
      anyNA(variable_names) || any(!nzchar(variable_names))) {
    stop("`variable_names` must contain one non-empty name per data column.",
         call. = FALSE)
  }

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
  missing    <- target_dates[is.na(target_obs)]
  if (length(missing) > 0) {
    stop("date(s) not found in `dates`: ",
         paste(format(missing, "%Y-%m-%d"), collapse = ", "),
         call. = FALSE)
  }
  # Sort by date so datestar is anchored at the earliest event
  target_obs   <- sort(target_obs)
  target_dates <- dates[target_obs]
  if (any(target_obs <= lags)) {
    stop("Target dates must occur after the first `lags` observations.",
         call. = FALSE)
  }
  datestar     <- length(dates) - target_obs[1]

  effective_dates <- dates[(lags + 1L):length(dates)]
  parse_date_bound <- function(value, argument, default) {
    if (is.null(value)) return(default)

    parsed <- tryCatch(as.Date(value), error = function(e) as.Date(NA))
    if (length(parsed) != 1L || is.na(parsed)) {
      stop("`", argument, "` must be one non-missing date or `NULL`.",
           call. = FALSE)
    }
    parsed
  }

  start_date <- parse_date_bound(start_date, "start_date",
                                 effective_dates[1L])
  end_date <- parse_date_bound(end_date, "end_date",
                               effective_dates[length(effective_dates)])
  if (start_date > end_date) {
    stop("`start_date` must be on or before `end_date`.", call. = FALSE)
  }

  estimation_rows <- which(
    effective_dates >= start_date & effective_dates <= end_date
  )
  if (length(estimation_rows) == 0L) {
    stop("No effective-sample observations fall between `start_date` and ",
         "`end_date`.", call. = FALSE)
  }
  estimation_dates <- effective_dates[estimation_rows]

  # --- 1. Estimate reduced-form VAR ---
  if (verbose) {
    cat("Estimating VAR(", lags, ") from ",
        format(estimation_dates[1L], "%Y-%m-%d"), " to ",
        format(estimation_dates[length(estimation_dates)], "%Y-%m-%d"),
        " (", length(estimation_dates), " observations) ...\n", sep = "")
  }
  var_est <- estimate_var(data, lags = lags, constant = constant,
                          trend = trend, nimp_wold = nimp,
                          estimation_rows = estimation_rows)
  var_est$estimation_dates <- estimation_dates
  var_est$estimation_start_date <- estimation_dates[1L]
  var_est$estimation_end_date <- estimation_dates[length(estimation_dates)]

  # --- 2. Optimise rotation angles ---
  # target_obs are in raw-data indices; adjust to effective-sample indices
  target_eff <- target_obs - lags

  if (optimize) {
    if (verbose) cat("Optimising rotation angles ...\n")
    opt <- etvar_optimize(
      start        = start_angles,
      var_est      = var_est,
      target_obs   = target_eff,
      shock_index  = shock_index,
      n_starts     = n_starts,
      seed         = seed,
      control_nm   = control_nm,
      control_bfgs = control_bfgs,
      n_cores      = n_cores,
      verbose      = verbose
    )
    theta_opt <- opt$par
    obj_value <- opt$value
  } else {
    if (is.null(start_angles)) start_angles <- rep(0, n * (n - 1) / 2)
    theta_opt <- start_angles
    obj_value <- etvar_objective(start_angles, var_est, target_eff, shock_index)
    if (verbose) cat("Using supplied angles. Objective:", obj_value, "\n")
  }

  # --- 3. Structural impact matrix ---
  Q  <- euler_rotation(theta_opt, n)
  A0 <- -t(chol(var_est$SIGMA)) %*% Q

  # --- 4. Impulse response functions ---
  if (verbose) cat("Computing IRFs ...\n")
  irfs <- compute_structural_irfs(var_est, A0, nimp)

  # --- 5. Historical decomposition ---
  if (verbose) cat("Running historical decomposition ...\n")
  hd <- hist_decomposition_3d(
    A0       = A0,
    ALPHA    = var_est$ALPHA,
    Y        = var_est$Y,
    X        = var_est$X,
    n        = n,
    lags     = lags,
    trend    = trend,
    constant = as.integer(constant),
    datestar = datestar
  )

  # --- 6. Rearrange shock contributions ---
  T_obs  <- nrow(data)
  shocks <- rearrange_shocks(
    hd$hist_matrix, hd$hist_matrix_stop, T_obs, n, shock_index
  )

  # --- Assemble output ---
  structure(list(
    var            = var_est,
    theta          = theta_opt,
    A0             = A0,
    irfs           = irfs,
    hd             = hd,
    shocks         = shocks,
    dates          = dates,
    data           = data,
    target_dates   = target_dates,
    target_date    = target_dates[1],   # backward-compat alias
    target_obs     = target_obs,
    start_date     = start_date,
    end_date       = end_date,
    estimation_dates = estimation_dates,
    datestar       = datestar,
    shock_index    = shock_index,
    variable_names = variable_names,
    objective      = obj_value
  ), class = "etvar")
}


#' Print an etvar object
#' @param x An `etvar` object.
#' @param ... Ignored.
#' @export
print.etvar <- function(x, ...) {
  cat("Event-Targeted VAR (ETVAR)\n")
  cat("  Variables:    ", x$var$n, "\n")
  cat("  Lags:         ", x$var$lags, "\n")
  if (length(x$target_dates) == 1L) {
    cat("  Target date:  ", format(x$target_date, "%Y-%m-%d"), "\n")
  } else {
    cat("  Target window:", format(min(x$target_dates), "%Y-%m-%d"), "to",
        format(max(x$target_dates), "%Y-%m-%d"),
        paste0("(", length(x$target_dates), " dates)"), "\n")
  }
  if (!is.null(x$var$estimation_start_date)) {
    cat("  VAR sample:   ",
        format(x$var$estimation_start_date, "%Y-%m-%d"), "to",
        format(x$var$estimation_end_date, "%Y-%m-%d"),
        paste0("(", x$var$nobs_estimation, " observations)"), "\n")
  }
  cat("  Target shock: ", x$shock_index, "\n")
  cat("  Objective:    ", round(x$objective, 6), "\n")
  cat("  IRF horizons: ", dim(x$irfs)[3], "\n")
  invisible(x)
}


#' Summarise an etvar object
#' @param object An `etvar` object.
#' @param ... Ignored.
#' @export
summary.etvar <- function(object, ...) {
  cat("========================================\n")
  cat("Event-Targeted VAR - Summary\n")
  cat("========================================\n\n")

  cat("Model specification:\n")
  cat("  Variables (n):", object$var$n, "\n")
  cat("  Lags (p):     ", object$var$lags, "\n")
  cat("  Constant:     ", ifelse(object$var$constant, "Yes", "No"), "\n")
  trend_label <- c("None", "Linear", "Quadratic")[object$var$trend + 1L]
  cat("  Trend:        ", trend_label, "\n")
  cat("  T (raw):      ", nrow(object$data), "\n")
  cat("  T (effective):", nrow(object$var$Y), "\n")
  if (!is.null(object$var$estimation_start_date)) {
    cat("  VAR fit sample:",
        format(object$var$estimation_start_date, "%Y-%m-%d"), "to",
        format(object$var$estimation_end_date, "%Y-%m-%d"),
        paste0("(", object$var$nobs_estimation, " observations)"), "\n")
  }
  cat("\n")

  cat("Identification:\n")
  if (length(object$target_dates) == 1L) {
    cat("  Target date:  ", format(object$target_date, "%Y-%m-%d"), "\n")
  } else {
    cat("  Target window:", format(min(object$target_dates), "%Y-%m-%d"), "to",
        format(max(object$target_dates), "%Y-%m-%d"), "\n")
    cat("  Event dates:  ",
        paste(format(object$target_dates, "%Y-%m-%d"), collapse = ", "), "\n")
  }
  cat("  Target shock: ", object$shock_index, "\n")
  cat("  N angles:     ", length(object$theta), "\n")
  cat("  Objective:    ", round(object$objective, 6), "\n\n")

  cat("Variable names:\n")
  for (i in seq_along(object$variable_names)) {
    cat("  ", i, ". ", object$variable_names[i], "\n", sep = "")
  }

  max_eigenvalue <- max(abs(
    eigen(object$var$companion, only.values = TRUE)$values
  ))
  cat("\nMax companion eigenvalue:", round(max_eigenvalue, 4), "\n")

  invisible(object)
}


#' Plot an etvar object
#'
#' S3 plotting method that dispatches to specialised plotting functions based
#' on the `type` argument.
#'
#' @param x    An `etvar` object.
#' @param type Character. Type of plot:
#'   \describe{
#'     \item{`"irfs"`}{Structural impulse response functions.}
#'     \item{`"hd"`}{Historical decomposition stacked bar charts.}
#'     \item{`"counterfactual"`}{Counterfactual vs realised path (specify
#'       `variable` argument).}
#'     \item{`"errors"`}{Structural error time series (specify `shock`
#'       argument).}
#'   }
#' @param variable Integer. For `type = "counterfactual"`, which variable
#'   to display. Default 1.
#' @param shock    Integer. Shock to display for IRF, counterfactual, and error
#'   plots. Default: the target shock.
#' @param scaling_index Numeric vector for `type = "irfs"`. See
#'   [plot_etvar_irfs()].
#' @param scaling_multiplier Numeric scalar for `type = "irfs"`. Applied
#'   after normalisation; e.g., `0.5` scales a unit-normalised shock to 0.5.
#'   Default 1.
#' @param common_ylim Logical for `type = "irfs"`. If `TRUE`, use a shared
#'   y-axis across all IRF panels. If `FALSE`, each panel uses its own
#'   y-range. Default `FALSE`.
#' @param event_dates  Date vector for `type = "errors"`. Vertical event
#'   markers.
#' @param event_labels Character vector for `type = "errors"`. Event labels.
#' @param n_tail  Integer. Number of trailing observations for counterfactual
#'   and error plots. Default 40.
#' @param x_date_breaks Character. Date break spacing for ggplot-based
#'   date axes (`"counterfactual"`, `"errors"`). Passed to
#'   [ggplot2::scale_x_date()] (e.g. `"1 week"`, `"2 weeks"`, `"1 month"`).
#'   Default `"1 week"`.
#' @param x_break_every Integer. For `type = "hd"` (base-R), show every
#'   `x_break_every`-th date label on the x-axis. Default 2.
#' @param merge_target_components Logical. For `type = "hd"`, merge
#'   the event-date target-shock component and continued target-shock
#'   component into one bar segment. Default `FALSE`.
#' @param title Character, `NULL`, or `NA`.
#'   - `NULL`: use default per-plot title.
#'   - `NA`: remove title.
#'   - Character scalar: custom title.
#'   - Character vector: for multi-panel variable plots, custom titles by
#'     variable/panel. Named vectors are matched to variable names.
#' @param conf_int Logical. For `type = "irfs"`, draw bootstrap confidence
#'   intervals around the IRFs. Default `FALSE`.
#' @param conf_level Numeric. Confidence level used when `conf_int = TRUE`.
#'   Default `0.95`.
#' @param n_boot Integer. Number of residual-bootstrap draws used when
#'   `conf_int = TRUE`. Default `500`.
#' @param ci_seed Integer or `NULL`. Optional random seed for bootstrap
#'   confidence intervals.
#' @param irf_ci Optional pre-computed IRF confidence interval list from
#'   [compute_irf_confint()]. Supplying this avoids recomputing the bootstrap.
#' @param ci_col Fill colour for IRF confidence bands.
#' @param ci_border Border colour for IRF confidence bands. Default `NA`.
#' @param ...     Additional arguments (currently unused).
#' @return A ggplot object (for `"counterfactual"` and `"errors"`) or `NULL`
#'   invisibly (for base-R plots).
#' @export
plot.etvar <- function(x, type = c("irfs", "hd", "counterfactual", "errors"),
                       variable = 1L, shock = NULL,
                       scaling_index = NULL,
                       scaling_multiplier = 1,
                       common_ylim = FALSE,
                       event_dates = NULL, event_labels = NULL,
                       n_tail = 40L,
                       x_date_breaks = "1 week",
                       x_break_every = 2L,
                       merge_target_components = FALSE,
                       title = NULL,
                       conf_int = FALSE,
                       conf_level = 0.95,
                       n_boot = 500L,
                       ci_seed = NULL,
                       irf_ci = NULL,
                       ci_col = grDevices::adjustcolor("gray70",
                                                       alpha.f = 0.45),
                       ci_border = NA,
                       ...) {
  if (!inherits(x, "etvar")) {
    stop("`x` must be an `etvar` object.", call. = FALSE)
  }
  type <- match.arg(type)
  if (is.null(shock)) shock <- x$shock_index
  shock <- .etvar_integer_scalar(shock, "shock", min = 1L, max = x$var$n)
  variable <- .etvar_integer_scalar(variable, "variable",
                                     min = 1L, max = x$var$n)
  conf_int <- .etvar_flag(conf_int, "conf_int")

  switch(type,
    irfs = {
      if (isTRUE(conf_int) && is.null(irf_ci)) {
        irf_ci <- compute_irf_confint(
          object = x,
          level  = conf_level,
          n_boot = n_boot,
          seed   = ci_seed,
          nimp   = dim(x$irfs)[3]
        )
      }
      plot_etvar_irfs(
        irfs           = x$irfs,
        variable_names = x$variable_names,
        shock_index    = shock,
        scaling_index  = scaling_index,
        scaling_multiplier = scaling_multiplier,
        common_ylim    = common_ylim,
        title          = title,
        irf_ci         = irf_ci,
        ci_col         = ci_col,
        ci_border      = ci_border
      )
    },
    hd = {
      plot_etvar_hd(
        date_vec       = x$dates,
        data           = x$data,
        deterministic  = x$hd$deterministic,
        target_stop    = x$shocks$target_stop,
        target_rest    = x$shocks$target_rest,
        others         = x$shocks$others,
        datestar       = x$datestar,
        variable_names = x$variable_names,
        n_tail         = n_tail,
        x_break_every  = x_break_every,
        merge_target_components = merge_target_components,
        title          = title
      )
    },
    counterfactual = {
      cf_title <- etvar_select_panel_title(
        title          = title,
        default_titles = x$variable_names,
        variable_index = variable,
        panel_index    = 1L,
        n_panels       = 1L
      )
      plot_etvar_counterfactual(
        date_vec       = x$dates,
        realized       = x$data[, variable],
        counterfactual = x$hd$all_data[, variable, shock],
        n_tail         = n_tail,
        x_date_breaks  = x_date_breaks,
        title          = cf_title
      )
    },
    errors = {
      err_dates <- x$dates[(x$var$lags + 1):length(x$dates)]
      err_title <- etvar_select_scalar_title(
        title         = title,
        default_title = paste("Structural error - shock", shock)
      )
      plot_etvar_errors(
        date_vec     = err_dates,
        errors       = x$hd$errors_str[, shock],
        n_tail       = n_tail,
        event_dates  = event_dates,
        event_labels = event_labels,
        x_date_breaks = x_date_breaks,
        title        = err_title
      )
    }
  )
}
