# ==============================================================================
# Main user-facing function and S3 methods
# ==============================================================================

#' Estimate an Event-Targeted VAR model
#'
#' This is the main entry point for the package. It runs the full ETVAR
#' pipeline in a single call:
#'
#' 1. Estimate a reduced-form VAR.
#' 2. Optimise the target-shock direction to target the event date.
#' 3. Build the structural impact matrix A0.
#' 4. Compute impulse response functions.
#' 5. Perform historical decomposition (with intervention path).
#' 6. Rearrange shock contributions for easy plotting.
#'
#' The event-targeting criterion identifies one unit direction: the selected
#' target shock. Consequently, estimation uses `n - 1` parameters rather than
#' a full `n*(n-1)/2` rotation. The remaining columns of the returned rotation
#' are a deterministic orthogonal completion. Individual non-target shocks
#' should not be interpreted as separately identified; their aggregate
#' contribution is invariant to the chosen completion.
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
#' @param start_angles Numeric vector or `NULL`. Target-direction angles of
#'   length `n - 1` can be supplied to warm-start the optimisation. Legacy
#'   full Euler angles of length `n*(n-1)/2` are also accepted and converted
#'   to their target-shock direction. When `NULL`, the Cholesky target
#'   direction is used as the first start.
#' @param start_parameterization Character. Interpretation of `start_angles`:
#'   `"target"`, `"full"`, or `"auto"` to infer from its length. In the
#'   ambiguous two-variable case, `"auto"` selects `"target"`.
#' @param optimize       Logical. Run the three-stage optimisation? Default
#'   `TRUE`. Set to `FALSE` when supplying final `start_angles` directly.
#' @param n_starts       Integer. Total number of optimiser starting points.
#'   Only used when `start_angles = NULL`; the Cholesky target direction is
#'   always included. A higher number gives a better chance of finding the
#'   global minimum at the cost of runtime. Default 1.
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
#' @param verbose        Logical. Print progress messages (via [message()])?
#'   Default `FALSE`. When `TRUE`, sequential optimisation also displays the
#'   live [optim()] objective trace.
#' @return An S3 object of class `"etvar"` containing:
#'   \describe{
#'     \item{var}{`etvar_var` object. Reduced-form VAR estimates.}
#'     \item{theta}{Numeric vector of length `n - 1`. Optimal target-direction
#'       angles.}
#'     \item{target_direction}{Numeric vector of length n. Identified unit
#'       direction of the target shock.}
#'     \item{rotation}{Matrix (n x n). Deterministic orthogonal completion of
#'       `target_direction`.}
#'     \item{parameterization}{The string `"target_direction"`.}
#'     \item{A0}{Matrix (n x n). Structural impact matrix.}
#'     \item{irfs}{Array (n x n x nimp). Structural IRFs.}
#'     \item{hd}{Named list from [hist_decomposition_3d()].}
#'     \item{shock_contributions}{Named list of historical shock-contribution
#'       matrices from [rearrange_shocks()].}
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
#'     \item{optimization}{Named list with optimiser diagnostics: convergence
#'       code and message, selected start and stage, objective values by start
#'       and stage, evaluation counts, elapsed time, controls, random seed, and
#'       starting values.}
#'   }
#' @references
#' Pinter, G., Smets, F. and Üslü, S. (2025). Market whiplash after the
#' 2025 tariff shock: an event-targeted VAR approach. *BIS Working Papers*,
#' No. 1282. Bank for International Settlements.
#' <https://www.bis.org/publications/working-paper-1282-market-whiplash-after-2025-tariff-shock-event-targeted-var-approach>
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
                  verbose = FALSE,
                  start_date = NULL, end_date = NULL,
                  start_parameterization = c("auto", "target", "full")) {

  # --- Input validation ---
  data         <- .etvar_numeric_matrix(data)
  dates        <- .etvar_as_date(dates, "dates")
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

  target_start <- if (is.null(start_angles)) {
    NULL
  } else {
    .etvar_prepare_target_start(
      start_angles, n, shock_index, start_parameterization
    )
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

  # --- 2. Optimise the target-shock direction ---
  # target_obs are in raw-data indices; adjust to effective-sample indices
  target_eff <- target_obs - lags

  if (optimize) {
    if (verbose) message("Optimising target-shock direction ...")
    opt <- etvar_optimize(
      start        = target_start,
      var_est      = var_est,
      target_obs   = target_eff,
      shock_index  = shock_index,
      start_parameterization = "target",
      n_starts     = n_starts,
      seed         = seed,
      control_nm   = control_nm,
      control_bfgs = control_bfgs,
      n_cores      = n_cores,
      verbose      = verbose
    )
    theta_opt <- opt$par
    Q <- opt$rotation
    obj_value <- opt$value
    optimization <- opt$optimization
  } else {
    if (is.null(target_start)) {
      cholesky_direction <- numeric(n)
      cholesky_direction[shock_index] <- 1
      target_start <- .etvar_angles_from_direction(
        cholesky_direction, shock_index
      )
    }
    theta_opt <- target_start
    Q <- .etvar_target_rotation(theta_opt, n, shock_index)
    objective_start_time <- proc.time()[["elapsed"]]
    obj_value <- etvar_objective(
      theta_opt, var_est, target_eff, shock_index,
      parameterization = "target"
    )
    objective_elapsed <- proc.time()[["elapsed"]] - objective_start_time
    start_values <- matrix(theta_opt, nrow = 1L)
    colnames(start_values) <- paste0("angle", seq_len(ncol(start_values)))
    optimization <- list(
      convergence = 0L,
      message = "Optimization skipped; objective evaluated at supplied target-direction angles.",
      selected_stage = "not_optimized",
      selected_start = NA_integer_,
      objective_by_start_stage = data.frame(
        start = NA_integer_,
        stage = "not_optimized",
        method = NA_character_,
        objective = as.numeric(obj_value),
        convergence = 0L,
        message = "Optimization skipped.",
        fn_evaluations = 1L,
        gr_evaluations = 0L,
        elapsed_seconds = as.numeric(objective_elapsed),
        stringsAsFactors = FALSE
      ),
      n_evaluations = 1L,
      n_gradient_evaluations = 0L,
      elapsed_seconds = as.numeric(objective_elapsed),
      controls = list(nelder_mead = control_nm, bfgs = control_bfgs,
                      n_cores = n_cores, parallel = FALSE),
      seed = seed,
      n_starts = 0L,
      start_values = start_values
    )
    if (verbose) message("Using supplied angles. Objective: ", obj_value)
  }

  # --- 3. Structural impact matrix ---
  A0 <- -t(chol(var_est$SIGMA)) %*% Q

  # --- 4. Impulse response functions ---
  if (verbose) message("Computing IRFs ...")
  irfs <- compute_structural_irfs(var_est, A0, nimp)

  # --- 5. Historical decomposition ---
  if (verbose) message("Running historical decomposition ...")
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
  shock_contributions <- rearrange_shocks(
    hd$hist_matrix, hd$hist_matrix_stop, T_obs, n, shock_index
  )

  # --- Assemble output ---
  structure(list(
    var            = var_est,
    theta          = theta_opt,
    target_direction = drop(Q[, shock_index]),
    rotation       = Q,
    parameterization = "target_direction",
    A0             = A0,
    irfs           = irfs,
    hd             = hd,
    shock_contributions = shock_contributions,
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
    objective      = obj_value,
    optimization   = optimization
  ), class = "etvar")
}


#' Print an etvar object
#'
#' Prints a short description of a fitted ETVAR model.
#'
#' @param x An `etvar` object.
#' @param ... Ignored.
#' @return `x`, invisibly. Called for its side effect of printing to the
#'   console.
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
#'
#' Prints the model specification, identification settings, and stability
#' diagnostics of a fitted ETVAR model.
#'
#' @param object An `etvar` object.
#' @param ... Ignored.
#' @return `object`, invisibly. Called for its side effect of printing to
#'   the console.
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
  cat("  Target-direction parameters:", length(object$theta), "\n")
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
#'     \item{`"fevd"`}{Structural forecast error variance decomposition,
#'       grouping the selected target shock against all other shocks.}
#'     \item{`"hd"`}{Historical decomposition stacked bar charts.}
#'     \item{`"counterfactual"`}{Counterfactual vs realised paths.}
#'     \item{`"errors"`}{Target structural error time series.}
#'   }
#' @param variable Integer or `NULL`. Backward-compatible alias for selecting
#'   one variable. Prefer `variables`. If both are supplied, an error is
#'   raised.
#' @param variables Integer vector or `NULL`. Variables to display for IRF,
#'   FEVD, counterfactual, and historical-decomposition plots. `NULL` displays
#'   all variables.
#' @param shock Integer or `NULL`. Shock to display for IRF, FEVD,
#'   counterfactual, and error plots. The default is the target shock.
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
#' @param n_tail  Integer. Number of trailing observations for counterfactual,
#'   historical-decomposition, and error plots. Default 40.
#' @param x_date_breaks Character. Date break spacing for ggplot-based
#'   date axes (`"counterfactual"`, `"errors"`). Passed to
#'   [ggplot2::scale_x_date()] (e.g. `"1 week"`, `"2 weeks"`, `"1 month"`).
#'   Default `"1 week"`.
#' @param x_break_every Integer. For `type = "hd"`, show every
#'   `x_break_every`-th date label on the x-axis. Default 2.
#' @param merge_target_components Logical. For `type = "hd"`, merge
#'   the event-date target-shock component and continued target-shock
#'   component into one bar segment. Default `FALSE`.
#' @param facet_ncol Integer. Number of columns in faceted variable plots.
#'   Default 2.
#' @param title Character scalar, `NULL`, or `NA`. `NULL` uses the default plot
#'   title, `NA` removes it, and a character scalar supplies a custom title.
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
#' @return A [ggplot2::ggplot] object.
#' @export
plot.etvar <- function(
    x,
    type = c("irfs", "fevd", "hd", "counterfactual", "errors"),
    variable = NULL,
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
  if (!inherits(x, "etvar")) {
    stop("`x` must be an `etvar` object.", call. = FALSE)
  }
  type <- match.arg(type)
  if (!is.null(variable)) {
    variable <- .etvar_integer_scalar(variable, "variable",
                                       min = 1L, max = x$var$n)
    if (!is.null(variables)) {
      stop("Use only one of `variable` and `variables`.", call. = FALSE)
    }
    variables <- variable
  }
  if (identical(type, "errors") && !is.null(variables)) {
    stop("`variables` does not apply to structural-error plots; use `shock`.",
         call. = FALSE)
  }
  if (!is.null(variables)) {
    variables <- .etvar_index_vector(variables, "variables", x$var$n)
  }
  facet_ncol <- .etvar_integer_scalar(facet_ncol, "facet_ncol", min = 1L)
  conf_int <- .etvar_flag(conf_int, "conf_int")
  if (is.null(shock)) shock <- x$shock_index
  shock <- .etvar_integer_scalar(shock, "shock",
                                 min = 1L, max = x$var$n)

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
        variables      = variables,
        facet_ncol     = facet_ncol,
        title          = title,
        irf_ci         = irf_ci,
        ci_col         = ci_col,
        ci_border      = ci_border
      )
    },
    fevd = {
      fevd_result <- compute_fevd(x, target_shock = shock)
      plot_etvar_fevd(
        fevd = fevd_result,
        variables = variables,
        title = title,
        facet_ncol = facet_ncol
      )
    },
    hd = {
      plot_etvar_hd(
        date_vec       = x$dates,
        data           = x$data,
        deterministic  = x$hd$deterministic,
        target_stop    = x$shock_contributions$target_stop,
        target_rest    = x$shock_contributions$target_rest,
        others         = x$shock_contributions$others,
        datestar       = x$datestar,
        variable_names = x$variable_names,
        variables      = variables,
        n_tail         = n_tail,
        x_break_every  = x_break_every,
        merge_target_components = merge_target_components,
        facet_ncol     = facet_ncol,
        title          = title
      )
    },
    counterfactual = {
      plot_etvar_counterfactual(
        date_vec       = x$dates,
        realized       = x$data,
        counterfactual = x$hd$all_data[, , shock],
        variable_names = x$variable_names,
        variables      = variables,
        n_tail         = n_tail,
        x_date_breaks  = x_date_breaks,
        facet_ncol     = facet_ncol,
        title          = title
      )
    },
    errors = {
      err_dates <- x$dates[(x$var$lags + 1):length(x$dates)]
      panel_title <- if (shock == x$shock_index) {
        paste0("Target shock (", x$shock_index, ")")
      } else {
        paste("Shock", shock)
      }
      plot_etvar_errors(
        date_vec     = err_dates,
        errors       = x$hd$errors_str[, shock],
        panel_title  = panel_title,
        line_color   = "black",
        n_tail       = n_tail,
        event_dates  = event_dates,
        event_labels = event_labels,
        x_date_breaks = x_date_breaks,
        title        = title
      )
    }
  )
}
