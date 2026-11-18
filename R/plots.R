# Column names used inside aes(); declared so R CMD check does not flag them.
utils::globalVariables(c(
  "Value", "Series", "Structural_error", "dates",
  "Horizon", "Share", "Component", "Variable",
  "Response", "Lower", "Upper", "Shock"
))

.etvar_plot_theme <- function(legend_position = "bottom") {
  ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = legend_position,
      legend.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title.position = "plot"
    )
}

# Facet labels for the selected variables; make.unique() keeps panels apart
# when two variables share a name.
panel_names <- function(variable_names, variables) {
  make.unique(variable_names[variables], sep = " #")
}

rotated_dates <- function() {
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


#' Plot an ETVAR forecast error variance decomposition
#'
#' Stacked bars of the forecast error variance shares of the target shock and
#' of all other shocks combined, by horizon, for each selected variable. The
#' two shares add up to one.
#'
#' @param fevd An `"etvar_fevd"` object from [compute_fevd()].
#' @param variables Integer vector or `NULL`. Variables to show; `NULL` shows
#'   all of them.
#' @param title Character, `NULL` or `NA`. `NULL` gives a default title and
#'   `NA` removes it.
#' @param component_labels Character vector of length two, or `NULL`. Legend
#'   labels for the target shock and the other shocks.
#' @param colors Character vector of length two. Fill colours for the target
#'   shock and the other shocks.
#' @param facet_ncol Integer. Number of panel columns. Default 2.
#' @return A [ggplot2::ggplot] object.
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3,
#'             dimnames = list(NULL, c("x1", "x2", "x3")))
#' dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 100)
#' fit <- etvar(y, dates, target_date = dates[80], lags = 2, nimp = 10)
#' plot_etvar_fevd(compute_fevd(fit))
#' @export
plot_etvar_fevd <- function(fevd, variables = NULL, title = NULL,
                            component_labels = NULL,
                            colors = c("#D95319", "#4DBEEE"),
                            facet_ncol = 2L) {
  if (!inherits(fevd, "etvar_fevd")) {
    stop("`fevd` must be an `etvar_fevd` object from `compute_fevd()`.",
         call. = FALSE)
  }
  target <- .etvar_numeric_matrix(fevd$target, "fevd$target")
  others <- .etvar_numeric_matrix(fevd$others, "fevd$others")
  if (!identical(dim(target), dim(others))) {
    stop("`fevd$target` and `fevd$others` must have the same dimensions.",
         call. = FALSE)
  }

  n <- nrow(target)
  nimp <- ncol(target)
  if (is.null(variables)) variables <- seq_len(n)
  variables <- .etvar_index_vector(variables, "variables", n)
  facet_ncol <- .etvar_integer_scalar(facet_ncol, "facet_ncol", min = 1L)

  variable_names <- fevd$variable_names
  if (!is.character(variable_names) || length(variable_names) != n ||
      anyNA(variable_names) || any(!nzchar(variable_names))) {
    variable_names <- paste0("V", seq_len(n))
  }
  if (is.null(component_labels)) {
    first <- if (isFALSE(fevd$is_identified_shock)) {
      paste("Shock", fevd$target_shock)
    } else {
      paste0("Target shock (", fevd$target_shock, ")")
    }
    component_labels <- c(first, "Other shocks")
  }
  if (!is.character(component_labels) || length(component_labels) != 2L ||
      anyNA(component_labels) || any(!nzchar(component_labels)) ||
      anyDuplicated(component_labels)) {
    stop("`component_labels` must contain two distinct, non-empty labels.",
         call. = FALSE)
  }
  if (!is.character(colors) || length(colors) != 2L || anyNA(colors) ||
      any(!nzchar(colors))) {
    stop("`colors` must contain two non-empty colour values.", call. = FALSE)
  }

  labels <- panel_names(variable_names, variables)
  panel <- factor(rep(labels, each = nimp), levels = labels)
  horizon <- rep(seq_len(nimp), times = length(variables))
  plot_data <- rbind(
    data.frame(Horizon = horizon,
               Share = as.vector(t(target[variables, , drop = FALSE])),
               Component = component_labels[1L], Variable = panel),
    data.frame(Horizon = horizon,
               Share = as.vector(t(others[variables, , drop = FALSE])),
               Component = component_labels[2L], Variable = panel)
  )
  plot_data$Component <- factor(plot_data$Component, levels = component_labels)

  ggplot2::ggplot(plot_data,
                  ggplot2::aes(x = Horizon, y = Share, fill = Component)) +
    ggplot2::geom_col(width = 0.9, color = "white", linewidth = 0.1) +
    ggplot2::facet_wrap(~Variable, ncol = facet_ncol) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(colors, component_labels),
      breaks = component_labels
    ) +
    ggplot2::scale_x_continuous(
      breaks = function(limits) unique(pmax(1, round(pretty(limits))))
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, by = 0.25),
      labels = function(x) paste0(round(100 * x), "%"),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Forecast horizon",
      y = "Share of forecast error variance",
      fill = NULL,
      title = .etvar_title(title,
                           "Structural forecast error variance decomposition")
    ) +
    .etvar_plot_theme()
}


#' Plot observed data against a shock-only path
#'
#' Compares the observed data with a model-implied path, typically the path
#' the model implies when only the target shock operates
#' (`fit$hd$all_data[, , fit$shock_index]`).
#'
#' @param date_vec Vector of dates, one per row of `realized`.
#' @param realized Numeric vector or matrix of observed data, variables in
#'   columns.
#' @param counterfactual Numeric vector or matrix of the same shape: the
#'   model-implied paths.
#' @param n_tail Integer. Number of trailing observations to show. Default 40.
#' @param x_date_breaks Character. Spacing of the date labels, passed to
#'   [ggplot2::scale_x_date()], e.g. `"2 weeks"`. Default `"1 week"`.
#' @param title Character, `NULL` or `NA`. `NULL` gives a default title and
#'   `NA` removes it.
#' @param ylab Character. Y-axis label. Default `""`.
#' @param legend_position Legend position for [ggplot2::theme()]. Default
#'   `"bottom"`.
#' @param variable_names Character vector, one name per variable, or `NULL`
#'   for the column names of `realized` (or `V1`, `V2`, ...).
#' @param variables Integer vector or `NULL`. Variables to show; `NULL` shows
#'   all of them.
#' @param facet_ncol Integer. Number of panel columns. Default 2.
#' @return A [ggplot2::ggplot] object.
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3,
#'             dimnames = list(NULL, c("x1", "x2", "x3")))
#' dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 100)
#' fit <- etvar(y, dates, target_date = dates[80], lags = 2, nimp = 10)
#' plot_etvar_counterfactual(
#'   date_vec = fit$dates,
#'   realized = fit$data,
#'   counterfactual = fit$hd$all_data[, , fit$shock_index],
#'   variable_names = fit$variable_names,
#'   n_tail = 30
#' )
#' @export
plot_etvar_counterfactual <- function(date_vec, realized, counterfactual,
                                      n_tail = 40L,
                                      x_date_breaks = "1 week",
                                      title = NULL, ylab = "",
                                      legend_position = "bottom",
                                      variable_names = NULL,
                                      variables = NULL,
                                      facet_ncol = 2L) {
  date_vec <- .etvar_as_date(date_vec, "date_vec")
  n_tail <- .etvar_integer_scalar(n_tail, "n_tail", min = 1L)
  as_columns <- function(x, arg) {
    if (is.numeric(x) && is.null(dim(x))) {
      return(matrix(as.double(x), ncol = 1L))
    }
    .etvar_numeric_matrix(x, arg)
  }
  realized <- as_columns(realized, "realized")
  counterfactual <- as_columns(counterfactual, "counterfactual")
  if (!identical(dim(realized), dim(counterfactual)) ||
      nrow(realized) != length(date_vec) || any(!is.finite(realized)) ||
      any(!is.finite(counterfactual))) {
    stop("`realized` and `counterfactual` must be finite numeric vectors or ",
         "matrices with one row per date and matching dimensions.",
         call. = FALSE)
  }
  .etvar_string(x_date_breaks, "x_date_breaks")
  n <- ncol(realized)
  if (is.null(variable_names)) {
    variable_names <- colnames(realized)
    if (is.null(variable_names)) variable_names <- paste0("V", seq_len(n))
  }
  .etvar_check_labels(variable_names, n, "variable_names")
  if (is.null(variables)) variables <- seq_len(n)
  variables <- .etvar_index_vector(variables, "variables", n)
  facet_ncol <- .etvar_integer_scalar(facet_ncol, "facet_ncol", min = 1L)

  rows <- utils::tail(seq_along(date_vec), n_tail)
  labels <- panel_names(variable_names, variables)
  plot_data <- do.call(rbind, lapply(seq_along(variables), function(i) {
    v <- variables[i]
    data.frame(
      dates = rep(date_vec[rows], times = 2L),
      Series = rep(c("Counterfactual", "Realized"), each = length(rows)),
      Value = c(counterfactual[rows, v], realized[rows, v]),
      Variable = labels[i]
    )
  }))
  plot_data$Variable <- factor(plot_data$Variable, levels = labels)

  ggplot2::ggplot(plot_data, ggplot2::aes(x = dates, y = Value,
                                          linetype = Series, color = Series)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~Variable, ncol = facet_ncol, scales = "free_y") +
    ggplot2::scale_color_manual(
      values = c(Realized = "black", Counterfactual = "#0072B2")
    ) +
    ggplot2::scale_linetype_manual(
      values = c(Realized = "solid", Counterfactual = "dashed")
    ) +
    ggplot2::scale_x_date(date_labels = "%b %d, %y",
                          date_breaks = x_date_breaks) +
    ggplot2::labs(
      title = .etvar_title(title, "Counterfactual versus realised paths"),
      x = NULL,
      y = .etvar_title(ylab, "", arg = "ylab")
    ) +
    .etvar_plot_theme(legend_position) +
    rotated_dates()
}


#' Plot a structural error series
#'
#' Line plot of one structural error series, optionally with vertical
#' markers and labels at event dates.
#'
#' @param date_vec Vector of dates, one per error (the effective sample).
#' @param errors Numeric vector of structural errors.
#' @param n_tail Integer. Number of trailing observations to show. Default
#'   100.
#' @param event_dates Dates to mark with dashed red lines, or `NULL`.
#' @param event_labels Character vector of labels, one per event date, or
#'   `NULL`.
#' @param x_date_breaks Character. Spacing of the date labels, passed to
#'   [ggplot2::scale_x_date()], e.g. `"2 weeks"`. Default `"1 week"`.
#' @param title Character, `NULL` or `NA`. `NULL` gives a default title and
#'   `NA` removes it.
#' @param panel_title Character or `NULL`. Label of the panel strip; `NULL`
#'   uses `"Target shock"`.
#' @param line_color Character. Colour of the line. Default `"black"`.
#' @return A [ggplot2::ggplot] object.
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3,
#'             dimnames = list(NULL, c("x1", "x2", "x3")))
#' dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 100)
#' fit <- etvar(y, dates, target_date = dates[80], lags = 2, nimp = 10)
#' plot_etvar_errors(
#'   date_vec = fit$dates[-(1:2)],
#'   errors = fit$hd$errors_str[, fit$shock_index],
#'   n_tail = 40,
#'   event_dates = fit$target_date,
#'   event_labels = "Event"
#' )
#' @export
plot_etvar_errors <- function(date_vec, errors, n_tail = 100L,
                              event_dates = NULL, event_labels = NULL,
                              x_date_breaks = "1 week",
                              title = NULL,
                              panel_title = NULL,
                              line_color = "black") {
  date_vec <- .etvar_as_date(date_vec, "date_vec")
  n_tail <- .etvar_integer_scalar(n_tail, "n_tail", min = 1L)
  if (!is.numeric(errors) || !is.null(dim(errors)) ||
      length(errors) != length(date_vec) || any(!is.finite(errors))) {
    stop("`errors` must be a finite numeric vector with one value per date.",
         call. = FALSE)
  }
  .etvar_string(x_date_breaks, "x_date_breaks")
  .etvar_string(line_color, "line_color")
  panel_title <- .etvar_string(
    .etvar_title(panel_title, "Target shock", arg = "panel_title"),
    "panel_title"
  )
  if (!is.null(event_dates)) {
    event_dates <- .etvar_as_date(event_dates, "event_dates")
  }
  if (!is.null(event_labels) &&
      (!is.character(event_labels) || anyNA(event_labels) ||
       length(event_labels) != length(event_dates))) {
    stop("`event_labels` must contain one non-missing label per event date.",
         call. = FALSE)
  }

  plot_data <- data.frame(
    dates = utils::tail(date_vec, n_tail),
    Structural_error = utils::tail(errors, n_tail),
    Shock = factor(panel_title)
  )

  pl <- ggplot2::ggplot(plot_data,
                        ggplot2::aes(x = dates, y = Structural_error)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "gray55", linewidth = 0.5) +
    ggplot2::geom_line(colour = line_color, linewidth = 0.8) +
    ggplot2::facet_wrap(~Shock, ncol = 1L) +
    ggplot2::scale_x_date(date_labels = "%b %d, %y",
                          date_breaks = x_date_breaks) +
    ggplot2::labs(title = .etvar_title(title, "Structural error"), x = NULL,
                  y = "Structural error") +
    .etvar_plot_theme() +
    rotated_dates()

  if (is.null(event_dates)) {
    return(pl)
  }
  pl <- pl + ggplot2::geom_vline(xintercept = event_dates, linetype = "dashed",
                                 linewidth = 0.6, color = "red")
  # Labels of the first half of the events go at the bottom, the rest at the
  # top, so that neighbouring labels do not overlap.
  for (k in seq_along(event_labels)) {
    bottom <- k <= length(event_dates) / 2
    pl <- pl + ggplot2::annotate(
      "text", x = event_dates[k], y = if (bottom) -Inf else Inf,
      label = event_labels[k], hjust = 1.05, vjust = if (bottom) -0.5 else 1.5,
      size = 3
    )
  }
  pl
}


#' Plot structural impulse responses
#'
#' One panel per selected variable showing its response to a chosen
#' structural shock, optionally with confidence bands.
#'
#' @param irfs Array (n x n x nimp) of impulse responses, as in
#'   [compute_structural_irfs()].
#' @param variable_names Character vector of n variable names.
#' @param shock_index Integer. The shock whose responses are shown. Default 1.
#' @param scaling_index Integer vector `c(variable, shock, position)` or
#'   `NULL`. If given, every response (and every bootstrap draw) is divided
#'   by the element `irfs[variable, shock, position]`. `position` counts
#'   along the horizon dimension from 1, so 1 is the impact response, shown
#'   as horizon 0 in the plot. For example `c(2, 1, 1)` expresses the
#'   responses relative to the impact response of variable 2.
#' @param scaling_multiplier Number. Multiplies the responses after the
#'   normalisation, e.g. `0.5` for a shock half the normalised size.
#'   Default 1.
#' @param common_ylim Logical. Use the same y axis in every panel? Default
#'   `FALSE`.
#' @param variables Integer vector or `NULL`. Variables to show; `NULL` shows
#'   all of them.
#' @param title Character, `NULL` or `NA`. `NULL` gives a default title and
#'   `NA` removes it.
#' @param irf_ci Optional list with `lower` and `upper` arrays shaped like
#'   `irfs`, typically from [compute_irf_confint()]. If given, the bands are
#'   drawn. When it also holds the bootstrap `draws` and the responses are
#'   normalised with `scaling_index`, the bands are recomputed from the
#'   normalised draws.
#' @param ci_col,ci_border Fill and border colour of the bands.
#' @param facet_ncol Integer. Number of panel columns. Default 2.
#' @return A [ggplot2::ggplot] object.
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3,
#'             dimnames = list(NULL, c("x1", "x2", "x3")))
#' dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 100)
#' fit <- etvar(y, dates, target_date = dates[80], lags = 2, nimp = 10)
#' plot_etvar_irfs(fit$irfs, fit$variable_names, shock_index = 1)
#' @export
plot_etvar_irfs <- function(irfs, variable_names, shock_index = 1L,
                            scaling_index = NULL,
                            scaling_multiplier = 1,
                            common_ylim = FALSE,
                            variables = NULL,
                            title = NULL,
                            irf_ci = NULL,
                            ci_col = grDevices::adjustcolor("gray70",
                                                            alpha.f = 0.45),
                            ci_border = NA,
                            facet_ncol = 2L) {
  irfs <- .etvar_numeric_array(irfs, "irfs")
  if (length(dim(irfs)) != 3L || dim(irfs)[1L] != dim(irfs)[2L]) {
    stop("`irfs` must be an n by n by h numeric array.", call. = FALSE)
  }
  n <- dim(irfs)[1L]
  .etvar_check_labels(variable_names, n, "variable_names")
  shock_index <- .etvar_integer_scalar(shock_index, "shock_index",
                                       min = 1L, max = n)
  if (is.null(variables)) variables <- seq_len(n)
  variables <- .etvar_index_vector(variables, "variables", n)
  facet_ncol <- .etvar_integer_scalar(facet_ncol, "facet_ncol", min = 1L)
  if (!is.numeric(scaling_multiplier) || length(scaling_multiplier) != 1L ||
      !is.finite(scaling_multiplier)) {
    stop("`scaling_multiplier` must be one finite number.", call. = FALSE)
  }
  common_ylim <- .etvar_flag(common_ylim, "common_ylim")

  if (!is.null(irf_ci)) {
    if (!is.list(irf_ci) || is.null(irf_ci$lower) || is.null(irf_ci$upper)) {
      stop("`irf_ci` must be a list with `lower` and `upper` arrays.",
           call. = FALSE)
    }
    .etvar_numeric_array(irf_ci$lower, "irf_ci$lower", dim(irfs))
    .etvar_numeric_array(irf_ci$upper, "irf_ci$upper", dim(irfs))
  }

  scale <- 1
  if (!is.null(scaling_index)) {
    scaling_index <- .etvar_index_vector(scaling_index, "scaling_index",
                                         max(dim(irfs)))
    if (length(scaling_index) != 3L || any(scaling_index > dim(irfs))) {
      stop("`scaling_index` must contain valid variable, shock, and horizon ",
           "indices.", call. = FALSE)
    }
    scale <- irfs[scaling_index[1L], scaling_index[2L], scaling_index[3L]]
    if (abs(scale) <= sqrt(.Machine$double.eps)) {
      stop("The IRF selected by `scaling_index` is zero or too close to zero.",
           call. = FALSE)
    }
  }

  if (!is.null(irf_ci)) {
    lower <- irf_ci$lower / scale
    upper <- irf_ci$upper / scale
    draws <- irf_ci$draws
    if (!is.null(scaling_index) && !is.null(draws) &&
        is.numeric(irf_ci$level) && identical(dim(draws)[1:3], dim(irfs))) {
      # Bands of the normalised responses: divide every draw by its own
      # scaling element before taking the quantiles.
      alpha <- 1 - irf_ci$level
      ratios <- sweep(draws, 4L, draws[scaling_index[1L], scaling_index[2L],
                                       scaling_index[3L], ], "/")
      ratios <- matrix(ratios, ncol = dim(draws)[4L])
      lower <- array(row_quantiles(ratios, alpha / 2), dim = dim(irfs))
      upper <- array(row_quantiles(ratios, 1 - alpha / 2), dim = dim(irfs))
    }
    # A negative scale or multiplier swaps the bounds.
    lower <- lower * scaling_multiplier
    upper <- upper * scaling_multiplier
    band_low <- pmin(lower, upper)
    upper <- pmax(lower, upper)
    lower <- band_low
  }

  labels <- panel_names(variable_names, variables)
  plot_data <- do.call(rbind, lapply(seq_along(variables), function(i) {
    v <- variables[i]
    response <- irfs[v, shock_index, ] / scale * scaling_multiplier
    panel <- data.frame(
      Horizon = seq_along(response) - 1L,
      Response = response,
      Variable = labels[i]
    )
    if (!is.null(irf_ci)) {
      panel$Lower <- lower[v, shock_index, ]
      panel$Upper <- upper[v, shock_index, ]
    }
    panel
  }))
  plot_data$Variable <- factor(plot_data$Variable, levels = labels)

  pl <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Horizon, y = Response))
  if (!is.null(irf_ci)) {
    pl <- pl + ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper),
                                    fill = ci_col, colour = ci_border)
  }
  pl +
    ggplot2::geom_hline(yintercept = 0, color = "gray55", linetype = "dashed",
                        linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(shape = 1, size = 1.8) +
    ggplot2::facet_wrap(~Variable, ncol = facet_ncol,
                        scales = if (common_ylim) "fixed" else "free_y") +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(
      x = "Forecast horizon",
      y = "Response",
      title = .etvar_title(title, "Structural impulse response functions")
    ) +
    .etvar_plot_theme()
}


#' Plot the historical decomposition around the event
#'
#' Stacked bars of three contributions -- the target shock on the event date,
#' target shocks after it, and everything else -- with the data (net of the
#' deterministic component) as a black line and a red line at the event date.
#' The three contributions add up to the line. Positive and negative
#' contributions are stacked separately, each starting from zero with the
#' event shock.
#'
#' [plot.etvar()] counts the target shocks *before* the event among the other
#' shocks. To show the target shocks up to and including the event instead,
#' pass `fit$shock_contributions$target_stop` as `target_event` and
#' `fit$shock_contributions$others` as `others`.
#'
#' @param date_vec Vector of dates, one per row of `data`.
#' @param data Matrix (T x n) of observed data.
#' @param deterministic Matrix (T x n). Deterministic component.
#' @param target_event Matrix (T x n). Contribution of the event-date target
#'   shock, e.g. `fit$shock_contributions$target_event`.
#' @param target_rest Matrix (T x n). Contribution of target shocks after the
#'   event date.
#' @param others Matrix (T x n). Contribution of everything else.
#' @param datestar Integer. Number of observations after the event date.
#' @param variable_names Character vector of n variable names.
#' @param n_tail Integer. Number of trailing observations to show. Default 41.
#' @param x_break_every Integer. Label every `x_break_every`-th date. Default
#'   2.
#' @param title Character, `NULL` or `NA`. `NULL` gives a default title and
#'   `NA` removes it.
#' @param merge_target_components Logical. Show `target_event + target_rest`
#'   as a single target-shock component? Default `FALSE`.
#' @param shock_labels Character vector of three legend labels for the
#'   components.
#' @param variables Integer vector or `NULL`. Variables to show; `NULL` shows
#'   all of them.
#' @param facet_ncol Integer. Number of panel columns. Default 2.
#' @return A [ggplot2::ggplot] object.
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3,
#'             dimnames = list(NULL, c("x1", "x2", "x3")))
#' dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 100)
#' fit <- etvar(y, dates, target_date = dates[80], lags = 2, nimp = 10)
#' parts <- fit$shock_contributions
#' plot_etvar_hd(
#'   date_vec = fit$dates,
#'   data = fit$data,
#'   deterministic = fit$hd$deterministic,
#'   target_event = parts$target_event,
#'   target_rest = parts$target_rest,
#'   others = parts$others + parts$target_before,
#'   datestar = fit$datestar,
#'   variable_names = fit$variable_names,
#'   n_tail = 30
#' )
#' @export
plot_etvar_hd <- function(date_vec, data, deterministic, target_event,
                          target_rest, others, datestar, variable_names,
                          n_tail = 41L,
                          x_break_every = 2L,
                          title = NULL,
                          merge_target_components = FALSE,
                          shock_labels = c("Event shock", "Continued shocks",
                                           "Other shocks"),
                          variables = NULL,
                          facet_ncol = 2L) {
  date_vec <- .etvar_as_date(date_vec, "date_vec")
  data <- .etvar_numeric_matrix(data)
  T_obs <- nrow(data)
  n <- ncol(data)
  if (length(date_vec) != T_obs) {
    stop("`date_vec` must contain one date per row of `data`.", call. = FALSE)
  }
  parts <- list(deterministic = deterministic, target_event = target_event,
                target_rest = target_rest, others = others)
  parts <- Map(.etvar_numeric_matrix, parts, names(parts))
  if (!all(vapply(parts, function(x) identical(dim(x), dim(data)),
                  logical(1)))) {
    stop("Decomposition matrices must have the same dimensions as `data`.",
         call. = FALSE)
  }
  .etvar_check_labels(variable_names, n, "variable_names")
  datestar <- .etvar_integer_scalar(datestar, "datestar", min = 0L,
                                    max = T_obs - 1L)
  if (is.null(variables)) variables <- seq_len(n)
  variables <- .etvar_index_vector(variables, "variables", n)
  n_tail <- .etvar_integer_scalar(n_tail, "n_tail", min = 1L)
  x_break_every <- .etvar_integer_scalar(x_break_every, "x_break_every",
                                         min = 1L)
  facet_ncol <- .etvar_integer_scalar(facet_ncol, "facet_ncol", min = 1L)
  merge_target_components <- .etvar_flag(merge_target_components,
                                         "merge_target_components")
  if (!is.character(shock_labels) || length(shock_labels) != 3L ||
      anyNA(shock_labels) || any(!nzchar(shock_labels))) {
    stop("`shock_labels` must contain three non-empty labels.", call. = FALSE)
  }

  rows <- seq.int(max(1L, T_obs - n_tail + 1L), T_obs)
  if (merge_target_components) {
    components <- list(parts$target_event + parts$target_rest, parts$others)
    component_labels <- shock_labels[c(1L, 3L)]
    component_colors <- c("#D95319", "#4DBEEE")
  } else {
    components <- parts[c("target_event", "target_rest", "others")]
    component_labels <- shock_labels
    component_colors <- c("#D95319", "#EDB120", "#4DBEEE")
  }

  labels <- panel_names(variable_names, variables)
  bars <- do.call(rbind, lapply(seq_along(variables), function(i) {
    do.call(rbind, lapply(seq_along(components), function(k) {
      data.frame(dates = date_vec[rows],
                 Value = components[[k]][rows, variables[i]],
                 Component = component_labels[k],
                 Variable = labels[i])
    }))
  }))
  bars$Component <- factor(bars$Component, levels = component_labels)
  bars$Variable <- factor(bars$Variable, levels = labels)

  detrended <- data - parts$deterministic
  line <- do.call(rbind, lapply(seq_along(variables), function(i) {
    data.frame(dates = date_vec[rows],
               Value = detrended[rows, variables[i]],
               Variable = labels[i])
  }))
  line$Variable <- factor(line$Variable, levels = labels)

  breaks <- date_vec[rows][seq.int(1L, length(rows), by = x_break_every)]
  pl <- ggplot2::ggplot(bars, ggplot2::aes(x = dates, y = Value,
                                           fill = Component)) +
    ggplot2::geom_col(width = 0.8, colour = NA,
                      position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
    ggplot2::geom_line(data = line,
                       ggplot2::aes(x = dates, y = Value, colour = "Data"),
                       inherit.aes = FALSE, linewidth = 0.8) +
    ggplot2::facet_wrap(~Variable, ncol = facet_ncol, scales = "free_y") +
    ggplot2::scale_fill_manual(
      values = stats::setNames(component_colors, component_labels)
    ) +
    ggplot2::scale_colour_manual(values = c(Data = "black")) +
    ggplot2::scale_x_date(breaks = breaks, date_labels = "%b %d, %Y") +
    ggplot2::labs(x = NULL, y = "Detrended values",
                  title = .etvar_title(title, "Historical decomposition")) +
    .etvar_plot_theme() +
    rotated_dates()

  event_row <- T_obs - datestar
  if (event_row %in% rows) {
    pl <- pl + ggplot2::geom_vline(xintercept = date_vec[event_row],
                                   colour = "red", linewidth = 0.8)
  }
  pl
}
