# ==============================================================================
# Plotting functions for ETVAR results
# ==============================================================================

# Suppress R CMD check NOTEs for ggplot2 NSE columns
utils::globalVariables(c(
  "Value", "Series", "Structural_error", "dates",
  "Horizon", "Share", "Component", "Variable",
  "Response", "Lower", "Upper", "Shock"
))

etvar_select_scalar_title <- function(title, default_title, arg = "title") {
  if (is.null(title)) {
    return(default_title)
  }

  if (length(title) != 1L) {
    stop("`", arg, "` must be NULL, NA, or a character scalar for this plot.")
  }

  title <- as.character(title)
  if (is.na(title)) "" else title
}


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


#' Plot an ETVAR forecast error variance decomposition
#'
#' Draws a stacked-bar FEVD for each selected variable. The two components are
#' the event-targeted structural shock and the sum of all other structural
#' shocks. At every variable and horizon, the two plotted shares sum to one.
#'
#' @param fevd An `"etvar_fevd"` object from [compute_fevd()].
#' @param variables Integer vector or `NULL`. Variables to display. The default
#'   displays every variable.
#' @param title Character, `NULL`, or `NA`. A scalar plot title. `NULL` uses
#'   `"Structural forecast error variance decomposition"` and `NA` removes the
#'   title.
#' @param component_labels Character vector of length two. Labels for the target
#'   shock and all other shocks, respectively.
#' @param colors Character vector of length two. Fill colors for the target
#'   shock and all other shocks, respectively.
#' @param facet_ncol Integer. Number of columns in the faceted display.
#'   Default 2.
#' @return A [ggplot2::ggplot] object.
#' @export
plot_etvar_fevd <- function(
    fevd,
    variables = NULL,
    title = NULL,
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
    component_labels <- c(
      paste0("Target shock (", fevd$target_shock, ")"),
      "Other shocks"
    )
  }
  if (!is.character(component_labels) || length(component_labels) != 2L ||
      anyNA(component_labels) || any(!nzchar(component_labels)) ||
      anyDuplicated(component_labels)) {
    stop("`component_labels` must contain two distinct, non-empty labels.",
         call. = FALSE)
  }
  if (!is.character(colors) || length(colors) != 2L ||
      anyNA(colors) || any(!nzchar(colors))) {
    stop("`colors` must contain two non-empty color values.", call. = FALSE)
  }

  plot_title <- etvar_select_scalar_title(
    title,
    "Structural forecast error variance decomposition"
  )
  selected_names <- make.unique(variable_names[variables], sep = " #")
  horizons <- seq_len(nimp)

  target_values <- as.vector(t(target[variables, , drop = FALSE]))
  other_values <- as.vector(t(others[variables, , drop = FALSE]))
  variable_factor <- factor(
    rep(selected_names, each = nimp),
    levels = selected_names
  )

  plot_data <- rbind(
    data.frame(
      Horizon = rep(horizons, times = length(variables)),
      Share = target_values,
      Component = component_labels[1L],
      Variable = variable_factor
    ),
    data.frame(
      Horizon = rep(horizons, times = length(variables)),
      Share = other_values,
      Component = component_labels[2L],
      Variable = variable_factor
    )
  )
  plot_data$Component <- factor(
    plot_data$Component,
    levels = component_labels
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Horizon, y = Share, fill = Component)
  ) +
    ggplot2::geom_col(
      width = 0.9,
      color = "white",
      linewidth = 0.1
    ) +
    ggplot2::facet_wrap(~Variable, ncol = facet_ncol) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(colors, component_labels),
      breaks = component_labels
    ) +
    ggplot2::scale_x_continuous(
      breaks = function(limits) {
        unique(pmax(1, round(pretty(limits))))
      }
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
      title = plot_title
    ) +
    .etvar_plot_theme()
}


#' Plot counterfactual versus realised paths
#'
#' Creates a faceted ggplot comparing the model-implied counterfactual paths
#' (driven by the selected structural shock only) against the observed data.
#'
#' @param date_vec       Date vector. Full date index.
#' @param realized       Numeric vector or matrix. Observed data, with variables
#'   in columns.
#' @param counterfactual Numeric vector or matrix with the same dimensions as
#'   `realized`. Model-implied paths.
#' @param n_tail   Integer. Number of trailing observations to display.
#'   Default 40.
#' @param x_date_breaks Character. Date break spacing passed to
#'   [ggplot2::scale_x_date()] (e.g. `"1 week"`, `"2 weeks"`, `"1 month"`).
#'   Default `"1 week"`.
#' @param title Character scalar, `NULL`, or `NA`. `NULL` uses
#'   `"Counterfactual versus realised paths"`, and `NA` removes the title.
#' @param ylab     Character. Y-axis label. Default `""`.
#' @param legend_position Legend position passed to [ggplot2::theme()]. Default
#'   `"bottom"`.
#' @param variable_names Character vector with one label per variable. If
#'   `NULL`, column names are used when available, otherwise `V1`, `V2`, ...
#' @param variables Integer vector or `NULL`. Variables to display. `NULL`
#'   displays all variables.
#' @param facet_ncol Integer. Number of facet columns. Default 2.
#' @return A [ggplot2::ggplot] object.
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
  if (is.numeric(realized) && is.null(dim(realized))) {
    realized <- matrix(as.double(realized), ncol = 1L)
  } else {
    realized <- .etvar_numeric_matrix(realized, "realized")
  }
  if (is.numeric(counterfactual) && is.null(dim(counterfactual))) {
    counterfactual <- matrix(as.double(counterfactual), ncol = 1L)
  } else {
    counterfactual <- .etvar_numeric_matrix(counterfactual, "counterfactual")
  }
  if (!identical(dim(realized), dim(counterfactual)) ||
      nrow(realized) != length(date_vec) || anyNA(realized) ||
      anyNA(counterfactual) || any(!is.finite(realized)) ||
      any(!is.finite(counterfactual))) {
    stop("`realized` and `counterfactual` must be finite numeric vectors or ",
         "matrices with one row per date and matching dimensions.",
         call. = FALSE)
  }
  if (!is.character(x_date_breaks) || length(x_date_breaks) != 1L ||
      is.na(x_date_breaks) || !nzchar(x_date_breaks)) {
    stop("`x_date_breaks` must be a non-empty character scalar.",
         call. = FALSE)
  }
  n <- ncol(realized)
  if (is.null(variable_names)) {
    variable_names <- colnames(realized)
    if (is.null(variable_names)) variable_names <- paste0("V", seq_len(n))
  }
  if (!is.character(variable_names) || length(variable_names) != n ||
      anyNA(variable_names) || any(!nzchar(variable_names))) {
    stop("`variable_names` must contain one non-empty name per variable.",
         call. = FALSE)
  }
  if (is.null(variables)) variables <- seq_len(n)
  variables <- .etvar_index_vector(variables, "variables", n)
  facet_ncol <- .etvar_integer_scalar(facet_ncol, "facet_ncol", min = 1L)
  title <- etvar_select_scalar_title(
    title, "Counterfactual versus realised paths"
  )
  ylab <- etvar_select_scalar_title(ylab, "", arg = "ylab")

  selected_dates <- utils::tail(date_vec, n_tail)
  selected_rows <- utils::tail(
    seq_len(nrow(realized)), length(selected_dates)
  )
  selected_names <- make.unique(variable_names[variables], sep = " #")
  plot_long <- do.call(
    rbind,
    lapply(seq_along(variables), function(panel_index) {
      variable_index <- variables[panel_index]
      data.frame(
        dates = rep(selected_dates, times = 2L),
        Series = rep(c("Counterfactual", "Realized"),
                     each = length(selected_dates)),
        Value = c(
          counterfactual[selected_rows, variable_index],
          realized[selected_rows, variable_index]
        ),
        Variable = selected_names[panel_index]
      )
    })
  )
  plot_long$Variable <- factor(plot_long$Variable, levels = selected_names)

  ggplot2::ggplot(plot_long,
                  ggplot2::aes(x = dates, y = Value,
                               linetype = Series, color = Series)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(
      ~Variable, ncol = facet_ncol, scales = "free_y"
    ) +
    ggplot2::scale_color_manual(
      values = c("Realized" = "black", "Counterfactual" = "#0072B2")) +
    ggplot2::scale_linetype_manual(
      values = c("Realized" = "solid", "Counterfactual" = "dashed")) +
    ggplot2::scale_x_date(date_labels = "%b %d, %y",
                          date_breaks = x_date_breaks) +
    ggplot2::labs(title = title, x = NULL, y = ylab) +
    .etvar_plot_theme(legend_position) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}


#' Plot structural error time series
#'
#' Creates a ggplot of one structural error series with optional vertical event
#' lines and text annotations. A facet strip identifies the plotted shock so
#' that the visual style matches the other ETVAR plots.
#'
#' @param date_vec      Date vector. Full date index for the structural errors.
#' @param errors Numeric vector. Structural error series.
#' @param n_tail        Integer. Number of trailing observations. Default 100.
#' @param event_dates   Date vector (or `NULL`). Dates at which to draw
#'   vertical red dashed lines.
#' @param event_labels  Character vector. Labels for the event dates.
#' @param x_date_breaks Character. Date break spacing passed to
#'   [ggplot2::scale_x_date()] (e.g. `"1 week"`, `"2 weeks"`, `"1 month"`).
#'   Default `"1 week"`.
#' @param title Character scalar, `NULL`, or `NA`. `NULL` uses
#'   `"Structural error"`, and `NA` removes the title.
#' @param panel_title Character scalar or `NULL`. Label displayed in the facet
#'   strip. `NULL` uses `"Target shock"`.
#' @param line_color Character scalar. Colour of the structural-error line.
#'   Default `"black"`.
#' @return A [ggplot2::ggplot] object.
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
      length(errors) != length(date_vec) || anyNA(errors) ||
      any(!is.finite(errors))) {
    stop("`errors` must be a finite numeric vector with one value per date.",
         call. = FALSE)
  }
  if (!is.character(x_date_breaks) || length(x_date_breaks) != 1L ||
      is.na(x_date_breaks) || !nzchar(x_date_breaks)) {
    stop("`x_date_breaks` must be a non-empty character scalar.",
         call. = FALSE)
  }
  panel_title <- etvar_select_scalar_title(
    panel_title, "Target shock", arg = "panel_title"
  )
  if (!nzchar(panel_title)) {
    stop("`panel_title` must be a non-empty character scalar or `NULL`.",
         call. = FALSE)
  }
  if (!is.character(line_color) || length(line_color) != 1L ||
      is.na(line_color) || !nzchar(line_color)) {
    stop("`line_color` must be a non-empty character scalar.", call. = FALSE)
  }
  title <- etvar_select_scalar_title(title, "Structural error")
  if (!is.null(event_dates)) {
    event_dates <- .etvar_as_date(event_dates, "event_dates")
  }
  if (!is.null(event_labels) &&
      (!is.character(event_labels) || length(event_labels) != length(event_dates) ||
       anyNA(event_labels))) {
    stop("`event_labels` must contain one non-missing label per event date.",
         call. = FALSE)
  }

  plot_df <- data.frame(
    dates = utils::tail(date_vec, n_tail),
    Structural_error = utils::tail(errors, n_tail),
    Shock = panel_title
  )
  plot_df$Shock <- factor(plot_df$Shock, levels = panel_title)

  pl <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = dates, y = Structural_error)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "gray55", linewidth = 0.5) +
    ggplot2::geom_line(colour = line_color, linewidth = 0.8) +
    ggplot2::facet_wrap(~Shock, ncol = 1L) +
    ggplot2::scale_x_date(date_labels = "%b %d, %y",
                          date_breaks = x_date_breaks) +
    ggplot2::labs(title = title, x = NULL, y = "Structural error") +
    .etvar_plot_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  if (!is.null(event_dates)) {
    pl <- pl + ggplot2::geom_vline(
      xintercept = event_dates, linetype = "dashed",
      linewidth = 0.6, color = "red"
    )
    if (!is.null(event_labels)) {
      n_evt <- length(event_dates)
      for (k in seq_len(n_evt)) {
        vj <- if (k <= n_evt / 2) -0.5 else 1.5
        yi <- if (k <= n_evt / 2) -Inf else Inf
        pl <- pl + ggplot2::annotate(
          "text", x = event_dates[k], y = yi,
          label = event_labels[k], hjust = 1.05, vjust = vj, size = 3
        )
      }
    }
  }

  pl
}


#' Plot structural impulse response functions
#'
#' Produces a faceted ggplot showing the normalised response of each selected
#' variable to a chosen structural shock.
#'
#' @param irfs     Array (n x n x nimp). Structural IRFs from
#'   [compute_structural_irfs()].
#' @param variable_names Character vector of length n. Variable labels.
#' @param shock_index    Integer. Which shock to display. Default 1.
#' @param scaling_index  Numeric vector `c(var, shock, horizon)`. The IRF
#'   element used to normalise all responses. Default `NULL` (no normalisation).
#' @param scaling_multiplier Numeric scalar. Multiplier applied after
#'   normalisation. Use e.g. `0.5` to scale a unit-normalised shock to 0.5.
#'   Default 1.
#' @param common_ylim Logical. If `TRUE`, use a common y-axis range across
#'   all displayed variables. If `FALSE`, each panel uses its own y-range.
#'   Default `FALSE`.
#' @param variables Integer vector. Which variables to plot (indices). Default:
#'   all variables.
#' @param title Character scalar, `NULL`, or `NA`. `NULL` uses
#'   `"Structural impulse response functions"`, and `NA` removes the title.
#' @param irf_ci Optional list with `lower` and `upper` arrays matching `irfs`,
#'   typically from [compute_irf_confint()]. When supplied, confidence bands are
#'   drawn around each IRF.
#' @param conf_int Logical or `NULL`. If `TRUE`, draw confidence bands from
#'   `irf_ci`; if `FALSE`, suppress them. The default `NULL` draws bands when
#'   `irf_ci` is supplied. The standalone plotter cannot compute confidence
#'   intervals from `irfs` alone.
#' @param ci_col Fill colour for the confidence band.
#' @param ci_border Border colour for the confidence band. Default `NA`.
#' @param facet_ncol Integer. Number of facet columns. Default 2.
#' @return A [ggplot2::ggplot] object.
#' @export
plot_etvar_irfs <- function(irfs, variable_names, shock_index = 1L,
                             scaling_index = NULL,
                             scaling_multiplier = 1,
                             common_ylim = FALSE,
                             variables = NULL,
                             title = NULL,
                             irf_ci = NULL,
                             conf_int = NULL,
                             ci_col = grDevices::adjustcolor("gray70",
                                                             alpha.f = 0.45),
                             ci_border = NA,
                             facet_ncol = 2L) {
  irfs <- .etvar_numeric_array(irfs, "irfs")
  if (length(dim(irfs)) != 3L || dim(irfs)[1L] != dim(irfs)[2L]) {
    stop("`irfs` must be an n by n by h numeric array.", call. = FALSE)
  }
  n <- length(variable_names)
  if (!is.character(variable_names) || n != dim(irfs)[1L] ||
      anyNA(variable_names) || any(!nzchar(variable_names))) {
    stop("`variable_names` must contain one non-empty name per variable.",
         call. = FALSE)
  }
  shock_index <- .etvar_integer_scalar(shock_index, "shock_index",
                                        min = 1L, max = n)
  if (is.null(variables)) variables <- seq_len(n)
  variables <- .etvar_index_vector(variables, "variables", n)
  facet_ncol <- .etvar_integer_scalar(facet_ncol, "facet_ncol", min = 1L)
  if (!is.numeric(scaling_multiplier) || length(scaling_multiplier) != 1L ||
      is.na(scaling_multiplier) || !is.finite(scaling_multiplier)) {
    stop("`scaling_multiplier` must be one finite number.", call. = FALSE)
  }
  common_ylim <- .etvar_flag(common_ylim, "common_ylim")

  if (!is.null(conf_int) &&
      (!is.logical(conf_int) || length(conf_int) != 1L || is.na(conf_int))) {
    stop("`conf_int` must be TRUE, FALSE, or NULL.", call. = FALSE)
  }

  draw_ci <- if (is.null(conf_int)) !is.null(irf_ci) else isTRUE(conf_int)

  if (draw_ci && is.null(irf_ci)) {
    stop("`plot_etvar_irfs()` cannot compute confidence intervals from ",
         "`irfs` alone. Use `compute_irf_confint(fit)` and pass the result ",
         "as `irf_ci`, or call `plot(fit, type = \"irfs\", conf_int = TRUE)`.")
  }

  if (draw_ci) {
    if (!is.list(irf_ci) || is.null(irf_ci$lower) || is.null(irf_ci$upper)) {
      stop("`irf_ci` must be a list with `lower` and `upper` arrays.",
           call. = FALSE)
    }
    if (!identical(dim(irf_ci$lower), dim(irfs)) ||
        !identical(dim(irf_ci$upper), dim(irfs))) {
      stop("`irf_ci$lower` and `irf_ci$upper` must have the same dimensions ",
           "as `irfs`.", call. = FALSE)
    }
    .etvar_numeric_array(irf_ci$lower, "irf_ci$lower", dim(irfs))
    .etvar_numeric_array(irf_ci$upper, "irf_ci$upper", dim(irfs))
  }

  scaling_factor <- if (!is.null(scaling_index)) {
    scaling_index <- .etvar_index_vector(
      scaling_index, "scaling_index", max(dim(irfs))
    )
    if (length(scaling_index) != 3L ||
        any(scaling_index > dim(irfs))) {
      stop("`scaling_index` must contain valid variable, shock, and horizon ",
           "indices.", call. = FALSE)
    }
    irfs[scaling_index[1], scaling_index[2], scaling_index[3]]
  } else {
    1
  }
  if (!is.finite(scaling_factor) ||
      abs(scaling_factor) <= sqrt(.Machine$double.eps)) {
    stop("The IRF selected by `scaling_index` is zero or too close to zero.",
         call. = FALSE)
  }

  selected_names <- make.unique(variable_names[variables], sep = " #")
  plot_data <- do.call(
    rbind,
    lapply(seq_along(variables), function(panel_index) {
      variable_index <- variables[panel_index]
      response <- (irfs[variable_index, shock_index, ] / scaling_factor) *
        scaling_multiplier
      panel_data <- data.frame(
        Horizon = seq_along(response) - 1L,
        Response = as.double(response),
        Variable = selected_names[panel_index]
      )
      if (draw_ci) {
        lower <- (irf_ci$lower[variable_index, shock_index, ] /
                    scaling_factor) * scaling_multiplier
        upper <- (irf_ci$upper[variable_index, shock_index, ] /
                    scaling_factor) * scaling_multiplier
        panel_data$Lower <- as.double(pmin(lower, upper))
        panel_data$Upper <- as.double(pmax(lower, upper))
      }
      panel_data
    })
  )
  plot_data$Variable <- factor(plot_data$Variable, levels = selected_names)
  plot_title <- etvar_select_scalar_title(
    title, "Structural impulse response functions"
  )

  pl <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Horizon, y = Response)
  )

  if (draw_ci) {
    pl <- pl + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = Lower, ymax = Upper),
      fill = ci_col,
      colour = ci_border
    )
  }

  pl +
    ggplot2::geom_hline(
      yintercept = 0, color = "gray55", linetype = "dashed", linewidth = 0.5
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(shape = 1, size = 1.8) +
    ggplot2::facet_wrap(
      ~Variable,
      ncol = facet_ncol,
      scales = if (common_ylim) "fixed" else "free_y"
    ) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(
      x = "Forecast horizon",
      y = "Response",
      title = plot_title
    ) +
    .etvar_plot_theme()
}


#' Plot historical decomposition as stacked bar charts
#'
#' For each variable, produces a stacked bar chart showing the contributions of
#' (1) the event-date target shock, (2) continued target shocks after the event,
#' and (3) all other shocks. The actual detrended data is overlaid as a black
#' line, and a vertical red line marks the event date. Within the positive and
#' negative parts of each bar, the event shock is stacked from zero, followed
#' by continued target shocks and then all other shocks.
#'
#' @param date_vec         Date vector of length T.
#' @param data             Matrix (T x n). Observed data.
#' @param deterministic    Matrix (T x n). Deterministic trend component.
#' @param target_stop      Matrix (T x n). Event-date target shock
#'   (intervention path).
#' @param target_rest      Matrix (T x n). Continued target shocks after the
#'   event.
#' @param others           Matrix (T x n). Combined non-target shocks.
#' @param datestar         Integer. Offset from end marking the event date.
#' @param variable_names   Character vector of length n.
#' @param variables Integer vector or `NULL`. Variables to display. `NULL`
#'   displays all variables.
#' @param n_tail           Integer. Number of trailing observations. Default 41.
#' @param x_break_every Integer. Show every k-th date label on the x-axis
#'   Default 2.
#' @param title Character scalar, `NULL`, or `NA`. `NULL` uses
#'   `"Historical decomposition"`, and `NA` removes the title.
#' @param merge_target_components Logical. If `TRUE`, merge `target_stop`
#'   (event shock) and `target_rest` (continued target shocks) into a single
#'   target-shock bar component. Default `FALSE`.
#' @param shock_labels     Character vector of length 3. Labels for the three
#'   bar components. Default: `c("Event shock", "Continued shocks",
#'   "Other shocks")`.
#' @param open_device Deprecated compatibility argument. It is validated but
#'   otherwise ignored because the function now returns a ggplot object.
#' @param ask Deprecated compatibility argument. It is validated but otherwise
#'   ignored; faceted plots never prompt between variables.
#' @param facet_ncol Integer. Number of facet columns. Default 2.
#' @return A [ggplot2::ggplot] object.
#' @export
plot_etvar_hd <- function(date_vec, data, deterministic, target_stop,
                           target_rest, others, datestar, variable_names,
                           n_tail = 41L,
                           x_break_every = 2L,
                           title = NULL,
                           merge_target_components = FALSE,
                           shock_labels = c("Event shock",
                                            "Continued shocks",
                                            "Other shocks"),
                           open_device = FALSE,
                           ask = interactive(),
                           variables = NULL,
                           facet_ncol = 2L) {
  date_vec <- .etvar_as_date(date_vec, "date_vec")
  data <- .etvar_numeric_matrix(data)
  T_obs <- nrow(data)
  n <- ncol(data)
  if (length(date_vec) != T_obs) {
    stop("`date_vec` must contain one date per row of `data`.", call. = FALSE)
  }
  expected_dim <- c(T_obs, n)
  deterministic <- .etvar_numeric_matrix(deterministic, "deterministic")
  target_stop <- .etvar_numeric_matrix(target_stop, "target_stop")
  target_rest <- .etvar_numeric_matrix(target_rest, "target_rest")
  others <- .etvar_numeric_matrix(others, "others")
  matrices <- list(deterministic = deterministic, target_stop = target_stop,
                   target_rest = target_rest, others = others)
  if (any(!vapply(matrices, function(x) identical(dim(x), expected_dim),
                  logical(1)))) {
    stop("Decomposition matrices must have the same dimensions as `data`.",
         call. = FALSE)
  }
  if (!is.character(variable_names) || length(variable_names) != n ||
      anyNA(variable_names) || any(!nzchar(variable_names))) {
    stop("`variable_names` must contain one non-empty name per variable.",
         call. = FALSE)
  }
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
  open_device <- .etvar_flag(open_device, "open_device")
  ask <- .etvar_flag(ask, "ask")
  if (!is.character(shock_labels) || length(shock_labels) != 3L ||
      anyNA(shock_labels) || any(!nzchar(shock_labels))) {
    stop("`shock_labels` must contain three non-empty labels.", call. = FALSE)
  }

  start_idx  <- max(1, T_obs - n_tail + 1)
  plot_idx   <- start_idx:T_obs
  selected_dates <- date_vec[plot_idx]
  selected_names <- make.unique(variable_names[variables], sep = " #")

  if (merge_target_components) {
    components <- list(
      target_stop[plot_idx, , drop = FALSE] +
        target_rest[plot_idx, , drop = FALSE],
      others[plot_idx, , drop = FALSE]
    )
    component_labels <- c(shock_labels[1L], shock_labels[3L])
    component_colors <- c("#D95319", "#4DBEEE")
  } else {
    components <- list(
      target_stop[plot_idx, , drop = FALSE],
      target_rest[plot_idx, , drop = FALSE],
      others[plot_idx, , drop = FALSE]
    )
    component_labels <- shock_labels
    component_colors <- c("#D95319", "#EDB120", "#4DBEEE")
  }

  component_data <- do.call(
    rbind,
    lapply(seq_along(variables), function(panel_index) {
      variable_index <- variables[panel_index]
      do.call(
        rbind,
        lapply(seq_along(components), function(component_index) {
          data.frame(
            dates = selected_dates,
            Value = components[[component_index]][, variable_index],
            Component = component_labels[component_index],
            Variable = selected_names[panel_index]
          )
        })
      )
    })
  )
  component_data$Component <- factor(
    component_data$Component, levels = component_labels
  )
  component_data$Variable <- factor(
    component_data$Variable, levels = selected_names
  )

  line_data <- do.call(
    rbind,
    lapply(seq_along(variables), function(panel_index) {
      variable_index <- variables[panel_index]
      data.frame(
        dates = selected_dates,
        Value = data[plot_idx, variable_index] -
          deterministic[plot_idx, variable_index],
        Variable = selected_names[panel_index]
      )
    })
  )
  line_data$Variable <- factor(line_data$Variable, levels = selected_names)

  plot_title <- etvar_select_scalar_title(title, "Historical decomposition")
  date_breaks <- selected_dates[
    seq.int(1L, length(selected_dates), by = x_break_every)
  ]
  pl <- ggplot2::ggplot(
    component_data,
    ggplot2::aes(x = dates, y = Value, fill = Component)
  ) +
    ggplot2::geom_col(
      width = 0.8,
      colour = NA,
      position = ggplot2::position_stack(reverse = TRUE)
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
    ggplot2::geom_line(
      data = line_data,
      ggplot2::aes(x = dates, y = Value, colour = "Data"),
      inherit.aes = FALSE,
      linewidth = 0.8
    ) +
    ggplot2::facet_wrap(~Variable, ncol = facet_ncol, scales = "free_y") +
    ggplot2::scale_fill_manual(
      values = stats::setNames(component_colors, component_labels)
    ) +
    ggplot2::scale_colour_manual(values = c("Data" = "black")) +
    ggplot2::scale_x_date(
      breaks = date_breaks,
      date_labels = "%b %d, %Y"
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Detrended values",
      title = plot_title
    ) +
    .etvar_plot_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  target_index <- T_obs - datestar
  if (target_index %in% plot_idx) {
    pl <- pl + ggplot2::geom_vline(
      xintercept = date_vec[target_index],
      colour = "red",
      linewidth = 0.8
    )
  }

  pl
}
