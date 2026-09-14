# ==============================================================================
# Plotting functions for ETVAR results
# ==============================================================================

# Suppress R CMD check NOTEs for ggplot2 NSE columns
utils::globalVariables(c("Value", "Series", "Tariff_err", "dates"))

etvar_select_panel_title <- function(title, default_titles, variable_index,
                                     panel_index = variable_index,
                                     n_panels = length(default_titles),
                                     arg = "title") {
  if (is.null(title)) {
    return(default_titles[variable_index])
  }

  if (length(title) == 0L) {
    stop("`", arg, "` must be NULL, NA, or a non-empty character vector.")
  }

  title <- as.character(title)

  if (length(title) == 1L) {
    if (is.na(title)) return("")
    return(title)
  }

  title_names <- names(title)
  default_title <- default_titles[variable_index]
  if (!is.null(title_names) && default_title %in% title_names) {
    selected <- title[[default_title]]
    if (is.na(selected)) return("")
    return(selected)
  }

  if (length(title) == length(default_titles)) {
    selected <- title[variable_index]
    if (is.na(selected)) return("")
    return(selected)
  }

  if (!is.null(n_panels) && length(title) == n_panels) {
    selected <- title[panel_index]
    if (is.na(selected)) return("")
    return(selected)
  }

  stop("`", arg, "` must be NULL, NA, a character scalar, or a character ",
       "vector with one entry per variable/panel.")
}


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


#' Plot counterfactual versus realised paths
#'
#' Creates a ggplot comparing the model-implied counterfactual path (driven by
#' the target structural shock only) against the observed data for a single
#' variable.
#'
#' @param date_vec       Date vector. Full date index.
#' @param realized       Numeric vector. Observed data for the variable.
#' @param counterfactual Numeric vector. Model-implied path (shock only).
#' @param n_tail   Integer. Number of trailing observations to display.
#'   Default 40.
#' @param x_date_breaks Character. Date break spacing passed to
#'   [ggplot2::scale_x_date()] (e.g. `"1 week"`, `"2 weeks"`, `"1 month"`).
#'   Default `"1 week"`.
#' @param title    Character. Plot title. Default `""`.
#' @param ylab     Character. Y-axis label. Default `""`.
#' @return A [ggplot2::ggplot] object.
#' @export
plot_etvar_counterfactual <- function(date_vec, realized, counterfactual,
                                      n_tail = 40L,
                                      x_date_breaks = "1 week",
                                      title = "", ylab = "",
                                      legend_position = c(0.75, 0.2)) {
  plot_df <- data.frame(
    dates          = utils::tail(date_vec, n_tail),
    Counterfactual = as.double(utils::tail(counterfactual, n_tail)),
    Realized       = as.double(utils::tail(realized, n_tail))
  )

  plot_long <- tidyr::pivot_longer(
    plot_df, cols = c("Counterfactual", "Realized"),
    names_to = "Series", values_to = "Value"
  )

  ggplot2::ggplot(plot_long,
                  ggplot2::aes(x = dates, y = Value,
                               linetype = Series, color = Series)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(
      values = c("Realized" = "black", "Counterfactual" = "blue")) +
    ggplot2::scale_linetype_manual(
      values = c("Realized" = "solid", "Counterfactual" = "dashed")) +
    ggplot2::scale_x_date(date_labels = "%b %d, %y",
                          date_breaks = x_date_breaks) +
    ggplot2::labs(title = title, x = NULL, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position  = legend_position,
      legend.title     = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )
}


#' Plot structural error time series
#'
#' Creates a ggplot of a single structural error series with optional vertical
#' event lines and text annotations.
#'
#' @param date_vec      Date vector. Full date index for the structural errors.
#' @param errors        Numeric vector. The structural error series to plot.
#' @param n_tail        Integer. Number of trailing observations. Default 100.
#' @param event_dates   Date vector (or `NULL`). Dates at which to draw
#'   vertical red dashed lines.
#' @param event_labels  Character vector. Labels for the event dates.
#' @param x_date_breaks Character. Date break spacing passed to
#'   [ggplot2::scale_x_date()] (e.g. `"1 week"`, `"2 weeks"`, `"1 month"`).
#'   Default `"1 week"`.
#' @param title         Character. Plot title. Default `"Structural error"`.
#' @return A [ggplot2::ggplot] object.
#' @export
plot_etvar_errors <- function(date_vec, errors, n_tail = 100L,
                               event_dates = NULL, event_labels = NULL,
                               x_date_breaks = "1 week",
                               title = "Structural error") {
  plot_df <- data.frame(
    dates      = utils::tail(date_vec, n_tail),
    Tariff_err = utils::tail(errors, n_tail)
  )

  pl <- ggplot2::ggplot(plot_df, ggplot2::aes(x = dates, y = Tariff_err)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "gray55", linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_x_date(date_labels = "%b %d, %y",
                          date_breaks = x_date_breaks) +
    ggplot2::labs(title = title, x = NULL, y = "%") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.title     = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(event_dates)) {
    pl <- pl + ggplot2::geom_vline(
      xintercept = event_dates, linetype = "dashed",
      linewidth = 0.6, color = "red"
    )
    n_evt <- min(length(event_dates), length(event_labels))
    for (k in seq_len(n_evt)) {
      vj <- if (k <= n_evt / 2) -0.5 else 1.5
      yi <- if (k <= n_evt / 2) -Inf else Inf
      pl <- pl + ggplot2::annotate(
        "text", x = event_dates[k], y = yi,
        label = event_labels[k], hjust = 1.05, vjust = vj, size = 3
      )
    }
  }

  pl
}


#' Plot structural impulse response functions
#'
#' Produces a multi-panel base-R figure showing the normalised response of
#' each variable to a chosen structural shock.
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
#' @param title Character, `NULL`, or `NA`. `NULL` uses `variable_names`, `NA`
#'   removes titles, a scalar applies to every panel, and a vector supplies
#'   per-variable/panel titles. Named vectors are matched to `variable_names`.
#' @param irf_ci Optional list with `lower` and `upper` arrays matching `irfs`,
#'   typically from [compute_irf_confint()]. When supplied, confidence bands are
#'   drawn around each IRF.
#' @param conf_int Logical or `NULL`. If `TRUE`, draw confidence bands from
#'   `irf_ci`; if `FALSE`, suppress them. The default `NULL` draws bands when
#'   `irf_ci` is supplied. The standalone plotter cannot compute confidence
#'   intervals from `irfs` alone.
#' @param ci_col Fill colour for the confidence band.
#' @param ci_border Border colour for the confidence band. Default `NA`.
#' @return `NULL` (invisibly). Produces a plot on the current device.
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
                             ci_border = NA) {
  n <- length(variable_names)
  if (is.null(variables)) variables <- seq_len(n)
  n_plot <- length(variables)

  if (!is.null(conf_int) &&
      (!is.logical(conf_int) || length(conf_int) != 1L || is.na(conf_int))) {
    stop("`conf_int` must be TRUE, FALSE, or NULL.")
  }

  draw_ci <- if (is.null(conf_int)) !is.null(irf_ci) else isTRUE(conf_int)

  if (draw_ci && is.null(irf_ci)) {
    stop("`plot_etvar_irfs()` cannot compute confidence intervals from ",
         "`irfs` alone. Use `compute_irf_confint(fit)` and pass the result ",
         "as `irf_ci`, or call `plot(fit, type = \"irfs\", conf_int = TRUE)`.")
  }

  if (draw_ci) {
    if (!is.list(irf_ci) || is.null(irf_ci$lower) || is.null(irf_ci$upper)) {
      stop("`irf_ci` must be a list with `lower` and `upper` arrays.")
    }
    if (!identical(dim(irf_ci$lower), dim(irfs)) ||
        !identical(dim(irf_ci$upper), dim(irfs))) {
      stop("`irf_ci$lower` and `irf_ci$upper` must have the same dimensions ",
           "as `irfs`.")
    }
  }

  scaling_factor <- if (!is.null(scaling_index)) {
    irfs[scaling_index[1], scaling_index[2], scaling_index[3]]
  } else {
    1
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  nr <- ceiling(n_plot / 2)
  graphics::par(mfrow = c(nr, 2), mar = c(4, 4, 3, 1))

  y_limits_common <- NULL
  if (isTRUE(common_ylim)) {
    scaled_series <- lapply(variables, function(idx) {
      med_vals <- (irfs[idx, shock_index, ] / scaling_factor) *
        scaling_multiplier
      if (!draw_ci) return(med_vals)

      lower_vals <- (irf_ci$lower[idx, shock_index, ] / scaling_factor) *
        scaling_multiplier
      upper_vals <- (irf_ci$upper[idx, shock_index, ] / scaling_factor) *
        scaling_multiplier
      c(med_vals, pmin(lower_vals, upper_vals), pmax(lower_vals, upper_vals))
    })
    all_vals <- unlist(scaled_series, use.names = FALSE)
    y_min <- min(0, min(all_vals, na.rm = TRUE))
    y_max <- max(0, max(all_vals, na.rm = TRUE))
    y_span <- y_max - y_min
    y_pad <- if (isTRUE(all.equal(y_span, 0))) max(abs(y_max), 1) * 0.1 else y_span * 0.1
    y_limits_common <- c(y_min - y_pad, y_max + y_pad)
  }

  for (panel_pos in seq_along(variables)) {
    idx <- variables[panel_pos]
    series   <- irfs[idx, shock_index, ]
    med_vals <- (series / scaling_factor) * scaling_multiplier
    lower_vals <- upper_vals <- NULL

    if (draw_ci) {
      lower_raw <- (irf_ci$lower[idx, shock_index, ] / scaling_factor) *
        scaling_multiplier
      upper_raw <- (irf_ci$upper[idx, shock_index, ] / scaling_factor) *
        scaling_multiplier
      lower_vals <- pmin(lower_raw, upper_raw)
      upper_vals <- pmax(lower_raw, upper_raw)
    }

    if (isTRUE(common_ylim)) {
      y_limits <- y_limits_common
    } else {
      panel_vals <- if (!draw_ci) {
        med_vals
      } else {
        c(med_vals, lower_vals, upper_vals)
      }
      y_min_i <- min(0, min(panel_vals, na.rm = TRUE))
      y_max_i <- max(0, max(panel_vals, na.rm = TRUE))
      y_span_i <- y_max_i - y_min_i
      y_pad_i <- if (isTRUE(all.equal(y_span_i, 0))) max(abs(y_max_i), 1) * 0.1 else y_span_i * 0.1
      y_limits <- c(y_min_i - y_pad_i, y_max_i + y_pad_i)
    }

    panel_title <- etvar_select_panel_title(
      title          = title,
      default_titles = variable_names,
      variable_index = idx,
      panel_index    = panel_pos,
      n_panels       = n_plot
    )

    horizons <- seq_along(med_vals)
    graphics::plot(horizons, med_vals, type = "n",
                   ylim = y_limits, xlab = "Horizon", ylab = "",
                   main = panel_title,
                   cex.lab = 1.2, cex.main = 1.2, cex.axis = 1.1)
    if (draw_ci) {
      graphics::polygon(
        x = c(horizons, rev(horizons)),
        y = c(lower_vals, rev(upper_vals)),
        col = ci_col,
        border = ci_border
      )
    }
    graphics::abline(h = 0, col = "gray55", lty = 2, lwd = 1)
    graphics::lines(horizons, med_vals, type = "o", pch = 1,
                    col = "black", lwd = 1.5)
    graphics::box()
  }

  invisible(NULL)
}


#' Plot historical decomposition as stacked bar charts
#'
#' For each variable, produces a stacked bar chart showing the contributions of
#' (1) the event-date target shock, (2) continued target shocks after the event,
#' and (3) all other shocks. The actual detrended data is overlaid as a black
#' line, and a vertical red line marks the event date.
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
#' @param n_tail           Integer. Number of trailing observations. Default 41.
#' @param x_break_every Integer. Show every k-th date label on the x-axis
#'   (base-R historical decomposition plot). Default 2.
#' @param title Character, `NULL`, or `NA`.
#'   - Character scalar: custom title for every panel.
#'   - Character vector: per-variable/panel titles. Named vectors are matched
#'     to `variable_names`.
#'   - `NA`: no panel title.
#'   - `NULL`: use `variable_names`.
#' @param merge_target_components Logical. If `TRUE`, merge `target_stop`
#'   (event shock) and `target_rest` (continued target shocks) into a single
#'   target-shock bar component. Default `FALSE`.
#' @param shock_labels     Character vector of length 3. Labels for the three
#'   bar components. Default: `c("Event shock", "Continued shocks",
#'   "Other shocks")`.
#' @param open_device Logical. Open a new graphics device per variable?
#'   Default `FALSE`. Set to `TRUE` to open separate windows.
#' @param ask Logical. If `TRUE` and multiple variables are plotted in
#'   the current device, prompt the user before advancing to the next plot.
#'   Default `interactive()`.
#' @return `NULL` (invisibly). Produces plots (one per variable).
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
                           ask = interactive()) {
  T_obs      <- nrow(data)
  n          <- ncol(data)
  bar_colors <- c("#4DBEEE", "#D95319", "#EDB120")

  start_idx  <- max(1, T_obs - n_tail + 1)
  plot_idx   <- start_idx:T_obs
  n_obs_plot <- length(plot_idx)

  if (open_device) {
    # Open a single new device (not per-variable)
    grDevices::dev.new(width = 10, height = 6, noRStudioGD = TRUE)
  }

  if (!open_device && n > 1L && ask) {
    old_ask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(old_ask), add = TRUE)
  }

  for (i in seq_len(n)) {
    graphics::par(mar = c(7, 5, 4, 2) + 0.1)

    if (merge_target_components) {
      raw_matrix <- t(cbind(target_stop[plot_idx, i] + target_rest[plot_idx, i],
                            others[plot_idx, i]))
      bar_labels <- c("Target shock", "Other shocks")
      if (length(shock_labels) >= 2L) {
        bar_labels <- c(shock_labels[1], shock_labels[length(shock_labels)])
      }
      bar_colors_i <- bar_colors[c(2, 3)]
    } else {
      raw_matrix <- t(cbind(target_stop[plot_idx, i],
                            target_rest[plot_idx, i],
                            others[plot_idx, i]))
      bar_labels <- shock_labels
      bar_colors_i <- bar_colors
    }

    data_pos <- raw_matrix; data_pos[data_pos < 0] <- 0
    data_neg <- raw_matrix; data_neg[data_neg > 0] <- 0

    line_values <- data[plot_idx, i] - deterministic[plot_idx, i]

    y_max   <- max(c(colSums(data_pos), line_values), na.rm = TRUE)
    y_min   <- min(c(colSums(data_neg), line_values), na.rm = TRUE)
    my_ylim <- c(y_min, y_max) * 1.2

    bp <- graphics::barplot(data_pos, col = bar_colors_i, border = NA,
                            space = 0.2, ylab = "", xlab = "", main = "",
                            ylim = my_ylim)
    graphics::barplot(data_neg, col = bar_colors_i, border = NA, space = 0.2,
                      axes = FALSE, add = TRUE)

    graphics::grid(nx = NA, ny = NULL)
    graphics::abline(h = 0, col = "black", lwd = 1)
    graphics::lines(x = bp, y = line_values, col = "black", lwd = 2.5)

    target_in_subset <- n_obs_plot - datestar
    if (target_in_subset > 0 && target_in_subset <= n_obs_plot) {
      graphics::abline(v = bp[target_in_subset], col = "red", lwd = 2.5)
    }

    graphics::legend("topright",
                     legend = c(bar_labels, "Data"),
                     fill   = c(bar_colors_i, NA),
                     border = c(rep("black", length(bar_labels)), NA),
                     lty    = c(rep(NA, length(bar_labels)), 1),
                     lwd    = c(rep(NA, length(bar_labels)), 2.5),
                     col    = c(rep(NA, length(bar_labels)), "black"),
                     cex = if (open_device) 0.9 else 0.55, bty = "n")

    panel_title <- etvar_select_panel_title(
      title          = title,
      default_titles = variable_names,
      variable_index = i,
      panel_index    = i,
      n_panels       = n
    )
    graphics::title(main = panel_title,
                    cex.main = if (open_device) 2.0 else 1.2)
    graphics::title(ylab = "Detrended values",
                    cex.lab = if (open_device) 1.2 else 0.8)

    date_labels <- format(date_vec[plot_idx], "%b %d, %Y")
    if (x_break_every > 1L) {
      keep_idx <- seq(1L, length(date_labels), by = x_break_every)
      date_labels[-keep_idx] <- NA
    }
    if (open_device) {
      graphics::axis(1, at = bp, labels = FALSE)
      graphics::text(x = bp,
                     y = graphics::par("usr")[3] -
                         0.05 * diff(graphics::par("usr")[3:4]),
                     labels = date_labels, srt = 45, adj = 1, xpd = TRUE,
                     cex = 0.8)
    } else {
      non_na <- which(!is.na(date_labels))
      graphics::axis(1, at = bp[non_na], labels = date_labels[non_na],
                     las = 2, cex.axis = 0.45)
    }
  }

  invisible(NULL)
}
