# ==============================================================================
# window_replication.R
#
# Demonstrates the ETVAR target-window interface, where the structural shock
# is identified by targeting movements across *multiple* event dates
# simultaneously rather than a single day.
#
# The objective function sums the scaled deviations over all event dates, so
# the optimiser finds the rotation that best explains movements on every
# date in the window at once.
#
# Requirement: install the etvar package before running this script.
# ==============================================================================

library(etvar)

# --- Load data ---
# Daily US market data, 2021-01-04 to 2025-05-15 (see ?etvar_example).
data("etvar_example", package = "etvar")

data   <- as.matrix(etvar_example[, -1])
dates  <- etvar_example$date

variable_names <- colnames(data)


# ==============================================================================
# Target window: three major tariff-related event dates in 2025
#   2025-04-03  "Liberation Day"   — broad tariffs announced
#   2025-04-09  "90-day pause"     — temporary tariff relief
#   2025-05-12  "China-US truce"   — tariff reduction agreement
# ==============================================================================

target_window <- as.Date(c("2025-04-03", "2025-04-09", "2025-05-12"))


# =========================================================================
# A) ONE-CALL INTERFACE
# =========================================================================
# Pass a vector to `target_date` to activate the target-window mode.

fit_window <- etvar(
  data           = data,
  dates          = dates,
  target_date    = target_window,      # <-- vector of dates
  lags           = 8,
  constant       = TRUE,
  trend          = 0,
  shock_index    = 1,
  n_starts       = 4,
  seed           = 1,
  nimp           = 30,
  variable_names = variable_names,
  verbose        = TRUE
)

# The print and summary methods show the full window
print(fit_window)
summary(fit_window)

# --- Plots (identical API to single-date fit) ---

plot(fit_window, type = "irfs", scaling_index = c(2, 1, 1))

plot(fit_window, type = "counterfactual", variables = 1, n_tail = 40)

plot(fit_window, type = "errors", n_tail = 100,
     event_dates  = target_window,
     event_labels = c("Liberation Day", "90-day pause", "China-US truce"))

plot(fit_window, type = "hd", n_tail = 41)


# =========================================================================
# B) STEP-BY-STEP INTERFACE
# =========================================================================

# 1. Estimate the reduced-form VAR
var_est <- estimate_var(data, lags = 8,
                        constant = TRUE,
                        trend    = 0)

# 2. Map all event dates to effective-sample indices
target_obs_raw <- match(target_window, dates)    # indices in raw data
target_eff     <- sort(target_obs_raw) - var_est$lags   # adjust for lags

# 3. Optimise over the full window (objective sums across all dates)
opt_window <- etvar_optimize(
  start      = fit_window$theta,
  var_est    = var_est,
  target_obs = target_eff,             # <-- vector of indices
  verbose    = TRUE
)

cat("Final objective (window):", opt_window$value, "\n")

# For comparison: single-date objective on Liberation Day only
opt_single <- etvar_optimize(
  start      = fit_window$theta,
  var_est    = var_est,
  target_obs = target_eff[1],          # single date
  verbose    = FALSE
)
cat("Final objective (single date):", opt_single$value, "\n")

# 4. Build structural impact matrix from window solution
n  <- var_est$n
Q  <- opt_window$rotation
A0 <- -t(chol(var_est$SIGMA)) %*% Q

# 5. Impulse response functions
irfs_window <- compute_structural_irfs(var_est, A0, nimp = 30)

# 6. Historical decomposition anchored at the *first* (earliest) event date
datestar <- length(dates) - target_obs_raw[1]
hd_window <- hist_decomposition_3d(
  A0       = A0,
  ALPHA    = var_est$ALPHA,
  Y        = var_est$Y,
  X        = var_est$X,
  n        = n,
  lags     = 8,
  trend    = 0,
  constant = example$constant,
  datestar = datestar
)

# 7. Rearrange shocks
shock_contributions_window <- rearrange_shocks(
  hd_window$hist_matrix, hd_window$hist_matrix_stop,
  nrow(data), n, target_shock = 1
)

# 8. Standalone plots
plot_etvar_irfs(irfs_window, variable_names, shock_index = 1,
                scaling_index = c(2, 1, 1))

plot_etvar_counterfactual(
  date_vec       = dates,
  realized       = data[, 1],
  counterfactual = hd_window$all_data[, 1, 1],
  n_tail         = 40,
  title          = variable_names[1]
)

plot_etvar_hd(
  date_vec       = dates,
  data           = data,
  deterministic  = hd_window$deterministic,
  target_stop    = shock_contributions_window$target_stop,
  target_rest    = shock_contributions_window$target_rest,
  others         = shock_contributions_window$others,
  datestar       = datestar,
  variable_names = variable_names,
  n_tail         = 41,
  shock_labels   = c("Liberation Day", "New tariff shocks", "Other shocks")
)
