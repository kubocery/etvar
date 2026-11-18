# ==============================================================================
# replication.R
#
# Full replication of the ETVAR Tariff Shock analysis using the etvar package.
#
# This script demonstrates both:
#   A) The one-call interface: etvar() + plot()
#   B) The step-by-step interface for fine-grained control
#
# Requirement: install the etvar package before running this script.
# ==============================================================================
library(etvar)

# --- Load data ---
# Daily US market data, 2021-01-04 to 2025-05-15 (see ?etvar_example).
data("etvar_example", package = "etvar")

data    <- as.matrix(etvar_example[, -1])
dates   <- data.frame(dates = etvar_example$date)

variable_names <- colnames(data)


# =========================================================================
# A) ONE-CALL INTERFACE
# =========================================================================
# The etvar() function runs the entire pipeline in one call.
# Supply `start_angles` with `optimize = FALSE` to skip the optimisation.

fit <- etvar(
  data           = data,
  dates          = dates$dates,
  target_date    = "2025-04-03",
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

# Print and summarise
print(fit)
summary(fit)

# --- Plots via the S3 method ---

# Impulse responses relative to the impact response of the VIX (variable 2),
# with 80% bootstrap bands
plot(fit, type = "irfs", scaling_index = c(2, 1, 1),
     conf_int = TRUE, conf_level = 0.8)

# Forecast error variance shares: target shock versus all other shocks
plot(fit, type = "fevd")

# Counterfactual vs realised for SP500 (variable 1). Omit `variables` to show
# every variable in one faceted plot.
plot(fit, type = "counterfactual", variables = 1, n_tail = 40)

# Target structural error with event markers.
plot(fit, type = "errors", n_tail = 100,
     event_dates  = as.Date(c("2025-04-03", "2025-04-09", "2025-05-12")),
     event_labels = c("Liberation Day", "90-day pause", "China-US truce"))

# Historical decomposition stacked bar charts (all variables; no paging)
plot(fit, type = "hd", n_tail = 41)


# =========================================================================
# B) STEP-BY-STEP INTERFACE
# =========================================================================
# For users who need fine-grained control over each step.

# 1. Estimate the reduced-form VAR
var_est <- estimate_var(data, lags = 8,
                        constant = TRUE,
                        trend = 0)

# 2. Find the target date in the effective sample
target_eff <- which(dates$dates == as.Date("2025-04-03")) - var_est$lags

# 3. Optimise the target-shock direction (3-stage: NM -> BFGS -> NM)
opt <- etvar_optimize(
  start      = fit$theta,   # warm start from the one-call fit
  var_est    = var_est,
  target_obs = target_eff,
  verbose    = TRUE
)

# 4. Build structural impact matrix. The optimiser returns a deterministic
#    orthogonal completion whose target column is the identified direction.
n  <- var_est$n
Q  <- opt$rotation
A0 <- -t(chol(var_est$SIGMA)) %*% Q

# 5. Impulse response functions
irfs <- compute_structural_irfs(var_est, A0, nimp = 30)

# 6. Historical decomposition with intervention path
datestar <- length(dates$dates) - which(dates$dates == as.Date("2025-04-03"))
hd <- hist_decomposition_3d(
  A0 = A0, ALPHA = var_est$ALPHA, Y = var_est$Y, X = var_est$X,
  n = n, lags = 8, trend = 0, constant = 1L,
  datestar = datestar
)

# 7. Rearrange shocks
shock_contributions <- rearrange_shocks(
  hd$hist_matrix, hd$hist_matrix_stop,
  nrow(data), n, target_shock = 1
)

# target_stop holds all target shocks up to and including the event day.
# Cutting the path one day earlier gives the target shocks before the event,
# and the difference is the event-day shock on its own.
hd_before <- hist_decomposition_3d(
  A0 = A0, ALPHA = var_est$ALPHA, Y = var_est$Y, X = var_est$X,
  n = n, lags = 8, trend = 0, constant = 1L,
  datestar = datestar + 1
)
target_before <- hd_before$hist_matrix_stop[, , 1]
target_event <- shock_contributions$target_stop - target_before

# 8. Standalone plotting functions
plot_etvar_irfs(irfs, variable_names, shock_index = 1,
                scaling_index = c(2, 1, 1))

plot_etvar_counterfactual(
  date_vec = dates$dates, realized = data[, 6],
  counterfactual = hd$all_data[, 6, 1], n_tail = 40,
  title = variable_names[6]
)

plot_etvar_hd(
  date_vec = dates$dates, data = data, deterministic = hd$deterministic,
  target_event = target_event,
  target_rest = shock_contributions$target_rest,
  others = shock_contributions$others + target_before,
  datestar = datestar, variable_names = variable_names, n_tail = 41,
  shock_labels = c("Liberation Day", "New tariff shocks", "Other shocks")
)
