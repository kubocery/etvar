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
setwd(r"(D:\PhD\clanky\cran\etvar)")
if (requireNamespace("devtools", quietly = TRUE)) devtools::load_all(".")
library(etvar)

# --- Load data ---
data_dir    <- system.file("extdata", package = "etvar")
example     <- readRDS(file.path(data_dir, "etvar_example.rds"))

data    <- example$data
dates   <- data.frame(dates = example$dates)

variable_names <- colnames(data)

# Pre-computed starting angles from the original MATLAB solution
theta_start <- readRDS(file.path(data_dir, "start_angles.rds"))


# =========================================================================
# A) ONE-CALL INTERFACE
# =========================================================================
# The etvar() function runs the entire pipeline in one call.
# Set optimize = FALSE to use the pre-computed angles directly.

fit <- etvar(
  data           = data,
  dates          = dates$dates,
  target_date    = "2025-04-03",
  lags           = 1,
  constant       = as.logical(example$constant),
  trend          = example$trend,
  shock_index    = 1,
  n_starts = 1,
 # start_angles   = theta_start,
  optimize       = TRUE,    # Set TRUE to re-run optimisation
  nimp           = 30,
  variable_names = variable_names,
  verbose        = TRUE
)

# Print and summarise
print(fit)
summary(fit)

# --- Plots via the S3 method ---

# Impulse response functions (normalised by VIX response at horizon 1)
plot(fit, type = "irfs", scaling_index = c(2, 1, 1), conf_int = TRUE, conf_level = 0.8)

# Forecast error variance shares: target shock versus all other shocks
plot(fit, type = "fevd")

# Counterfactual vs realised for 1Y Inflation Swap Rate (variable 6)
plot(fit, type = "counterfactual", variable = 1, n_tail = 40)

# Structural error series for shock 1 with event markers
plot(fit, type = "errors", shock = 1, n_tail = 100,
     event_dates  = as.Date(c("2025-04-03", "2025-04-09", "2025-05-12")),
     event_labels = c("Liberation Day", "90-day pause", "China-US truce"))

# Historical decomposition stacked bar charts (all variables)
plot(fit, type = "hd", n_tail = 41)


# =========================================================================
# B) STEP-BY-STEP INTERFACE
# =========================================================================
# For users who need fine-grained control over each step.

# 1. Estimate the reduced-form VAR
var_est <- estimate_var(data, lags = 8,
                        constant = as.logical(example$constant),
                        trend = example$trend)

# 2. Find the target date in the effective sample
target_eff <- which(dates$dates == as.Date("2025-04-03")) - var_est$lags

# 3. Optimise rotation angles (3-stage: NM -> BFGS -> NM)
opt <- etvar_optimize(
  start      = theta_start,
  var_est    = var_est,
  target_obs = target_eff,
  verbose    = TRUE
)

# 4. Build structural impact matrix
n  <- var_est$n
Q  <- euler_rotation(opt$par, n)
A0 <- -t(chol(var_est$SIGMA)) %*% Q

# 5. Impulse response functions
irfs <- compute_structural_irfs(var_est, A0, nimp = 30)

# 6. Historical decomposition with intervention path
datestar <- length(dates$dates) - which(dates$dates == as.Date("2025-04-03"))
hd <- hist_decomposition_3d(
  A0 = A0, ALPHA = var_est$ALPHA, Y = var_est$Y, X = var_est$X,
  n = n, lags = 8, trend = example$trend, constant = example$constant,
  datestar = datestar
)

# 7. Rearrange shocks
shocks <- rearrange_shocks(hd$hist_matrix, hd$hist_matrix_stop,
                           nrow(data), n, target_shock = 1)

# 8. Standalone plotting functions
plot_etvar_irfs(irfs, variable_names, shock_index = 1,
                scaling_index = c(2, 1, 1))

plot_etvar_counterfactual(
  date_vec = dates$dates, realized = data[, 6],
  counterfactual = hd$all_data[, 6, 1], n_tail = 40,
  title = "1Y Inflation Swap Rate"
)

plot_etvar_hd(
  date_vec = dates$dates, data = data, deterministic = hd$deterministic,
  target_stop = shocks$target_stop, target_rest = shocks$target_rest,
  others = shocks$others, datestar = datestar,
  variable_names = variable_names, n_tail = 41,
  shock_labels = c("Liberation Day", "New tariff shocks", "Other shocks")
)
