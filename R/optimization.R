# ==============================================================================
# Three-stage optimisation routine
# ==============================================================================

#' Optimise the ETVAR rotation angles
#'
#' Runs a three-stage optimisation alternating between Nelder-Mead (simplex)
#' and BFGS (gradient-based) to find the rotation angles that minimise the
#' [etvar_objective()].
#'
#' When `start = NULL` and `n_starts > 1`, a multi-start strategy is used:
#' angles are drawn uniformly from \eqn{[0, 2\pi)} and the optimiser is run
#' from each starting point (plus the zero vector). The solution with the
#' lowest objective value is returned. This is the recommended approach when
#' no pre-computed angles are available.
#'
#' @param start   Numeric vector of length n*(n-1)/2, or `NULL`. Initial
#'   rotation angles. When `NULL`, defaults to zeros (Cholesky identification)
#'   unless `n_starts > 1`, in which case random multi-start is used.
#' @param var_est An `etvar_var` object from [estimate_var()].
#' @param target_obs Integer or integer vector. Row index (indices) of the
#'   event date(s) in the effective sample. See [etvar_objective()].
#' @param shock_index Integer. Which structural shock to target. Default 1.
#' @param n_starts Integer. Number of random starting points. Only used when
#'   `start = NULL`. The zero vector is always included as one of the starts.
#'   Default 1 (no random starts, just zeros).
#' @param seed Integer or `NULL`. Random seed for reproducible multi-start
#'   draws. Default `NULL` (no seed set).
#' @param control_nm   List. Control parameters for Nelder-Mead. See [optim()].
#' @param control_bfgs List. Control parameters for BFGS. See [optim()].
#' @param n_cores Integer. Number of CPU cores for parallel execution. Only
#'   used when multiple starts are generated (`start = NULL` and
#'   `n_starts > 1`). Default 1 (sequential). On Windows a PSOCK cluster is
#'   used; on Unix-like systems forking via [parallel::mclapply()] is used.
#' @param verbose Logical. Print progress messages? Default `TRUE`.
#' @return Named list with `par` (optimal angles), `value` (final objective),
#'   and `convergence` (0 = success).
#' @importFrom parallel makeCluster stopCluster clusterExport
#'   parLapply mclapply
#' @export
etvar_optimize <- function(start = NULL, var_est, target_obs,
                           shock_index = 1L,
                           n_starts = 1L, seed = NULL,
                           control_nm   = list(reltol = 1e-9, maxit = 20000,
                                               trace = 0),
                           control_bfgs = list(reltol = 1e-9, maxit = 4000,
                                               trace = 0, REPORT = 1),
                           n_cores = 1L,
                           verbose = TRUE) {
  n      <- var_est$n
  n_ang  <- n * (n - 1) / 2

  # --- Generate starting point(s) ---
  if (is.null(start) && n_starts > 1L) {
    if (!is.null(seed)) set.seed(seed)
    starts <- vector("list", n_starts)
    starts[[1]] <- rep(0, n_ang)   # always include zeros (Cholesky)
    for (i in seq(2, n_starts)) {
      starts[[i]] <- stats::runif(n_ang, 0, 2 * pi)
    }
  } else {
    if (is.null(start)) start <- rep(0, n_ang)
    starts <- list(start)
  }

  # --- Worker function: three-stage optimisation for a single start ---
  run_three_stage <- function(init, var_est, target_obs, shock_index,
                              control_nm, control_bfgs) {
    obj_fn <- function(theta) {
      etvar_objective(theta, var_est, target_obs, shock_index)
    }
    opt1 <- stats::optim(init, obj_fn, method = "Nelder-Mead",
                         control = control_nm)
    opt2 <- stats::optim(opt1$par, obj_fn, method = "BFGS",
                         control = control_bfgs)
    opt3 <- stats::optim(opt2$par, obj_fn, method = "Nelder-Mead",
                         control = control_nm)
    list(par = opt3$par, value = opt3$value, convergence = opt3$convergence)
  }

  # --- Parallel multi-start ---
  use_parallel <- (length(starts) > 1L) && (n_cores > 1L)

  if (use_parallel) {
    n_cores <- min(n_cores, length(starts))
    if (verbose) {
      cat(sprintf("Running %d starts in parallel on %d cores ...\n",
                  length(starts), n_cores))
    }

    # Suppress tracing in parallel workers (output would be garbled)
    ctrl_nm_quiet   <- control_nm
    ctrl_bfgs_quiet <- control_bfgs
    ctrl_nm_quiet$trace   <- 0
    ctrl_bfgs_quiet$trace <- 0

    if (.Platform$OS.type == "unix") {
      # Fork-based parallelism (not available on Windows)
      results <- parallel::mclapply(starts, function(init) {
        run_three_stage(init, var_est, target_obs, shock_index,
                        ctrl_nm_quiet, ctrl_bfgs_quiet)
      }, mc.cores = n_cores, mc.set.seed = TRUE)
    } else {
      # PSOCK cluster (Windows-compatible)
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)

      # Load the full etvar namespace on each worker. This works whether
      # the package is formally installed or loaded via devtools::load_all().
      ns <- asNamespace("etvar")
      all_fns <- ls(ns, all.names = TRUE)
      parallel::clusterExport(cl, varlist = all_fns, envir = ns)

      # Export objects needed by the worker function
      parallel::clusterExport(cl, varlist = c(
        "run_three_stage", "var_est", "target_obs", "shock_index",
        "ctrl_nm_quiet", "ctrl_bfgs_quiet"
      ), envir = environment())

      results <- parallel::parLapply(cl, starts, function(init) {
        run_three_stage(init, var_est, target_obs, shock_index,
                        ctrl_nm_quiet, ctrl_bfgs_quiet)
      })
    }

    # Find the best result
    values <- vapply(results, function(r) r$value, numeric(1))
    best_idx <- which.min(values)
    best <- results[[best_idx]]

    if (verbose) {
      cat(sprintf("Best objective across %d starts: %g (start #%d)\n",
                  length(starts), best$value, best_idx))
    }

  } else {
    # --- Sequential execution (original behaviour) ---

    # When verbose, enable optim() tracing unless the user explicitly set trace
    if (verbose) {
      if (is.null(control_nm$trace)   || control_nm$trace   == 0)
        control_nm$trace   <- 1
      if (is.null(control_bfgs$trace) || control_bfgs$trace == 0)
        control_bfgs$trace <- 1
      if (is.null(control_bfgs$REPORT))
        control_bfgs$REPORT <- 1
    }

    obj_fn <- function(theta) {
      etvar_objective(theta, var_est, target_obs, shock_index)
    }

    best <- list(par = NULL, value = Inf, convergence = 1L)

    for (s in seq_along(starts)) {
      init <- starts[[s]]
      if (verbose && length(starts) > 1L) {
        cat(sprintf("=== Start %d/%d ===\n", s, length(starts)))
      }

      if (verbose) cat("Stage 1/3: Nelder-Mead ...\n")
      opt1 <- stats::optim(init, obj_fn, method = "Nelder-Mead",
                            control = control_nm)
      if (verbose) cat("  Objective:", opt1$value, "\n")

      if (verbose) cat("Stage 2/3: BFGS ...\n")
      opt2 <- stats::optim(opt1$par, obj_fn, method = "BFGS",
                            control = control_bfgs)
      if (verbose) cat("  Objective:", opt2$value, "\n")

      if (verbose) cat("Stage 3/3: Nelder-Mead refinement ...\n")
      opt3 <- stats::optim(opt2$par, obj_fn, method = "Nelder-Mead",
                            control = control_nm)
      if (verbose) {
        cat("  Final objective:", opt3$value, "\n")
        cat("  Convergence:",
            ifelse(opt3$convergence == 0, "Success", "Warning"), "\n")
      }

      if (opt3$value < best$value) {
        best <- list(par = opt3$par, value = opt3$value,
                     convergence = opt3$convergence)
      }
    }

    if (verbose && length(starts) > 1L) {
      cat(sprintf("Best objective across %d starts: %g\n",
                  length(starts), best$value))
    }
  }

  best
}
