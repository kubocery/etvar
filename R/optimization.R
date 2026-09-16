# ==============================================================================
# Three-stage optimisation routine
# ==============================================================================

#' Optimise the ETVAR target-shock direction
#'
#' Runs a three-stage optimisation alternating between Nelder-Mead (simplex)
#' and BFGS (gradient-based) to find the target-shock direction that minimises
#' the [etvar_objective()]. Only `n - 1` target-direction angles are estimated;
#' rotations among the non-target shocks are not identified by the criterion
#' and are therefore not optimised.
#'
#' When `start = NULL` and `n_starts > 1`, a multi-start strategy is used:
#' target directions are drawn uniformly on the unit sphere and the optimiser
#' is run from each starting point (plus the Cholesky target direction). The
#' solution with the lowest objective value is returned. This is the
#' recommended approach when no pre-computed direction is available.
#'
#' The returned non-target columns are a deterministic orthogonal completion.
#' They are required for a square structural system, but they are not
#' separately identified by the event-targeting objective.
#'
#' @param start Numeric vector or `NULL`. Initial target-direction angles of
#'   length `n - 1`. Legacy full Euler angles of length `n*(n-1)/2` are also
#'   accepted; their target column is extracted before optimisation. When
#'   `NULL`, the Cholesky target direction is used as the first start.
#' @param start_parameterization Character. Interpretation of `start`:
#'   `"target"`, `"full"`, or `"auto"` to infer from its length. In the
#'   ambiguous two-variable case, `"auto"` selects `"target"`.
#' @param var_est An `etvar_var` object from [estimate_var()].
#' @param target_obs Integer or integer vector. Row index (indices) of the
#'   event date(s) in the effective sample. See [etvar_objective()].
#' @param shock_index Integer. Which structural shock to target. Default 1.
#' @param n_starts Integer. Total number of starting points. Only used when
#'   `start = NULL`. The Cholesky target direction is included, with the
#'   remaining directions drawn uniformly on the unit sphere. Default 1.
#' @param seed Integer or `NULL`. Random seed for reproducible multi-start
#'   draws. Default `NULL` (no seed set).
#' @param control_nm   List. Control parameters for Nelder-Mead. See [optim()].
#' @param control_bfgs List. Control parameters for BFGS. See [optim()].
#' @param n_cores Integer. Number of CPU cores for parallel execution. Only
#'   used when multiple starts are generated (`start = NULL` and
#'   `n_starts > 1`). Default 1 (sequential); the maximum is 2 to respect CRAN
#'   resource limits. On Windows a PSOCK cluster is used; on Unix-like systems
#'   forking via [parallel::mclapply()] is used.
#' @param verbose Logical. Print progress messages? Default `TRUE`. For
#'   sequential runs, this also enables [optim()] tracing when `trace` is
#'   absent or zero. Set `verbose = FALSE` to suppress optimiser output.
#' @return Named list with `par` (the `n - 1` optimal target-direction angles),
#'   `target_direction` (the corresponding unit vector), `rotation` (a
#'   deterministic orthogonal completion), `value` (best objective),
#'   `convergence` (0 indicates successful convergence), `stage` (the
#'   optimisation stage that produced the best value), and `parameterization`.
#' @export
etvar_optimize <- function(start = NULL, var_est, target_obs,
                           shock_index = 1L,
                           n_starts = 1L, seed = NULL,
                           control_nm   = list(reltol = 1e-9, maxit = 20000,
                                               trace = 0),
                           control_bfgs = list(reltol = 1e-9, maxit = 4000,
                                               trace = 0, REPORT = 1),
                           n_cores = 1L,
                           verbose = TRUE,
                           start_parameterization = c("auto", "target", "full")) {
  if (!inherits(var_est, "etvar_var")) {
    stop("`var_est` must be an `etvar_var` object from `estimate_var()`.",
         call. = FALSE)
  }
  n      <- var_est$n
  n_par  <- n - 1L
  shock_index <- .etvar_integer_scalar(shock_index, "shock_index",
                                        min = 1L, max = n)
  start_parameterization <- match.arg(start_parameterization)
  target_obs <- .etvar_index_vector(target_obs, "target_obs", nrow(var_est$Y))
  n_starts <- .etvar_integer_scalar(n_starts, "n_starts", min = 1L)
  n_cores <- .etvar_integer_scalar(n_cores, "n_cores", min = 1L, max = 2L)
  verbose <- .etvar_flag(verbose, "verbose")
  if (!is.list(control_nm) || !is.list(control_bfgs)) {
    stop("`control_nm` and `control_bfgs` must be lists.", call. = FALSE)
  }
  if (!is.null(seed)) {
    seed <- .etvar_integer_scalar(seed, "seed", min = 0L,
                                  max = .Machine$integer.max)
  }
  if (!is.null(start)) {
    start <- .etvar_prepare_target_start(
      start, n, shock_index, start_parameterization
    )
  }

  # --- Generate starting point(s) ---
  cholesky_direction <- numeric(n)
  cholesky_direction[shock_index] <- 1
  cholesky_start <- .etvar_angles_from_direction(
    cholesky_direction, shock_index
  )

  if (is.null(start) && n_starts > 1L) {
    if (!is.null(seed)) set.seed(seed)
    starts <- vector("list", n_starts)
    starts[[1]] <- cholesky_start
    for (i in seq(2, n_starts)) {
      random_direction <- stats::rnorm(n)
      random_direction <- random_direction / sqrt(sum(random_direction^2))
      starts[[i]] <- .etvar_angles_from_direction(
        random_direction, shock_index
      )
    }
  } else {
    if (is.null(start)) start <- cholesky_start
    starts <- list(start)
  }

  if (verbose) cat("Preparing event-targeting objective ...\n")
  objective_context <- .etvar_prepare_objective(
    var_est, target_obs, shock_index
  )

  # --- Worker function: three-stage optimisation for a single start ---
  run_three_stage <- function(init, objective_context,
                              control_nm, control_bfgs) {
    obj_fn <- function(theta) {
      .etvar_objective_prepared(
        theta, objective_context, parameterization = "target"
      )
    }
    opt1 <- stats::optim(init, obj_fn, method = "Nelder-Mead",
                         control = control_nm)
    opt2 <- stats::optim(opt1$par, obj_fn, method = "BFGS",
                         control = control_bfgs)
    opt3 <- stats::optim(opt2$par, obj_fn, method = "Nelder-Mead",
                         control = control_nm)
    candidates <- list(stage1 = opt1, stage2 = opt2, stage3 = opt3)
    values <- vapply(candidates, function(x) x$value, numeric(1))
    best_stage <- which.min(values)
    best <- candidates[[best_stage]]
    list(par = best$par, value = best$value,
         convergence = best$convergence, stage = names(candidates)[best_stage])
  }

  if (n_par == 0L) {
    rotation <- matrix(1, 1L, 1L)
    return(list(
      par = numeric(),
      target_direction = 1,
      rotation = rotation,
      value = .etvar_objective_prepared(
        numeric(), objective_context, parameterization = "target"
      ),
      convergence = 0L,
      stage = "none",
      parameterization = "target_direction"
    ))
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
        run_three_stage(init, objective_context,
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
        "run_three_stage", "objective_context",
        "ctrl_nm_quiet", "ctrl_bfgs_quiet"
      ), envir = environment())

      results <- parallel::parLapply(cl, starts, function(init) {
        run_three_stage(init, objective_context,
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
    # --- Sequential execution ---

    # Preserve the original package behaviour: verbose sequential runs show
    # the objective values emitted by optim() while each stage is running.
    if (verbose) {
      if (is.null(control_nm$trace) || control_nm$trace == 0) {
        control_nm$trace <- 1
      }
      if (is.null(control_bfgs$trace) || control_bfgs$trace == 0) {
        control_bfgs$trace <- 1
      }
      if (is.null(control_bfgs$REPORT)) {
        control_bfgs$REPORT <- 1
      }
    }

    obj_fn <- function(theta) {
      .etvar_objective_prepared(
        theta, objective_context, parameterization = "target"
      )
    }

    best <- list(par = NULL, value = Inf, convergence = 1L, stage = NA_character_)

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
      }

      candidates <- list(stage1 = opt1, stage2 = opt2, stage3 = opt3)
      values <- vapply(candidates, function(x) x$value, numeric(1))
      best_stage <- which.min(values)
      candidate <- candidates[[best_stage]]

      if (verbose) {
        cat("  Selected stage:", best_stage, "\n")
        cat("  Selected convergence:",
            ifelse(candidate$convergence == 0, "Success", "Warning"), "\n")
      }

      if (candidate$value < best$value) {
        best <- list(par = candidate$par, value = candidate$value,
                     convergence = candidate$convergence,
                     stage = names(candidates)[best_stage])
      }
    }

    if (verbose && length(starts) > 1L) {
      cat(sprintf("Best objective across %d starts: %g\n",
                  length(starts), best$value))
    }
  }

  best$target_direction <- .etvar_direction_from_angles(
    best$par, n, shock_index
  )
  best$rotation <- .etvar_complete_target_rotation(
    best$target_direction, shock_index
  )
  best$parameterization <- "target_direction"
  best
}
