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
#'   draws. Default `NULL` (no seed set). The caller's random-number state is
#'   restored on exit.
#' @param control_nm   List. Control parameters for Nelder-Mead. See [optim()].
#' @param control_bfgs List. Control parameters for BFGS. See [optim()].
#' @param n_cores Integer. Number of CPU cores for parallel execution. Only
#'   used when multiple starts are generated (`start = NULL` and
#'   `n_starts > 1`). Default 1 (sequential); the maximum is 2 to respect CRAN
#'   resource limits. On Windows a PSOCK cluster is used; on Unix-like systems
#'   forking via [parallel::mclapply()] is used.
#' @param verbose Logical. Print progress messages (via [message()])?
#'   Default `FALSE`. When `TRUE`, sequential runs also enable [optim()]
#'   tracing when `trace` is absent or zero.
#' @return Named list with `par` (the `n - 1` optimal target-direction angles),
#'   `target_direction` (the corresponding unit vector), `rotation` (a
#'   deterministic orthogonal completion), `value` (best objective),
#'   `convergence` (0 indicates successful convergence), `stage` (the
#'   optimisation stage that produced the best value), `optimization` (full
#'   optimiser diagnostics), and `parameterization`.
#' @export
etvar_optimize <- function(start = NULL, var_est, target_obs,
                           shock_index = 1L,
                           n_starts = 1L, seed = NULL,
                           control_nm   = list(reltol = 1e-9, maxit = 20000,
                                               trace = 0),
                           control_bfgs = list(reltol = 1e-9, maxit = 4000,
                                               trace = 0, REPORT = 1),
                           n_cores = 1L,
                           verbose = FALSE,
                           start_parameterization = c("auto", "target", "full")) {
  optimization_start_time <- proc.time()[["elapsed"]]
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
    if (!is.null(seed)) {
      old_rng <- .etvar_get_rng()
      on.exit(.etvar_restore_rng(old_rng), add = TRUE)
      set.seed(seed)
    }
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

  if (verbose) message("Preparing event-targeting objective ...")
  objective_context <- .etvar_prepare_objective(
    var_est, target_obs, shock_index
  )

  convergence_message <- function(code, message = NULL) {
    if (!is.null(message) && length(message) > 0L && !is.na(message[1L])) {
      return(as.character(message[1L]))
    }
    switch(as.character(code),
      "0" = "Successful convergence.",
      "1" = "Iteration limit reached.",
      "10" = "Degeneracy of the Nelder-Mead simplex.",
      "51" = "Warning from L-BFGS-B.",
      "52" = "Error from L-BFGS-B.",
      paste0("optim returned convergence code ", code, ".")
    )
  }

  count_value <- function(counts, name) {
    if (is.null(counts) || is.null(names(counts)) ||
        !(name %in% names(counts))) {
      return(NA_integer_)
    }
    as.integer(unname(counts[[name]]))
  }

  stage_record <- function(start_index, stage, method, opt, elapsed_seconds) {
    data.frame(
      start = as.integer(start_index),
      stage = as.character(stage),
      method = as.character(method),
      objective = as.numeric(opt$value),
      convergence = as.integer(opt$convergence),
      message = convergence_message(opt$convergence, opt$message),
      fn_evaluations = count_value(opt$counts, "function"),
      gr_evaluations = count_value(opt$counts, "gradient"),
      elapsed_seconds = as.numeric(elapsed_seconds),
      stringsAsFactors = FALSE
    )
  }

  timed_optim <- function(par, fn, method, control) {
    stage_start_time <- proc.time()[["elapsed"]]
    opt <- stats::optim(par, fn, method = method, control = control)
    list(
      opt = opt,
      elapsed_seconds = proc.time()[["elapsed"]] - stage_start_time
    )
  }

  # --- Worker function: three-stage optimisation for a single start ---
  run_three_stage <- function(init, objective_context,
                              control_nm, control_bfgs,
                              start_index = 1L) {
    start_time <- proc.time()[["elapsed"]]
    obj_fn <- function(theta) {
      .etvar_objective_prepared(
        theta, objective_context, parameterization = "target"
      )
    }
    opt1_timed <- timed_optim(init, obj_fn, "Nelder-Mead", control_nm)
    opt1 <- opt1_timed$opt
    opt2_timed <- timed_optim(opt1$par, obj_fn, "BFGS", control_bfgs)
    opt2 <- opt2_timed$opt
    opt3_timed <- timed_optim(opt2$par, obj_fn, "Nelder-Mead", control_nm)
    opt3 <- opt3_timed$opt
    candidates <- list(stage1 = opt1, stage2 = opt2, stage3 = opt3)
    records <- rbind(
      stage_record(start_index, "stage1", "Nelder-Mead", opt1,
                   opt1_timed$elapsed_seconds),
      stage_record(start_index, "stage2", "BFGS", opt2,
                   opt2_timed$elapsed_seconds),
      stage_record(start_index, "stage3", "Nelder-Mead", opt3,
                   opt3_timed$elapsed_seconds)
    )
    values <- vapply(candidates, function(x) x$value, numeric(1))
    best_stage <- which.min(values)
    best <- candidates[[best_stage]]
    list(par = best$par, value = best$value,
         convergence = best$convergence,
         message = convergence_message(best$convergence, best$message),
         stage = names(candidates)[best_stage],
         start = as.integer(start_index),
         records = records,
         elapsed_seconds = proc.time()[["elapsed"]] - start_time)
  }

  if (n_par == 0L) {
    objective_start_time <- proc.time()[["elapsed"]]
    objective_value <- .etvar_objective_prepared(
      numeric(), objective_context, parameterization = "target"
    )
    elapsed_seconds <- proc.time()[["elapsed"]] - objective_start_time
    objective_by_start_stage <- data.frame(
      start = 1L,
      stage = "none",
      method = NA_character_,
      objective = as.numeric(objective_value),
      convergence = 0L,
      message = "No free rotation parameters.",
      fn_evaluations = 1L,
      gr_evaluations = 0L,
      elapsed_seconds = as.numeric(elapsed_seconds),
      stringsAsFactors = FALSE
    )
    total_elapsed <- proc.time()[["elapsed"]] - optimization_start_time
    rotation <- matrix(1, 1L, 1L)
    return(list(
      par = numeric(),
      target_direction = 1,
      rotation = rotation,
      value = objective_value,
      convergence = 0L,
      message = "No free rotation parameters.",
      stage = "none",
      start = 1L,
      optimization = list(
        convergence = 0L,
        message = "No free rotation parameters.",
        selected_stage = "none",
        selected_start = 1L,
        objective_by_start_stage = objective_by_start_stage,
        n_evaluations = 1L,
        n_gradient_evaluations = 0L,
        elapsed_seconds = as.numeric(total_elapsed),
        controls = list(nelder_mead = control_nm, bfgs = control_bfgs,
                        n_cores = n_cores, parallel = FALSE),
        seed = seed,
        n_starts = length(starts),
        start_values = matrix(numeric(), nrow = 1L, ncol = 0L)
      ),
      parameterization = "target_direction"
    ))
  }

  # --- Parallel multi-start ---
  use_parallel <- (length(starts) > 1L) && (n_cores > 1L)

  if (use_parallel) {
    n_cores <- min(n_cores, length(starts))
    if (verbose) {
      message(sprintf("Running %d starts in parallel on %d cores ...",
                  length(starts), n_cores))
    }

    # Suppress tracing in parallel workers (output would be garbled)
    ctrl_nm_quiet   <- control_nm
    ctrl_bfgs_quiet <- control_bfgs
    ctrl_nm_quiet$trace   <- 0
    ctrl_bfgs_quiet$trace <- 0

    if (.Platform$OS.type == "unix") {
      # Fork-based parallelism (not available on Windows)
      results <- parallel::mclapply(seq_along(starts), function(i) {
        init <- starts[[i]]
        run_three_stage(init, objective_context,
                        ctrl_nm_quiet, ctrl_bfgs_quiet, i)
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
        "ctrl_nm_quiet", "ctrl_bfgs_quiet", "starts"
      ), envir = environment())

      results <- parallel::parLapply(cl, seq_along(starts), function(i) {
        init <- starts[[i]]
        run_three_stage(init, objective_context,
                        ctrl_nm_quiet, ctrl_bfgs_quiet, i)
      })
    }

    # Find the best result
    values <- vapply(results, function(r) r$value, numeric(1))
    best_idx <- which.min(values)
    best <- results[[best_idx]]
    objective_by_start_stage <- do.call(rbind, lapply(results, `[[`, "records"))

    if (verbose) {
      message(sprintf("Best objective across %d starts: %g (start #%d)",
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

    best <- list(par = NULL, value = Inf, convergence = 1L,
                 message = NA_character_, stage = NA_character_,
                 start = NA_integer_)
    records <- vector("list", length(starts))

    for (s in seq_along(starts)) {
      init <- starts[[s]]
      if (verbose && length(starts) > 1L) {
        message(sprintf("=== Start %d/%d ===", s, length(starts)))
      }

      if (verbose) message("Stage 1/3: Nelder-Mead ...")
      opt1_timed <- timed_optim(init, obj_fn, "Nelder-Mead", control_nm)
      opt1 <- opt1_timed$opt
      if (verbose) message("  Objective: ", opt1$value)

      if (verbose) message("Stage 2/3: BFGS ...")
      opt2_timed <- timed_optim(opt1$par, obj_fn, "BFGS", control_bfgs)
      opt2 <- opt2_timed$opt
      if (verbose) message("  Objective: ", opt2$value)

      if (verbose) message("Stage 3/3: Nelder-Mead refinement ...")
      opt3_timed <- timed_optim(opt2$par, obj_fn, "Nelder-Mead", control_nm)
      opt3 <- opt3_timed$opt
      if (verbose) {
        message("  Final objective: ", opt3$value)
      }

      candidates <- list(stage1 = opt1, stage2 = opt2, stage3 = opt3)
      records[[s]] <- rbind(
        stage_record(s, "stage1", "Nelder-Mead", opt1,
                     opt1_timed$elapsed_seconds),
        stage_record(s, "stage2", "BFGS", opt2,
                     opt2_timed$elapsed_seconds),
        stage_record(s, "stage3", "Nelder-Mead", opt3,
                     opt3_timed$elapsed_seconds)
      )
      values <- vapply(candidates, function(x) x$value, numeric(1))
      best_stage <- which.min(values)
      candidate <- candidates[[best_stage]]

      if (verbose) {
        message("  Selected stage: ", best_stage)
        message("  Selected convergence: ",
            ifelse(candidate$convergence == 0, "Success", "Warning"))
      }

      if (candidate$value < best$value) {
        best <- list(par = candidate$par, value = candidate$value,
                     convergence = candidate$convergence,
                     message = convergence_message(
                       candidate$convergence, candidate$message
                     ),
                     stage = names(candidates)[best_stage],
                     start = as.integer(s))
      }
    }
    objective_by_start_stage <- do.call(rbind, records)

    if (verbose && length(starts) > 1L) {
      message(sprintf("Best objective across %d starts: %g",
                  length(starts), best$value))
    }
  }

  fn_evaluations <- objective_by_start_stage$fn_evaluations
  gr_evaluations <- objective_by_start_stage$gr_evaluations
  total_elapsed <- proc.time()[["elapsed"]] - optimization_start_time
  start_values <- do.call(rbind, starts)
  if (is.null(dim(start_values))) {
    start_values <- matrix(start_values, nrow = length(starts), byrow = TRUE)
  }
  colnames(start_values) <- paste0("angle", seq_len(ncol(start_values)))

  best$optimization <- list(
    convergence = as.integer(best$convergence),
    message = as.character(best$message),
    selected_stage = as.character(best$stage),
    selected_start = as.integer(best$start),
    objective_by_start_stage = objective_by_start_stage,
    n_evaluations = as.integer(sum(fn_evaluations, na.rm = TRUE)),
    n_gradient_evaluations = as.integer(sum(gr_evaluations, na.rm = TRUE)),
    elapsed_seconds = as.numeric(total_elapsed),
    controls = list(nelder_mead = control_nm, bfgs = control_bfgs,
                    n_cores = n_cores, parallel = use_parallel),
    seed = seed,
    n_starts = length(starts),
    start_values = start_values
  )

  if (!identical(as.integer(best$convergence), 0L)) {
    warning("ETVAR optimisation returned nonzero convergence code ",
            best$convergence, " at ", best$stage, " of start ",
            best$start, ": ", best$message, call. = FALSE)
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
