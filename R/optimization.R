#' Optimise the ETVAR target-shock direction
#'
#' Minimises [etvar_objective()] over the `n - 1` angles of the target
#' direction. Rotations among the other shocks do not change the objective,
#' so they are not estimated; the returned `rotation` completes the target
#' direction in a fixed, deterministic way.
#'
#' With two variables there is a single angle and the objective repeats
#' every pi radians, so it is minimised globally: a one-degree grid over all
#' directions, refined with Brent's method. With more variables, each starting
#' point goes through Nelder-Mead, then BFGS from that solution, then
#' Nelder-Mead again, and the best of the three stages is kept. When
#' `start = NULL` and `n_starts > 1`, the Cholesky direction is used together
#' with directions drawn uniformly on the unit sphere, and the best start
#' wins.
#'
#' The objective cannot tell a direction from its negative, which only flips
#' the sign of the target shock. The result is signed so that the target
#' shock is positive on the first event date.
#'
#' @param start Numeric vector or `NULL`. Starting angles: `n - 1`
#'   target-direction angles or the `n * (n - 1) / 2` angles of a full
#'   rotation (see `start_parameterization`). `NULL` starts from the Cholesky
#'   direction.
#' @param var_est An `etvar_var` object from [estimate_var()].
#' @param target_obs Integer vector. Effective-sample row(s) of the event
#'   date(s); see [etvar_objective()].
#' @param shock_index Integer. Which structural shock to target. Default 1.
#' @param n_starts Integer. Number of starting points when `start = NULL`.
#'   Default 1.
#' @param seed Integer or `NULL`. Seed for the random starting points. The
#'   caller's random-number state is restored afterwards.
#' @param control_nm,control_bfgs Lists of [optim()] controls for the
#'   Nelder-Mead and BFGS stages.
#' @param n_cores Integer. Number of cores for running several starts in
#'   parallel. Windows uses a PSOCK cluster, other systems fork with
#'   [parallel::mclapply()]. Default 1.
#' @param verbose Logical. Print progress messages? Sequential runs then also
#'   show the [optim()] trace, unless the controls set `trace` themselves.
#'   Default `FALSE`.
#' @param start_parameterization Character. How to read `start`: `"target"`,
#'   `"full"`, or `"auto"` to decide from its length (with two variables
#'   `"auto"` reads it as `"target"`).
#' @return A list with
#'   \describe{
#'     \item{par}{The `n - 1` optimal angles.}
#'     \item{target_direction}{The corresponding unit vector.}
#'     \item{rotation}{An n x n rotation whose target column is
#'       `target_direction`.}
#'     \item{value}{The minimised objective.}
#'     \item{convergence, message}{[optim()] convergence code and message of
#'       the selected stage (0 means success).}
#'     \item{stage, start}{The stage and the starting point that produced
#'       the solution.}
#'     \item{optimization}{Diagnostics: objective values, convergence codes
#'       and evaluation counts by start and stage, elapsed time, controls,
#'       seed and starting values.}
#'   }
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3)
#' var_est <- estimate_var(y, lags = 2, nimp_wold = 10)
#' opt <- etvar_optimize(var_est = var_est, target_obs = 78,
#'                       n_starts = 3, seed = 1)
#' opt$value
#' opt$target_direction
#' @export
etvar_optimize <- function(start = NULL, var_est, target_obs,
                           shock_index = 1L,
                           n_starts = 1L, seed = NULL,
                           control_nm = list(reltol = 1e-9, maxit = 20000,
                                             trace = 0),
                           control_bfgs = list(reltol = 1e-9, maxit = 4000,
                                               trace = 0, REPORT = 1),
                           n_cores = 1L,
                           verbose = FALSE,
                           start_parameterization = c("auto", "target",
                                                      "full")) {
  started <- proc.time()[["elapsed"]]
  if (!inherits(var_est, "etvar_var")) {
    stop("`var_est` must be an `etvar_var` object from `estimate_var()`.",
         call. = FALSE)
  }
  n <- var_est$n
  shock_index <- .etvar_integer_scalar(shock_index, "shock_index",
                                       min = 1L, max = n)
  start_parameterization <- match.arg(start_parameterization)
  target_obs <- .etvar_index_vector(target_obs, "target_obs", nrow(var_est$Y))
  n_starts <- .etvar_integer_scalar(n_starts, "n_starts", min = 1L)
  n_cores <- .etvar_integer_scalar(n_cores, "n_cores", min = 1L)
  verbose <- .etvar_flag(verbose, "verbose")
  if (!is.list(control_nm) || !is.list(control_bfgs)) {
    stop("`control_nm` and `control_bfgs` must be lists.", call. = FALSE)
  }
  if (!is.null(seed)) {
    seed <- .etvar_integer_scalar(seed, "seed", min = 0L,
                                  max = .Machine$integer.max)
  }

  if (!is.null(start)) {
    starts <- list(.etvar_prepare_target_start(
      start, n, shock_index, start_parameterization
    ))
  } else {
    starts <- list(cholesky_angles(n, shock_index))
    if (n_starts > 1L) {
      if (!is.null(seed)) {
        old_rng <- .etvar_get_rng()
        on.exit(.etvar_restore_rng(old_rng), add = TRUE)
        set.seed(seed)
      }
      for (i in seq.int(2L, n_starts)) {
        starts[[i]] <- .etvar_angles_from_direction(stats::rnorm(n),
                                                    shock_index)
      }
    }
  }

  if (verbose) message("Preparing event-targeting objective ...")
  context <- .etvar_prepare_objective(var_est, target_obs, shock_index)
  objective <- function(theta) {
    .etvar_objective_prepared(theta, context, parameterization = "target")
  }
  controls <- list(nelder_mead = control_nm, bfgs = control_bfgs,
                   n_cores = n_cores, parallel = FALSE)

  if (n == 1L) {
    # Nothing to rotate: the only direction is the Cholesky shock itself.
    value <- objective(numeric())
    record <- stage_record(1L, "none", NA_character_,
                           list(value = value, convergence = 0L,
                                message = "No free rotation parameters.",
                                counts = c(`function` = 1L)),
                           proc.time()[["elapsed"]] - started)
    return(list(
      par = numeric(),
      target_direction = 1,
      rotation = matrix(1, 1L, 1L),
      value = value,
      convergence = 0L,
      message = "No free rotation parameters.",
      stage = "none",
      start = 1L,
      optimization = optimization_summary(
        record, record, controls, seed, starts, started
      )
    ))
  }

  use_parallel <- length(starts) > 1L && n_cores > 1L
  if (use_parallel) {
    control_nm$trace <- 0
    control_bfgs$trace <- 0
    controls$parallel <- TRUE
    if (verbose) {
      message(sprintf("Running %d starts in parallel on %d cores ...",
                      length(starts), min(n_cores, length(starts))))
    }
    worker <- function(i) {
      optimise_from(starts[[i]], i, objective, control_nm, control_bfgs)
    }
    results <- run_in_parallel(seq_along(starts), worker, n_cores)
  } else {
    if (verbose) {
      # Show the optim() trace unless the user configured it.
      if (is.null(control_nm$trace) || control_nm$trace == 0) {
        control_nm$trace <- 1
      }
      if (is.null(control_bfgs$trace) || control_bfgs$trace == 0) {
        control_bfgs$trace <- 1
      }
      if (is.null(control_bfgs$REPORT)) control_bfgs$REPORT <- 1
    }
    results <- lapply(seq_along(starts), function(i) {
      if (verbose && length(starts) > 1L) {
        message(sprintf("=== Start %d/%d ===", i, length(starts)))
      }
      optimise_from(starts[[i]], i, objective, control_nm, control_bfgs,
                    verbose = verbose)
    })
  }
  controls$nelder_mead <- control_nm
  controls$bfgs <- control_bfgs

  best <- results[[which.min(vapply(results, `[[`, numeric(1), "value"))]]
  if (verbose && length(starts) > 1L) {
    message(sprintf("Best objective across %d starts: %g (start %d)",
                    length(starts), best$value, best$start))
  }

  direction <- .etvar_direction_from_angles(best$par, n, shock_index)
  if (target_shock_sign(direction, var_est, min(target_obs)) < 0) {
    direction <- -direction
  }
  # optim() may return any of the equivalent angle vectors (e.g. shifted by
  # 2 pi). Recomputing them from the direction makes the completed rotation,
  # and so the non-target shocks, depend on the direction alone.
  par <- .etvar_angles_from_direction(direction, shock_index)

  records <- do.call(rbind, lapply(results, `[[`, "records"))
  selected <- records$start == best$start & records$stage == best$stage
  if (!identical(as.integer(best$convergence), 0L)) {
    warning("ETVAR optimisation returned nonzero convergence code ",
            best$convergence, " at ", best$stage, " of start ", best$start,
            ": ", best$message, call. = FALSE)
  }

  list(
    par = par,
    target_direction = direction,
    rotation = .etvar_target_rotation(par, n, shock_index),
    value = best$value,
    convergence = as.integer(best$convergence),
    message = best$message,
    stage = best$stage,
    start = best$start,
    optimization = optimization_summary(
      records, records[selected, , drop = FALSE], controls, seed, starts,
      started
    )
  )
}


# Optimisation helpers ---------------------------------------------------------

cholesky_angles <- function(n, shock_index) {
  direction <- numeric(n)
  direction[shock_index] <- 1
  .etvar_angles_from_direction(direction, shock_index)
}

# Sign of the target shock on an event date. With A0 = -L Q the structural
# error of the target shock is -q'z, where z = L^{-1} u is the whitened
# residual.
target_shock_sign <- function(direction, var_est, event_row) {
  z <- forwardsolve(t(chol(var_est$SIGMA)), var_est$residuals[event_row, ])
  -sign(sum(direction * z))
}

# Optimises from one starting point and returns the best stage together with
# a diagnostics row for every stage.
optimise_from <- function(init, start_index, objective, control_nm,
                          control_bfgs, verbose = FALSE) {
  stages <- if (length(init) == 1L) {
    one_angle_search(init, objective, verbose)
  } else {
    three_stage_search(init, objective, control_nm, control_bfgs, verbose)
  }

  values <- vapply(stages, function(s) s$opt$value, numeric(1))
  best <- stages[[which.min(values)]]
  if (verbose) {
    message("  Selected stage: ", best$stage, " (",
            if (best$opt$convergence == 0) "converged" else "not converged",
            ")")
  }
  list(
    par = best$opt$par,
    value = best$opt$value,
    convergence = best$opt$convergence,
    message = convergence_message(best$opt$convergence, best$opt$message),
    stage = best$stage,
    start = start_index,
    records = do.call(rbind, lapply(stages, function(s) {
      stage_record(start_index, s$stage, s$method, s$opt, s$seconds)
    }))
  )
}

three_stage_search <- function(init, objective, control_nm, control_bfgs,
                               verbose) {
  run <- function(stage, label, par, method, control) {
    if (verbose) message(label)
    result <- timed_optim(stage, method, par, objective, method = method,
                          control = control)
    if (verbose) message("  Objective: ", result$opt$value)
    result
  }
  nm1 <- run("stage1", "Stage 1/3: Nelder-Mead ...", init, "Nelder-Mead",
             control_nm)
  bfgs <- run("stage2", "Stage 2/3: BFGS ...", nm1$opt$par, "BFGS",
              control_bfgs)
  nm2 <- run("stage3", "Stage 3/3: Nelder-Mead refinement ...",
             bfgs$opt$par, "Nelder-Mead", control_nm)
  list(nm1, bfgs, nm2)
}

# One angle (two variables). The objective depends on the direction only
# through q q', so it has period pi in the angle: a one-degree grid over one
# period cannot miss the global basin, and Brent's method polishes the best
# grid point. Nelder-Mead is unreliable in one dimension and gets stuck in
# local minima here.
one_angle_search <- function(init, objective, verbose) {
  if (verbose) message("Grid search over all directions ...")
  started <- proc.time()[["elapsed"]]
  grid <- init + (-90:89) * pi / 180
  values <- vapply(grid, objective, numeric(1))
  best <- which.min(values)
  grid_stage <- list(
    stage = "grid",
    method = "grid search",
    opt = list(par = grid[best], value = values[best], convergence = 0L,
               message = NULL, counts = c(`function` = length(grid))),
    seconds = proc.time()[["elapsed"]] - started
  )

  if (verbose) message("Brent refinement ...")
  step <- pi / 180
  brent_stage <- timed_optim("brent", "Brent", grid[best], objective,
                             method = "Brent", lower = grid[best] - step,
                             upper = grid[best] + step)
  if (verbose) message("  Objective: ", brent_stage$opt$value)
  list(grid_stage, brent_stage)
}

# Arguments after `algorithm` go to optim(). Avoid formals that a name such as
# `method` or `lower` could partially match.
timed_optim <- function(stage, algorithm, ...) {
  started <- proc.time()[["elapsed"]]
  opt <- stats::optim(...)
  list(stage = stage, method = algorithm, opt = opt,
       seconds = proc.time()[["elapsed"]] - started)
}

run_in_parallel <- function(indices, worker, n_cores) {
  n_cores <- min(n_cores, length(indices))
  if (.Platform$OS.type == "unix") {
    results <- parallel::mclapply(indices, worker, mc.cores = n_cores,
                                  mc.set.seed = TRUE)
  } else {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # The workers need the package's internal functions. Exporting the whole
    # namespace also works under devtools::load_all(), when there is no
    # installed copy of the package for the workers to load.
    ns <- asNamespace("etvar")
    helpers <- Filter(function(name) is.function(get(name, envir = ns)),
                      ls(ns, all.names = TRUE))
    parallel::clusterExport(cl, helpers, envir = ns)
    results <- parallel::parLapply(cl, indices, worker)
  }
  failed <- vapply(results, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop("A parallel optimisation worker failed: ",
         conditionMessage(attr(results[[which(failed)[1L]]], "condition")),
         call. = FALSE)
  }
  results
}

convergence_message <- function(code, message = NULL) {
  if (length(message) > 0L && !is.na(message[1L])) {
    return(as.character(message[1L]))
  }
  switch(as.character(code),
    "0" = "Successful convergence.",
    "1" = "Iteration limit reached.",
    "10" = "Degeneracy of the Nelder-Mead simplex.",
    paste0("optim() returned convergence code ", code, ".")
  )
}

stage_record <- function(start, stage, method, opt, seconds) {
  count <- function(name) {
    if (name %in% names(opt$counts)) {
      as.integer(opt$counts[[name]])
    } else {
      NA_integer_
    }
  }
  data.frame(
    start = as.integer(start),
    stage = stage,
    method = method,
    objective = as.numeric(opt$value),
    convergence = as.integer(opt$convergence),
    message = convergence_message(opt$convergence, opt$message),
    fn_evaluations = count("function"),
    gr_evaluations = count("gradient"),
    elapsed_seconds = as.numeric(seconds),
    stringsAsFactors = FALSE
  )
}

optimization_summary <- function(records, selected, controls, seed, starts,
                                 started) {
  start_values <- do.call(rbind, lapply(starts, function(s) {
    matrix(s, nrow = 1L)
  }))
  # sprintf(), unlike paste0(), gives character(0) when there are no angles.
  colnames(start_values) <- sprintf("angle%d", seq_len(ncol(start_values)))
  list(
    convergence = selected$convergence[1L],
    message = selected$message[1L],
    selected_stage = selected$stage[1L],
    selected_start = selected$start[1L],
    objective_by_start_stage = records,
    n_evaluations = as.integer(sum(records$fn_evaluations, na.rm = TRUE)),
    n_gradient_evaluations = as.integer(sum(records$gr_evaluations,
                                            na.rm = TRUE)),
    elapsed_seconds = as.numeric(proc.time()[["elapsed"]] - started),
    controls = controls,
    seed = seed,
    n_starts = length(starts),
    start_values = start_values
  )
}
