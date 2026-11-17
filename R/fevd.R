# ==============================================================================
# Structural forecast error variance decomposition
# ==============================================================================

#' Compute an ETVAR forecast error variance decomposition
#'
#' Computes the structural forecast error variance decomposition (FEVD) implied
#' by an estimated ETVAR model. Because the ETVAR impact matrix identifies
#' mutually orthogonal, unit-variance structural shocks, the contribution of
#' shock \eqn{s} to variable \eqn{i} at forecast horizon \eqn{H} is
#'
#' \deqn{
#'   FEVD_{i,s}(H) =
#'   \frac{\sum_{h=0}^{H-1} IRF_{i,s}(h)^2}
#'        {\sum_{j=1}^{n}\sum_{h=0}^{H-1} IRF_{i,j}(h)^2}.
#' }
#'
#' The function returns both the individual structural-shock shares and a
#' two-component aggregation consisting of the selected target shock and all
#' other structural shocks together. Shares are expressed as proportions
#' between zero and one.
#'
#' This is a structural FEVD, not a Pesaran--Shin generalized FEVD. A
#' generalized FEVD is based on correlated reduced-form innovations, whose
#' indices do not correspond to the event-targeted structural shock.
#'
#' ETVAR identifies the target-shock direction but not a unique basis within
#' the orthogonal non-target subspace. Thus the target share and the aggregate
#' `others` share are invariant and identified by the event criterion, whereas
#' individual non-target entries in `shares` depend on the deterministic
#' orthogonal completion stored in the fitted object.
#'
#' @param object An `etvar` object.
#' @param nimp Integer or `NULL`. Number of forecast horizons. Horizon 1 uses
#'   the impact response only; horizon \eqn{H} accumulates IRFs from 0 through
#'   \eqn{H-1}. The default uses every horizon stored in `object$irfs`.
#' @param target_shock Integer or `NULL`. Structural shock to report separately.
#'   The default uses `object$shock_index`.
#' @return An object of class `"etvar_fevd"` containing:
#'   \describe{
#'     \item{shares}{Array (variable x shock x horizon) with the FEVD share of
#'       every structural shock.}
#'     \item{target}{Matrix (variable x horizon) containing the selected
#'       target-shock share.}
#'     \item{others}{Matrix (variable x horizon) containing the sum of all
#'       non-target structural-shock shares.}
#'     \item{horizons}{Integer forecast horizons, beginning at 1.}
#'     \item{target_shock}{The separately reported structural shock.}
#'     \item{variable_names}{Variable labels inherited from `object`.}
#'     \item{method}{The string `"structural"`.}
#'   }
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(300), ncol = 3,
#'             dimnames = list(NULL, c("x1", "x2", "x3")))
#' dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 100)
#' fit <- etvar(y, dates, target_date = dates[80], lags = 2, nimp = 10)
#' fevd <- compute_fevd(fit)
#' round(fevd$target[, c(1, 5, 10)], 3)  # target-shock share by horizon
#' @export
compute_fevd <- function(object, nimp = NULL, target_shock = NULL) {
  if (!inherits(object, "etvar")) {
    stop("`object` must be an `etvar` object.", call. = FALSE)
  }

  irfs <- .etvar_numeric_array(object$irfs, "object$irfs")
  if (length(dim(irfs)) != 3L || dim(irfs)[1L] != dim(irfs)[2L]) {
    stop("`object$irfs` must be an n by n by h numeric array.",
         call. = FALSE)
  }

  n <- dim(irfs)[1L]
  n_available <- dim(irfs)[3L]

  if (is.null(nimp)) nimp <- n_available
  nimp <- .etvar_integer_scalar(nimp, "nimp", min = 1L,
                                 max = n_available)

  if (is.null(target_shock)) target_shock <- object$shock_index
  target_shock <- .etvar_integer_scalar(target_shock, "target_shock",
                                         min = 1L, max = n)

  variable_names <- object$variable_names
  if (!is.character(variable_names) || length(variable_names) != n ||
      anyNA(variable_names) || any(!nzchar(variable_names))) {
    variable_names <- paste0("V", seq_len(n))
  }

  shock_names <- paste("Shock", seq_len(n))
  horizon_names <- as.character(seq_len(nimp))
  shares <- array(
    0,
    dim = c(n, n, nimp),
    dimnames = list(
      variable = variable_names,
      shock = shock_names,
      horizon = horizon_names
    )
  )

  cumulative_squares <- matrix(0, nrow = n, ncol = n)

  for (h in seq_len(nimp)) {
    cumulative_squares <- cumulative_squares + irfs[, , h]^2
    total_variance <- rowSums(cumulative_squares)

    if (any(!is.finite(total_variance)) ||
        any(total_variance <= sqrt(.Machine$double.eps))) {
      stop("Forecast error variance is zero or non-finite at horizon ", h,
           ".", call. = FALSE)
    }

    shares[, , h] <- sweep(cumulative_squares, 1L, total_variance, "/")
  }

  target <- matrix(
    shares[, target_shock, , drop = FALSE],
    nrow = n,
    ncol = nimp,
    dimnames = list(variable = variable_names, horizon = horizon_names)
  )

  if (n == 1L) {
    others <- matrix(
      0,
      nrow = 1L,
      ncol = nimp,
      dimnames = dimnames(target)
    )
  } else {
    other_shocks <- setdiff(seq_len(n), target_shock)
    others <- apply(
      shares[, other_shocks, , drop = FALSE],
      c(1L, 3L),
      sum
    )
    dimnames(others) <- dimnames(target)
  }

  structure(
    list(
      shares = shares,
      target = target,
      others = others,
      horizons = seq_len(nimp),
      target_shock = target_shock,
      variable_names = variable_names,
      method = "structural"
    ),
    class = "etvar_fevd"
  )
}
