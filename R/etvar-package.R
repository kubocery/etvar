#' etvar: Event-Targeted Vector Autoregression
#'
#' @description
#' Structural VAR models identified by one or more event dates: the target
#' shock is the rotation of the Cholesky shocks that best explains how all
#' variables moved on those dates.
#'
#' Start with [etvar()], which fits the whole model, and look at the result
#' with [plot.etvar()], [compute_fevd()] and [compute_irf_confint()]. The
#' building blocks are available on their own: [estimate_var()] for the
#' reduced-form VAR, [etvar_objective()] and [etvar_optimize()] for the
#' identification, [compute_structural_irfs()], [hist_decomposition()] and
#' [hist_decomposition_3d()] for the structural analysis.
#'
#' @references
#' Pinter, G., Smets, F. and Üslü, S. (2025). Market whiplash after the
#' 2025 tariff shock: an event-targeted VAR approach. *BIS Working Papers*,
#' No. 1282. Bank for International Settlements.
#' <https://www.bis.org/publications/working-paper-1282-market-whiplash-after-2025-tariff-shock-event-targeted-var-approach>
"_PACKAGE"
