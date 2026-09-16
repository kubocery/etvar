#' etvar: Event-Targeted Vector Autoregression
#'
#' @description
#' Estimate structural VAR models identified by targeting one or more event
#' dates.
#' The package finds the orthogonal rotation of the Cholesky factor that best
#' explains observed event-date movements using a single structural shock. It
#' provides impulse responses, structural forecast error variance
#' decompositions, historical decompositions, and diagnostic plots.
#'
#' @section Main functions:
#' \describe{
#'   \item{[etvar()]}{One-call estimation: VAR + optimisation + decomposition.}
#'   \item{[estimate_var()]}{Reduced-form VAR estimation.}
#'   \item{[etvar_optimize()]}{Three-stage optimisation of the target-shock
#'     direction.}
#'   \item{[compute_fevd()]}{Structural forecast error variance decomposition.}
#'   \item{[hist_decomposition()]}{Historical decomposition.}
#'   \item{[hist_decomposition_3d()]}{Historical decomposition with
#'     intervention path.}
#' }
#'
#' @docType package
#' @name etvar-package
"_PACKAGE"
