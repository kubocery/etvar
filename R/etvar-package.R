#' etvar: Event-Targeted Vector Autoregression
#'
#' @description
#' Estimate structural VAR models identified by targeting a specific event date.
#' The package finds the orthogonal rotation of the Cholesky factor that best
#' explains observed movements on the event date using a single structural
#' shock.
#'
#' @section Main functions:
#' \describe{
#'   \item{[etvar()]}{One-call estimation: VAR + optimisation + decomposition.}
#'   \item{[estimate_var()]}{Reduced-form VAR estimation.}
#'   \item{[etvar_optimize()]}{Three-stage optimisation of rotation angles.}
#'   \item{[hist_decomposition()]}{Historical decomposition.}
#'   \item{[hist_decomposition_3d()]}{Historical decomposition with
#'     intervention path.}
#' }
#'
#' @docType package
#' @name etvar-package
"_PACKAGE"
