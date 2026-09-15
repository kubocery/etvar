# Internal validation helpers -------------------------------------------------

.etvar_integer_scalar <- function(x, arg, min = NULL, max = NULL) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    x == floor(x)

  if (!valid || (!is.null(min) && x < min) || (!is.null(max) && x > max)) {
    bounds <- if (!is.null(min) && !is.null(max)) {
      paste0(" between ", min, " and ", max)
    } else if (!is.null(min)) {
      paste0(" greater than or equal to ", min)
    } else if (!is.null(max)) {
      paste0(" less than or equal to ", max)
    } else {
      ""
    }
    stop("`", arg, "` must be one integer", bounds, ".", call. = FALSE)
  }

  as.integer(x)
}

.etvar_flag <- function(x, arg) {
  if (length(x) != 1L || is.na(x) ||
      !(is.logical(x) || (is.numeric(x) && is.finite(x) && x %in% c(0, 1)))) {
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  }
  isTRUE(as.logical(x))
}

.etvar_numeric_matrix <- function(x, arg = "data") {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("`", arg, "` must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (is.data.frame(x) && !all(vapply(x, is.numeric, logical(1)))) {
    stop("Every column of `", arg, "` must be numeric.", call. = FALSE)
  }

  x <- as.matrix(x)
  if (!is.numeric(x) || length(dim(x)) != 2L || ncol(x) < 1L) {
    stop("`", arg, "` must be a numeric matrix with at least one column.",
         call. = FALSE)
  }
  if (nrow(x) < 1L || anyNA(x) || any(!is.finite(x))) {
    stop("`", arg, "` must contain only finite, non-missing values.",
         call. = FALSE)
  }

  storage.mode(x) <- "double"
  x
}

.etvar_numeric_array <- function(x, arg, dimensions = NULL) {
  if (!is.numeric(x) || is.null(dim(x)) || anyNA(x) || any(!is.finite(x))) {
    stop("`", arg, "` must be a finite numeric array.", call. = FALSE)
  }
  if (!is.null(dimensions) && !identical(as.integer(dim(x)),
                                          as.integer(dimensions))) {
    stop("`", arg, "` must have dimensions ",
         paste(dimensions, collapse = " x "), ".", call. = FALSE)
  }
  x
}

.etvar_as_date <- function(x, arg) {
  parsed <- tryCatch(
    as.Date(x),
    error = function(e) {
      stop("`", arg, "` must be coercible to `Date`.", call. = FALSE)
    }
  )
  if (length(parsed) == 0L || anyNA(parsed)) {
    stop("`", arg, "` must contain at least one non-missing date.",
         call. = FALSE)
  }
  parsed
}

.etvar_index_vector <- function(x, arg, upper, allow_empty = FALSE) {
  valid <- is.numeric(x) && !anyNA(x) && all(is.finite(x)) &&
    all(x == floor(x))
  if (!valid || (!allow_empty && length(x) == 0L) ||
      any(x < 1L | x > upper)) {
    stop("`", arg, "` must contain integer indices between 1 and ", upper,
         ".", call. = FALSE)
  }
  as.integer(x)
}

.etvar_validate_decomposition <- function(A0, ALPHA, Y, X, n, lags,
                                           trend, constant) {
  n <- .etvar_integer_scalar(n, "n", min = 1L)
  lags <- .etvar_integer_scalar(lags, "lags", min = 1L)
  trend <- .etvar_integer_scalar(trend, "trend", min = 0L, max = 2L)
  constant <- as.integer(.etvar_flag(constant, "constant"))

  Y <- .etvar_numeric_matrix(Y, "Y")
  X <- .etvar_numeric_matrix(X, "X")
  ALPHA <- .etvar_numeric_matrix(ALPHA, "ALPHA")
  A0 <- .etvar_numeric_matrix(A0, "A0")

  expected_k <- constant + trend + n * lags
  if (!identical(dim(A0), c(n, n))) {
    stop("`A0` must be an `n` by `n` matrix.", call. = FALSE)
  }
  if (qr(A0)$rank < n) {
    stop("`A0` must be nonsingular.", call. = FALSE)
  }
  if (ncol(Y) != n || nrow(X) != nrow(Y) || ncol(X) != expected_k) {
    stop("`Y` and `X` dimensions are inconsistent with `n`, `lags`, ",
         "`trend`, and `constant`.", call. = FALSE)
  }
  if (!identical(dim(ALPHA), c(expected_k, n))) {
    stop("`ALPHA` dimensions are inconsistent with the VAR specification.",
         call. = FALSE)
  }

  list(A0 = A0, ALPHA = ALPHA, Y = Y, X = X, n = n, lags = lags,
       trend = trend, constant = constant)
}
