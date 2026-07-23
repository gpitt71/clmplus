#' Predict reserves
#'
#' Re-export of [ReSurv::predictReserve()].
#' @param object A fitted model or prediction object.
#' @param ... Additional arguments passed to a method.
#' @export
predictReserve <- ReSurv::predictReserve

#' Convert clmplus predictions to the ReSurv reserving interface
#'
#' @param object A `clmplusmodel` or `clmpluspredictions` object.
#' @param granularity Ignored; retained for ReSurv compatibility.
#' @param lower_triangle_only Return only forecast lower-triangle cells.
#' @param ... Arguments passed to [predict.clmplusmodel()].
#' @return A data frame containing `AP`, `DP`, `CP`, and `IBNR`.
#' @importFrom ReSurv predictReserve
#' @export
#' @method predictReserve clmplusmodel
predictReserve.clmplusmodel <- function(object, granularity = NULL,
                                        lower_triangle_only = TRUE, ...) {
  if (!inherits(object, "clmplusmodel")) {
    stop("`object` must be a `clmplusmodel` object.", call. = FALSE)
  }
  predictReserve.clmpluspredictions(
    stats::predict(object, ...), granularity = granularity,
    lower_triangle_only = lower_triangle_only
  )
}

#' @rdname predictReserve.clmplusmodel
#' @export
#' @method predictReserve clmpluspredictions
predictReserve.clmpluspredictions <- function(object, granularity = NULL,
                                              lower_triangle_only = TRUE, ...) {
  if (!inherits(object, "clmpluspredictions")) {
    stop("`object` must be a `clmpluspredictions` object.", call. = FALSE)
  }
  full <- as.matrix(object$full_triangle)
  lower <- as.matrix(object$lower_triangle)
  if (!identical(dim(full), dim(lower)) || nrow(full) != ncol(full)) {
    stop("Prediction triangles must be square matrices of equal dimensions.", call. = FALSE)
  }
  incremental <- full
  incremental[, -1L] <- full[, -1L, drop = FALSE] - full[, -ncol(full), drop = FALSE]
  index <- which(if (isTRUE(lower_triangle_only)) !is.na(lower) else !is.na(incremental),
                 arr.ind = TRUE)
  out <- data.frame(
    AP = as.integer(index[, 1L]),
    DP = as.integer(index[, 2L]),
    CP = as.integer(rowSums(index) - 1L),
    IBNR = as.numeric(incremental[index])
  )
  out <- out[is.finite(out$IBNR), , drop = FALSE]
  out <- out[order(out$AP, out$DP), , drop = FALSE]
  rownames(out) <- NULL
  out
}
