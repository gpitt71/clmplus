#' Predict reserve cells from a clmplus model
#'
#' Converts a fitted clmplus model, or its prediction object, into the
#' ReSurv reserving interface: AP, DP, CP, IBNR.
#'
#' @param object A clmplusmodel or clmpluspredictions object.
#' @param granularity Ignored. Included for compatibility with ReSurv.
#' @param lower_triangle_only Logical. If TRUE, return only predicted lower-triangle cells.
#' @param ... Passed to predict.clmplusmodel().
#'
#' @return A data.frame with AP, DP, CP, IBNR.
#'
#' @importFrom ReSurv predictReserve
#' @export
#' @method predictReserve clmplusmodel
predictReserve.clmplusmodel <- function(object,
                                        granularity = NULL,
                                        lower_triangle_only = TRUE,
                                        ...) {
  
  if (!inherits(object, "clmplusmodel")) {
    stop("`object` must be a `clmplusmodel` object.", call. = FALSE)
  }
  
  pred <- predict(object, ...)
  
  predictReserve.clmpluspredictions(
    object = pred,
    granularity = granularity,
    lower_triangle_only = lower_triangle_only
  )
}


#' @rdname predictReserve.clmplusmodel
#' @export
#' @method predictReserve clmpluspredictions
predictReserve.clmpluspredictions <- function(object,
                                              granularity = NULL,
                                              lower_triangle_only = TRUE,
                                              ...) {
  
  if (!inherits(object, "clmpluspredictions")) {
    stop("`object` must be a `clmpluspredictions` object.", call. = FALSE)
  }
  
  if (is.null(object$full_triangle)) {
    stop("`object$full_triangle` is missing.", call. = FALSE)
  }
  
  if (is.null(object$lower_triangle)) {
    stop("`object$lower_triangle` is missing.", call. = FALSE)
  }
  
  full_triangle <- as.matrix(object$full_triangle)
  lower_triangle <- as.matrix(object$lower_triangle)
  
  if (!all(dim(full_triangle) == dim(lower_triangle))) {
    stop(
      "`full_triangle` and `lower_triangle` must have the same dimensions.",
      call. = FALSE
    )
  }
  
  J <- nrow(full_triangle)
  
  if (ncol(full_triangle) != J) {
    stop("`full_triangle` must be square.", call. = FALSE)
  }
  
  incremental_full <- full_triangle
  
  incremental_full[, 1L] <- full_triangle[, 1L]
  
  if (J >= 2L) {
    for (jj in 2:J) {
      incremental_full[, jj] <- full_triangle[, jj] - full_triangle[, jj - 1L]
    }
  }
  
  lower_mask <- !is.na(lower_triangle)
  
  if (!isTRUE(lower_triangle_only)) {
    lower_mask <- matrix(TRUE, nrow = J, ncol = J)
  }
  
  AP <- rep(seq_len(J), times = J)
  DP <- rep(seq_len(J), each = J)
  
  out <- data.frame(
    AP = AP,
    DP = DP,
    CP = AP + DP - 1L,
    IBNR = as.numeric(incremental_full[cbind(AP, DP)])
  )
  
  out <- out[as.vector(lower_mask), , drop = FALSE]
  
  out$AP <- as.integer(out$AP)
  out$DP <- as.integer(out$DP)
  out$CP <- as.integer(out$CP)
  out$IBNR <- as.numeric(out$IBNR)
  
  out <- out[is.finite(out$IBNR), , drop = FALSE]
  out <- out[order(out$AP, out$DP), , drop = FALSE]
  
  rownames(out) <- NULL
  
  out
}