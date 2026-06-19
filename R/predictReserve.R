#' Predict reserve cells
#'
#' Generic for reserve-cell predictions.
#'
#' @param object A fitted reserving model.
#' @param ... Additional arguments.
#'
#' @export
predictReserve <- function(object, ...) {
  UseMethod("predictReserve")
}


#' Predict reserve cells from a clmplus model
#'
#' Converts a fitted \code{clmplusmodel} or a \code{clmpluspredictions}
#' object into the reserving interface used by ReSurv:
#' \code{AP}, \code{DP}, \code{CP}, \code{IBNR}.
#'
#' The returned \code{IBNR} values are incremental lower-triangle predictions.
#'
#' @param object A \code{clmplusmodel} or \code{clmpluspredictions} object.
#' @param granularity Ignored. Included for compatibility with ReSurv's
#'   \code{predictReserve()} interface.
#' @param lower_triangle_only Logical. If \code{TRUE}, return only cells in the
#'   predicted lower triangle.
#' @param ... Arguments passed to \code{predict.clmplusmodel()} when
#'   \code{object} is a fitted \code{clmplusmodel}.
#'
#' @return A data frame with columns \code{AP}, \code{DP}, \code{CP}, \code{IBNR}.
#'
#' @export
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
    stop(
      "`full_triangle` must be a square cumulative triangle.",
      call. = FALSE
    )
  }
  
  ## Convert cumulative full triangle to incremental full triangle.
  incremental_full <- full_triangle
  
  if (J >= 1L) {
    incremental_full[, 1L] <- full_triangle[, 1L]
  }
  
  if (J >= 2L) {
    for (jj in 2:J) {
      incremental_full[, jj] <- full_triangle[, jj] - full_triangle[, jj - 1L]
    }
  }
  
  ## The clmplus lower triangle is cumulative. Its non-NA cells identify
  ## the predicted lower-triangle cells, also when forecasting_horizon is used.
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
  
  out <- out[order(out$AP, out$DP), , drop = FALSE]
  
  rownames(out) <- NULL
  
  out$AP <- as.integer(out$AP)
  out$DP <- as.integer(out$DP)
  out$CP <- as.integer(out$CP)
  out$IBNR <- as.numeric(out$IBNR)
  
  out
}