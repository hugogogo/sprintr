#' Calculate prediction from a \code{cv.sprinter} object.
#'
#' @param object a fitted \code{cv.sprinter} object.
#' @param newdata a design matrix of all the \code{p} main effects of some new observations of which predictions are to be made.
#' @param ... additional argument (not used here, only for S3 generic/method consistency)
#' @return The prediction of \code{newdata} by the cv.sprinter fit \code{object}.
#' @examples
#' n <- 100
#' p <- 200
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] + 2 * x[, 2] - 3 * x[, 1] * x[, 2] + rnorm(n)
#' mod <- cv.sprinter(x = x, y = y)
#' fitted <- predict(mod, newdata = x)
#'
#' @export
predict.cv.sprinter <- function(object, newdata, ...) {
  # input check
  stopifnot(ncol(newdata) == object$p)
  # if in Step 3, the residual is fitted
  # then in prediction, we add the prediciton in step 1
  if(object$type == 2){
    if(object$square){
      x_step1 <- cbind(newdata, myscale(newdata)^2)
    }
    else{
      x_step1 <- newdata
    }
    fitted_step1 <- as.numeric(object$step1$a0 + x_step1 %*% object$step1$beta)
  }
  else{
    fitted_step1 <- rep(0, nrow(newdata))
  }
  idx <- object$compact[, 1:2, drop = FALSE]
  # selected indices for main effects
  idxm <- idx[idx[, 1] == 0, 2]
  # selected index pairs for interactions
  idxi <- idx[idx[, 1] != 0, , drop = FALSE]

  # need to standardize the main effects to construct interactions
  xm <- myscale(newdata)
  xint <- xm[, idxi[, 1]] * xm[, idxi[, 2]]

  fitted_step3 <- as.numeric(object$a0 + cbind(newdata[, idxm], xint) %*% object$compact[, 3])
  return(fitted_step1 + fitted_step3)
}
