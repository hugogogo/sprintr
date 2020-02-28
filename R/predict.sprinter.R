#' Calculate prediction from a \code{sprinter} object.
#'
#' @param object a fitted \code{sprinter} object.
#' @param newdata a design matrix of all the \code{p} main effects of some new observations of which predictions are to be made.
#' @param ... additional argument (not used here, only for S3 generic/method consistency)
#' @return The prediction of \code{newdata} by the sprinter fit \code{object}.
#' @examples
#' n <- 100
#' p <- 200
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] + 2 * x[, 2] - 3 * x[, 1] * x[, 2] + rnorm(n)
#' mod <- sprinter(x = x, y = y)
#' fitted <- predict(mod, newdata = x)
#'
#' @export
predict.sprinter <- function(object, newdata, ...) {
  # input check
  stopifnot(ncol(newdata) == object$p)
  n_num_keep <- length(object$step3)
  fitted <- vector("list", length = n_num_keep)
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
  for(i in seq(n_num_keep)){
    obj_curr <- object$step3[[i]]
    idx <- obj_curr$idx[, , drop = FALSE]
    # selected indices for main effects
    idxm <- idx[idx[, 1] == 0, 2]
    # selected index pairs for interactions
    idxi <- idx[idx[, 1] != 0, , drop = FALSE]

    # need to standardize the main effects to construct interactions
    xm <- myscale(newdata)
    xint <- xm[, idxi[, 1]] * xm[, idxi[, 2]]

    fitted[[i]]$fitted <- fitted_step1 + as.matrix(cbind(newdata[, idxm], xint) %*% obj_curr$coef)
    colnames(fitted[[i]]$fitted) <- NULL
    fitted[[i]]$fitted <- t(obj_curr$a0 + t(fitted[[i]]$fitted))
  }
  return(fitted)
}
