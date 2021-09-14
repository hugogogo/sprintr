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
  nlam1 <- length(object$lambda1)

  out <- vector("list", length = nlam1)

  if(object$square)
    x_step1 <- cbind(newdata, myscale(newdata)^2)
  else
    x_step1 <- newdata

  # need to standardize the main effects to construct interactions
  xm <- myscale(newdata)

  # we add the prediction in Step1 (for each lambda1 value) and in Step3 (for a path of lambda3 values)
  for(k in seq(nlam1)){
    fitted_step1 <- as.numeric(object$step1$a0[k] + x_step1 %*% object$step1$beta[, k])

    idx <- object$step2[[k]]
    xint <- xm[, idx[, 1]] * xm[, idx[, 2]]
    design <- cbind(x_step1, xint)

    out[[k]] <- fitted_step1 + as.matrix(design %*% object$step3[[k]]$coef)
    colnames(out[[k]]) <- NULL
    out[[k]] <- t(object$step3[[k]]$a0 + t(out[[k]]))
  }
  return(out)
}
