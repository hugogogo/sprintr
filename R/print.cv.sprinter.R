#' Print the cross validation information of cv.sprinter
#'
#' Print a summary of the cross-validation information for running cv.sprinter.
#'
#' This function takes in a \code{cv.sprinter} object and produces summary of the cross-validation informationabout the tuning parameters (in Step 3) selected by \code{lambda.min} and \code{lambda.1se}.
#' Adopted from the function \code{print.cv.rgam} from package \code{relgam} by Kenneth Tay and Robert Tibshirani.
#'
#' @param fit A fitted \code{cv.sprinter} object.
#' @param digits Significant digits in printout.
#'
#' @seealso \code{\link{cv.sprinter}}, \code{\link{print.printer}}.
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 100
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)
#'
#' fit.cv <- cv.sprinter(x = x, y = y)
#' print(fit.cv)
#'
#' @export
print.cv.sprinter <- function(fit, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(fit$call), "\n\n")

  # get indices for optimal lambda values
  which <- c(fit$i_lambda1_best, fit$i_lambda3_best)
  lams <- c(fit$fit$lambda1[which[1]], fit$fit$lambda3[which[2], which[1]])
  out <- cbind(signif(lams[1], digits),
               signif(lams[2], digits),
               signif(fit$cvm[which[1], which[2]], digits),
               signif(fit$cvse[which[1], which[2]], digits),
               fit$fit$step3[[which[1]]]$nzm[which[2]],
               fit$fit$step3[[which[1]]]$nzi[which[2]])
  dimnames(out) = list(NULL, c("lambda1", "lambda3", "mean(vali-err)", "se(vali-err)", "#nonzero-main", "#nonzero-inter"))
  print(out)
}
