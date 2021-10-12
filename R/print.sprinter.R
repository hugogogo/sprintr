#' Print a summary of the sprinter fit
#'
#' Print a summary of the sprinter fit at each step along the path of tuning parameters used in Step 3, for any given tuning parameter in Step 1.
#'
#' The function produces a three-column matrix with tuning parameter values (in Step 3), number of nonzero main effects, and the number of nonzero interactions.
#' Adopted from the function \code{print.rgam} from package \code{relgam} by Kenneth Tay and Robert Tibshirani.
#'
#' @param fit A \code{sprinter} object.
#' @param which Which tuning parameter of Step 1 to print. Default is 1.
#' @param digits Significant digits in printout.
#' @param ... Additional print arguments.
#'
#' @seealso \code{\link{sprinter}}.
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 100
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)
#' fit <- sprinter(x = x, y = y)
#'
#' print(fit, which = 3)
#'
#' @export
print.sprinter <- function(fit, which = 1, digits = max(3, getOption("digits") - 3), ...) {
  stopifnot(which <= length(fit$lambda1))
  cat("\nCall: ", deparse(fit$call), "\n\n")

  this <- fit$step3[[which]]

  out <- cbind(signif(fit$lambda3[, which], digits), this$nzm, this$nzi)
  colnames(out) = c("lambda", "#nz main", "#nz interaction")
  print(out)
  invisible(out)
}
