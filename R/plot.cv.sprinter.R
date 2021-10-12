#' Plot function of cv.sprinter fit
#'
#' This function produces plots of cross-validation for cv.sprinter.
#'
#' The orange pairs on the top of the plot shows the number of non-zero (main effects, interactions) selected by each value of lambda. Adopted from the function \code{plot.cv.rgam} from package \code{relgam} by Kenneth Tay and Robert Tibshirani.
#'
#' @param fit A "\code{cv.sprinter}" object.
#'
#' @seealso \code{\link{cv.sprinter}}.
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 200
#' # dense input
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)
#' mod <- cv.sprinter(x = x, y = y)
#'
#' plot(mod)
#' @export
plot.cv.sprinter <- function(fit) {
  xlab <- expression(paste("-log(", lambda, ")"))

  lambda <- fit$fit$lambda3[, fit$i_lambda1_best]
  lam <- -log(lambda)

  cvm <- fit$cvm[fit$i_lambda1_best, ]
  cvse <- fit$cvse[fit$i_lambda1_best, ]
  cvup <- cvm + cvse
  cvlo <- cvm - cvse

  plot.args <- list(x = lam,
                    y = cvm,
                    ylim = range(cvlo, cvup),
                    xlab = xlab, ylab = "validation set error", type = "n")
  do.call("plot", plot.args)

  # started plotting error bars
  width <- 0.01
  barw <- diff(range(lam)) * width
  segments(lam, cvup, lam, cvlo, "grey")
  segments(lam - barw, cvup, lam + barw, cvup)
  segments(lam - barw, cvlo, lam + barw, cvlo)

  points(lam, cvm, pch = 20, col = "red")

  nzm <- fit$fit$step3[[fit$i_lambda1_best]]$nzm
  nzi <- fit$fit$step3[[fit$i_lambda1_best]]$nzi

  axis(side = 3, at = lam, labels = paste(paste(nzm), "/", paste(nzi)),
       tick = FALSE, line = 0, col.axis = "Orange")
  legend("topright", paste0("With lambda_1 = ", signif(fit$lambda1_best, 3), " selected for Step 1."), bty = "n", pch = NA, cex = 1.2)

  lam_min <- lambda[fit$i_lambda3_best]
  abline(v = -log(lam_min), lty = 2, col = "Blue")
}
