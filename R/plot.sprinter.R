#' Plot function of sprinter fit
#'
#' Produces a two-panel plot of the sprinter object showing coefficient paths for both main effects and interactions.
#'
#' A two panel plot is produced, that summarizes the main effects (left) and interaction (right) coefficients, as a function of lambda. Adopted from the function \code{summary.rgam} from package \code{relgam} by Kenneth Tay and Robert Tibshirani.
#'
#' @param fit Fitted \code{sprinter} object.
#' @param which The tuning parameter considered in Step 2.
#' @param label If \code{TRUE} (default), annotate the plot with variable labels.
#' @param index Lambda indices to plot
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 100
#' # dense input
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)
#' fit <- sprinter(x = x, y = y)
#'
#' plot(fit)
#'
#' @export
plot.sprinter <- function (fit, which = 1, label = TRUE, index = NULL)
{
  p <- fit$p
  pp <- nrow(fit$step3[[which]]$coef)
  colours_m <- rainbow(n = p)
  colours_i <- rainbow(n = pp - p)

  # get the relevant lambda
  if (is.null(index)) {
    lambda <- fit$lambda3[, which]
    index <- 1:length(lambda)
  } else {
    lambda <- fit$lambda3[index, which]
  }

  # internal function for drawing line for one feature
  nzlines <- function(lambda, alpha, ...) {
    if (any(abs(alpha) > 0)) {
      num_lambda <- length(lambda)
      start = max(1, min(seq(num_lambda)[abs(alpha) > 0]) - 1)
      whichnz = seq(from = start, to = num_lambda)
      if (length(whichnz) > 1) {
        lines(lambda[whichnz], alpha[whichnz], ...)
      }
    }
    invisible()
  }

  par(mfrow = c(1, 2))
  # plot for linear coefficients
  lin_beta <- fit$step3[[which]]$coef[1:p, index, drop = FALSE]
  plot(0, type = "n", xlab = expression(lambda), ylab = "Coefficient",
       xlim = c(max(lambda) * 1.1, min(lambda) * (0.8 - 0.2 * label)),
       ylim = range(lin_beta), main = "Main Effects", log = "x")
  abline(h = 0, lty = 3)
  for (j in 1:p) {
    nzlines(lambda, lin_beta[j, ], col = colours_m[j], lwd = 2,
            type = "l", pch = "", cex = 0)
  }
  if (label) {
    text(rep(min(lambda) * 0.85, p), lin_beta[1:p, length(lambda)],
         labels = 1:p, col = colours_m, cex = 0.6)
  }

  # plot for interaction coefficients
  inter_beta <- fit$step3[[which]]$coef[(p + 1): pp, index, drop = FALSE]
  plot(0, type = "n", xlab = expression(lambda), ylab = "Coefficient",
       xlim = c(max(lambda) * 1.1, min(lambda) * (0.8 - 0.2 * label)),
       ylim = range(inter_beta), main = "Interactions", log = "x")
  abline(h = 0, lty = 3)
  for (j in 1:nrow(inter_beta)) {
    nzlines(lambda, inter_beta[j, ], col = colours_i[j], lwd = 2,
            type = "l", pch = "", cex = 0)
  }
  label_i <- fit$step2[[which]][, 1:2]
  label_i <- paste0("(", label_i[, 1], ", ", label_i[, 2], ")")
  if (label) {
    text(rep(min(lambda) * 0.85, p), inter_beta[1:(pp - p), length(lambda)],
         labels = label_i, col = colours_i, cex = 0.6)
  }

  par(mfrow = c(1, 1))
}
