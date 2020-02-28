#' Plot function of sprinter fit
#'
#' This function produces plots of the predicted values for specified main effects at a given tuning parameter values pair.
#'
#' For each given query of interaction (j, k), this function produces four plots. The top two panels are the dependence of the predicted value only on the two main effects. The bottom left panel shows the dependence of the predicted value on the interaction x_j * x_k. And the bottom right panel shows the dependence of the predicted value on both the main effects and their interaction. Adopted from the function \code{plot.rgam} from package \code{relgam} by Kenneth Tay and Robert Tibshirani.
#'
#' @param fit The \code{sprinter} fit object.
#' @param newdata A new design matrix of all the main effects.
#' @param index Indices pair of the two tuning parameters at which the plot is produced. The default value is \code{c(1, length(fit$lambda))}.
#' @param which Indices pair of the interaction to plot. The default value is (1, 2).
#' @param rugplot If \code{TRUE} (default), adds a rugplot showing the values of x at the bottom of each fitted function plot.
#' @param grid_length The number of points to evaluate the estimated function at.
#' @param ... Optional graphical parameters in the function \code{plot}.
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
#' plot(fit, newdata = x, index = c(1, 50), which = c(4, 5))
#'
#' # plot the fit from cv.sprinter
#' fit.cv <- cv.sprinter(x = x, y = y)
#' plot(fit.cv$fit, newdata = x, index = c(fit.cv$num_keep_best, fit.cv$lam_best))
#' @import plot3D
#'
#' @export
plot.sprinter <- function(fit, newdata, index = NULL, which = NULL,
                          rugplot = TRUE, grid_length = 100, ...) {
  x = newdata
  p = ncol(x)
  if (is.null(index)) {
    warning(paste("Tuning parameter pair indices not provided. Plot the pair (", 1, ",", length(fit$lambda[[1]]), ") by default.", sep = ""))
    index = c(1, length(fit$lambda))
  }

  if (is.null(which)) {
    warning(paste("Interaction indices not provided. Plot the interaction that has the highest score in Step 2.", sep = ""))
    which <- fit$step3[[index[1]]]$idx[p + 1, 1:2]
  }
  if(which[1] > which[2])
    which <- which[c(2, 1)]

  x1lab <- bquote("x"[.(which[1])])
  x2lab <- bquote("x"[.(which[2])])

  ylab_m1 <- bquote(hat(beta)[.(which[1])] * "x"[.(which[1])])
  ylab_m2 <- bquote(hat(beta)[.(which[2])] * "x"[.(which[2])])

  # get xrange for the plot
  x1range <- seq(min(x[, which[1]]), max(x[, which[1]]), length.out = grid_length)
  x2range <- seq(min(x[, which[2]]), max(x[, which[2]]), length.out = grid_length)

  # extract model from the sprinter fit
  mod <- fit$step3[[index[1]]]
  # start with the main effects part
  beta1 <- mod$coef[which[1], index[2]]
  beta2 <- mod$coef[which[2], index[2]]

  # predicted value from main effects
  a0 <- mod$a0[index[2]]
  y_x1 <- x1range * beta1
  y_x2 <- x2range * beta2

  # main effects plots
  par(mfrow = c(2, 2))
  par(mar=c(2,2,2,2))
  layout(matrix(c(1, 3, 2, 4), 2), widths = c(1, 1), height = c(1, 1.5), respect = FALSE)

  # first main effects
  plot(x = x1range, y = y_x1, type = "l", col = "Blue", lwd = 2, xlab = x1lab, ylab = "", main = ylab_m1)
  if (rugplot) {
    rug(x[, which[1]])
  }
  # second main effects
  plot(x = x2range, y = y_x2, type = "l", col = "Blue", lwd = 2, xlab = x2lab, ylab = "", main = ylab_m2)
  if (rugplot) {
    rug(x[, which[2]])
  }

  # get the interaction part
  # check if the interaction is selected
  int_check <- rowSums(cbind(mod$idx[, 1] == which[1], mod$idx[, 2] == which[2]))
  gamma <- mod$coef[which(int_check == 2), index[2]]
  if (length(gamma) == 0 || gamma == 0){
    warning(paste("Coefficient estimate of the interaction between ", x1lab, " and ", x2lab, " is zero!", sep = ""))
    gamma <- 0
  }

  # plot pure interaction effects
  ff <- gamma * tcrossprod(x1range, x2range)

  # unfortunately, persp3D does not work with either expression or bquote
  title_inter <- paste0(toString(round(gamma, 2)), " x", which[1], " * x", which[2])
  plot3D::persp3D(x = x1range, y = x2range, z = ff, colkey = TRUE, col = gg2.col(200), phi = 25, theta = -25, xlab = paste0("x", which[1]), ylab = paste0("x", which[2]), zlab = "", zlim = range(ff), main = title_inter)

  # plot the full effects
  getz <- function(x, y){
    result <- matrix(NA, length(x), length(y))
    for(i in seq(length(x))){
      for(j in seq(length(y))){
        result[i, j] <- a0 + beta1 * x[i] + beta2 * y[j] + gamma * x[i] * y[j]
      }
    }
    return(result)
  }
  fff <- getz(x1range, x2range)
  #mm <- paste0("f (x", which[1], ", x", which[2], ")")
  if(beta2 >= 0){
    title_main <- paste0(toString(round(beta1, 2)), " x", which[1], " + ", toString(round(beta2, 2)), " x", which[2])
  }
  else{
    title_main <- paste0(toString(round(beta1, 2)), " x", which[1], " ", toString(round(beta2, 2)), " x", which[2])
  }

  if(gamma >= 0){
    title <- paste0(title_main, " + ", title_inter)
  }
  else{
    title <- paste0(title_main, " ", title_inter)
  }

  plot3D::persp3D(x = x1range, y = x2range, z = fff, colkey = TRUE, col = gg2.col(200), phi = 25, theta = -25, xlab = paste0("x", which[1]), ylab = paste0("x", which[2]), zlab = "", zlim = range(ff), main = title)

  # set back to original plotting par
  par(mfrow = c(1, 1))
}
