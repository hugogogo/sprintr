#' Sure independence screening followed by lasso
#'
#'
#' @param x An \code{n} by \code{p} design matrix of main effects. Each row is an observation of \code{p} main effects.
#' @param y A response vector of size \code{n}.
#' @param num_keep Number of variables to keep in the screening phase
#' @param ... other arguments to be passed to the \code{glmnet} calls, such as \code{alpha} or \code{penalty.factor}
#'
#' @return An object of S3 class "\code{cv.hier}".
#'  \describe{
#'   \item{\code{n}}{The sample size.}
#'   \item{\code{p}}{The number of main effects.}
#'   \item{\code{fit}}{The whole \code{cv.glmnet} fit object.}
#'   \item{\code{compact}}{A compact representation of the selected variables. \code{compact} has three columns, with the first two columns representing the indices of a selected variable (main effects with first index = 0), and the last column representing the estimate of coefficients.}
#'  }
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 200
#' # dense input
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)
#' mod <- hier_lasso(x = x, y = y)
#'
#' @import glmnet
#' @export
sis_lasso <- function(x, y,
                      num_keep = NULL,
                      lam_min_ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04),
                      nfold = 5, foldid = NULL, ...){
  n <- nrow(x)
  p <- ncol(x)
  stopifnot(n == length(y))

  if(is.null(num_keep))
    num_keep <- ceiling(n / log(n))

  if(inherits(x, "sparseMatrix")){
    idx <- screen_sparse_cpp(x = x, y = y, num_keep = num_keep, square = FALSE, main_effect = FALSE)
  }
  else{
    idx <- screen_cpp(x = x, y = y, num_keep = num_keep, square = FALSE, main_effect = FALSE)
  }

  # step 2: fit cv.glment with all main effects and constructed (hierarchical) interactions
  xx <- myscale(cbind(x, x[, idx[, 1]] * x[, idx[, 2]]))
  mean_y <- mean(y)

  fit <- cv.glmnet(x = xx, y = y,
                   lambda = get_lambda(x = xx, y = y - mean_y),
                   intercept = FALSE,
                   standardize = FALSE, ...)


  coef <- fit$glmnet.fit$beta[, which.min(fit$cvm)]
  # in-sample fitted value
  fitted <- as.numeric(mean_y + xx %*% coef)
  # scale estimates back to the original scale of x
  coef <- coef / attr(xx, which = "scaled:scale")
  a0 <- as.numeric(mean_y - crossprod(attr(xx, which = "scaled:center"), coef))

  idx_all <- rbind(cbind(rep(0, p), seq(p)), idx[, 1:2])
  compact <- cbind(idx_all[which(coef != 0), , drop = FALSE], coef[coef != 0])
  rownames(compact) <- NULL
  colnames(compact) <- c("index_1", "index_2", "coefficient")

  # finally return the best lambda
  out <- list(n = n,
              p = p,
              type = 1,
              fit = fit,
              a0 = a0,
              fitted = fitted,
              compact = compact,
              call = match.call())
  class(out) <- "other"
  return(out)
}
