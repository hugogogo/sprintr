#' Sure independence screening followed by lasso
#'
#'
#' @param x An \code{n} by \code{p} design matrix of main effects. Each row is an observation of \code{p} main effects.
#' @param y A response vector of size \code{n}.
#' @param num_keep Number of variables to keep in the screening phase
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
                      nfold = 5, foldid = NULL){
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
  xx <- x[, idx[, 1]] * x[, idx[, 2]]
   # note that we use penalty.factor to make sure that main effects are not penalized
  fit <- cv.glmnet(x = cbind(x, xx), y = y)

  coef <- fit$glmnet.fit$beta[, which.min(fit$cvm)]
  idx_all <- rbind(cbind(rep(0, p), seq(p)), idx[, 1:2])
  compact <- cbind(idx_all[which(coef != 0), , drop = FALSE], coef[coef != 0])
  rownames(compact) <- NULL
  colnames(compact) <- c("index_1", "index_2", "coefficient")

  # finally return the best lambda
  out <- list(n = n,
              p = p,
              fit = fit,
              compact = compact,
              call = match.call())
  class(out) <- "cv.hier"
  return(out)
}
