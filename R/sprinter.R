#' Sparse Reluctant Interaction Modeling
#'
#' This is the main function that fits interaction models with a path of tuning parameters (for Step 3).
#'
#' @param x An \code{n} by \code{p} design matrix of main effects. Each row is an observation of \code{p} main effects.
#' @param y A response vector of size \code{n}.
#' @param square Indicator of whether squared effects should be fitted in Step 1. Default to be FALSE.
#' @param type An indicator: if type == 1, in Step 3 we fit y on (main effects, selected interactions); if type == 2, in Step 3 we fit the residual from step 1 on only the selected interactions
#' @param num_keep A user specified list of number of candidate interactions to keep in Step 2. If \code{num_keep} is not specified (as default), it will be set to a sequence from \code{[n / log n]} to 0, where the length of the sequence is set to the default value of \code{n_num_keep}.
#' @param n_num_keep The number of \code{num_keep} values. Default to be \code{5}.
#' @param lambda A user specified list of tuning parameter. \code{lambda} is a list object of length \code{n_num_keep}, and an element \code{lambda[[i]]} is a vector of length \code{nlam}. Default to be NULL, and the program will compute its own \code{lambda} paths based on \code{n_num_keep}, \code{nlam} and \code{lam_min_ratio}.
#' @param nlam the number of values in each \code{lambda[[i]]} path. If not specified, they will be all set to \code{100}.
#' @param lam_min_ratio The ratio of the smallest and the largest values in each \code{lambda[[i]]}. The largest value in \code{lambda} is usually the smallest value for which all coefficients are set to zero. Default to be \code{1e-2} in the \code{n} < \code{p} setting.
#' @param verbose If \code{TRUE} (default), a progress bar shows the progress of the fitting.
#'
#' @return An object of S3 class "\code{sprinter}".
#'  \describe{
#'   \item{\code{type}}{The \code{type} parameter passed into sprinter}
#'   \item{\code{square}}{The \code{square} parameter passed into sprinter}
#'   \item{\code{step1}}{The output from fitting Step 1}
#'   \item{\code{lambda}}{The path of tuning parameters passed into / computed for fitting Step 3}
#'   \item{\code{num_keep}}{The path of tuning parameters for Step 2}
#'   \item{\code{n}}{The sample size.}
#'   \item{\code{p}}{The number of main effects.}
#'   \item{\code{step3}}{The output from fitting Step 3}
#'   \item{\code{call}}{Function call.}
#'  }
#' @seealso
#'   \code{\link{cv.sprinter}}
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 100
#' # dense input
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)
#' mod <- sprinter(x = x, y = y)
#'
#' # sparse input
#' library(Matrix)
#' x <- Matrix::Matrix(0, n, p)
#' idx <- cbind(sample(seq(n), size = 10, replace = TRUE), sample(seq(p), size = 10, replace = TRUE))
#' x[idx] <- 1
#' y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)
#' mod <- sprinter(x = x, y = y)
#'
#' @import glmnet
#' @export
sprinter <- function(x, y, square = FALSE, type = 1,
                     num_keep = NULL, n_num_keep = 5,
                     lambda = NULL, nlam = 100,
                     lam_min_ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04),
                     verbose = TRUE){
  n <- nrow(x)
  p <- ncol(x)

  # num_keep = ceiling(seq(ceiling(n / log(n)), 1, length.out = n_num_keep))
  # by default, max(num_keep) is set to be [n / log(n)]
  # and num_keep is decreasing
  if(is.null(num_keep)){
    num_keep = ceiling(seq(ceiling(n / log(n)), 1, length.out = n_num_keep))
  }
  else{
    stopifnot(max(num_keep) <= ifelse(square, p*(p - 1) / 2, p * (p - 1)/ 2 + p))
  }
  n_num_keep <- length(num_keep)

  # we always standardize the design matrix to get main effects
  # squared effects and interactions are built upon standardized main effects
  x <- myscale(x)
  # xm is the (standardized) design matrix of main effects
  col_mean <- attr(x = x, which = "scaled:center")
  col_sd <- attr(x = x, which = "scaled:scale")
  xm <- x
  mean_y <- mean(y)

  # The First Step
  # run CV-lasso on
  # (1) main effects M (square == FALSE)
  # (2) main effects M + squared main effects M^2 (square == TRUE)
  # return the fitted value of response
  if(square){
    x_sq <- myscale(x^2)

    col_mean <- c(col_mean, attr(x_sq, which = "scaled:center"))
    col_sd <- c(col_sd, attr(x_sq, which = "scaled:scale"))
    x <- cbind(x, x_sq)
  }

  # progress bar initialized
  if (verbose)
    pb <- utils::txtProgressBar(min = 0, max = 3, style = 3)

  # Step 1
  fit <- glmnet::cv.glmnet(x = x, y = y - mean_y,
                           lambda = get_lambda(x = x, y = y - mean_y),
                           intercept = FALSE,
                           standardize = FALSE)
  # residual
  fitted_first <- as.numeric(mean_y + x %*% as.numeric(fit$glmnet.fit$beta[, which.min(fit$cvm)]))
  r <- as.numeric(y - fitted_first)

  # output from step 1
  step1 <- list()
  step1$fitted <- fitted_first
  step1$r <- r
  step1$beta <- as.numeric(fit$glmnet.fit$beta[, which.min(fit$cvm)]) / col_sd
  step1$a0 <- as.numeric(mean_y - crossprod(col_mean, as.matrix(step1$beta)))

  # update progress bar after Step 1
  if (verbose) utils::setTxtProgressBar(pb, 1)

  # Step 2
  # find num_keep higher order terms from
  # (1) squared main effects M^2 + Interaction effects I
  #     (square == FALSE)
  # (2) Interaction effects I (square == TRUE)
  # with largest absolute correlation with the residuals r from first step
  # return the selected variables set B
  if(inherits(x, "sparseMatrix")){
    idx <- screen_sparse_cpp(x = xm, y = r, num_keep = max(num_keep), square = square)
  }
  else{
    idx <- screen_cpp(x = xm, y = r, num_keep = max(num_keep), square = square)
  }

  # update progress bar after Step 2
  if (verbose) utils::setTxtProgressBar(pb, 2)

  # preparing for Step 3
  idx <- idx[order(idx[, 3], decreasing = TRUE), ]

  # if lambda matrix is given
  if(is.null(lambda)){
    lambda <- vector("list", n_num_keep)
  }
  else if (any(unlist(lapply(lambda, is.null)) == TRUE)){
    stop("lambda needs to be a list of tuning parameter values")
  }

  # pre-specify the returns
  out <- vector("list", n_num_keep)

  for(i in seq(n_num_keep)){
    # construct design matrix of selected interactions
    idx_i <- idx[1:num_keep[i], , drop = FALSE]
    # idx_out has three columns
    # which are the j,k indices of nonzero elements
    # main effect index is of form (0, k)
    # squared effect index is of form (k, k)
    # interaction effect index is of form (j, k) for j < k
    # and the third column is the score
    if(square){
      idx_out <- rbind(cbind(rep(0, p), seq(p), rep(0, p)), cbind(seq(p), seq(p), rep(0, p)), idx_i)
    }
    else{
      idx_out <- rbind(cbind(rep(0, p), seq(p), rep(0, p)), idx_i)
    }
    idx_out <- as.matrix(idx_out)
    colnames(idx_out) <- c("index_1", "index_2", "score")
    out[[i]]$idx <- idx_out

    design_i <- myscale(xm[, idx_i[, 1]] * xm[, idx_i[, 2]])

    col_mean_i <- c(col_mean,
                    attr(design_i, which = "scaled:center"))
    col_sd_i <- c(col_sd,
                  attr(design_i, which = "scaled:scale"))
    # the total design matrix
    design_i <- cbind(x, design_i)

    # if type == 1, in Step 3 we fit y on design_i
    # else if type == 2, in Step 3 we fit the residual with design_i
    if(type == 1){
      mean_response <- mean_y
      response <- y - mean_response
    }
    else{
      mean_response <- mean(r)
      response <- r - mean_response
    }

    if(is.null(lambda) | is.null(lambda[[i]])){
      out[[i]]$lambda <- get_lambda(x = design_i, y = response,
                                    nlam = nlam)
      lambda[[i]] <- out[[i]]$lambda
    }
    else{
      out[[i]]$lambda <- lambda[[i]]
    }

    # The Third Step:
    #     run lasso of response y on A + B
    #     corresponding to the best lambda
    fit <- glmnet::glmnet(x = design_i, y = response,
                          lambda = out[[i]]$lambda,
                          intercept = FALSE,
                          standardize = FALSE)
    coef_i <- fit$beta
    colnames(coef_i) <- NULL
    rownames(coef_i) <- NULL
    # drop the names of the matrix object returned by glmnet
    out[[i]]$fitted <- mean_response + design_i %*% coef_i
    # scale estimates back to the original scale of x
    out[[i]]$coef <- coef_i / col_sd_i
    out[[i]]$a0 <- as.numeric(mean_response - crossprod(col_mean_i, as.matrix(out[[i]]$coef)))

    out[[i]]$nzm <- colSums(as.matrix(coef_i[1:p, ] != 0))
    out[[i]]$nzi <- colSums(as.matrix(coef_i[(p + 1): nrow(coef_i), ] != 0))

    # update progress bar after Step 2
    if (verbose)
      utils::setTxtProgressBar(pb, 2 + i)
  }

  cat("\n")
  result <- list(type = type, square = square, step1 = step1, lambda = lambda, num_keep = num_keep, n = n, p = p, step3 = out, call = match.call())
  class(result) <- "sprinter"
  return(result)
}
