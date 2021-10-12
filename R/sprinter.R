#' Reluctant Interaction Modeling
#'
#' This is the main function that fits interaction models with a path of tuning parameters (for Step 3).
#' @param x An \code{n} by \code{p} design matrix of main effects. Each row is an observation of \code{p} main effects.
#' @param y A response vector of size \code{n}.
#' @param square Indicator of whether squared effects should be fitted in Step 1. Default to be FALSE.
#' @param num_keep A user specified number of candidate interactions to keep in Step 2. If \code{num_keep} is not specified (as default), it will be set to \code{round[n / log n]}.
#' @param lambda1 Tuning parameter values for Step 1. \code{lambda1} is a vector. Default to be NULL, and the program will compute its own \code{lambda1} based on \code{nlam1} and \code{lam_min_ratio}.
#' @param lambda3 Tuning parameter values for Step 3. \code{lambda3} is a matrix, where the k-th column is the list of tuning parameter in Step 3 corresponding to Step 1 using \code{lambda1[k]}. Default to be NULL, and the program will compute its own \code{lambda3} based on \code{nlam3} and \code{lam_min_ratio}.
#' @param cv_step1 Indicator of whether cross-validation of \code{lambda1} should be carried out in Step 1 before subsequent steps. Default is \code{FALSE}.
#' @param nlam1 the number of values in \code{lambda1}. If not specified, they will be all set to \code{10}.
#' @param nlam3 the number of values in each column of \code{lambda3}. If not specified, they will be all set to \code{100}.
#' @param lam_min_ratio The ratio of the smallest and the largest values in \code{lambda1} and each column of \code{lambda2}. The largest value is usually the smallest value for which all coefficients are set to zero. Default to be \code{1e-2} in the \code{n} < \code{p} setting.
#' @param ... other arguments to be passed to the \code{glmnet} calls, such as \code{alpha} or \code{penalty.factor}
#'
#' @return An object of S3 class "\code{sprinter}".
#'  \describe{
#'   \item{\code{square}}{The \code{square} parameter passed into sprinter}
#'   \item{\code{n}}{The number of observations in the dataset}
#'   \item{\code{p}}{The number of main effects}
#'   \item{\code{step1}}{The output from fitting Step 1}
#'   \item{\code{lambda1}}{The path of tuning parameters passed into / computed for fitting Step 1}
#'   \item{\code{step2}}{The output from the screening Step 2}
#'   \item{\code{num_keep}}{The path of tuning parameters for Step 2}
#'   \item{\code{step3}}{The output from fitting Step 3}
#'   \item{\code{lambda3}}{The path of tuning parameters passed into / computed for fitting Step 3}
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
sprinter <- function(x, y, square = FALSE, num_keep = NULL,
                     lambda1 = NULL, lambda3 = NULL,
                     cv_step1 = FALSE,
                     nlam1 = 10, nlam3 = 100,
                     lam_min_ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04), ...){
  n <- nrow(x)
  p <- ncol(x)

  q <- ifelse(square, p * (p - 1) / 2, p * (p - 1) / 2 + p)
  if(is.null(num_keep))
    num_keep <- ceiling(n / log(n))
  else
    stopifnot(num_keep > 0 & num_keep <= q)

  # we always standardize the design matrix to get main effects
  # squared effects and interactions are built upon standardized main effects
  x <- myscale(x)
  # xm is the (standardized) design matrix of main effects
  col_mean <- attr(x = x, which = "scaled:center")
  col_sd <- attr(x = x, which = "scaled:scale")
  xm <- x
  mean_y <- mean(y)

  # The First Step
  # run lasso on
  # (1) main effects M (square == FALSE)
  # (2) main effects M + squared main effects M^2 (square == TRUE)
  # return the fitted value of response
  if(square){
    x_sq <- myscale(x^2)

    col_mean <- c(col_mean, attr(x_sq, which = "scaled:center"))
    col_sd <- c(col_sd, attr(x_sq, which = "scaled:scale"))
    x <- cbind(x, x_sq)
  }

  # initiate lambda 1
  if(is.null(lambda1))
    lambda1 <- get_lambda(x = x, y = y - mean_y, nlam = nlam1, lam_min_ratio = lam_min_ratio)

  # Step 1
  if(cv_step1){
    fit <- glmnet::cv.glmnet(x = x, y = y - mean_y,
                             lambda = lambda1,
                             intercept = FALSE,
                             standardize = FALSE, ...)

    # grab coefficient estimate
    theta <- matrix(fit$glmnet.fit$beta[, which.min(fit$cvm)], ncol = 1)
    colnames(theta) <- NULL
    rownames(theta) <- NULL

    # fitted values
    fitted <- matrix(mean_y + x %*% theta, ncol = 1)
    r <- matrix(y - fitted, ncol = 1)

    # update lambda1
    lambda1 <- fit$lambda[which.min(fit$cvm)]
  }
  else{
    fit <- glmnet::glmnet(x = x, y = y - mean_y,
                          lambda = lambda1,
                          intercept = FALSE,
                          standardize = FALSE, ...)

    # grab coefficient estimate
    theta <- fit$beta
    colnames(theta) <- NULL
    rownames(theta) <- NULL

    # fitted values
    fitted <- as.matrix(mean_y + x %*% theta)
    # residual is a matrix
    r <- y - fitted
  }

  # update lambda1 and initiate lambda3
  nlam1 <- length(lambda1)
  if(is.null(lambda3))
    lambda3 <- matrix(NA, nlam3, nlam1)

  stopifnot(is.matrix(lambda3))
  nlam3 <- nrow(lambda3)

  # output from step 1
  step1 <- list()
  step1$fitted <- fitted
  step1$r <- r
  step1$beta <- as.matrix(theta / col_sd)
  step1$a0 <- as.numeric(mean_y - crossprod(col_mean, step1$beta))

  # pre-specify the returns
  step2 <- vector("list", nlam1)
  step3 <- vector("list", nlam1)

  for(k in seq(nlam1)){
    # Step 2
    # find num_keep higher order terms from
    # (1) squared main effects M^2 + Interaction effects I
    #     (square == FALSE)
    # (2) Interaction effects I (square == TRUE)
    # with largest absolute correlation with the residuals r from first step
    # return the selected variables set B
    if(inherits(x, "sparseMatrix")){
      idx <- screen_sparse_cpp(x = xm, y = r[, k],
                               num_keep = num_keep, square = square)
    }
    else{
      idx <- screen_cpp(x = xm, y = r[, k],
                        num_keep = num_keep, square = square)
    }

    # preparing for Step 3
    idx <- idx[order(idx[, 3], decreasing = TRUE), , drop = FALSE]
    colnames(idx) <- c("index_1", "index_2", "score")
    step2[[k]] <- idx

    if(nrow(idx) == 1)
      design <- myscale(matrix(xm[, idx[, 1]] * xm[, idx[, 2]], ncol = 1))
    else
      design <- myscale(xm[, idx[, 1]] * xm[, idx[, 2]])

    col_mean3 <- c(col_mean,
                  attr(design, which = "scaled:center"))
    col_sd3 <- c(col_sd,
                attr(design, which = "scaled:scale"))
    # the total design matrix
    design <- cbind(x, design)

    # in Step3, we fit the residual from Step 1
    mean_response <- mean(r[, k])
    response <- r[, k] - mean_response

    if(any(is.na(lambda3[, k])))
      lambda3[, k] <- get_lambda(x = design, y = response,
                                 nlam = nlam3,
                                 lam_min_ratio = lam_min_ratio)

    # The Third Step:
    #     run lasso of response y on A + B
    #     corresponding to the best lambda
    fit <- glmnet::glmnet(x = design, y = response,
                          lambda = lambda3[, k],
                          intercept = FALSE,
                          standardize = FALSE, ...)
    coef <- fit$beta
    colnames(coef) <- NULL
    rownames(coef) <- NULL
    # drop the names of the matrix object returned by glmnet
    step3[[k]]$fitted <- mean_response + design %*% coef
    # scale estimates back to the original scale of design
    step3[[k]]$coef <- as.matrix(coef / col_sd3)
    step3[[k]]$a0 <- as.numeric(mean_response - crossprod(col_mean3, step3[[k]]$coef))

    # number of non-zero main effects & interactions
    step3[[k]]$nzm <- colSums(as.matrix(coef[1:p, ] != 0))
    step3[[k]]$nzi <- colSums(as.matrix(coef[(p + 1): nrow(coef), ] != 0))
  }

  result <- list(square = square, n = n, p = p, cv_step1 = cv_step1,
                 step1 = step1, lambda1 = lambda1,
                 step2 = step2, num_keep = num_keep,
                 step3 = step3, lambda3 = lambda3,
                 call = match.call())
  class(result) <- "sprinter"
  return(result)
}
