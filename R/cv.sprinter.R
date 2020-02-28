#' Running sprinter with cross-validation
#'
#' The main cross-validation function to select the best sprinter fit for a path of tuning parameters.
#'
#' @param x An \code{n} by \code{p} design matrix of main effects. Each row is an observation of \code{p} main effects.
#' @param y A response vector of size \code{n}.
#' @param num_keep Number of candidate interactions to keep in Step 2. If \code{num_keep} is not specified (as default), it will be set to \code{[n / log n]}.
#' @param square Indicator of whether squared effects should be fitted in Step 1. Default to be FALSE.
#' @param num_keep A user specified list of number of candidate interactions to keep in Step 2. If \code{num_keep} is not specified (as default), it will be set to a sequence from \code{[n / log n]} to 0, where the length of the sequence is set to the default value of \code{n_num_keep}.
#' @param n_num_keep The number of \code{num_keep} values. Default to be \code{5}.
#' @param lambda A user specified list of tuning parameter. \code{lambda} is a list object of length \code{n_num_keep}, and an element \code{lambda[[i]]} is a vector of length \code{nlam[i]}. Default to be NULL, and the program will compute its own \code{lambda} paths based on \code{n_num_keep}, \code{nlam} and \code{lam_min_ratio}.
#' @param nlam A vector of length \code{n_num_keep}, where \code{nlam[i]} is the number of values in \code{lambda[[i]]}. If not specified, they will be all set to \code{100}.
#' @param lam_min_ratio The ratio of the smallest and the largest values in each \code{lambda[[i]]}. The largest value in \code{lambda} is usually the smallest value for which all coefficients are set to zero. Default to be \code{1e-2} in the \code{n} < \code{p} setting.
#' @param nfold Number of folds in cross-validation. Default value is 5. If each fold gets too view observation, a warning is thrown and the minimal \code{nfold = 3} is used.
#' @param foldid A vector of length \code{n} representing which fold each observation belongs to. Default to be \code{NULL}, and the program will generate its own randomly.
#' @return An object of S3 class "\code{sprinter}".
#'  \describe{
#'   \item{\code{n}}{The sample size.}
#'   \item{\code{p}}{The number of main effects.}
#'   \item{\code{a0}}{estimate of intercept corresponding to the CV-selected model.}
#'   \item{\code{compact}}{A compact representation of the selected variables. \code{compact} has three columns, with the first two columns representing the indices of a selected variable (main effects with first index = 0), and the last column representing the estimate of coefficients.}
#'   \item{\code{fit}}{The whole \code{glmnet} fit object in Step 3.}
#'   \item{\code{fitted}}{fitted value of response corresponding to the CV-selected model.}
#'   \item{\code{lambda}}{The sequence of \code{lambda} values used.}
#'   \item{\code{cvm}}{The averaged estimated prediction error on the test sets over K folds.}
#'   \item{\code{cvse}}{The standard error of the estimated prediction error on the test sets over K folds.}
#'   \item{\code{foldid}}{Fold assignment. A vector of length \code{n}.}
#'   \item{\code{ibest}}{The index in \code{lambda} that is chosen by CV.}
#'   \item{\code{call}}{Function call.}
#'  }
#' @seealso
#'   \code{\link{predict.cv.sprinter}}
#' @examples
#' n <- 100
#' p <- 100
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)
#' mod <- cv.sprinter(x = x, y = y)
#'
#' @import glmnet
#' @export
cv.sprinter <- function(x, y, square = FALSE, type = 1,
                        num_keep = NULL, n_num_keep = 5,
                        lambda = NULL, nlam = 100,
                        lam_min_ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04),
                        nfold = 5, foldid = NULL){
  n <- nrow(x)
  p <- ncol(x)
  stopifnot(n == length(y))

  # by default, the maximum value of num_keep is set to be [n / log(n)]
  if(is.null(num_keep)){
    num_keep = ceiling(seq(ceiling(n / log(n)), 1, length.out = n_num_keep))
  }
  else{
    stopifnot(max(num_keep) <= ifelse(square, p*(p - 1) / 2, p * (p - 1)/ 2 + p))
  }
  n_num_keep <- length(num_keep)

  cat("cv initial:", fill = TRUE)
  # first fit the sprinter using all data with all lambdas
  fit <- sprinter(x = x, y = y, square = square, type = type,
                  num_keep = num_keep, n_num_keep = n_num_keep,
                  nlam = nlam, lam_min_ratio = lam_min_ratio)

#  if(is.null(lambda)){
#    # if lambda is not provided
#    lambda <- fit$lambda
#  }
#  nlam <- length(lambda)

  # use cross-validation to select the best lambda
  if (is.null(foldid)){
    # foldid is a vector of values between 1 and nfold
    # identifying what fold each observation is in.
    # If supplied, nfold can be missing.
    foldid <- sample(seq(nfold), size = n, replace = TRUE)
  }

  # mse of lasso estimate of coef
  err <- vector("list", nfold)
  err_mat <- matrix(0, nrow = n_num_keep, ncol = nlam)
  for (i in seq(nfold)){
    cat(paste("cv fold ", i, ":"), fill = TRUE)
    # train on all but i-th fold
    id_tr <- (foldid != i)
    id_te <- (foldid == i)
    # training/testing data set
    # standardize the training data
    x_tr <- myscale(x[id_tr, ])
    # use the scaling for training data
    x_te <- myscale(x[id_te, ],
                  center = attr(x = x_tr, which = "scaled:center"),
                  scale = attr(x = x_tr, which = "scaled:scale"))
    y_tr <- y[id_tr]
    y_te <- y[id_te]

    # get the fit using lasso on training data
    # and fit/obj on the test data
    fit_tr <- sprinter(x = x_tr, y = y_tr,
                       square = square, type = type,
                       num_keep = fit$num_keep,
                       lambda = fit$lambda)
    err[[i]] <- matrix(NA, nrow = n_num_keep, ncol = nlam)
    pred_te <- predict.sprinter(object = fit_tr, newdata = x_te)
    for(k in seq(length(fit_tr$step3))){
      err[[i]][k, ] <- sqrt(as.numeric(colMeans((y_te - pred_te[[k]]$fitted)^2)))
    }
    err_mat <- err_mat + err[[i]]
  }

  # extract information from CV
  # the mean cross-validated error, a n_num_keep-by-nlam matrix
  err_mat <- err_mat / nfold

  # also compute the se_mat
  se_mat <- matrix(0, nrow = n_num_keep, ncol = nlam)
  for (i in seq(nfold)){
    se_mat <- se_mat + (err[[i]] - err_mat)^2
  }
  se_mat <- sqrt(se_mat / (nfold - 1))

  # the index of best lambda
  ibest <- tail(which(err_mat == min(err_mat), arr.ind = TRUE), 1)
  colnames(ibest) <- NULL

  num_keep_best <- ibest[1]
  lam_best <- ibest[2]
  lam_1se <- which(err_mat[num_keep_best, ] < err_mat[num_keep_best, lam_best] + se_mat[num_keep_best, lam_best])[1]

  idx <- fit$step3[[num_keep_best]]$idx
  coef <- fit$step3[[num_keep_best]]$coef[, lam_best]
  compact <- cbind(idx[which(coef != 0), 1:2, drop = FALSE], coef[coef != 0])
  colnames(compact) <- c("index_1", "index_2", "coefficient")

  fitted <- fit$step3[[num_keep_best]]$fitted[, lam_best]
  if(type == 2)
    fitted <- fit$step1$fitted + fitted

  # finally return the best lambda
  out <- list(n = n,
              p = p,
              type = type,
              square = square,
              step1 = fit$step1,
              a0 = fit$step3[[num_keep_best]]$a0[lam_best],
              compact = compact,
              fit = fit,
              fitted = fitted,
              lambda = fit$lambda,
              num_keep = num_keep,
              cvm = err_mat,
              cvse = se_mat,
              foldid = foldid,
              num_keep_best = num_keep_best,
              lam_best = lam_best,
              lam_1se = lam_1se,
              call = match.call())
  class(out) <- "cv.sprinter"
  return(out)
}
