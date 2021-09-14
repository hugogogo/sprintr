#' Two-stage hierarchical lasso
#'
#' An implementation of the two-stage lasso studied in Hao et, al (2018).
#'
#' @param x An \code{n} by \code{p} design matrix of main effects. Each row is an observation of \code{p} main effects.
#' @param y A response vector of size \code{n}.
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
hier_lasso <- function(x, y,
                       lambda = NULL, nlam = 100,
                       lam_choice = "min",
                       lam_min_ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04),
                       nfold = 5, foldid = NULL){
  n <- nrow(x)
  p <- ncol(x)
  stopifnot(n == length(y))
  stopifnot(lam_choice == "min" | lam_choice == "1se")

  x <- myscale(x)
  # first fit the sprinter using all data with all lambdas
  fit <- cv.glmnet(x = x, y = y, standardize = FALSE)

  # construct all pairs interactions from selected main effects by CV-lasso
  if (lam_choice == "min"){
    best <- fit$glmnet.fit$beta[, which.min(fit$cvm)]
  }
  else if(lam_choice == "1se"){
    best <- fit$glmnet.fit$beta[, which(fit$lambda == fit$lambda.1se)]
  }
  names(best) <- NULL

  #main_idx <- which(best != 0)
  main_idx <- head(sort(abs(best), decreasing = TRUE, index.return = TRUE)$ix, 500)
  #cat(paste("number of selected main effects = ", length(main_idx)), fill = TRUE)
  if(length(main_idx) != 0){
    int_idx <- t(combn(main_idx, 2))
    #cat(paste("number of constructed interactions = ", nrow(int_idx)), fill = TRUE)

    # step 2: fit cv.glmnet with all main effects and constructed (hierarchical) interactions
    xx <- x[, int_idx[, 1]] * x[, int_idx[, 2]]
    # note that we use penalty.factor to make sure that main effects are not penalized
    fit <- cv.glmnet(x = cbind(x[, main_idx], xx), y = y, penalty.factor = c(rep(0, length(main_idx)), rep(1, nrow(int_idx))), standardize = FALSE)

    coef <- fit$glmnet.fit$beta[, which.min(fit$cvm)]
    a0 <- fit$glmnet.fit$a0[which.min(fit$cvm)]
    idx <- rbind(cbind(rep(0, length(main_idx)), main_idx), int_idx)
    compact <- cbind(idx[which(coef != 0), , drop = FALSE], coef[coef != 0])
    rownames(compact) <- NULL
    colnames(compact) <- c("index_1", "index_2", "coefficient")
  }
  else{
    a0 <- fit$glmnet.fit$a0[which.min(fit$cvm)]
    coef <- fit$glmnet.fit$beta[, which.min(fit$cvm)]
    idx <- rbind(cbind(rep(0, length(main_idx)), main_idx))
    compact <- NULL
  }

  # finally return the best lambda
  out <- list(n = n,
              p = p,
              fit = fit,
              a0 = as.numeric(a0),
              compact = compact,
              call = match.call())
  class(out) <- "other"
  return(out)
}
