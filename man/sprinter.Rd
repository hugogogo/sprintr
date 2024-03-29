% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sprinter.R
\name{sprinter}
\alias{sprinter}
\title{Reluctant Interaction Modeling}
\usage{
sprinter(
  x,
  y,
  square = FALSE,
  num_keep = NULL,
  lambda1 = NULL,
  lambda3 = NULL,
  cv_step1 = FALSE,
  nlam1 = 10,
  nlam3 = 100,
  lam_min_ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04),
  ...
)
}
\arguments{
\item{x}{An \code{n} by \code{p} design matrix of main effects. Each row is an observation of \code{p} main effects.}

\item{y}{A response vector of size \code{n}.}

\item{square}{Indicator of whether squared effects should be fitted in Step 1. Default to be FALSE.}

\item{num_keep}{A user specified number of candidate interactions to keep in Step 2. If \code{num_keep} is not specified (as default), it will be set to \code{round[n / log n]}.}

\item{lambda1}{Tuning parameter values for Step 1. \code{lambda1} is a vector. Default to be NULL, and the program will compute its own \code{lambda1} based on \code{nlam1} and \code{lam_min_ratio}.}

\item{lambda3}{Tuning parameter values for Step 3. \code{lambda3} is a matrix, where the k-th column is the list of tuning parameter in Step 3 corresponding to Step 1 using \code{lambda1[k]}. Default to be NULL, and the program will compute its own \code{lambda3} based on \code{nlam3} and \code{lam_min_ratio}.}

\item{cv_step1}{Indicator of whether cross-validation of \code{lambda1} should be carried out in Step 1 before subsequent steps. Default is \code{FALSE}.}

\item{nlam1}{the number of values in \code{lambda1}. If not specified, they will be all set to \code{10}.}

\item{nlam3}{the number of values in each column of \code{lambda3}. If not specified, they will be all set to \code{100}.}

\item{lam_min_ratio}{The ratio of the smallest and the largest values in \code{lambda1} and each column of \code{lambda2}. The largest value is usually the smallest value for which all coefficients are set to zero. Default to be \code{1e-2} in the \code{n} < \code{p} setting.}

\item{...}{other arguments to be passed to the \code{glmnet} calls, such as \code{alpha} or \code{penalty.factor}}
}
\value{
An object of S3 class "\code{sprinter}".
 \describe{
  \item{\code{square}}{The \code{square} parameter passed into sprinter}
  \item{\code{n}}{The number of observations in the dataset}
  \item{\code{p}}{The number of main effects}
  \item{\code{step1}}{The output from fitting Step 1}
  \item{\code{lambda1}}{The path of tuning parameters passed into / computed for fitting Step 1}
  \item{\code{step2}}{The output from the screening Step 2}
  \item{\code{num_keep}}{The path of tuning parameters for Step 2}
  \item{\code{step3}}{The output from fitting Step 3}
  \item{\code{lambda3}}{The path of tuning parameters passed into / computed for fitting Step 3}
  \item{\code{main_center}}{Column centers of the input main effects}
  \item{\code{main_scale}}{Column scales of the input main effects}
  \item{\code{call}}{Function call.}
 }
}
\description{
This is the main function that fits interaction models with a path of tuning parameters (for Step 3).
}
\examples{
set.seed(123)
n <- 100
p <- 100
# dense input
x <- matrix(rnorm(n * p), n, p)
y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)
mod <- sprinter(x = x, y = y)

# sparse input
library(Matrix)
x <- Matrix::Matrix(0, n, p)
idx <- cbind(sample(seq(n), size = 10, replace = TRUE), sample(seq(p), size = 10, replace = TRUE))
x[idx] <- 1
y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)
mod <- sprinter(x = x, y = y)

}
\seealso{
\code{\link{cv.sprinter}}
}
