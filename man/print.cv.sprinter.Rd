% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.cv.sprinter.R
\name{print.cv.sprinter}
\alias{print.cv.sprinter}
\title{Print the cross validation information of cv.sprinter}
\usage{
\method{print}{cv.sprinter}(fit, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
\item{fit}{A fitted \code{cv.sprinter} object.}

\item{digits}{Significant digits in printout.}
}
\description{
Print a summary of the cross-validation information for running cv.sprinter.
}
\details{
This function takes in a \code{cv.sprinter} object and produces summary of the cross-validation informationabout the tuning parameters (in Step 3) selected by \code{lambda.min} and \code{lambda.1se}.
Adopted from the function \code{print.cv.rgam} from package \code{relgam} by Kenneth Tay and Robert Tibshirani.
}
\examples{
set.seed(123)
n <- 100
p <- 100
x <- matrix(rnorm(n * p), n, p)
y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)

fit.cv <- cv.sprinter(x = x, y = y)
print(fit.cv)

}
\seealso{
\code{\link{cv.sprinter}}, \code{\link{print.printer}}.
}
