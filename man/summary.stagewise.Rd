\name{summary.stagewise}
\alias{summary.stagewise}

\title{
  Summary for SDR
}

\description{
  Summary table of stagewise model.
}

\usage{
\method{summary}{stagewise}(object,
  digits = max(3, getOption("digits") - 3),
  mstart = round(0.5 * length(object$logLik)),
  mstop = length(object$logLik), ...)
}

\arguments{
  \item{object}{Object returned from function \code{\link{sdr}}.}
  \item{digits}{Integer, number of displayed digits.}
  \item{mstart}{Integer, first iteration to be considered.
    Coefficients estimates are based on the range from \code{mstart} up to \code{mstop}. }
  \item{mstop}{Integer, last iteration to be considered.
    Coefficients estimates are based on the range of \code{mstart} and \code{mstop}.}
  \item{\dots}{Not used.}
}

\details{
  Computes summary statistics of estimated stagewise model. Coefficient estimates are computed with 
  respect to the specified range of mstart and mstop. Returned are both mean and quantiles of coefficients.
  Additionally log likelihood of mstart and mstop is also returned.
  
}

\seealso{
  \code{\link{sdr}}
}

\examples{
\donttest{
set.seed(123)

## Draw 100 uniform distributed variables.
nobs <- 1000
p <- 100
d <- matrix(runif(nobs * p, -1, 1), ncol = p)

colnames(d) <- paste("x", 1:p, sep = "")
d <- as.data.frame(d)

## Create additive predictors.
d$eta.mu <- d$x1 + 2 * d$x2 + 0.5 * d$x3 - 1*d$x4 
d$eta.sigma <- 0.5 * d$x3 + 0.25 * d$x4 -
  0.25 * d$x5  - 0.5 * d$x6 

## Draw normal distributed
## response for given predictors.
d$y <- rNO(nobs, mu = d$eta.mu,
           sigma = exp(d$eta.sigma))

## Model formula.
f <- as.formula(paste("y ~ ",
                      paste0("x", 1:p, collapse = "+")))
f <- list(f, update(f, sigma ~ .))

## Variable selection with correlation filtering
## and best subset updating.
b <- sdr(formula = f, data = d, CF = TRUE, updating = "bs",
         family = NO, maxit = 300)

## Determine early stopping through BIC.
bic <- BIC(b)
itermax_bic <- which.min(bic)

## Determine starting values for final boosting and
## new model formula corresponding to early stopping.
cb <- coef(b, mstop = itermax_bic, refit = TRUE)
nf <- newformula(b, mstop = itermax_bic, name = "y")

## Refit until convergence.
b_final <- sdr(formula = nf, data = d, updating = "bs",
               family = NO, CF = FALSE, maxit = 300,
               coef_start = cb)


# summary
summary(b_final)
}
}

\keyword{regression}
\keyword{models}
