\name{AIC.stagewise}
\alias{AIC.stagewise}
\alias{BIC.stagewise}
\alias{logLik.stagewise}



\title{
  AIC, BIC and logLik Method for Stagewise Model.
}

\description{
  AIC and BIC calculate information criterion for each iteration of stagewise model with their respective penalization.
  By default \code{AIC = -2*loglik + 2*df}, where \code{df} is the number of
  non-zero coefficient in model, is computed. BIC uses \code{BIC = -2*loglik + K*df}, where K is the number of the sample size. 
  In the batchwise updating case the number corresponds to the batchsize.
  \code{logLik.stagewise} computes the log likelihood of the model. 
}

\usage{
## AIC/BIC extractor function.
\method{AIC}{stagewise}(object, K = 2, ...)
\method{BIC}{stagewise}(object, ...)

## logLik extractor function.
\method{logLik}{stagewise}(object, mstart = 1, mstop = object$maxit,
  all = TRUE,...)

}

\arguments{
  \item{object}{Object returned from function \code{\link{sdr}}.}
  \item{K}{Numeric, penalization factor in \code{IC = -2*loglik + K*df}.}
  \item{mstart}{Integer, first iteration to be considered.}
  \item{mstop}{Integer, last iteration to be considered.}
  \item{all}{Logical, \code{all = TRUE} returns sequence of loglikelihoods/BIC/AIC.
    \code{all = FALSE} returns mean of the specified range.}
  \item{\dots}{Not used.}
  
}

\details{
  When comparing models to the same data, the smaller the AIC or BIC, the better the fit.
  For the log-likelihood, higher values are better. 
}

\value{
  A vector with class "stagewise_AIC" for BIC and AIC function and a vector of class "logLik.stagewise" for function logLik.
}

\seealso{
  \code{\link{sdr}}.
}

\examples{\dontrun{ 
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
b <- sdr(formula = f, data = d, CF = TRUE,
         family = NO, maxit = 300)

bic <- BIC(b)
aic <- AIC(b)
ll <- logLik(b)


}
}

\keyword{regression}
\keyword{models}

