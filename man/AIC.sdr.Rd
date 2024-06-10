\name{AIC.sdr}
\alias{AIC.sdr}
\alias{BIC.sdr}
\alias{logLik.sdr}

\title{
  AIC, BIC and logLik Method.
}

\description{
  AIC and BIC calculate the information criterion for each iteration of
  stagewise boosted model with the respective penalization.
  By default, \code{AIC = -2 * logLik + 2 * df} is calculated,
  where \code{df} is the number of non-zero coefficients in the model.
  BIC uses \code{BIC = -2 * logLik + K * df}, where K is the logarithm of
  the number of obseravtions in the data used for estimation. 
  In the batchwise updating case the number corresponds to the batchsize.
  \code{logLik.sdr} computes the log likelihood of the model. 
}

\usage{
## AIC/BIC extractor function.
\method{AIC}{sdr}(object, K = 2, ...)
\method{BIC}{sdr}(object, ...)

## logLik extractor function.
\method{logLik}{sdr}(object, mstart = 1, mstop = length(object$logLik),
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

\value{
  A vector of class "sdr_AIC" for object retuened by \code{AIC()} and \code{BIC()}, for
  objects returned ny \code{logLik()} an object of class "logLik.sdr".
}

\seealso{
  \code{\link{sdr}}.
}

\examples{\dontrun{## Set the seed for reproducibility.
set.seed(123)

## Draw 100 uniform distributed variables.
nobs <- 1000
p <- 100
d <- matrix(runif(nobs * p, -1, 1), ncol = p)

colnames(d) <- paste("x", 1:p, sep = "")
d <- as.data.frame(d)

## Create additive predictors.
d$mu <- d$x1 + 2 * d$x2 + 0.5 * d$x3 - 1*d$x4 
d$sigma <- 0.5 * d$x3 + 0.25 * d$x4 -
  0.25 * d$x5  - 0.5 * d$x6 

## Draw normal distributed
## response for given predictors.
d$y <- rNO(nobs, mu = d$mu,
           sigma = exp(d$sigma))

## Model formula.
f <- as.formula(paste("y ~ ",
                      paste0("x", 1:p, collapse = "+")))
f <- list(f, update(f, sigma ~ .))

## Variable selection with correlation filtering
## and best subset updating.
b <- sdr(formula = f, data = d, CF = TRUE, updating = "bs",
         family = NO, maxit = 300)

bic <- BIC(b)
aic <- AIC(b)
ll <- logLik(b)
}}

\keyword{regression}
\keyword{models}

