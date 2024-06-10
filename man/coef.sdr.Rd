\name{coef.sdr}
\alias{coef.sdr}

\title{
  Coefficients Extractor Function
}

\description{
  Method to extract coefficients of fitted \code{\link{sdr}} object.
  Either all coefficients are extracted or only non-zero coefficients
  needed for refitting routine.
}

\usage{
\method{coef}{sdr}(object, model = NULL, refit = FALSE,
  mstop = NULL, mstart = mstop,...)
}

\arguments{
  \item{object}{Object returned from function \code{\link{sdr}}}
  \item{model}{Selects distributional parameter to be returned. Only works if
    \code{refit = FALSE}. Default is \code{model = NULL}, i.e., all parameters are returned.}
  \item{refit}{If \code{refit = TRUE}, only scaled non-zero coefficients for the
    refitting process are returned. Otherwise all unscaled coefficients are returned.}
  \item{mstop}{Integer, stopping iteration. If not provided, then it is
    set as last iteration.}
  \item{mstart}{Coefficients from \code{mstart} up to \code{mstop} are averaged
    and returned. Default value for \code{mstart = mstop}, i.e., only coefficients
    from iteration \code{mstop} are returned.}
  \item{\dots}{Not used.}
}

\details{
  Our stagewise framework supports a two step modelfitting procedure. The first step
  is to identify the important variables with the standard \code{\link{sdr}} function,
  with the argument \code{updating} any of \code{"bs"}, \code{"cyclic"} or \code{"noncyclic")
  and correlation filtering turned on, \code{CF = TRUE}. Then, extract the 
  important variables with \code{coef(..., refit = TRUE)} and the updated model
  formula with \code{newformula(...)}, (see function \code{\link{newformula}}).
  Finally, refit without correlation filtering and with the updated formula and
  coefficients. See the examples below.
}

\seealso{
  \code{\link{sdr}}.
}

\examples{
\dontrun{set.seed(123)

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

# Final coefficients
cb_final <- coef(b_final)
}}

\keyword{regression}
\keyword{models}

