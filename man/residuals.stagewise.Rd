\name{residuals.stagewise}
\alias{residuals.stagewise}
\alias{plot.stagewise_residuals}

\title{
  Stagewise Model Residuals
}
\description{
  Function to compute quantile and response residuals.
}

\usage{
## Extract residuals.
\method{residuals}{stagewise}(object, type = c("quantile", "response"), nsamps = NULL, ...)

## Plotting method.
\method{plot}{stagewise_residuals}(x, which = c("hist-resid", "qq-resid", "wp"),
  spar = TRUE, ...)
}

\arguments{
  \item{object}{Object returned from function \code{\link{sdr}}.}
  \item{type}{Character, \code{"quantile"} computes median residuals and \code{"response"}
    computes mean residuals with respect to the mu-predictor. The latter only makes sense if
    the mu-predictor corresponds to the mean as for a normal distribution.}
  \item{nsamps}{To get a first feeling it can
    be useful to compute residuals only based on \code{nsamps} samples, i.e., \code{nsamps}
    specifies the number of samples which are extracted on equidistant intervals.}
  \item{x}{Object returned from function \code{residuals.stagewise()}.}
  \item{which}{Should a histogram with kernel density estimates be plotted, a qq-plot or a worm plot? 
  Any subset of c("hist-resid", "qq-resid", "wp") can be specified. Per default all three are plotted.}
  \item{spar}{Should graphical parameters be set by the plotting function?}
  \item{\dots}{Not used.}
}

\details{
  Response residuals are the raw residuals, i.e., the response data minus the fitted distributional
  mean. If the distributional family object contains a function \code{$mu(par, \dots)}, then
  raw residuals are computed with \code{y - mu(par)} where \code{par} is the named list of fitted
  values of distributional parameters. If \code{$mu(par, ...)} is missing, then the fitted values
  of the first distributional parameter are used.

  Randomized quantile residuals are based on the cumulative distribution function of the
  distributional family object, i.e., the \code{$p(y, par, ...)} function.
}

\references{
  Dunn P. K., and Smyth G. K. (1996). Randomized Quantile Residuals.
    \emph{Journal of Computational and Graphical Statistics} \bold{5}, 236--244.

  van Buuren S., and Fredriks M. (2001) Worm Plot: Simple Diagnostic Device for Modelling Growth
    Reference Curves. \emph{Statistics in Medicine}, \bold{20}, 1259--1277
}

\seealso{
  \code{\link{sdr}}
}

\examples{
\dontrun{
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


# Residual diagnostics
res <- residuals(b_final)
plot(res)
}
}

\keyword{regression}
\keyword{models}

