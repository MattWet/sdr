\name{plot.stagewise}
\alias{plot.stagewise}

\title{
  Stagewise Object Plots
}

\description{
  A function that shows the evolution of the estimated coefficients or the corresponding
  AIC.
}

\usage{
\method{plot}{stagewise}(x, which = c("all", "coefficients", "AIC"), K = 2,
                           bw = 0, spar = TRUE, ...)
}

\arguments{
  \item{x}{Object returned from function \code{\link{sdr}}.}
  \item{which}{Specifies which plots are displayed. "coefficients", "AIC" or both with "all".}
  \item{K}{Specifies the penalty parameter in the information criterion. The default K = 2 corresponds to AIC. }
  \item{bw}{Specifies the moving average window considered for computing the information criterion. The default, 
  bw = 0, means no moving average is considered. Instead the iteration corresponding to the minimum over 
  all information criterion values is used. The chosen iteration is displayed with a vertical line.}
  \item{spar}{Should graphical parameters be set by the plotting function?}
  \item{\dots}{Not used.}
}

\details{
  The AIC incorporates the number of non-zero coefficients as degree of freedom. 
  
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
b <- sdr(formula = f, data = d, CF = TRUE,
         family = NO, maxit = 300)
plot(b)

## Determine early stopping through BIC.
bic <- BIC(b)
itermax_bic <- which.min(bic)

## Determine starting values for final boosting and
## new model formula corresponding to early stopping.
cb <- coef(b, mstop = itermax_bic, refit = TRUE)
nf <- newformula(b, mstop = itermax_bic, name = "y")

## Refit until convergence.
b_final <- sdr(formula = nf, data = d,
               family = NO, CF = FALSE, maxit = 300,
               coef_start = cb)

# Plot
plot(b_final)
}
}

\keyword{regression}
\keyword{models}
