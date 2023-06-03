\name{alpha2cap2}
\alias{alpha2cap2}

\title{
  Cap Value for Correlation Filtering
}

\description{
 The function computes the cap value for correlation filtering. 
}
\usage{
alpha2cap2(alpha = 0.01, nnobs = NULL,
  nvars = NULL, mean = 0)
}

\arguments{
  \item{alpha}{Numeric, significance level.}
  \item{nnobs}{Integer, number of observations.}
  \item{nvars}{Integer, number of variables.}
  \item{mean}{Numeric, mean of distribution. Use mean = 0.}
}

\value{
  Returns the optimal cap value. 
}

\seealso{
  \code{\link{sdr}}.
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

## Compute cap value for correlation filtering. 
cap <- alpha2cap2(alpha = 0.1, nnobs = nobs,
                  nvars = p, mean = 0)

## Variable selection with correlation filtering
## and best subset updating.
b <- sdr(formula = f, data = d, CF = TRUE, cap = cap,
         family = NO, maxit = 300)


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

# Final coefficients
cb_final <- coef(b_final)
}
}

\keyword{regression}
\keyword{models}
\keyword{htest}

