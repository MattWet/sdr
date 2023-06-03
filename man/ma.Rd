\name{ma}
\alias{ma}

\title{
  Moving Average
}

\description{
  Computes the moving average for every value of a given vector.
}

\usage{
ma(x, order = 20)
}

\arguments{
  \item{x}{A time series.}
  \item{order}{Integer, order of moving average.}
}

\details{
 Calculates for every entry i of a vector a moving average using the following rule:
 If i < order, then compute the average of i up to i + order. Else, compute the average of i - order up to i. 
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

## Batchwise best subset with correlation filtering.
b <- sdr(formula = f, data = d, family = NO, CF = TRUE,
         maxit = 300, batch_ids = 500)
plot(b)

bic <- BIC(b)
itermax_bic <- which.min(ma(bic, order = 20))

cb <- coef(b, mstop = itermax_bic, refit = TRUE)
nf <- newformula(b, mstop = itermax_bic, name = "y")

b_final <- sdr(formula = nf, data = d,
               family = NO, CF = FALSE, maxit = 300,
               coef_start = cb, batch_ids = 500)
plot(b_final)
}
}

\keyword{regression}
\keyword{models}
