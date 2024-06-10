\name{newformula}
\alias{newformula}

\title{
  Newformula function
}

\description{
  Updating of model formula after variable selection.
}

\usage{
newformula(object, mstop = NULL, name = NULL)
}

\arguments{
  \item{object}{Object returned from function \code{\link{sdr}}.}
  \item{mstop}{Integer, stopping iteration. If not provided, then it is
    set as last iteration.}
  \item{name}{Character, provide name of response variable,
    otherwise \code{name = "y"} is used.}
}

\details{
  Our stagewise framework supports a two step modelfitting procedure. The first step
  is to identify the important variabels with the standard \link{\link{sdr}} function,
  with the argument \code{updating} any of \code{"bs"}, \code{"cyclic"} or \code{"noncyclic"}
  and correlation filtering turned on, \code{CF = TRUE}. Then, extract the important 
  variables with \code{coef(..., refit = TRUE)} and the updated model formula with
  \code{newformula(...)}. Finally, refit without correlation filtering and with
  the updated formula and coefficients. See example below.
}

\seealso{
  \code{\link{sdr}}.
}

\examples{\dontrun{set.seed(123)

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
}}

\keyword{regression}
\keyword{models}

