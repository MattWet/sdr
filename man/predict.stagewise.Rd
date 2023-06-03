\name{predict.stagewise}
\alias{predict.stagewise}

\title{
  Compute Predictions
}

\description{
  Function to compute predictions given new data from fitted stagewise object. 
}

\usage{
\method{predict}{stagewise}(object, newdata = NULL, type = c("link", "parameter"), 
drop = TRUE, mstart = NULL, mstop = NULL, model = NULL,
   ...)
}

\arguments{
  \item{object}{Object returned from function \code{\link{sdr}}.}
  \item{newdata}{A data frame or list containing the values of the model
    covariates at which predictions are required. Note that depending on argument term,
    only covariates that are needed by the corresponding model terms need to be supplied.}
  \item{type}{Character, prediction for predictor (\code{"link"}) or
    distributional parameter (\code{"parameter"}) separately.}
  \item{drop}{Logical, if predictions for only one model are returned,
    the list structure is dropped.}
  \item{mstart}{Integer, prediction is based on mean of coefficients from
    iteration \code{mstart} up to \code{mstop}. If no mstart is provided, mstart is set to median of iterations.}
  \item{mstop}{Integer, prediction is based on mean of coefficients from
    iteration \code{mstart} up to \code{mstop}. If no mstop is provided, mstop is set to last iteration.}
  \item{model}{Character or integer, specifies the model for which predictions should be computed.}
  \item{\dots}{Not used.}
}

\value{
  The predictions are either at the linear predictor scale or transformed to the distributional
  parameter scale. 
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


# Prediction
pred <- predict(b_final)
mse1 <- mean((d$y - pred$mu)^2)
}
}

\keyword{regression}
\keyword{models}

