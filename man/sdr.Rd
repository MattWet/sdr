\name{sdr}
\alias{sdr}


\title{
  Stagewise Distributional Regression
}

\description{
  With this function distributional regression models using an efficient
  stagewise algorithm can be estimated. The function can also be applied
  using very large data sets in combination with a batched variant of the
  algorithm. Please see the examples.
}

\usage{
sdr(formula, family = NULL, data = NULL, batch_ids = NULL,
  updating = c("thresdesc", "bs", "cyclic", "noncyclic"), light = FALSE,
  CF = TRUE, cap = NULL,  caps = NULL,  scalex = TRUE, refitting = TRUE,
  eps = 0.01, nu = 0.1, ncaps = 10, quick_ffdf = FALSE, ...)
}


\arguments{
  \item{formula}{A formula or extended formula, i.e., the formula can be a
    list of formulas where each list entry specifies the details of one parameter of
    the modeled response distribution.}
  \item{family}{A family object (based on the family objects of the \pkg{bamlss} package),
    specifying the details of the modeled distribution such as the parameter names,
    the density function, link functions, etc.}
  \item{data}{A data.frame containing the model response variable(s) and covariates
    specified in the formula. By default the variables are taken from environment(formula):
    typically the environment from which \code{sdr} is called. Full support for the \pkg{ff} package. 
    The package \pkg{bigmemory} is only half supported, meaning if a big.matrix is used as input, then 
    big.matrix must contain all final covariates (e.g. interactions) of the formula. Also, big.matrix requires
    pre scaled covariates.}
  \item{batch_ids}{Either a number for the batchsize for automatically generating 
    random batches or a list of vectors, where each vector corresponds to sample ids.}
  \item{updating}{What updating should be used? Grid search over correlation filtering thresholds with best subset updating 
    "thresdesc", best subset updating "bs", cyclical updating "cyclic" or noncyclical updating "noncyclic".}
  \item{light}{Logical, if set to \code{light = FALSE}, the model frame used for modeling is
    not part of the return value.}
  \item{CF}{Logical for correlation filtering on or off. }
  \item{cap}{Threshold value for correlation filtering. If \code{CF = TRUE} and \code{cap = NULL} 
    leads to automatically setting the threshold value. Setting the threshold value manually 
    overrides the automatic selection for cap.}
  \item{caps}{Only relevant if \code{updating = "threshdesc"}. Explicitly defines different cap values in grid search.
    \code{caps = NULL} then the range is automatically set.}
  \item{scalex}{Logical, standardization of covariates.}
  \item{refitting}{Logical, only relevant if \code{updating = "threshdesc"}. Should a refitting be carried out after selection 
    step? If FALSE, models are used directly at variable selection step.}
 \item{eps}{Steplength parameter for coefficients.}
 \item{nu}{Shrinkage parameter to define the range for the clipping for coefficients.}
 \item{ncaps}{Only relevant if \code{updating = "threshdesc"}. How many different cap values should be used in threshold descent? }
 \item{\dots}{Arguments passed to the fitting function.}
 \item{quick_ffdf}{Logical, if FALSE and ff data.frame is provided, then a ff model.matrix is created. Otherwise ff data.frame 
   must be a model.matrix containing all covariates. }
}

\author{
  Mattias Wetscher \email{Mattias.Wetscher@uibk.ac.at},
  Nikolaus Umlauf \email{Nikolaus.Umlauf@uibk.ac.at}
}

\references{
  Mayr A., Fenske N, Hofner B, Kneib T, Schmid M (2012).
  Generalized additive models for location, scale and shape for
  high-dimensional data - a flexible approach based on boosting.
  \emph{Journal of the Royal Statistical Society, Series C (Applied Statistics)}
  \bold{61}(3), 403--427. \doi{10.1111/j.1467-9876.2011.01033.x}

  Rigby, R.A. and Stasinopoulos, D.M. (2005)
  Generalized Additive Models for Location, Scale and Shape.
  \emph{Journal of the Royal Statistical Society, Series C (Applied Statistics)}
  \bold{54}(3), 507--554. \doi{10.1111/j.1467-9876.2005.00510.x}
  
  Tibshirani RJ (2015).
  A General Framework for Fast Stagewise Algorithms.
  Journal of Machine Learning Research 
  \bold{16}(78), 2543--2588. \url{http://jmlr.org/papers/v16/tibshirani15a.html}

  Umlauf N, Klein N, Zeileis A (2019). BAMLSS: Bayesian Additive Models for Location,
  Scale and Shape (and Beyond). \emph{Journal of Computational and Graphical Statistics},
  \bold{27}(3), 612--627. \doi{10.1080/10618600.2017.1407325}
}

\seealso{
  \code{\link{plot.sdr}}, \code{\link{AIC.sdr}},
  \code{\link{predict.sdr}}, \code{\link{residuals.sdr}},
  \code{\link{plot.sdr_residuals}}, \code{\link{coef.sdr}},
  \code{\link{newformula}}
}

\details{
  Distributional regression models are characterized by the dependency of all
  distributional parameters on covariates in an additive manner. To facilitate model fitting,
  we use stagewise boosting optimization with a semi-constant step length.
}

\value{
  An object of class \code{"sdr"}. Contains, among other, the data used for modelfitting, 
  the evolution of coefficient paths and the evolution of the log-likelihood.
}


\examples{
\donttest{
set.seed(123)

## Draw 100 uniform distributed variables.
nobs <- 2000
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

# testing data
dt <- d[1001:2000,]

# training data
d <- d[1:1000,]


## Model formula.
f <- as.formula(paste("y ~ ",
                      paste0("x", 1:p, collapse = "+")))
f <- list(f, update(f, sigma ~ .))

## Variable selection with correlation filtering
## and best subset updating.
b <- sdr(formula = f, data = d, CF = TRUE, updating = "bs",
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
b_final <- sdr(formula = nf, data = d, updating = "bs",
               family = NO, CF = FALSE, maxit = 300,
               coef_start = cb)

# Plot and summary
plot(b_final)
summary(b_final, mstart = 100, mstop = 300)

# Residual diagnostics
res <- residuals(b_final)
plot(res)

# Prediction and logLik on new testing data
pred <- predict(b_final, newdata = dt, type = "parameter")
l1 <- sum(dnorm(dt$y, mean = pred$mu, sd = pred$sigma, log = TRUE))

## Other flavours.
## Batchwise noncyclical updating with correlation filtering.
b <- sdr(formula = f, data = d, family = NO, CF = TRUE, updating = "bs",
         maxit = 300, batch_ids = 500)

plot(b)

bic <- BIC(b)
itermax_bic <- which.min(ma(bic, order = 20))

cb <- coef(b, mstop = itermax_bic, refit = TRUE)
nf <- newformula(b, mstop = itermax_bic, name = "y")

b_final <- sdr(formula = nf, data = d, updating = "bs",
               family = NO, CF = FALSE, maxit = 300,
               coef_start = cb, batch_ids = 500)

plot(b_final)

# Prediction and logLik on new testing data
pred <- predict(b_final, newdata = dt, type = "parameter")
l2 <- sum(dnorm(dt$y, mean = pred$mu, sd = pred$sigma, log = TRUE))

## Batchwise non cyclic updating with correlation filtering.
b <- sdr(formula = f, data = d, updating = "noncyclic",
         family = NO, CF = TRUE, maxit = 300,
         batch_ids = 500)

plot(b)

bic <- BIC(b)
itermax_bic <- which.min(ma(bic, order = 20))

cb <- coef(b, mstop = itermax_bic, refit = TRUE)
nf <- newformula(b, mstop = itermax_bic, name = "y") 

b_final <- sdr(formula = nf, data = d, updating = "noncyclic",
               family = NO, CF = FALSE, maxit = 300,
               coef_start = cb, batch_ids = 500)

plot(b_final)

# Prediction and logLik on new testing data
pred <- predict(b_final, newdata = dt, type = "parameter")
l3 <- sum(dnorm(dt$y, mean = pred$mu, sd = pred$sigma, log = TRUE))

## Batchwise threshold descent with best subset updating.
## K is penalty parameter in BIC, i.e., BIC = -2*logLik + df*K.
batch_ids = 500
b <- sdr(formula = f, data = d, updating = "thresdesc",
         family = NO, CF = TRUE, maxit = 200, ncaps = 10,
         batch_ids = batch_ids, K = log(batch_ids))

## Coefficient paths from refitting of best model.
plot(b, which = "coefficients")

# Prediction and logLik on new testing data
pred <- predict(b, newdata = dt, type = "parameter")
l4 <- sum(dnorm(dt$y, mean = pred$mu, sd = pred$sigma, log = TRUE))

## Comparison of out of sample logLik
c(l1, l2, l3, l4)
}
}

\keyword{models}
\keyword{regression}

