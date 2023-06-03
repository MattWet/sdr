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
## Model fitting function.
sdr(formula, family = NULL, data = NULL, weights = NULL, batch_ids = NULL,
  cyclic = FALSE, subset = NULL, offset = NULL, contrasts = NULL, model = TRUE, 
  CF = FALSE, cap = NULL,  scalex = TRUE, boost = TRUE,
  eps = 0.01, nu = 0.1, ...)

}

\arguments{
  \item{formula}{A formula or extended formula, i.e., the formula can be a
    list of formulas where each list entry specifies the details of one parameter of
    the modeled response distribution.}
  \item{family}{A family object (based on the family objects of the \pkg{bamlss} package),
    specifying the details of the modeled distribution such as the parameter names,
    the density function, link functions, etc.}
  \item{data}{A data.frame or list containing the model response variable(s) and covariates
    specified in the formula. By default the variables are taken from environment(formula):
    typically the environment from which \code{sdr} is called.}
  \item{weights}{Prior weights on the data.}
  \item{batch_ids}{Either a number for the batchsize for automatically generating 
  random batches or a list of vectors, where each vector corresponds to sample ids.}
  \item{cyclic}{Logical, if \code{cyclic = TRUE}, cyclic boosting is used for updating.}
  \item{subset}{An optional vector specifying a subset of observations to be used
    in the fitting process.}
  \item{offset}{Can be used to supply model offsets for use in fitting.}
  \item{contrasts}{An optional list. See the contrasts.arg of model.matrix.default.}
  \item{model}{Logical, if set to \code{model = FALSE}, the model frame used for modeling is
    not part of the return value.}
  \item{CF}{Logical for correlation filtering on or off. Works only with \code{boost = T}}
  \item{cap}{Threshold value for correlation filtering. If \code{CF = TRUE}, then \code{cap = NULL} 
  leads to automatically setting the threshold value. Setting the threshold value manually 
  overrides the automatic selection for cap.}
  \item{scalex}{Logical, standardize covariates.}
  \item{boost}{Logical, specify \code{boost = TRUE} and \code{cyclic = FALSE}
    for best subset boosting, \code{boost = TRUE} and \code{cyclic = TRUE} for cyclical
    boosting or \code{boost = FALSE} for backfitting updating (updating of every parameter in 
    every iteration).}
 \item{eps}{Steplength parameter for coefficients.}
 \item{nu}{Shrinkage parameter to define the range for the clipping for coefficients.}
 \item{\dots}{Arguments passed to the fitting function.}
}

\author{
  Mattias Wetscher \email{Mattias.Wetscher@uibk.ac.at},
  Nikolaus Umlauf \email{Nikolaus.Umlauf@uibk.ac.at}
}

\references{
  Umlauf N, Klein N, Zeileis A (2019). BAMLSS: Bayesian Additive Models for Location,
  Scale and Shape (and Beyond). \emph{Journal of Computational and Graphical Statistics},
  \bold{27}(3), 612--627. \doi{10.1080/10618600.2017.1407325}

  Umlauf N, Klein N, Simon T, Zeileis A (2021).
  bamlss: A Lego Toolbox for Flexible Bayesian Regression (and Beyond).
  \emph{Journal of Statistical Software},
  \bold{100}(4), 1--53. \doi{10.18637/jss.v100.i04}
  
  Mayr A., Fenske N, Hofner B, Kneib T, Schmid M (2012).
     Generalized additive models for location, scale and shape for
     high-dimensional data - a flexible approach based on boosting.
     \emph{Journal of the Royal Statistical Society, Series C (Applied Statistics)}
     \bold{61}(3), 403--427.

Rigby, R.A. and Stasinopoulos, D.M. (2005)
     Generalized Additive Models for Location, Scale and Shape.
     \emph{Journal of the Royal Statistical Society, Series C (Applied Statistics)}
     \bold{54}(3), 507--554.
  
  Tibshirani RJ (2015).
  A General Framework for Fast Stagewise Algorithms.
  Journal of Machine Learning Research 
  \bold{16}(78), 2543--2588.
  

  
}

\seealso{
  \code{\link{plot.stagewise}}, \code{\link{AIC.stagewise}},
  \code{\link{predict.stagewise}}, \code{\link{residuals.stagewise}},
  \code{\link{plot.stagewise_residuals}}, \code{\link{coef.stagewise}},
  \code{\link{newformula}}
}

\details{
  Distributional regression models are models in which all distributional parameters 
  can depend in an additive way on covariates. We incorporate stagewise 
  boosting optimization for modelfitting with semi-constant step length.
}

\value{
  An object of class "stagewise". Contains, among other, the data used for modelfitting, 
  the evolution of coefficient paths and the evolution of the log-likelihood.
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

# Plot and summary
plot(b_final)
summary(b_final, mstart = 100, mstop = 300)

# Residual diagnostics
res <- residuals(b_final)
plot(res)

# Prediction
pred <- predict(b_final)
mse1 <- mean((d$y - pred$mu)^2)

## Other flavours.
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
pred <- predict(b_final)
mse2 <- mean((d$y - pred$mu)^2)


## Batchwise non cyclic with correlation filtering.
## length_ps = 1 means that only subsets with one element
## are considered in the updating selection. 
## This is equivalent to the non cyclic updating.
b <- sdr(formula = f, data = d,
         family = NO, CF = TRUE, maxit = 300,
         batch_ids = 500, length_ps = 1)
plot(b)

bic <- BIC(b)
itermax_bic <- which.min(ma(bic, order = 20))

cb <- coef(b, mstop = itermax_bic, refit = TRUE)
nf <- newformula(b, mstop = itermax_bic, name = "y") 

b_final <- sdr(formula = nf, data = d,
               family = NO, CF = FALSE, maxit = 300,
               coef_start = cb, batch_ids = 500, length_ps = 1)
plot(b_final)
pred <- predict(b_final)
mse3 <- mean((d$y - pred$mu)^2)

c(mse1,mse2,mse3)
}
}

\keyword{models}
\keyword{regression}

