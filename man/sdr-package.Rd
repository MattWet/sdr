\name{sdr-package}
\alias{sdr-package}
\alias{sdr}
\docType{package}

\title{
  \packageTitle{sdr}
}

\description{
  \packageDescription{sdr}
}

\details{
  Estimation of distributional regression using an efficient stagewise algorithm is performed
  by function \code{\link{sdr}}.

  In addition, extractor functions like \code{\link[stats]{coef}}, \code{\link[stats]{residuals}},
  \code{\link[stats]{predict}}, \code{\link[stats]{logLik}}, \code{\link[graphics]{plot}},
  etc., are provided.

  \packageIndices{stagewise}
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


\keyword{ package }
\seealso{
  \code{\link{sdr}}, \code{\link{predict.stagewise}}.
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

# Prediction
pred <- predict(b_final)
mse1 <- mean((d$y - pred$mu)^2)

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
pred <- predict(b_final)
mse2 <- mean((d$y - pred$mu)^2)


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
pred <- predict(b_final)
mse3 <- mean((d$y - pred$mu)^2)

## Batchwise threshold descent with best subset updating.
## K is penalty parameter in BIC, i.e. BIC = -2*logLik + df*K
batch_ids = 500
b <- sdr(formula = f, data = d, updating = "thresdesc",
         family = NO, CF = TRUE, maxit = 100, ncaps = 10,
         batch_ids = batch_ids, K = log(batch_ids))
# Coef-paths from refitting of best model
plot(b, which = "coefficients")

pred <- predict(b)
mse4 <- mean((d$y - pred$mu)^2)

c(mse1,mse2,mse3,mse4)

}
}
