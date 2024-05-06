\name{stagewise-package}
\alias{stagewise-package}
\docType{package}

\title{
  \packageTitle{stagewise}
}

\description{
  \packageDescription{stagewise}
}

\details{
  Distributional regression models are estimated using an efficient stagewise algorithm
  implemented in the function \code{\link{sdr}}.

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

\keyword{package}

\seealso{
  \code{\link{sdr}}, \code{\link{predict.stagewise}}.
}

\examples{
## See the ?sdr manual page.
}

