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
  Tibshirani RJ (2015).
  A General Framework for Fast Stagewise Algorithms.
  Journal of Machine Learning Research 
  \bold{16}(78), 2543--2588. \url{http://jmlr.org/papers/v16/tibshirani15a.html}

  Wetscher M, Seiler J, Stauffer R, Umlauf N (2024).
  Stagewise Boosting Distributional Regression.
  \emph{arXiv}. \doi{10.48550/arXiv.2405.18288}
}

\keyword{package}

\seealso{
  \code{\link{sdr}}, \code{\link{predict.sdr}}.
}

\examples{
## See the ?sdr manual page.
}

