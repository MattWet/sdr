\name{binmm}
\alias{binmm}
\alias{read.binmm}
\alias{dim.binmm}
\alias{nrow.binmm}
\alias{ncol.binmm}
\alias{dimnames.binmm}
\alias{head.binmm}
\alias{tail.binmm}
\alias{print.binmm}
\alias{summary.binmm}
\alias{[.binmm}

\title{
    TODO(R) REDO MANUAL Class for Handling Large Model Matrices stored as CSV
}

\description{
    The stagewise algorithm \code{\link{sdr}} allows to estimate
    regression models using batches. In each iteration only a (random)
    subset of the data set is used.

    Thus, we never have to store the whole data set in memory at any
    given time. The \code{binmm} class handles a connection to a CSV
    file (see Details) containing (large) model matrices and to efficiently
    only read the segment of the data needed (C++).

    This massively reduces the memory footprint, however, as the data is
    stored on disc the computation time will in crease, depending on the
    reading speed of the hard drive where the file is located.
}

\usage{
## Reading file/open file connection
read.binmm(file, skip = 0, header = TRUE, sep = ",", verbose = FALSE)

## Dimension of the data set
\method{dim}{binmm}(x)
\method{nrow}{binmm}(x)
\method{ncol}{binmm}(x)

## Dimension names (no row names, but contains all available column names)
\method{dimnames}{binmm}(x)

## Return first/last few entries of the data set
\method{head}{binmm}(x, n = 6, standardize = FALSE, \dots)
\method{tail}{binmm}(x, n = 6, standardize = FALSE, \dots)

## Default print method; shows some meta information and the first n rows of data
\method{print}{binmm}(x, n = 6, \dots)

## Very brief summary of the object; print is more informative
\method{summary}{binmm}(object, \dots)

## Subsetting; will read data on demand from file
\method{[}{binmm}(x, i, j, standardize = FALSE, drop = FALSE)
}


\arguments{
  \item{file}{the name of the file which contains the data.
      See section 'Details' for more information about the file format.}
  \item{skip}{integer: the number of lines of the data file to skip before
          beginning to read data.}
  \item{header}{a logical value indicating whether the file contains the
          names of the variables as its first line. If \code{FALSE}
          the columns will be auto-named.}
  \item{sep}{he field separator character. Values on each line of the
          file are separated by this character. Must be a single character,
          default is \code{","}.}
  \item{verbose}{logical: if set \code{TRUE} the functions will provide some additional
      information on runtime.}
  \item{n}{positive integer which specifies the number of rows to be returned.}
  \item{standardize}{logical: if set \code{TRUE} standardized values will be returned.
      See section 'Standardization' for more details.}
  \item{\dots}{currently ignored.}
  \item{i}{integer, rows to be read and returned; must be provided (can't be missing).}
  \item{j}{integer or character vector specifying which columns to read; can 
    be missing in which case all columns are returned.}
  \item{drop}{logical: if \code{TRUE} the dimensions will be dropped if possible.}
  \item{object,x}{object of class \code{"binmm"}.}
}

\author{
  Mattias Wetscher \email{Mattias.Wetscher@uibk.ac.at},
  Nikolaus Umlauf \email{Nikolaus.Umlauf@uibk.ac.at},
  Reto Stauffer \email{Reto.Stauffer@uibk.ac.at}
}

\details{
    This object has been designed for a very specific use case: reading
    a model matrix in batches to allow to process vast amounts of data
    without the need of massive amounts of memory. Thus, the flexibility
    of this type of object is limited and the files to be processed must
    have specific (although well known) format:

    \itemize{
        \item must be a tabular text file with single-character separator
        \item must only contain numeric values
        \item must not have row names
        \item can have a header line with column names
        \item can contain a series of lines at the beginning containing
            non-data information (can be skipped)
    }

    The example section contains a few minimals. If the format is violated
    the C++ code may throw an error when calling \code{\link{binmm}}
    on the first pass over the file.
}

\section{Standardization}{
    When calling \code{\link{read.binmm}} the data file is read once
    sequentially (line by line) to identify its dimension. At the same time
    row-wise meas and standard deviation is calculated (entire data set).

    This information is stored on the \code{binmm} object (see \code{$scaling}).
    When subsetting the object the data is returned in its original form by default
    (not standardized). However, if \code{standardize} is set \code{TRUE} the
    data will be returned standardized (\code{(x - mean(x)) / sd(x)}) on a row-by-row
    basis.
}

\value{
  An object of class \code{"binmm"}. Contains the following information:

  \itemize{
      \item \code{file}: Name of the file containing the data.
      \item \code{header}, \code{sep}, \code{skip}: Information regarding the format
        of the file used to process the content.
      \item \code{dim}: List of length two with number of rows and columns found
        in the file.
      \item \code{colnames}: Character vector with column names (either from the header
        or generic \code{V1}, \code{V2}, \dots.
      \item \code{scale}: List of length two with column-wise means and
        standard deviation. Used for standardization (see usage).
  }

}


\examples{
## Creating numeric 10x3 matrix to be stored
set.seed(1)
d1 <- matrix(10 + runif(30, -10, 10), ncol = 3,
             dimnames = list(NULL, LETTERS[1:3]))
head(d1, n = 3)

## Saving the data set to a temporary CSV file
write.csv(d1, file = tmpfile <- tempfile(fileext = ".csv"),
          row.names = FALSE)

## Reading the data using reto.binmm
## Will call the C++ function which makes one full pass
## trought the file.
## Returns properties of the data set as well as row-wise
## means and standard deviation. However, no data is returned
## at this stage.
rmat <- read.binmm(tmpfile)
summary(rmat)
dim(rmat)
rmat$scale

## To get the actual data subset is used. For each call
## the data will be extracted on the fly from the original
## data file.
rmat[1:3, 2:3]
rmat[1:10, ]

## Get standardized results (row-wise)
(tmp <- rmat[1:10, , standardize = TRUE])
round(apply(tmp, 2, mean), 8) # Should be 0
round(apply(tmp, 2, sd),   8) # Should be 1

## Default print method; includes a rough guess
## of the amount of memory required when the data
## set would be fully loaded.
print(rmat)


## Second example: Creating a second temporary
## file with some non-data content at the beginning,
## no header (no column names) and a different separator.
(tmpfile <- tempfile(fileext = ".dat"))

meta <- c("# -----------------------------",
          "# Meta information coming alongside",
          "# the data; needs to be skipped",
          "# when reading the data set",
          "# -----------------------------")
writeLines(meta, con = tmpfile)

## Taking the mtcars data set for example:
data("mtcars")
mtcars <- transform(mtcars,
                    vs = factor(vs, 0:1, c("V-shaped", "straight")),
                    am = factor(am, 0:1, c("automatic", "manual")))
head(mtcars, n = 2)

## Imagine being interested mpg (miles per gallon). Currently the
## data set contains factor variables, tough binmm requires
## a numeric model matrix. We'll build it and save it to the file
## we started above.
mf <- model.frame(mpg ~ ., data = mtcars)
mm <- cbind(mpg = model.response(mf),
            model.matrix(mf, data = mtcars)[, -1])
head(mm)
write.table(mm, file = tmpfile, append = TRUE,
          row.names = FALSE, col.names = FALSE, sep = ";")

## File content ...
cat(readLines(tmpfile, 8), sep = "\n")

## Reading the data set with binmm
rmat2 <- read.binmm(tmpfile, skip = 5, header = FALSE, sep = ";")
rmat2

## We have (by purpose) lost the column names

mod <- sdr(formula = list(V0    ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10,
                   sigma ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10),
           data = rmat2[seq_len(rmat2$dim$nrow), ],
           updating = "bs",
           batch_ids = 25)
plot(mod)
}

\keyword{models}
\keyword{regression}

