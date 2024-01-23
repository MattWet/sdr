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
    Class for Handling Large Model Matrices via File
}

\description{
    The stagewise algorithm \code{\link{sdr}} allows to estimate
    regression models using batches. In each iteration only a (random)
    subset of the data set is used.

    Thus, we never have to store the whole data set in memory at any
    given time. \code{binmm} allows to take an existing tabular text file (CSV)
    containing a (possibly large) model matrix and converts it into a custom
    binary file format. This reduces the file size and allows fast access
    to specific rows/columns via disc I/O.

    This massively reduces the memory footprint, however, as the data is
    stored on disc the computation time will increase depending on the
    reading speed of the hard drive where the file is located.
}

\usage{
## Reading file/open file connection
read.binmm(file, format = c("auto", "csv", "binary"),
           type = c("float", "double"), binfile = NULL, skip = 0,
           header = TRUE, sep = ",", ntest = 10, verbose = FALSE)

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
\method{[}{binmm}(x, i, j, standardize = FALSE, drop = FALSE, verbose = FALSE)
}


\arguments{
  \item{file}{the name of the file which contains the data. Either the original
      tabular text file or the binary \code{binmm} file. See section
      'Details'/'Reading binmm' for more information about file format and
      processing.}
  \item{format}{if \code{"auto"} the function tries to auto-detect if the
      \code{file} is already a binary \code{binmm} file (\code{"binary"}) or a
      tabular text file (\code{"csv"}). Depending on this many of the additional
      arguments will be ignored (see 'Reading binmm' section).}
  \item{binfile}{if \code{NULL} and \code{file} is a tabular text file,
      the name of the binary file will be generated on \code{file}.}
  \item{type}{either \code{"float"} or \code{"double"}. If \code{"float"}
      data will be stored with single precision (4 bytes) the binary file,
      else double precision (8 bytes).}
  \item{skip}{integer: the number of lines of the data file to skip before
      beginning to read data.}
  \item{header}{a logical value indicating whether the file contains the
      names of the variables as its first line. If \code{FALSE}
      the columns will be auto-named.}
  \item{sep}{he field separator character. Values on each line of the
          file are separated by this character. Must be a single character,
          default is \code{","}.}
  \item{ntest}{used if \code{file} is a tabular text file; number of rows read
      using \code{read.csv} to check the file content.}
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
  Reto Stauffer \email{Reto.Stauffer@uibk.ac.at}
}

\details{
    The class is designed for a very specific use case: reading a model
    matrix from a tabular text file (CSV-alike file) and store it in a
    custom binary file format (\code{binmm}; big model matrix).
    This allows to read the model matrix in batches to process large amounts of
    data without the need of massive amounts of memory. Thus, the flexibility
    of this type of object is limited and the files to be processed must have
    specific (although well known) format:

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

\section{Reading binmm}{
    The function \code{\link{read.binmm}} is the class constructor function
    and allows to re-use existing binary \code{binmm} files.

    In the simplest case only the \code{file} argument is set. In this case
    \code{read.binmm} tries to identify if the file is already a binary \code{binmm}
    file or not. If not, it is assumed that this is a tabular text file. This is
    the default behavior (\code{format = "auto"}) which can be overruled if necessary.

    If \code{file} points to a file which is not a binary \code{binmm} file, the
    function will check if the file can be read using \code{read.csv} given the
    options specified (\code{skip}, \code{header}, \code{sep}) by reading the first
    \code{ntest} lines. If successful, and the data returned are all numeric, the
    content of this tabular text file will be converted into the binary file format.

    By default, the binary file name (\code{binfile}) is derived from the original
    \code{file} name. If \code{file = "my_data.CSV"}, the resulting binary file
    will be named \code{"my_data.binmm"} if not specified otherwise.

    The internal C++ function will make one single pass over the tabular text file
    and stores a series of information about the content of the file:

    \itemize{
        \item \code{original_file}/\code{binfile}: Will be preserved for later.
        \item \code{dim}: Number of rows/columns.
        \item \code{colnames}: Column names (if \code{header = TRUE}).
        \item \code{scale}: Row-wise mean and standard deviation for standartization
            (if required); calculated on-line on first pass.
        \item \code{bytes}: Bytes (the precision) the data are stored in the
            binary file.
        \item \code{data} ... followed by a binary stream of numeric values (the data).
    }

    The precision can be controlled via the \code{type} argument. If \code{"float"}
    is used, the data is stored in single precision (4 bytes; limits
    \code{1.175494351e-38 - 3.402823466e+38}), if set to \code{"double"} double
    precision is used (8 bytes; {2.2250738585072014e-308 - 1.7976931348623158e+308}).
    The C++ code keeps track of the largest absolute value and would throw an error
    if the data are outside range (hopefully). Storing the data in single precision
    can save a decent amount of disc space; on the cost of losing some precision
    (see limits given above).

    All this information is stored in the \code{binfile} which is used to retrieve
    the data. Once the \code{binfile} has been created, the binary file can be used
    directly by handing over the path/name to this binary file on the \code{file}
    argument. If \code{format = "auto"} it will be detected as a binary \code{binmm}
    file automatically. In this case all other arguments to \code{read.binmm}
    will be silently ignored (except \code{verbose}).
}

\section{Standardization}{
    When calling \code{\link{read.binmm}} the data file is read once
    sequentially (line by line) to identify its dimension. At the same time
    row-wise meas and standard deviation is calculated (entire data set).

    This information is stored on the \code{binmm} object (see \code{$scaling}).
    When subsetting the object the data is returned in its original form by default
    (not standardized). If \code{standardize} is set \code{TRUE} the
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

## Read tabular text file via read.binmm; creates new binary file.
rmat <- read.binmm(tmpfile)
rmat$binfile
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

## Default print method; includes a rough guess of the amount of memory
## required when the data set would be fully loaded.
print(rmat)

## Once the binary file has been produced it can be used directly.
rmat2 <- read.binmm(rmat$binfile)
rmat2[1:3, 2:3]


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

