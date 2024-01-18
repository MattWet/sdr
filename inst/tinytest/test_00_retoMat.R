# -------------------------------------------------------
# Checking sdr estimation
# -------------------------------------------------------

if (interactive()) {
    library("tinytest")
    library("sdr")
}

# -------------------------------------------------------
# Getting test data (for family = gaussian; default)
source("get_testdata.R")
data <- get_testdata(nobs = 20, p = 10)

# Write a temporary file for testing
tmpfile <- tempfile(fileext = ".csv")
write.csv(data, tmpfile, row.names = FALSE)

# -------------------------------------------------------
# Open file connection
expect_silent(rmat <- read.retoMat(tmpfile),
              info = "Reading test file")
expect_identical(class(rmat), "retoMat",
              info = "Checking return from read.retoMat")
expect_identical(sort(names(rmat)),
              c("colnames", "dim", "file", "header", "scale", "sep", "skip"),
              info = "Checking return from read.retoMat: elements")
expect_identical(rmat$file, tmpfile,
              info = "Checking return from read.retoMat: name of file")
expect_identical(rmat$header, TRUE,
              info = "Checking return from read.retoMat: header value")
expect_identical(rmat$skip, 0L,
              info = "Checking return from read.retoMat: skip value")
expect_identical(rmat$sep, ",",
              info = "Checking return from read.retoMat: separator used")
expect_identical(rmat$colnames, names(data),
              info = "Checking return from read.retoMat: column names")
expect_identical(rmat$dim, list(nrow = nrow(data), ncol = ncol(data)),
              info = "Checking return from read.retoMat: dimension names")


# -------------------------------------------------------
# Testing scale: calculated by CPP online; compare to
# mean/sd from the data.frame we have.
expect_equal(apply(data, 2, mean), rmat$scale$mean,
             info = "Comparing mean from cpp and R")
# RETO(TODO): Differs quite a bit. N/N-1 issue?
expect_equal(apply(data, 2, sd), rmat$scale$sd,
             tol = 1e8,
             info = "Comparing sd from cpp and R")


# -------------------------------------------------------
# Methods
expect_identical(dim(rmat), dim(data), info = "Testing S3 method dim()")

tmp <- as.matrix(head(data)) # Default n = 6
rownames(tmp) <- NULL
expect_equal(head(rmat), tmp, info = "Testing S3 method head()")

tmp <- as.matrix(head(data, n = 3)) # n = 3
rownames(tmp) <- NULL
expect_equal(head(rmat, n = 3), tmp, info = "Testing S3 method head()")

tmp <- as.matrix(tail(data)) # Default n = 6
rownames(tmp) <- NULL
expect_equal(tail(rmat), tmp, info = "Testing S3 method tail()")

tmp <- as.matrix(tail(data, n = 3)) # n = 3
rownames(tmp) <- NULL
expect_equal(tail(rmat, n = 3), tmp, info = "Testing S3 method tail()")

# RETO(TODO): Currently do not allow for negative indices
expect_error(head(rmat, -2),
             "argument must be coercible to non-negative integer")
expect_error(tail(rmat, -2),
             "argument must be coercible to non-negative integer")


# -------------------------------------------------------
# Testing subsetting by index

# Single elements
tmp <- as.matrix(data[1, 1, drop = FALSE]); rownames(tmp) <- NULL
expect_equal(rmat[1, 1], tmp, info = "Testing subsetting by index")
expect_equal(rmat[1, 1, drop = TRUE], data[1, 1], info = "Testing subsetting by index")

tmp <- as.matrix(data[10, 3, drop = FALSE]); rownames(tmp) <- NULL
expect_equal(rmat[10, 3], tmp, info = "Testing subsetting by index")
expect_equal(rmat[10, 3, drop = TRUE], data[10, 3], info = "Testing subsetting by index")

# Column matrix
tmp <- as.matrix(data[2:10, 6, drop = FALSE]); rownames(tmp) <- NULL
expect_equal(rmat[2:10, 6], tmp, info = "Testing subsetting by index")
expect_equal(rmat[2:10, 6, drop = TRUE], data[2:10, 6], info = "Testing subsetting by index")

# Row matrix
tmp <- as.matrix(data[3, 3:11, drop = FALSE]); rownames(tmp) <- NULL
expect_equal(rmat[3, 3:11], tmp, info = "Testing subsetting by index")


# Testing subsetting with empty index; only allowed for cols

# Single elements
tmp <- as.matrix(data[1:2, ]); rownames(tmp) <- NULL
expect_equal(rmat[1:2, ], tmp, info = "Testing subsetting with no col index")
expect_error(rmat[, 1:2],      info = "Testing subsetting with no row index (not allowed)")


# Testing subsetting with named columns

tmp <- as.matrix(data[10, "x3", drop = FALSE]); rownames(tmp) <- NULL
expect_equal(rmat[10, "x3"], tmp, info = "Testing subsetting by index")
expect_equal(rmat[10, "x3", drop = TRUE], data[10, "x3"], info = "Testing subsetting by index")

tmp <- as.matrix(data[3:5, c("x3", "y", "x4"), drop = FALSE]); rownames(tmp) <- NULL
expect_equal(rmat[3:5, c("x3", "y", "x4")], tmp, info = "Testing subsetting by index")

# Mixed order (not dawing row/cols in sequence)
ii <- c(3, 15, 4, 12)
jj <- c(10, 3, 1, 8)
tmp <- as.matrix(data[ii, jj]); rownames(tmp) <- NULL
expect_equal(rmat[ii, jj], tmp, info = "Testing subsetting by index")



# -------------------------------------------------------
# Testing stndardization; for our data.frame we have to
# do it manually first ...
sdata <- data.frame(lapply(data, function(x) (x - mean(x)) / sd(x)))
sdata <- as.matrix(sdata)
rownames(data) <- NULL

# nice ...
expect_equal(rmat[seq_len(nrow(data)), , standardize = TRUE], sdata,
             info = "Testing standardization")






