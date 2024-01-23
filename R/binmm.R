

read.binmm <- function(file, format = c("auto", "csv", "binary"),
                       type = c("float", "double"), binfile = NULL,
                       skip = 0, header = TRUE, sep = ",", ntest = 10, verbose = FALSE) {

    # Staying sane ...
    stopifnot("argument `file` must path/name to one file" =
              is.character(file) && length(file) == 1)
    stopifnot("file specified on argument `file` not found (does not exist)" =
              file.exists(file))
    stopifnot("argument `verbose` must be logical" = isTRUE(verbose) || isFALSE(verbose))
    stopifnot("argument `binfile` must be NULL or character of length 1" =
              is.null(binfile) || (is.character(binfile) && length(binfile) == 1L))

    # Auto-detect file type. If magic number matches it is one of our
    # binary files. If not, we assume it is a csv file which needs to be
    # processed first.
    format <- match.arg(format)
    if (format == "auto") {
        format <- if (readBin(file, "integer", 1) == 74143467L) "binary" else "csv"
        if (verbose) message("File format ", format, " (auto detect)")
    }

    # Testing if we can read the CSV file (first rows) and if what we get
    # is all numeric. If yes, we will proceed, else stop here already.
    if (format == "csv") {
        # Define binary file (if not specified by user)
        if (is.null(binfile)) {
            binfile <- paste(gsub("\\.\\w{0,3}$", "", file), "binmm", sep = ".")
        }
        type <- match.arg(type)
        
        # If that output file already exists: Throw an error!
        if (file.exists(binfile)) {
            stop("file \"", binfile, "\" already exists. Would be overwritten by ",
                 "this function. Delete the file manually and call this function again.")
        }

        # Sanity checks for arguments only needed when testing csv
        ntest <- as.integer(ntest)[1]
        stopifnot(ntest > 0)
        skip <- as.integer(skip)[1]
        stopifnot("argument `skip` must be integer >= 0" = is.integer(skip) && skip >= 0)
        stopifnot("argument `header` must be logical" = isTRUE(header) || isFALSE(header))
        stopifnot("argument `sep` must be a single character (`nchar(sep) == 1`)" =
              is.character(sep) && length(sep) == 1 && nchar(sep) == 1L)

        # Checking ...
        if (verbose) message("Checking reading CSV file (first ", ntest, " rows)")
        x <- tryCatch(read.csv(file, skip = skip, header = header, sep = sep, nrow = ntest,
                               quote = "\"", strip.white = TRUE),
                      error = function(e) stop(e),
                      warning = function(w) stop(w))
        if (!all(sapply(x, is.numeric))) {
            stop(paste("Error when checking CSV file format: Not all columns numeric.",
                       "Check arguments `skip`, `header`, `sep` and the content of the ",
                       "file (", file, ")."))
        }

        # Convert CSV to binary model matrix binmm
        unused <- create_binmm(file, binfile = binfile, type = type, skip = skip,
                               header = header, sep = sep, verbose = verbose)
    }

    # Reading meta information from new or existing binary file
    res <- meta_binmm(if (format == "binary") file else binfile)
    res$pointer_pos <- NULL # Not needed by end-user
    return(res)
}


`[.binmm` <- function (x, i, j, standardize = FALSE, drop = FALSE, verbose = FALSE) {

    i <- as.integer(i)
    if (missing(j)) j <- seq_len(x$dim$ncol)

    # If character, translate to integer
    if (is.character(j)) {
        idx <- which(!j %in% x$colnames)
        if (length(idx) > 0)
            stop("column names not in data set: ", paste(x$colnames[idx], collapse = ", "))
        j <- match(j, x$colnames)
    }

    # Out of range check
    if (any(i <= 0) | any(i > x$dim$nrow))
        stop("index `i` out of range (must be within {1, ", x$dim$nrow, "}")
    if (any(j <= 0) | any(j > x$dim$ncol))
        stop("index `j` out of range (must be within {1, ", x$dim$ncol, "}")

    # Calling cpp for getting the data; indices in cpp zero-based (-1)
    res <- subset_binmm(x$binfile,
                        sort(unique(i)) - 1, sort(unique(j)) - 1,
                        standardize = standardize, verbose = verbose)

    res <- res[match(i, sort(unique(i))), match(j, sort(unique(j))), drop = FALSE]

    if (length(j) == 1 & drop) {
        res <- as.vector(res)
    } else if (length(i) == 1 & drop) {
        res <- setNames(as.vector(res), colnames(res))
    }

    return(res)
}

# S3 methods related to dimension
dimnames.binmm <- function(x) list(NULL, x$colnames)
dim.binmm <- function(x) c(x$dim$nrow, x$dim$ncol)
nrow.binmm <- function(x) x$dim$nrow
ncol.binmm <- function(x) x$dim$nrow

# Head and tail
head.binmm <- function(x, n = 6, standardize = FALSE, ...) {
    i <- seq_len(n)
    x[i[i <= x$dim$nrow], , standardize = standardize]
}
tail.binmm <- function(x, n = 6, standardize = FALSE, ...) {
    i <- rev(x$dim$nrow - seq_len(n) + 1)
    x[i[i >= 1], , standardize = standardize]
}

# Default print
print.binmm <- function(x, n = 6, ...) {
    # Estimated size if fully loaded in MB, assuming
    # * 8 bytes for each value (double)
    # * 120 bytes for column names (char)
    # * 2 bytes for row names (int)
    # Rough guess, tough.
    mest <- ceiling((x$dim$nrow * x$dim$ncol * 8 +
                     x$dim$ncol * 100 + x$dim$nrow * 2) * 10) / 10

    cat("Warning: This is not advised.\n")
    cat("Here is some information and the head of the matrix.\n\n")

    cat("Object of class", paste(class(x), collapse = ", "), "\n")
    cat("    Original file:  ", x$original_file, "\n")
    cat("    Binary file:    ", x$binfile, "\n")
    cat("    Dimension:      ", x$dim$nrow, "x", x$dim$ncol, "\n")
    cat("    Data stored w/: ", x$bytes, "bytes (",
        if (x$bytes == 8) "double" else "float", ")\n")
    cat("    Estimated size when fully loaded:", sprintf("%.1f MB", mest / 1024^2), "\n")

    print(head(x, n = n))
    invisible(x)
}

# Brief summary
summary.binmm <- function(object, ...) {
    data.frame(Class = class(object),
               nrow = object$dim$nrow,
               ncol = object$dim$ncol,
               bytes = object$bytes,
               original_file = object$original_file,
               binary_file = object$binfile)
}


