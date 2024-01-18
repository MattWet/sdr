

read.retoMat <- function(file, skip = 0, header = TRUE, sep = ",", verbose = FALSE) {

    # Staying sane ...
    stopifnot("argument `file` must path/name to one file" =
              is.character(file) && length(file) == 1)
    stopifnot("file specified on argument `file` not found (does not exist)" =
              file.exists(file))

    skip <- as.integer(skip)[1]
    stopifnot("argument `skip` must be integer >= 0" = is.integer(skip) && skip >= 0)
    stopifnot("argument `header` must be logical" = isTRUE(header) || isFALSE(header))
    stopifnot("argument `sep` must be a single character (`nchar(sep) == 1`)" =
              is.character(sep) && length(sep) == 1 && nchar(sep) == 1L)
    stopifnot("argument `verbose` must be logical" = isTRUE(verbose) || isFALSE(verbose))

    # Calling cpp method, return result
    return(retoMat(file, skip = skip, header = header, sep = sep, verbose = verbose))
}


`[.retoMat` <- function (x, i, j, standardize = FALSE, drop = FALSE) {

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
    if (any(i < 0) | any(i > x$dim$nrow))
        stop("index `i` out of range (must be within {1, ", x$dim$nrow, "}")
    if (any(j < 0) | any(j > x$dim$ncol))
        stop("index `i` out of range (must be within {1, ", x$dim$ncol, "}")

    # Calling cpp for getting the data
    res <- retoMat_subset(x, sort(unique(i)), sort(unique(j)),
                          standardize = standardize, sep = x$sep)

    res <- res[match(i, sort(unique(i))), match(j, sort(unique(j))), drop = FALSE]

    if (length(j) == 1 & drop) {
        res <- as.vector(res)
    } else if (length(i) == 1 & drop) {
        res <- setNames(as.vector(res), colnames(res))
    }
    return(res)

}

dim.retoMat <- function(x) c(x$dim$nrow, x$dim$ncol)

head.retoMat <- function(x, n = 6, standardize = FALSE, ...) {
    i <- seq_len(n)
    x[i[i <= x$dim$nrow], , standardize = standardize]
}

tail.retoMat <- function(x, n = 6, standardize = FALSE, ...) {
    i <- rev(x$dim$nrow - seq_len(n) + 1)
    x[i[i >= 1], , standardize = standardize]
}

print.retoMat <- function(x, n = 6, ...) {
    # Estimated size if fully loaded in MB, assuming
    # * 8 bytes for each value (double)
    # * 120 bytes for column names (char)
    # * 2 bytes for row names (int)
    # Rough guess, tough.
    mest <- ceiling((x$dim$nrow * x$dim$ncol * 8 + x$dim$ncol * 100 + x$dim$nrow * 2) * 10) / 10

    cat("Warning: This is not advised.\n")
    cat("Here is some information and the head of the matrix.\n\n")

    cat("Object of class", paste(class(x), collapse = ", "), "\n")
    cat("    Dimension:", x$dim$nrow, "x", x$dim$ncol, "\n")
    cat("    Source file:", sprintf("\"%s\"", x$file), "\n")
    cat("    Estimated size when fully loaded:", sprintf("%.1f MB", mest / 1024^2), "\n")

    if (!file.exists(x$file)) {
        cat("\n[ERROR] File cannot be found\n")
    } else {
        print(head(x, n = n))
    }
    invisible(x)
}

summary.retoMat <- function(object, ...) {
    data.frame(Class = class(object),
               nrow = object$dim$nrow,
               ncol = object$dim$ncol,
               file = object$file)
}


