


# Returns a data.frame with test data according to
# the example in sdr(); the size/dimension can be
# controlled via nobs and p; seed can be changed,
# defaults to 1.
# p must be >= 6 as we rely on x1-x6 to create the response y.
get_testdata <- function(nobs, p, seed = 1) {
    nobs <- as.integer(nobs)[1]
    p    <- as.integer(p)[1]
    stopifnot(is.integer(nobs) && length(nobs) == 1 & nobs > 0)
    stopifnot(is.integer(p) && length(p) == 1 & p >= 6)
    stopifnot(length(seed) == 1 && !is.null(seed))

    set.seed(seed)
    # Draw random data
    d <- matrix(runif(nobs * p, -1, 1), ncol = p,
                dimnames = list(NULL, paste0("x", seq_len(p))))
    d <- as.data.frame(d)
    
    ## Create additive predictors.
    eta.mu    <- d$x1 + 2 * d$x2 + 0.5 * d$x3 - 1*d$x4 
    eta.sigma <- 0.5 * d$x3 + 0.25 * d$x4 - 0.25 * d$x5  - 0.5 * d$x6 
    
    ## Draw normal distributed
    ## response for given predictors.
    d$y <- rnorm(nobs, mean = eta.mu, sd = exp(eta.sigma))

    return(d[, c("y", names(d)[!names(d) == "y"])])

}
