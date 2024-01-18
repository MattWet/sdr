# -------------------------------------------------------
# Checking sdr calls
# -------------------------------------------------------

if (interactive()) {
    library("tinytest")
    library("sdr")
}


# -------------------------------------------------------
# Demo data
N   <- 20
d   <- as.data.frame(matrix(runif(N * 4, -1, 1), ncol = 4, dimnames = list(NULL, paste0("x", 1:4))))
d$y <- rnorm(nrow(d), mean = 5 + d$x1 + 2 * d$x2, sd = exp(0.25 + 2 * d$x2))


# -------------------------------------------------------
# Calling sdr with wrong input arguments to test the sanity checks

# Formula; unnamed arguments to also test the order
expect_error(sdr(10, NULL, data),
             pattern = "argument `formula` must be a formula or list of formulae objects",
             info    = "Error for invalid input argument `formula`")
expect_error(sdr(list(a = 10, b = 100), NULL, data),
             pattern = "argument `formula` must be a formula or list of formulae objects",
             info    = "Error for invalid input argument `formula`")
expect_error(sdr(NA, NULL, data),
             pattern = "argument `formula` must be a formula or list of formulae objects",
             info    = "Error for invalid input argument `formula`")

# batch_ids: Must be null, numeric, or list of numerics.
expect_error(sdr("y ~ x1 + x2 + x3", data = data, batch_ids = "foo"),
             pattern = "argument `batch_ids` incorrect",
             info    = "Error for invalid input argument `batch_ids`")
expect_error(sdr("y ~ x1 + x2 + x3", data = data, batch_ids = list("foo")),
             pattern = "argument `batch_ids` incorrect",
             info    = "Error for invalid input argument `batch_ids`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, batch_ids = NA),
             pattern = "argument `batch_ids` incorrect",
             info    = "Error for invalid input argument `batch_ids`")

# updating: must be one of the allowed options
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, updating = "foo"),
             pattern = "should be one of",
             info    = "Error for invalid input argument `updating`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, updating = NA),
             pattern = "must be NULL or a character vector",
             info    = "Error for invalid input argument `updating`, testing match.arg()")

# light, CF, scalex, refitting, quick_ffdf must be logical
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, light = "foo"),
             pattern = "argument `light` must be logical",
             info    = "Error for invalid input argument `light`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, light = c(TRUE, TRUE)),
             pattern = "argument `light` must be logical",
             info    = "Error for invalid input argument `light`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, light = NULL),
             pattern = "argument `light` must be logical",
             info    = "Error for invalid input argument `light`")

expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, CF = "foo"),
             pattern = "argument `CF` must be logical",
             info    = "Error for invalid input argument `CF`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, CF = c(TRUE, TRUE)),
             pattern = "argument `CF` must be logical",
             info    = "Error for invalid input argument `CF`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, CF = NULL),
             pattern = "argument `CF` must be logical",
             info    = "Error for invalid input argument `CF`")

expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, scalex = "foo"),
             pattern = "argument `scalex` must be logical",
             info    = "Error for invalid input argument `scalex`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, scalex = c(TRUE, TRUE)),
             pattern = "argument `scalex` must be logical",
             info    = "Error for invalid input argument `scalex`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, scalex = NULL),
             pattern = "argument `scalex` must be logical",
             info    = "Error for invalid input argument `scalex`")

expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, refitting = "foo"),
             pattern = "argument `refitting` must be logical",
             info    = "Error for invalid input argument `refitting`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, refitting = c(TRUE, TRUE)),
             pattern = "argument `refitting` must be logical",
             info    = "Error for invalid input argument `refitting`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, refitting = NULL),
             pattern = "argument `refitting` must be logical",
             info    = "Error for invalid input argument `refitting`")

expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, quick_ffdf = "foo"),
             pattern = "argument `quick_ffdf` must be logical",
             info    = "Error for invalid input argument `quick_ffdf`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, quick_ffdf = c(TRUE, TRUE)),
             pattern = "argument `quick_ffdf` must be logical",
             info    = "Error for invalid input argument `quick_ffdf`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, quick_ffdf = NULL),
             pattern = "argument `quick_ffdf` must be logical",
             info    = "Error for invalid input argument `quick_ffdf`")

# eps, nu, cap (could be NULL tough), ncaps: numeric length 1, positive
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, eps = "foo"),
             pattern = "argument `eps` must be positive numeric of length 1",
             info    = "Error for invalid input argument `eps`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, eps = NA),
             pattern = "argument `eps` must be positive numeric of length 1",
             info    = "Error for invalid input argument `eps`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, eps = c(1, 2)),
             pattern = "argument `eps` must be positive numeric of length 1",
             info    = "Error for invalid input argument `eps`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, eps = -0.01),
             pattern = "argument `eps` must be positive numeric of length 1",
             info    = "Error for invalid input argument `eps`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, eps = NULL),
             pattern = "argument `eps` must be positive numeric of length 1",
             info    = "Error for invalid input argument `eps`")


expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, ncaps = "foo"),
             pattern = "argument `ncaps` must be positive numeric of length 1",
             info    = "Error for invalid input argument `ncaps`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, ncaps = NA),
             pattern = "argument `ncaps` must be positive numeric of length 1",
             info    = "Error for invalid input argument `ncaps`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, ncaps = c(1, 2)),
             pattern = "argument `ncaps` must be positive numeric of length 1",
             info    = "Error for invalid input argument `ncaps`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, ncaps = -0.01),
             pattern = "argument `ncaps` must be positive numeric of length 1",
             info    = "Error for invalid input argument `ncaps`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, ncaps = NULL),
             pattern = "argument `ncaps` must be positive numeric of length 1",
             info    = "Error for invalid input argument `ncaps`")

expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, nu = "foo"),
             pattern = "argument `nu` must be positive numeric of length 1",
             info    = "Error for invalid input argument `nu`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, nu = NA),
             pattern = "argument `nu` must be positive numeric of length 1",
             info    = "Error for invalid input argument `nu`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, nu = c(1, 2)),
             pattern = "argument `nu` must be positive numeric of length 1",
             info    = "Error for invalid input argument `nu`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, nu = -0.01),
             pattern = "argument `nu` must be positive numeric of length 1",
             info    = "Error for invalid input argument `nu`")

expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, cap = "foo"),
             pattern = "argument `cap` must be NULL or single numeric",
             info    = "Error for invalid input argument `cap`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, cap = NA),
             pattern = "argument `cap` must be NULL or single numeric",
             info    = "Error for invalid input argument `cap`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, cap = c(1, 2)),
             pattern = "argument `cap` must be NULL or single numeric",
             info    = "Error for invalid input argument `cap`")
# NULL allowed thus not testing against cap = NULL
# Can be negative (!?)

# caps: NULL or numeric vector
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, caps = "foo"),
             pattern = "argument `caps` must be NULL or numeric",
             info    = "Error for invalid input argument `caps`")
expect_error(sdr(formula = "y ~ x1 + x2 + x3", data = data, caps = NA),
             pattern = "argument `caps` must be NULL or numeric",
             info    = "Error for invalid input argument `caps`")

