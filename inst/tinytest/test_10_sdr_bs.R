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
data <- get_testdata(nobs = 1000, p = 10)

# Building formula (for mu and sigma)
varnames <- names(data)[!names(data) == "y"]
f <- formula(paste("y ~", paste(varnames, collapse = "+")))
f <- list(f, update(f, "sigma ~ ."))

# -------------------------------------------------------
# Running estimation (batch size 100, maxit = 10)
expect_silent(mod <- sdr(f, data = data, updating = "bs",
                         batch_ids = 100, maxit = 10),
              info = "Running sdr with updating = 'bs' on test data")

# Testing class of returned object
expect_inherits(mod, "stagewise",   info = "Testing return class")
expect_true(is.list(mod),           info = "Testing return class")

# Testing content of the object (structure)
expect_identical(sort(names(mod)),
                 c("cap", "coefficients", "family", "formula", "logLik", "maxit", "nobs", "varnames", "X", "y"),
                 info = "Testing if the named list contains all expected elements")

# Testing result (content of returned object)
expect_identical(mod$maxit, list(var_selection = 10, refitting = 0),
                 info = "Testing content of $maxit")
expect_equal(mod$nobs, 100L,
                 info = "Testing content of $nobs; identical to batch size")
expect_identical(mod$cap, 0.175,
                 info = "Testing content of $cap; default upper limit set by sdr()")
expect_identical(mod$varnames, list(mu = varnames, sigma = varnames),
                 info = "Testing content of $varnames")
expect_inherits(mod$family, "family.bamlss",
                 info = "Testing $family (family.bamlss object; default)")
expect_identical(mod$family$family, "NO",
                 info = "Testing $family (NO; default)")

# Formula objects are not identical tough we can test if they are similar
expect_identical(gsub("\\s+", "", paste("y", format(mod$formula[[1]]))),
                 gsub("\\s+", "", format(f[[1]])),
                 info = "Testing $formula (mu), check if it matches our original formula")
expect_identical(gsub("\\s+", "", paste("sigma", format(mod$formula[[2]]))),
                 gsub("\\s+", "", format(f[[2]])),
                 info = "Testing $formula (sigma), check if it matches our original formula")

# X/y are returned as our input `data` was a simple data.frame; light = FALSE.
tmp <- structure(as.matrix(data[, varnames]), dimnames = list(seq_len(nrow(data)), varnames))
expect_identical(mod$X, tmp, info = "Testing $X returned by sdr()")
tmp <- structure(as.matrix(data[, "y", drop = FALSE]), dimnames = list(seq_len(nrow(data)), "y"))
expect_identical(mod$y, tmp, info = "Testing $y returned by sdr()")

# As we have a fixed setup and seeded data we can check the
# resulting logLik paths and the coefficients.
expect_true(is.numeric(mod$logLik) & length(mod$logLik) == 10, info = "Checking $logLik object")
expect_equal(mod$logLik,
             c(-203.2020, -199.7162, -203.7082, -206.2294, -203.4647,
               -200.5817, -201.0097, -198.7533, -196.8606, -192.2951),
             tol = 1e-3,
             info = "Testing content of $logLik")


# Coefficients; only testing structure and last iteration
expect_identical(names(mod$coefficients), c("mu", "sigma"),   info = "Testing object $coefficients")
expect_true(all(sapply(mod$coefficients, is.matrix)),         info = "Testing object $coefficients")
expect_true(all(sapply(mod$coefficients, nrow) == 10),
            info = "Testing dimension of $coefficients matrices")
expect_true(all(sapply(mod$coefficients, ncol) == 1 + length(varnames)),
            info = "Testing dimension of $coefficients matrices")
expect_identical(colnames(mod$coefficients$mu),   c("(Intercept)", varnames),
            info = "Testing colnames of $coefficients$mu matrices")
expect_identical(colnames(mod$coefficients$sigma), c("(Intercept)", varnames),
            info = "Testing colnames of $coefficients$sigma matrices")
# Testing last row of coefficients
expect_equal(as.vector(tail(mod$coefficients$mu, 1)),
             c(-0.1735, 0, 0.1403, 0, 0, 0, 0, 0, 0, 0, 0),
             tol = 1e-3,
             info = "Testing final coefficients (mu)")
expect_equal(as.vector(tail(mod$coefficients$sigma, 1)),
             c(0.6046, 0, 0, 0, 0.0249, 0, 0, 0, 0, 0, 0),
             tol = 1e-3,
             info = "Testing final coefficients (sigma)")
# Same but using the coef function
expect_equal(as.vector(coef(mod)$mu),
             c(-0.1735, 0, 0.1403, 0, 0, 0, 0, 0, 0, 0, 0),
             tol = 1e-3,
             info = "Testing final coefficients (mu)")
expect_equal(as.vector(coef(mod)$sigma),
             c(0.6046, 0, 0, 0, 0.0249, 0, 0, 0, 0, 0, 0),
             tol = 1e-3,
             info = "Testing final coefficients (sigma)")

# -------------------------------------------------------
# Testing some methods
# -------------------------------------------------------
# For this run also test AIC, BIC, ...
expect_inherits(AIC(mod), "stagewise_AIC", info = "Testing return of AIC()")
expect_equal(unclass(AIC(mod)),
             c(410.4040, 405.4324, 413.4163, 418.4587, 412.9293,
               409.1634, 410.0194, 405.5066, 401.7213, 392.5901),
             tol = 1e-3,
             info = "Testing values of AIC()")

expect_inherits(BIC(mod), "stagewise_BIC", info = "Testing return of BIC()")
expect_equal(unclass(BIC(mod)),
             c(415.6144, 413.2479, 421.2318, 426.2742, 420.7448,
               419.5841, 420.4401, 415.9273, 412.1420, 403.0108),
             tol = 1e-3,
             info = "Testing values of BIC()")

# Testing logLik method
expect_identical(class(logLik(mod)), c("logLik", "logLik.stagewise"),
                 info = "Testing return of logLik()")
expect_true("df" %in% names(attributes(logLik(mod))), info = "Testing df attribute on logLik() return")
expect_silent(df <- attr(logLik(mod), "df"),          info = "Testing df attribute on logLik() return")
expect_identical(df, setNames(as.numeric(rep(2:4, c(1, 4, 5))), 1:10),
                 info = "Testing df attribute on logLik() return")


# Predict method (works as we kept X in this case)
expect_silent(pred <- predict(mod),                             info = "Testing predict() method")
expect_identical(names(pred), c("mu", "sigma"),                 info = "Testing predict() method")
expect_true(all(sapply(pred, length) == 1000L),                 info = "Testing predict() method")
expect_true(all(sapply(pred, is.numeric)),                      info = "Testing predict() method")
expect_true(all(sapply(pred, function(x) sum(is.na(x)) == 0)),  info = "Testing predict() method")
# ... quick check - just calculating sum(abs(...)) to see if we got the correct result;
# not testing invididual predictions.
expect_equal(as.vector(sapply(pred, function(x) sum(abs(x)))),
             c(176.2504, 604.5879), tol = 1e-3,
             info = "Checking predictions (sum of absolute values of all predictions as a quick check)")
rm(pred)

# Residuals method
expect_silent(r <- residuals(mod),                             info = "Testing residuals() method")
expect_identical(attr(r, "type"), "Quantile",                  info = "Testing default return type of residuals()")
expect_identical(class(r), c("stagewise_residuals", "matrix", "array"),
                 info = "Testing residuals() return class")
expect_identical(dim(r), c(1000L, 1L),                         info = "Testing dimension of result of residuals()")
expect_identical(colnames(r), "y",                             info = "Testing dimnames of result of residuals()")
expect_identical(rownames(r), as.character(1:1000),            info = "Testing dimnames of result of residuals()")
expect_equal(sum(r^2), 949.1934, tol = 1e-3,                   info = "Checking residual sum of squares on residuals()")









# -------------------------------------------------------
# Running estimation again with light = TRUE; should
# result in no longer having X in the return (model matrix)
# -------------------------------------------------------
rm(mod)
expect_silent(mod <- sdr(f, data = data, updating = "bs",
                         batch_ids = 100, maxit = 10, light = TRUE),
              info = "Running sdr with updating = 'bs' and light = TRUE on test data")
expect_true(!"X" %in% names(mod), info = "Testing if $X is no longer present with light = TRUE")
expect_true("y" %in% names(mod),  info = "Testing that we still have $y when light = TRUE")




# -------------------------------------------------------
# Running estimation again with scalex = TRUE which should
# result in different results (scaling all covariates);
# Thus checking final coefficients and likelihood path once again
# -------------------------------------------------------
rm(mod)
expect_silent(mod <- sdr(f, data = data, updating = "bs",
                         batch_ids = 100, maxit = 10, scalex = TRUE),
              info = "Running sdr with updating = 'bs' and scalex = TRUE on test data")

expect_equal(mod$logLik,
             c(-201.8319, -203.6381, -191.1261, -213.4330, -193.6972,
               -200.4121, -202.1657, -200.9464, -209.3092, -214.2546),
             tol = 1e-3,
             info = "Testing content of $logLik")
expect_equal(as.vector(coef(mod)$mu),
             c(-0.02887, 0, 0.1149, 0, 0, 0, 0, 0, 0, 0, 0),
             tol = 1e-3,
             info = "Testing final coefficients (mu)")
expect_equal(as.vector(coef(mod)$sigma),
             c(0.6335, 0, -0.0111, 0, 0.02779, 0, -0.03116, 0, 0, 0, 0),
             tol = 1e-3,
             info = "Testing final coefficients (sigma)")
