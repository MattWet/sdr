# Libraries.
library("bamlss")
library("gamlss.dist")
library("stagewise")

setwd("data")

# model save directory
dir_save <- "bzanbi_train.rda" 
dir_save2 <- "bzanbi_refit.rda" 

# path data
dir_dat <- "training.csv" 
dt <- read.csv(dir_dat)
f <- setdiff(colnames(dt), "counts")

nvars <- length(f)
f <- paste0(f , collapse = "+")
f <- paste("counts~ ", f)
f <- as.formula(f)
f <- list(f, update(f, sigma ~ .), update(f, nu ~ .))
f

fam <- ZANBI
nobs_train <- nrow(dt) #906984
batch.size <- 10000
maxit <- 2000
maxit_refit <- 30000
ind <- 1:nrow(dt)
pos <- dt[, "counts"] > 0
ind0 <- ind[!pos]
ind1 <- ind[pos]

# Set upsampling proportion of nonzero lightning events
prop <- 0.5 
s1 <- round(prop*batch.size)
s0 <- round((1-prop)*batch.size)

# Selection model
set.seed(123456)
batch_ids <- lapply(1:(maxit), function(...) sample(c(sample(ind0, size = s0, replace = F),sample(ind1, size = s1, replace = F))))
ptm <- proc.time()
set.seed(123456)
b.zanbi1 <- sdr(formula = f, data = dt, family = fam, maxit = maxit, plot = F, eps = 0.1, nu = 0.01,
                batch_ids = batch_ids, scalex = FALSE, updating = "bs", quick_ffdf = TRUE,
                light = TRUE, CF = TRUE)

time.select <- c(proc.time() - ptm)[3]
time.select
gc()
saveRDS(b.zanbi1, file = dir_save)

bic <- AIC(b.zanbi1, K = log(batch.size))
itermax <- max(20,which.min(bic))
itermax

# Determine selected variables and new model formula consisting of selected variables only
f <- newformula(b.zanbi1, mstop = itermax, name = "counts")
c <- coef(b.zanbi1, mstop = itermax, refit = TRUE)

# Refitting for 30k iterations
set.seed(123456)
batch_ids <- lapply(1:(maxit_refit), function(...) sample(c(sample(ind0, size = s0, replace = F),sample(ind1, size = s1, replace = F))))
ptm <- proc.time()
set.seed(123456)
b.zanbi <- sdr(formula = f, data = dt, family = fam, maxit = maxit_refit, plot = FALSE, eps = 0.1, nu = 0.01,
               batch_ids = batch_ids, scalex = FALSE, updating = "bs", quick_ffdf = TRUE, CF = FALSE,
               light = TRUE, init = FALSE, coef_start = c)
time.refit <- c(proc.time() - ptm)[3]
time.refit

saveRDS(b.zanbi, file = dir_save2)


# update nu intercept to correctly represent lightning proportions
prop <- length(ind1)/nrow(dt) # = 0.0265
nu =  -log(prop/(1-prop))
nu1 = b.zanbi[["coefficients"]][["nu"]][,"(Intercept)"]
b.zanbi[["coefficients"]][["nu"]][,"(Intercept)"] <- nu1 + rep(nu, length(b.zanbi[["coefficients"]][["nu"]][,"(Intercept)"]))

saveRDS(b.zanbi, file = dir_save2)

gc()

time.select
time.refit
itermax