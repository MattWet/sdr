library("bamlss")
library("gamlss.dist")
library("gamboostLSS")
library("ASL")
library("ff")
library("stagewise")
source("~/SCRATCH/stagewise/prg_R/Var_deselection.R")

### Define data generating process
## Parameter for linear predictors in dgp
cm1 = 0.5; cm3 = -1; cm5 = 0.75; cm6 = 0.75
cs2 = 1; cs4 = -1.25; cs5 = 1
cn3 = 1; cn4 = -1; cn5 = -1

muc <- c(cm1,cm3, cm5, cm6)
sigmac <- c(cs2,cs4, cs5)
nuc <- c(cn3,cn4, cn5)


dgp <- function(nobs = 1000, nnoise = 10, rho = 0, seed = rnorm(1)) {
  
  p <- 6 + nnoise
  d <- matrix(runif(nobs * p, -1, 1), ncol = p)
  
  # Desired correlation in data
  S <- rho^as.matrix(dist(1:p))
  
  # Multiplication with cholesky matrix introduces correlation to data
  d <- matrix(d, ncol = p, nrow = nobs) %*% chol(S)
  set.seed(seed)
  d <- d[,sample(colnames(d))]
  
  colnames(d) <- paste("x", 1:p, sep = "")
  d <- as.data.frame(d)
  
  # Linear predictors and data generation
  d$eta.mu <- 0.5 + cm1 * d$x1 + cm3 * d$x3 + cm5 * d$x5 + cm6 * d$x6
  d$eta.sigma <--1 + cs2 * d$x2 + cs4 * d$x4 + cs5 * d$x5
  d$eta.nu <- -0.5 + cn3 * d$x3 + cn4 * d$x4 + cn5 * d$x5
  
  # transformation from linear predictor scale to dist parameter scale
  nu <- exp(d$eta.nu)
  nu <- nu/(1+nu)
  mu <- exp(d$eta.mu)
  sigma <- exp(d$eta.sigma)
  d$y <- rZANBI(nobs, mu = mu, sigma = sigma, nu = nu)
  
  return(d)
}

### Simulation function
sim <- function(nobs = 500, nnoise = 100, rho = 0){
  
  ## Generate training data
  seed = SGE_TASK_ID 
  d <- dgp(nobs = nobs, nnoise = nnoise, rho = rho/10, seed = seed)
  
  ## Generate validation data
  set.seed(123)
  nd <- dgp(nobs = 10000, nnoise = nnoise, rho = rho/10, seed = seed)
  
  # MCMC
  if(TRUE){
    maxit = NA
    
    fam <- ZANBI()
    
    ## (1) True model
    ftrue <- list(y ~ x1 + x3 + x5 + x6,
                  sigma ~ x2 + x4 + x5,
                  nu ~ x3 + x4 + x5)
    
    ptm <- proc.time()
    b <- bamlss(ftrue, data = d, family = fam, sampler = FALSE)
    time <- c(proc.time() - ptm)[3]
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    nu <- mean((nd$eta.nu - p$nu)^2)
    
    # transformation from linear predictor scale to dist parameter scale
    par.trafo <- function(list){
      list$mu <- exp(list$mu)
      list$sigma <- exp(list$sigma)
      list$nu <- 1/(1+exp(-list$nu))
      return(list)
    }
    p <- par.trafo(p)
    
    k = max(350, nd$y)
    mat <- fam$d(0, p)
    for(i in 1:k){
      pp <- fam$d(i, p)
      mat <- rbind(mat, pp)
      cat("\r");cat("i =",i, "/", k)
    }
    mat <- t(mat)
    colnames(mat) <- 0:k
    mean1 <- rowSums(mat %*% diag(c(0:k)))
    
    mse1 <- mean((nd$y - mean1)^2)
    crps1 <- try(rps(nd$y, mat),TRUE)
    if(!is.numeric(crps1)) crps1 <- NaN
    
    ll <- fam[["d"]]
    l <- try(sum(ll(nd$y, p, log = T)), T)
    if(!is.numeric(l)) l <- NaN
    
    # combine results
    m <- rbind(
      "likelihood" = c("mu" = mu, "sigma" = sigma, "nu" = nu, "mse" = mse1, "crps" = crps1, "logLiK" = l, 
                       "mstop" = maxit, time)
    )
    
    m <- as.data.frame(m)
    rownames(m) <- NULL
    
    m$method <- "likelihood"
    m$nobs <- nobs
    m$nnoise <- nnoise
    m$rho <- rho
    m$rep <- seed
    
    
    # Extract coefs and calc number of true and false positives
    c <- data.frame(t(coef(b)))
    
    nn <- colnames(c)
    nmu <- nn[grepl("mu",nn)]
    nsigma <- nn[grepl("sigma",nn)]
    nnu <- nn[grepl("nu",nn)]
    c.mu <- c[nmu]
    c.sigma <- c[nsigma]
    c.nu <- c[nnu]
    nmu <- gsub("mu.p.","",nmu)
    nsigma <- gsub("sigma.p.","",nsigma)
    nnu <- gsub("nu.p.","",nnu)
    colnames(c.mu) <- nmu
    colnames(c.sigma) <- nsigma
    colnames(c.nu) <- nnu
    colnames(c.mu)[1] <- "(Intercept)"
    colnames(c.sigma)[1] <- "(Intercept)"
    colnames(c.nu)[1] <- "(Intercept)"
    c <- list("mu" = unlist(c.mu), "sigma" = unlist(c.sigma), "nu" = unlist(c.nu))
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      c.nu <- c[["nu"]]
      cc.mu <- cc.sigma <- cc.nu <-  rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      cc.nu[1] <- c.nu["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
        cc.nu[i+1] <- ifelse(paste0("x",i) %in% names(c.nu) ,c.nu[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma, "nu" = cc.nu))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    nu.names <- c("nu.intercept", paste0("nu.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    names(c[["nu"]]) <- nu.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    coefs.nu <- rbind.data.frame(c[["nu"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- rownames(coefs.nu) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    colnames(coefs.nu) <- nu.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,c(2,4,6,7)] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,c(3,5,6)] != 0)
    resi$tp.nu <- rowSums(coefs.nu[,c(4,5,6)] != 0)
    
    resi$fp.mu <- rowSums(coefs.mu[,-c(1,2,4,6,7)] != 0)
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,3,5,6)] != 0)
    resi$fp.nu <- rowSums(coefs.nu[,-c(1,4,5,6)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2,4,6,7)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(3,5,6)] -  t(matrix(rep(sigmac, 1), nrow = 3)) ) ^2)
    resi$tp.nu.mse <- rowMeans((coefs.nu[,c(4,5,6)] - t(matrix(rep(nuc, 1), nrow = 3)) ) ^2)
    
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1,2,4,6,7)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,3,5,6)]^2)
    resi$fp.nu.mse <- rowMeans(coefs.nu[,-c(1,4,5,6)]^2)
    
    m <- cbind(m, coefs.mu[,c(2,4,6,7)], coefs.sigma[,c(3,5,6)],coefs.nu[,c(4,5,6)], resi)
    m1 <- m
  }
  
  # Thresdesc full batch
  if(nobs < 100000){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .),
      update(f, nu ~ .)
    )
    
    fam <- ZANBI()
    
    # Threshold descent
    ptm <- proc.time()
    maxit = 200
    maxit_refit = 400
    batch_ids = nobs
    b <- sdr(formula = f, family = fam, data = d, maxit = maxit, maxit_refit = maxit_refit, 
             batch_ids = batch_ids, K = log(batch_ids), updating = "thresdesc")
    time <- c(proc.time() - ptm)[3]
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    nu <- mean((nd$eta.nu - p$nu)^2)
    
    # transformation from linear predictor scale to dist parameter scale
    par.trafo <- function(list){
      list$mu <- exp(list$mu)
      list$sigma <- exp(list$sigma)
      list$nu <- 1/(1+exp(-list$nu))
      return(list)
    }
    p <- par.trafo(p)
    
    k = max(350, nd$y)
    mat <- fam$d(0, p)
    for(i in 1:k){
      pp <- fam$d(i, p)
      mat <- rbind(mat, pp)
      cat("\r");cat("i =",i, "/", k)
    }
    mat <- t(mat)
    colnames(mat) <- 0:k
    mean1 <- rowSums(mat %*% diag(c(0:k)))
    
    mse1 <- mean((nd$y - mean1)^2)
    crps1 <- try(rps(nd$y, mat),TRUE)
    if(!is.numeric(crps1)) crps1 <- NaN
    
    ll <- fam[["d"]]
    l <- try(sum(ll(nd$y, p, log = T)), T)
    if(!is.numeric(l)) l <- NaN
    
    
    # combine results
    m <- rbind(
      "thresdesc" = c("mu" = mu, "sigma" = sigma, "nu" = nu, "mse" = mse1, "crps" = crps1, "logLiK" = l, "mstop" = maxit_refit, time)
    )
    
    m <- as.data.frame(m)
    rownames(m) <- NULL
    
    m$method <- "thresdesc"
    m$nobs <- nobs
    m$nnoise <- nnoise
    m$rho <- rho
    m$rep <- seed
    
    
    # Extract coefs and calc number of true and false positives
    c <- coef(b)
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      c.nu <- c[["nu"]]
      cc.mu <- cc.sigma <- cc.nu <-  rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      cc.nu[1] <- c.nu["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
        cc.nu[i+1] <- ifelse(paste0("x",i) %in% names(c.nu) ,c.nu[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma, "nu" = cc.nu))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    nu.names <- c("nu.intercept", paste0("nu.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    names(c[["nu"]]) <- nu.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    coefs.nu <- rbind.data.frame(c[["nu"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- rownames(coefs.nu) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    colnames(coefs.nu) <- nu.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,c(2,4,6,7)] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,c(3,5,6)] != 0)
    resi$tp.nu <- rowSums(coefs.nu[,c(4,5,6)] != 0)
    
    resi$fp.mu <- rowSums(coefs.mu[,-c(1,2,4,6,7)] != 0)
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,3,5,6)] != 0)
    resi$fp.nu <- rowSums(coefs.nu[,-c(1,4,5,6)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2,4,6,7)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(3,5,6)] -  t(matrix(rep(sigmac, 1), nrow = 3)) ) ^2)
    resi$tp.nu.mse <- rowMeans((coefs.nu[,c(4,5,6)] - t(matrix(rep(nuc, 1), nrow = 3)) ) ^2)
    
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1,2,4,6,7)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,3,5,6)]^2)
    resi$fp.nu.mse <- rowMeans(coefs.nu[,-c(1,4,5,6)]^2)
    
    m <- cbind(m, coefs.mu[,c(2,4,6,7)], coefs.sigma[,c(3,5,6)],coefs.nu[,c(4,5,6)], resi)
    m1 <- rbind(m1, m)
  }
  
  # Thresdesc batchwise
  if(TRUE){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .),
      update(f, nu ~ .)
    )
    
    fam <- ZANBI()
    
    if(nobs > 10000) d <- as.ffdf(d)
    # Threshold descent bw
    ptm <- proc.time()
    maxit = 200
    maxit_refit = 400
    batch_ids <- min(round(2*nobs/3), 10000)
    b <- sdr(formula = f, family = fam, data = d, maxit = maxit, maxit_refit = maxit_refit, 
             batch_ids = batch_ids, K = log(batch_ids), updating = "thresdesc")
    time <- c(proc.time() - ptm)[3]
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    nu <- mean((nd$eta.nu - p$nu)^2)
    
    # transformation from linear predictor scale to dist parameter scale
    par.trafo <- function(list){
      list$mu <- exp(list$mu)
      list$sigma <- exp(list$sigma)
      list$nu <- 1/(1+exp(-list$nu))
      return(list)
    }
    p <- par.trafo(p)
    
    k = max(350, nd$y)
    mat <- fam$d(0, p)
    for(i in 1:k){
      pp <- fam$d(i, p)
      mat <- rbind(mat, pp)
      cat("\r");cat("i =",i, "/", k)
    }
    mat <- t(mat)
    colnames(mat) <- 0:k
    mean1 <- rowSums(mat %*% diag(c(0:k)))
    
    mse1 <- mean((nd$y - mean1)^2)
    crps1 <- try(rps(nd$y, mat),TRUE)
    if(!is.numeric(crps1)) crps1 <- NaN
    
    ll <- fam[["d"]]
    l <- try(sum(ll(nd$y, p, log = T)), T)
    if(!is.numeric(l)) l <- NaN
    
    
    # combine results
    m <- rbind(
      "thresdesc_bw" = c("mu" = mu, "sigma" = sigma, "nu" = nu, "mse" = mse1, "crps" = crps1, "logLiK" = l, "mstop" = maxit_refit, time)
    )
    
    m <- as.data.frame(m)
    rownames(m) <- NULL
    
    m$method <- "thresdesc_bw"
    m$nobs <- nobs
    m$nnoise <- nnoise
    m$rho <- rho
    m$rep <- seed
    
    
    # Extract coefs and calc number of true and false positives
    c <- coef(b)
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      c.nu <- c[["nu"]]
      cc.mu <- cc.sigma <- cc.nu <-  rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      cc.nu[1] <- c.nu["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
        cc.nu[i+1] <- ifelse(paste0("x",i) %in% names(c.nu) ,c.nu[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma, "nu" = cc.nu))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    nu.names <- c("nu.intercept", paste0("nu.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    names(c[["nu"]]) <- nu.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    coefs.nu <- rbind.data.frame(c[["nu"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- rownames(coefs.nu) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    colnames(coefs.nu) <- nu.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,c(2,4,6,7)] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,c(3,5,6)] != 0)
    resi$tp.nu <- rowSums(coefs.nu[,c(4,5,6)] != 0)
    
    resi$fp.mu <- rowSums(coefs.mu[,-c(1,2,4,6,7)] != 0)
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,3,5,6)] != 0)
    resi$fp.nu <- rowSums(coefs.nu[,-c(1,4,5,6)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2,4,6,7)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(3,5,6)] -  t(matrix(rep(sigmac, 1), nrow = 3)) ) ^2)
    resi$tp.nu.mse <- rowMeans((coefs.nu[,c(4,5,6)] - t(matrix(rep(nuc, 1), nrow = 3)) ) ^2)
    
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1,2,4,6,7)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,3,5,6)]^2)
    resi$fp.nu.mse <- rowMeans(coefs.nu[,-c(1,4,5,6)]^2)
    
    m <- cbind(m, coefs.mu[,c(2,4,6,7)], coefs.sigma[,c(3,5,6)],coefs.nu[,c(4,5,6)], resi)
    m1 <- rbind(m1,m)
  }
  
  # Gradient boosting with 10 fold cross validation (gamboostLSS)
  if(nobs < 100000){
    f <- grep("x", names(d), value = TRUE)
    f <- paste("y ~", paste(f, collapse = "+"))
    f <- as.formula(f)
    
    fam <- ZANBI()
    
    ptm <- proc.time()
    maxitlss <- ifelse(nobs >= 5000, 12000/2.5, 8000/2.5)
    b <- glmboostLSS(f,  families = as.families(ZANBI),
                     data = d, method = "noncyclic",
                     control = boost_control(mstop = maxitlss, nu = 0.25, trace = TRUE))
    set.seed(seed)
    folds <- mboost::cv(rep(1, nrow(d)), B = 10, type = "kfold")
    cvr <- cvrisk(b, folds = folds, papply = lapply)
    mstop(b) <- mstop(cvr)
    
    mstop <- mstop(b)[1]
    time <- c(proc.time() - ptm)[3]
    
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    nu <- mean((nd$eta.nu - p$nu)^2)
    
    # transformation from linear predictor scale to dist parameter scale
    par.trafo <- function(list){
      list$mu <- exp(list$mu)
      list$sigma <- exp(list$sigma)
      list$nu <- 1/(1+exp(-list$nu))
      return(list)
    }
    p <- par.trafo(p)
    
    k = max(350, nd$y)
    mat <- fam$d(0, p)
    for(i in 1:k){
      pp <- fam$d(i, p)
      mat <- rbind(mat, pp)
      cat("\r");cat("i =",i, "/", k)
    }
    mat <- t(mat)
    colnames(mat) <- 0:k
    mean1 <- rowSums(mat %*% diag(c(0:k)))
    
    mse1 <- mean((nd$y - mean1)^2)
    crps1 <- try(rps(nd$y, mat),TRUE)
    if(!is.numeric(crps1)) crps1 <- NaN
    
    ll <- fam[["d"]]
    l <- try(sum(ll(nd$y, p, log = T)), T)
    if(!is.numeric(l)) l <- NaN
    
    # combine results
    m <- rbind(
      "boosting" = c("mu" = mu, "sigma" = sigma, "nu" = nu, "mse" = mse1, "crps" = crps1, "logLiK" = l, "mstop" = mstop, time)
    )
    
    m <- as.data.frame(m)
    rownames(m) <- NULL
    
    m$method <- "boosting"
    m$nobs <- nobs
    m$nnoise <- nnoise
    m$rho <- rho
    m$rep <- seed
    
    
    # Extract coefs and calc number of true and false positives
    c <- data.frame(t(coef(b)))
    c <- coef(b, off2int = TRUE)
    for(i in names(c)) if(is.null(c[[i]])) c[[i]] <- c("(Intercept)" = 0)
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      c.nu <- c[["nu"]]
      cc.mu <- cc.sigma <- cc.nu <-  rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      cc.nu[1] <- c.nu["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
        cc.nu[i+1] <- ifelse(paste0("x",i) %in% names(c.nu) ,c.nu[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma, "nu" = cc.nu))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    nu.names <- c("nu.intercept", paste0("nu.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    names(c[["nu"]]) <- nu.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    coefs.nu <- rbind.data.frame(c[["nu"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- rownames(coefs.nu) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    colnames(coefs.nu) <- nu.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,c(2,4,6,7)] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,c(3,5,6)] != 0)
    resi$tp.nu <- rowSums(coefs.nu[,c(4,5,6)] != 0)
    
    resi$fp.mu <- rowSums(coefs.mu[,-c(1,2,4,6,7)] != 0)
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,3,5,6)] != 0)
    resi$fp.nu <- rowSums(coefs.nu[,-c(1,4,5,6)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2,4,6,7)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(3,5,6)] -  t(matrix(rep(sigmac, 1), nrow = 3)) ) ^2)
    resi$tp.nu.mse <- rowMeans((coefs.nu[,c(4,5,6)] - t(matrix(rep(nuc, 1), nrow = 3)) ) ^2)
    
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1,2,4,6,7)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,3,5,6)]^2)
    resi$fp.nu.mse <- rowMeans(coefs.nu[,-c(1,4,5,6)]^2)
    
    m <- cbind(m, coefs.mu[,c(2,4,6,7)], coefs.sigma[,c(3,5,6)],coefs.nu[,c(4,5,6)], resi)
    m1 <- rbind(m1, m)
  }
  # Var deselection
  # variables which contribute less than 1% of the total risk reduction are deselected
  if(nobs < 100000){
    f <- grep("x", names(d), value = TRUE)
    f <- paste("y ~", paste(f, collapse = "+"))
    f <- as.formula(f)
    
    fam <- ZANBI()
    
    # gradient boosting model b is input for vardeselect
    timegradboost <- time
    ptm <- proc.time()
    b <- DeselectBoost(b, fam = as.families(ZANBI), 
                       tau = 0.01)  
    
    time <- c(proc.time() - ptm)[3]
    time <- time + timegradboost
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- list()
    p$mu <- predict(b$mu,newdata=nd)
    p$sigma <- predict(b$sigma,newdata=nd)
    p$nu <- predict(b$nu,newdata=nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    nu <- mean((nd$eta.nu - p$nu)^2)
    
    # transformation from linear predictor scale to dist parameter scale
    par.trafo <- function(list){
      list$mu <- exp(list$mu)
      list$sigma <- exp(list$sigma)
      list$nu <- 1/(1+exp(-list$nu))
      return(list)
    }
    p <- par.trafo(p)
    
    k = max(350, nd$y)
    mat <- fam$d(0, p)
    for(i in 1:k){
      pp <- fam$d(i, p)
      mat <- rbind(mat, pp)
      cat("\r");cat("i =",i, "/", k)
    }
    mat <- t(mat)
    colnames(mat) <- 0:k
    mean1 <- rowSums(mat %*% diag(c(0:k)))
    
    mse1 <- mean((nd$y - mean1)^2)
    crps1 <- try(rps(nd$y, mat),TRUE)
    if(!is.numeric(crps1)) crps1 <- NaN
    
    ll <- fam[["d"]]
    l <- try(sum(ll(nd$y, p, log = T)), T)
    if(!is.numeric(l)) l <- NaN
    
    # combine results
    m <- rbind(
      "var_des" = c("mu" = mu, "sigma" = sigma, "nu" = nu, "mse" = mse1, "crps" = crps1, "logLiK" = l, 
                    "mstop" = NA, time)
    )
    
    m <- as.data.frame(m)
    rownames(m) <- NULL
    
    m$method <- "var_des"
    m$nobs <- nobs
    m$nnoise <- nnoise
    m$rho <- rho
    m$rep <- seed
    
    
    # Extract coefs and calc number of true and false positives
    c <- list()
    cc <- coef(b$mu, off2int = TRUE)
    if(is.null(cc)) c$mu <- unlist(list( "(Intercept)" = 0)) else c$mu <- cc
    cc <- coef(b$sigma, off2int = TRUE)
    if(is.null(cc)) c$sigma <- unlist(list( "(Intercept)" = 0)) else c$sigma <- cc
    cc <- coef(b$nu, off2int = TRUE)
    if(is.null(cc)) c$nu <- unlist(list( "(Intercept)" = 0)) else c$nu <- cc
    
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      c.nu <- c[["nu"]]
      cc.mu <- cc.sigma <- cc.nu <-  rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      cc.nu[1] <- c.nu["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
        cc.nu[i+1] <- ifelse(paste0("x",i) %in% names(c.nu) ,c.nu[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma, "nu" = cc.nu))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    nu.names <- c("nu.intercept", paste0("nu.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    names(c[["nu"]]) <- nu.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    coefs.nu <- rbind.data.frame(c[["nu"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- rownames(coefs.nu) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    colnames(coefs.nu) <- nu.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,c(2,4,6,7)] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,c(3,5,6)] != 0)
    resi$tp.nu <- rowSums(coefs.nu[,c(4,5,6)] != 0)
    
    resi$fp.mu <- rowSums(coefs.mu[,-c(1,2,4,6,7)] != 0)
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,3,5,6)] != 0)
    resi$fp.nu <- rowSums(coefs.nu[,-c(1,4,5,6)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2,4,6,7)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(3,5,6)] -  t(matrix(rep(sigmac, 1), nrow = 3)) ) ^2)
    resi$tp.nu.mse <- rowMeans((coefs.nu[,c(4,5,6)] - t(matrix(rep(nuc, 1), nrow = 3)) ) ^2)
    
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1,2,4,6,7)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,3,5,6)]^2)
    resi$fp.nu.mse <- rowMeans(coefs.nu[,-c(1,4,5,6)]^2)
    
    m <- cbind(m, coefs.mu[,c(2,4,6,7)], coefs.sigma[,c(3,5,6)],coefs.nu[,c(4,5,6)], resi)
    m1 <- rbind(m1, m)
  }
  
  # bamlss stability selection
  # Estimate selection frequencies of variables on B = 100 resampled data sets
  if(nobs < 100000){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .),
      update(f, nu ~ .)
    )
    
    fam <- ZANBI()
    
    # bamlss stability selection
    q = max(round(sqrt(6+nnoise)*3), 8)
    ptm <- proc.time()
    sel <- bamlss::stabsel(formula = f, data = d, q = q, B = 10, 
                           thr = 0.8, family = fam)
    stabf <- function(sel, thr = 0.8, B = 100, par = c("mu", "sigma")){
      tab <- sel[["table"]]
      tab <- tab[tab/B >= thr ]
      p <- names(tab)
      modelID <- substring(p, regexpr("(?<=\\.)[a-z0-9]+$", p, 
                                      perl = TRUE))
      p <- gsub("\\.[a-z0-9]+$", "", p)
      fn <- lapply(par, FUN = function(i){
        p2 <- p[modelID == i]
        p2 <- paste(p2, collapse = " + ", sep = "")
        if (p2 == "") p2 <- "1"
        lhs <- ifelse(i == "mu", "y", "")
        as.formula(paste0(lhs,"~",p2))
      }
      )
      
      return(fn)
    }
    nf <- stabf(sel, thr = 0.8, B = 10, par = c("mu", "sigma", "nu"))
    b <- bamlss(nf, data = d, family = fam)
    time <- c(proc.time() - ptm)[3]
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    nu <- mean((nd$eta.nu - p$nu)^2)
    
    # transformation from linear predictor scale to dist parameter scale
    par.trafo <- function(list){
      list$mu <- exp(list$mu)
      list$sigma <- exp(list$sigma)
      list$nu <- 1/(1+exp(-list$nu))
      return(list)
    }
    p <- par.trafo(p)
    
    k = max(350, nd$y)
    mat <- fam$d(0, p)
    for(i in 1:k){
      pp <- fam$d(i, p)
      mat <- rbind(mat, pp)
      cat("\r");cat("i =",i, "/", k)
    }
    mat <- t(mat)
    colnames(mat) <- 0:k
    mean1 <- rowSums(mat %*% diag(c(0:k)))
    
    mse1 <- mean((nd$y - mean1)^2)
    crps1 <- try(rps(nd$y, mat),TRUE)
    if(!is.numeric(crps1)) crps1 <- NaN
    
    ll <- fam[["d"]]
    l <- try(sum(ll(nd$y, p, log = T)), T)
    if(!is.numeric(l)) l <- NaN
    
    
    # combine results
    m <- rbind(
      "bamlss_stab" = c("mu" = mu, "sigma" = sigma, "nu" = nu, "mse" = mse1, "crps" = crps1, "logLiK" = l, "mstop" = NA, time)
    )
    
    m <- as.data.frame(m)
    rownames(m) <- NULL
    
    m$method <- "bamlss_stab"
    m$nobs <- nobs
    m$nnoise <- nnoise
    m$rho <- rho
    m$rep <- seed
    
    
    # Extract coefs and calc number of true and false positives
    c <- coef(b)[,1]
    
    nn <- names(c)
    nmu <- nn[grepl("mu",nn)]
    nsigma <- nn[grepl("sigma",nn)]
    nnu <- nn[grepl("nu",nn)]
    c.mu <- c[nmu]
    c.sigma <- c[nsigma]
    c.nu <- c[nnu]
    nmu <- gsub("mu.p.","",nmu)
    nsigma <- gsub("sigma.p.","",nsigma)
    nnu <- gsub("nu.p.","",nnu)
    names(c.mu) <- nmu
    names(c.sigma) <- nsigma
    names(c.nu) <- nnu
    names(c.mu)[1] <- "(Intercept)"
    names(c.sigma)[1] <- "(Intercept)"
    names(c.nu)[1] <- "(Intercept)"
    c.mu <- c.mu[-length(c.mu)] 
    c.sigma <- c.sigma[-length(c.sigma)] 
    c.nu <- c.nu[-length(c.nu)] 
    c <- list("mu" = c.mu,"sigma" = c.sigma, "nu" = c.nu)
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      c.nu <- c[["nu"]]
      cc.mu <- cc.sigma <- cc.nu <-  rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      cc.nu[1] <- c.nu["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
        cc.nu[i+1] <- ifelse(paste0("x",i) %in% names(c.nu) ,c.nu[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma, "nu" = cc.nu))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    nu.names <- c("nu.intercept", paste0("nu.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    names(c[["nu"]]) <- nu.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    coefs.nu <- rbind.data.frame(c[["nu"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- rownames(coefs.nu) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    colnames(coefs.nu) <- nu.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,c(2,4,6,7)] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,c(3,5,6)] != 0)
    resi$tp.nu <- rowSums(coefs.nu[,c(4,5,6)] != 0)
    
    resi$fp.mu <- rowSums(coefs.mu[,-c(1,2,4,6,7)] != 0)
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,3,5,6)] != 0)
    resi$fp.nu <- rowSums(coefs.nu[,-c(1,4,5,6)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2,4,6,7)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(3,5,6)] -  t(matrix(rep(sigmac, 1), nrow = 3)) ) ^2)
    resi$tp.nu.mse <- rowMeans((coefs.nu[,c(4,5,6)] - t(matrix(rep(nuc, 1), nrow = 3)) ) ^2)
    
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1,2,4,6,7)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,3,5,6)]^2)
    resi$fp.nu.mse <- rowMeans(coefs.nu[,-c(1,4,5,6)]^2)
    
    m <- cbind(m, coefs.mu[,c(2,4,6,7)], coefs.sigma[,c(3,5,6)],coefs.nu[,c(4,5,6)], resi)
    m1 <- rbind(m1, m)
  }
  
  # Full batch variant
  if(nobs < 100000){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .),
      update(f, nu ~ .)
    )
    
    fam <- ZANBI()
    
    # SDR1 noncyclic updating, cap = 0 and early stopping selection via bic
    ptm <- proc.time()
    maxit <- 1000
    batch_ids = nobs
    b <- sdr(formula = f, family = fam, data = d, maxit = maxit, updating = "noncyclic", 
             CF = FALSE, batch_ids = batch_ids)
    bic <- AIC(b, K = log(batch_ids))
    itermax <- max(20,which.min(bic))
    c <- coef(b, mstop = itermax, refit = TRUE)
    fu <- newformula(b, mstop = itermax)
    b <- sdr(formula = fu, family = fam, data = d, maxit = maxit, 
             coef_start = c, batch_ids = batch_ids, updating = "noncyclic",
             init = FALSE, CF = FALSE)
    time <- c(proc.time() - ptm)[3]
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    nu <- mean((nd$eta.nu - p$nu)^2)
    
    # transformation from linear predictor scale to dist parameter scale
    par.trafo <- function(list){
      list$mu <- exp(list$mu)
      list$sigma <- exp(list$sigma)
      list$nu <- 1/(1+exp(-list$nu))
      return(list)
    }
    p <- par.trafo(p)
    
    k = max(350, nd$y)
    mat <- fam$d(0, p)
    for(i in 1:k){
      pp <- fam$d(i, p)
      mat <- rbind(mat, pp)
      cat("\r");cat("i =",i, "/", k)
    }
    mat <- t(mat)
    colnames(mat) <- 0:k
    mean1 <- rowSums(mat %*% diag(c(0:k)))
    
    mse1 <- mean((nd$y - mean1)^2)
    crps1 <- try(rps(nd$y, mat),TRUE)
    if(!is.numeric(crps1)) crps1 <- NaN
    
    ll <- fam[["d"]]
    l <- try(sum(ll(nd$y, p, log = T)), T)
    if(!is.numeric(l)) l <- NaN
    
    
    # combine results
    m <- rbind(
      "sdr1" = c("mu" = mu, "sigma" = sigma, "nu" = nu, "mse" = mse1, "crps" = crps1, "logLiK" = l, "mstop" = maxit, time)
    )
    
    m <- as.data.frame(m)
    rownames(m) <- NULL
    
    m$method <- "sdr1"
    m$nobs <- nobs
    m$nnoise <- nnoise
    m$rho <- rho
    m$rep <- seed
    
    
    # Extract coefs and calc number of true and false positives
    c <- coef(b)
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      c.nu <- c[["nu"]]
      cc.mu <- cc.sigma <- cc.nu <-  rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      cc.nu[1] <- c.nu["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
        cc.nu[i+1] <- ifelse(paste0("x",i) %in% names(c.nu) ,c.nu[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma, "nu" = cc.nu))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    nu.names <- c("nu.intercept", paste0("nu.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    names(c[["nu"]]) <- nu.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    coefs.nu <- rbind.data.frame(c[["nu"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- rownames(coefs.nu) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    colnames(coefs.nu) <- nu.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,c(2,4,6,7)] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,c(3,5,6)] != 0)
    resi$tp.nu <- rowSums(coefs.nu[,c(4,5,6)] != 0)
    
    resi$fp.mu <- rowSums(coefs.mu[,-c(1,2,4,6,7)] != 0)
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,3,5,6)] != 0)
    resi$fp.nu <- rowSums(coefs.nu[,-c(1,4,5,6)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2,4,6,7)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(3,5,6)] -  t(matrix(rep(sigmac, 1), nrow = 3)) ) ^2)
    resi$tp.nu.mse <- rowMeans((coefs.nu[,c(4,5,6)] - t(matrix(rep(nuc, 1), nrow = 3)) ) ^2)
    
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1,2,4,6,7)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,3,5,6)]^2)
    resi$fp.nu.mse <- rowMeans(coefs.nu[,-c(1,4,5,6)]^2)
    
    m <- cbind(m, coefs.mu[,c(2,4,6,7)], coefs.sigma[,c(3,5,6)],coefs.nu[,c(4,5,6)], resi)
    m1 <- rbind(m1, m)
  }
  
  # SDR2 best subset updating, CF = FALSE and early stopping selection via bic
  if(nobs < 100000){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .),
      update(f, nu ~ .)
    )
    
    fam <- ZANBI()
    
    # SDR2 best subset updating, cap = 0 and early stopping selection via bic
    ptm <- proc.time()
    maxit <- 1000
    batch_ids = nobs
    b <- sdr(formula = f, family = fam, data = d, maxit = maxit, updating = "bs", 
             CF = FALSE, batch_ids = batch_ids)
    bic <- AIC(b, K = log(batch_ids))
    itermax <- max(20,which.min(bic))
    c <- coef(b, mstop = itermax, refit = TRUE)
    fu <- newformula(b, mstop = itermax)
    b <- sdr(formula = fu, family = fam, data = d, maxit = maxit, 
             coef_start = c, batch_ids = batch_ids, updating = "bs",
             init = FALSE, CF = FALSE)
    time <- c(proc.time() - ptm)[3]
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    nu <- mean((nd$eta.nu - p$nu)^2)
    
    # transformation from linear predictor scale to dist parameter scale
    par.trafo <- function(list){
      list$mu <- exp(list$mu)
      list$sigma <- exp(list$sigma)
      list$nu <- 1/(1+exp(-list$nu))
      return(list)
    }
    p <- par.trafo(p)
    
    k = max(350, nd$y)
    mat <- fam$d(0, p)
    for(i in 1:k){
      pp <- fam$d(i, p)
      mat <- rbind(mat, pp)
      cat("\r");cat("i =",i, "/", k)
    }
    mat <- t(mat)
    colnames(mat) <- 0:k
    mean1 <- rowSums(mat %*% diag(c(0:k)))
    
    mse1 <- mean((nd$y - mean1)^2)
    crps1 <- try(rps(nd$y, mat),TRUE)
    if(!is.numeric(crps1)) crps1 <- NaN
    
    ll <- fam[["d"]]
    l <- try(sum(ll(nd$y, p, log = T)), T)
    if(!is.numeric(l)) l <- NaN
    
    
    # combine results
    m <- rbind(
      "sdr2" = c("mu" = mu, "sigma" = sigma, "nu" = nu, "mse" = mse1, "crps" = crps1, "logLiK" = l, "mstop" = maxit, time)
    )
    
    m <- as.data.frame(m)
    rownames(m) <- NULL
    
    m$method <- "sdr2"
    m$nobs <- nobs
    m$nnoise <- nnoise
    m$rho <- rho
    m$rep <- seed
    
    
    # Extract coefs and calc number of true and false positives
    c <- coef(b)
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      c.nu <- c[["nu"]]
      cc.mu <- cc.sigma <- cc.nu <-  rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      cc.nu[1] <- c.nu["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
        cc.nu[i+1] <- ifelse(paste0("x",i) %in% names(c.nu) ,c.nu[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma, "nu" = cc.nu))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    nu.names <- c("nu.intercept", paste0("nu.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    names(c[["nu"]]) <- nu.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    coefs.nu <- rbind.data.frame(c[["nu"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- rownames(coefs.nu) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    colnames(coefs.nu) <- nu.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,c(2,4,6,7)] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,c(3,5,6)] != 0)
    resi$tp.nu <- rowSums(coefs.nu[,c(4,5,6)] != 0)
    
    resi$fp.mu <- rowSums(coefs.mu[,-c(1,2,4,6,7)] != 0)
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,3,5,6)] != 0)
    resi$fp.nu <- rowSums(coefs.nu[,-c(1,4,5,6)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2,4,6,7)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(3,5,6)] -  t(matrix(rep(sigmac, 1), nrow = 3)) ) ^2)
    resi$tp.nu.mse <- rowMeans((coefs.nu[,c(4,5,6)] - t(matrix(rep(nuc, 1), nrow = 3)) ) ^2)
    
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1,2,4,6,7)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,3,5,6)]^2)
    resi$fp.nu.mse <- rowMeans(coefs.nu[,-c(1,4,5,6)]^2)
    
    m <- cbind(m, coefs.mu[,c(2,4,6,7)], coefs.sigma[,c(3,5,6)],coefs.nu[,c(4,5,6)], resi)
    m1 <- rbind(m1, m)
  }
  
  # SDR3 noncyclic updating and CF = TRUE 
  if(nobs < 100000){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .),
      update(f, nu ~ .)
    )
    
    fam <- ZANBI()
    
    # SDR3 noncyclic updating, CF = TRUE and early stopping selection via bic
    ptm <- proc.time()
    maxit <- 1000
    batch_ids = nobs
    b <- sdr(formula = f, family = fam, data = d, maxit = maxit, updating = "noncyclic", 
             CF = TRUE, batch_ids = batch_ids)
    bic <- AIC(b, K = log(batch_ids))
    itermax <- max(20,which.min(bic))
    c <- coef(b, mstop = itermax, refit = TRUE)
    fu <- newformula(b, mstop = itermax)
    b <- sdr(formula = fu, family = fam, data = d, maxit = maxit, 
             coef_start = c, batch_ids = batch_ids, updating = "noncyclic",
             init = FALSE, CF = FALSE)
    time <- c(proc.time() - ptm)[3]
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    nu <- mean((nd$eta.nu - p$nu)^2)
    
    # transformation from linear predictor scale to dist parameter scale
    par.trafo <- function(list){
      list$mu <- exp(list$mu)
      list$sigma <- exp(list$sigma)
      list$nu <- 1/(1+exp(-list$nu))
      return(list)
    }
    p <- par.trafo(p)
    
    k = max(350, nd$y)
    mat <- fam$d(0, p)
    for(i in 1:k){
      pp <- fam$d(i, p)
      mat <- rbind(mat, pp)
      cat("\r");cat("i =",i, "/", k)
    }
    mat <- t(mat)
    colnames(mat) <- 0:k
    mean1 <- rowSums(mat %*% diag(c(0:k)))
    
    mse1 <- mean((nd$y - mean1)^2)
    crps1 <- try(rps(nd$y, mat),TRUE)
    if(!is.numeric(crps1)) crps1 <- NaN
    
    ll <- fam[["d"]]
    l <- try(sum(ll(nd$y, p, log = T)), T)
    if(!is.numeric(l)) l <- NaN
    
    
    # combine results
    m <- rbind(
      "sdr3" = c("mu" = mu, "sigma" = sigma, "nu" = nu, "mse" = mse1, "crps" = crps1, "logLiK" = l, "mstop" = maxit, time)
    )
    
    m <- as.data.frame(m)
    rownames(m) <- NULL
    
    m$method <- "sdr3"
    m$nobs <- nobs
    m$nnoise <- nnoise
    m$rho <- rho
    m$rep <- seed
    
    
    # Extract coefs and calc number of true and false positives
    c <- coef(b)
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      c.nu <- c[["nu"]]
      cc.mu <- cc.sigma <- cc.nu <-  rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      cc.nu[1] <- c.nu["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
        cc.nu[i+1] <- ifelse(paste0("x",i) %in% names(c.nu) ,c.nu[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma, "nu" = cc.nu))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    nu.names <- c("nu.intercept", paste0("nu.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    names(c[["nu"]]) <- nu.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    coefs.nu <- rbind.data.frame(c[["nu"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- rownames(coefs.nu) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    colnames(coefs.nu) <- nu.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,c(2,4,6,7)] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,c(3,5,6)] != 0)
    resi$tp.nu <- rowSums(coefs.nu[,c(4,5,6)] != 0)
    
    resi$fp.mu <- rowSums(coefs.mu[,-c(1,2,4,6,7)] != 0)
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,3,5,6)] != 0)
    resi$fp.nu <- rowSums(coefs.nu[,-c(1,4,5,6)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2,4,6,7)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(3,5,6)] -  t(matrix(rep(sigmac, 1), nrow = 3)) ) ^2)
    resi$tp.nu.mse <- rowMeans((coefs.nu[,c(4,5,6)] - t(matrix(rep(nuc, 1), nrow = 3)) ) ^2)
    
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1,2,4,6,7)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,3,5,6)]^2)
    resi$fp.nu.mse <- rowMeans(coefs.nu[,-c(1,4,5,6)]^2)
    
    m <- cbind(m, coefs.mu[,c(2,4,6,7)], coefs.sigma[,c(3,5,6)],coefs.nu[,c(4,5,6)], resi)
    m1 <- rbind(m1, m)
  }
  
  # SDR4 best subset updating and CF = TRUE
  if(nobs < 100000){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .),
      update(f, nu ~ .)
    )
    
    fam <- ZANBI()
    
    # SDR2 best subset updating, cap = 0 and early stopping selection via bic
    ptm <- proc.time()
    maxit <- 1000
    batch_ids = nobs
    b <- sdr(formula = f, family = fam, data = d, maxit = maxit, updating = "bs", 
             CF = TRUE, batch_ids = batch_ids)
    bic <- AIC(b, K = log(batch_ids))
    itermax <- max(20,which.min(bic))
    c <- coef(b, mstop = itermax, refit = TRUE)
    fu <- newformula(b, mstop = itermax)
    b <- sdr(formula = fu, family = fam, data = d, maxit = maxit, 
             coef_start = c, batch_ids = batch_ids, updating = "bs",
             init = FALSE, CF = FALSE)
    time <- c(proc.time() - ptm)[3]
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    nu <- mean((nd$eta.nu - p$nu)^2)
    
    # transformation from linear predictor scale to dist parameter scale
    par.trafo <- function(list){
      list$mu <- exp(list$mu)
      list$sigma <- exp(list$sigma)
      list$nu <- 1/(1+exp(-list$nu))
      return(list)
    }
    p <- par.trafo(p)
    
    k = max(350, nd$y)
    mat <- fam$d(0, p)
    for(i in 1:k){
      pp <- fam$d(i, p)
      mat <- rbind(mat, pp)
      cat("\r");cat("i =",i, "/", k)
    }
    mat <- t(mat)
    colnames(mat) <- 0:k
    mean1 <- rowSums(mat %*% diag(c(0:k)))
    
    mse1 <- mean((nd$y - mean1)^2)
    crps1 <- try(rps(nd$y, mat),TRUE)
    if(!is.numeric(crps1)) crps1 <- NaN
    
    ll <- fam[["d"]]
    l <- try(sum(ll(nd$y, p, log = T)), T)
    if(!is.numeric(l)) l <- NaN
    
    
    # combine results
    m <- rbind(
      "sdr4" = c("mu" = mu, "sigma" = sigma, "nu" = nu, "mse" = mse1, "crps" = crps1, "logLiK" = l, "mstop" = maxit, time)
    )
    
    m <- as.data.frame(m)
    rownames(m) <- NULL
    
    m$method <- "sdr4"
    m$nobs <- nobs
    m$nnoise <- nnoise
    m$rho <- rho
    m$rep <- seed
    
    
    # Extract coefs and calc number of true and false positives
    c <- coef(b)
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      c.nu <- c[["nu"]]
      cc.mu <- cc.sigma <- cc.nu <-  rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      cc.nu[1] <- c.nu["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
        cc.nu[i+1] <- ifelse(paste0("x",i) %in% names(c.nu) ,c.nu[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma, "nu" = cc.nu))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    nu.names <- c("nu.intercept", paste0("nu.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    names(c[["nu"]]) <- nu.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    coefs.nu <- rbind.data.frame(c[["nu"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- rownames(coefs.nu) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    colnames(coefs.nu) <- nu.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,c(2,4,6,7)] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,c(3,5,6)] != 0)
    resi$tp.nu <- rowSums(coefs.nu[,c(4,5,6)] != 0)
    
    resi$fp.mu <- rowSums(coefs.mu[,-c(1,2,4,6,7)] != 0)
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,3,5,6)] != 0)
    resi$fp.nu <- rowSums(coefs.nu[,-c(1,4,5,6)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2,4,6,7)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(3,5,6)] -  t(matrix(rep(sigmac, 1), nrow = 3)) ) ^2)
    resi$tp.nu.mse <- rowMeans((coefs.nu[,c(4,5,6)] - t(matrix(rep(nuc, 1), nrow = 3)) ) ^2)
    
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1,2,4,6,7)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,3,5,6)]^2)
    resi$fp.nu.mse <- rowMeans(coefs.nu[,-c(1,4,5,6)]^2)
    
    m <- cbind(m, coefs.mu[,c(2,4,6,7)], coefs.sigma[,c(3,5,6)],coefs.nu[,c(4,5,6)], resi)
    m1 <- rbind(m1, m)
  }
  
  # Batchwise variant
  # SDR1 bw noncyclic updating, CF = FALSE and early stopping selection via bic
  if(TRUE){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .),
      update(f, nu ~ .)
    )
    
    fam <- ZANBI()
    
    # SDR1 noncyclic updating, cap = 0 and early stopping selection via bic
    ptm <- proc.time()
    maxit <- 1000
    batch_ids <- min(round(2*nobs/3), 10000) 
    b <- sdr(formula = f, family = fam, data = d, maxit = maxit, updating = "noncyclic", 
             CF = FALSE, batch_ids = batch_ids)
    bic <- AIC(b, K = log(batch_ids))
    itermax <- max(20,which.min(bic))
    c <- coef(b, mstop = itermax, refit = TRUE)
    fu <- newformula(b, mstop = itermax)
    b <- sdr(formula = fu, family = fam, data = d, maxit = maxit, 
             coef_start = c, batch_ids = batch_ids, updating = "noncyclic",
             init = FALSE, CF = FALSE)
    time <- c(proc.time() - ptm)[3]
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    nu <- mean((nd$eta.nu - p$nu)^2)
    
    # transformation from linear predictor scale to dist parameter scale
    par.trafo <- function(list){
      list$mu <- exp(list$mu)
      list$sigma <- exp(list$sigma)
      list$nu <- 1/(1+exp(-list$nu))
      return(list)
    }
    p <- par.trafo(p)
    
    k = max(350, nd$y)
    mat <- fam$d(0, p)
    for(i in 1:k){
      pp <- fam$d(i, p)
      mat <- rbind(mat, pp)
      cat("\r");cat("i =",i, "/", k)
    }
    mat <- t(mat)
    colnames(mat) <- 0:k
    mean1 <- rowSums(mat %*% diag(c(0:k)))
    
    mse1 <- mean((nd$y - mean1)^2)
    crps1 <- try(rps(nd$y, mat),TRUE)
    if(!is.numeric(crps1)) crps1 <- NaN
    
    ll <- fam[["d"]]
    l <- try(sum(ll(nd$y, p, log = T)), T)
    if(!is.numeric(l)) l <- NaN
    
    
    # combine results
    m <- rbind(
      "sdr1_bw" = c("mu" = mu, "sigma" = sigma, "nu" = nu, "mse" = mse1, "crps" = crps1, "logLiK" = l, "mstop" = maxit, time)
    )
    
    m <- as.data.frame(m)
    rownames(m) <- NULL
    
    m$method <- "sdr1_bw"
    m$nobs <- nobs
    m$nnoise <- nnoise
    m$rho <- rho
    m$rep <- seed
    
    
    # Extract coefs and calc number of true and false positives
    c <- coef(b)
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      c.nu <- c[["nu"]]
      cc.mu <- cc.sigma <- cc.nu <-  rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      cc.nu[1] <- c.nu["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
        cc.nu[i+1] <- ifelse(paste0("x",i) %in% names(c.nu) ,c.nu[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma, "nu" = cc.nu))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    nu.names <- c("nu.intercept", paste0("nu.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    names(c[["nu"]]) <- nu.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    coefs.nu <- rbind.data.frame(c[["nu"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- rownames(coefs.nu) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    colnames(coefs.nu) <- nu.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,c(2,4,6,7)] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,c(3,5,6)] != 0)
    resi$tp.nu <- rowSums(coefs.nu[,c(4,5,6)] != 0)
    
    resi$fp.mu <- rowSums(coefs.mu[,-c(1,2,4,6,7)] != 0)
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,3,5,6)] != 0)
    resi$fp.nu <- rowSums(coefs.nu[,-c(1,4,5,6)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2,4,6,7)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(3,5,6)] -  t(matrix(rep(sigmac, 1), nrow = 3)) ) ^2)
    resi$tp.nu.mse <- rowMeans((coefs.nu[,c(4,5,6)] - t(matrix(rep(nuc, 1), nrow = 3)) ) ^2)
    
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1,2,4,6,7)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,3,5,6)]^2)
    resi$fp.nu.mse <- rowMeans(coefs.nu[,-c(1,4,5,6)]^2)
    
    m <- cbind(m, coefs.mu[,c(2,4,6,7)], coefs.sigma[,c(3,5,6)],coefs.nu[,c(4,5,6)], resi)
    m1 <- rbind(m1, m)
  }
  
  # SDR2 bw best subset updating, CF = FALSE and early stopping selection via bic
  if(TRUE){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .),
      update(f, nu ~ .)
    )
    
    fam <- ZANBI()
    
    # SDR2 best subset updating, cap = 0 and early stopping selection via bic
    ptm <- proc.time()
    maxit <- 1000
    batch_ids <- min(round(2*nobs/3), 10000) 
    b <- sdr(formula = f, family = fam, data = d, maxit = maxit, updating = "bs", 
             CF = FALSE, batch_ids = batch_ids)
    bic <- AIC(b, K = log(batch_ids))
    itermax <- max(20,which.min(bic))
    c <- coef(b, mstop = itermax, refit = TRUE)
    fu <- newformula(b, mstop = itermax)
    b <- sdr(formula = fu, family = fam, data = d, maxit = maxit, 
             coef_start = c, batch_ids = batch_ids, updating = "bs",
             init = FALSE, CF = FALSE)
    time <- c(proc.time() - ptm)[3]
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    nu <- mean((nd$eta.nu - p$nu)^2)
    
    # transformation from linear predictor scale to dist parameter scale
    par.trafo <- function(list){
      list$mu <- exp(list$mu)
      list$sigma <- exp(list$sigma)
      list$nu <- 1/(1+exp(-list$nu))
      return(list)
    }
    p <- par.trafo(p)
    
    k = max(350, nd$y)
    mat <- fam$d(0, p)
    for(i in 1:k){
      pp <- fam$d(i, p)
      mat <- rbind(mat, pp)
      cat("\r");cat("i =",i, "/", k)
    }
    mat <- t(mat)
    colnames(mat) <- 0:k
    mean1 <- rowSums(mat %*% diag(c(0:k)))
    
    mse1 <- mean((nd$y - mean1)^2)
    crps1 <- try(rps(nd$y, mat),TRUE)
    if(!is.numeric(crps1)) crps1 <- NaN
    
    ll <- fam[["d"]]
    l <- try(sum(ll(nd$y, p, log = T)), T)
    if(!is.numeric(l)) l <- NaN
    
    
    # combine results
    m <- rbind(
      "sdr2_bw" = c("mu" = mu, "sigma" = sigma, "nu" = nu, "mse" = mse1, "crps" = crps1, "logLiK" = l, "mstop" = maxit, time)
    )
    
    m <- as.data.frame(m)
    rownames(m) <- NULL
    
    m$method <- "sdr2_bw"
    m$nobs <- nobs
    m$nnoise <- nnoise
    m$rho <- rho
    m$rep <- seed
    
    
    # Extract coefs and calc number of true and false positives
    c <- coef(b)
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      c.nu <- c[["nu"]]
      cc.mu <- cc.sigma <- cc.nu <-  rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      cc.nu[1] <- c.nu["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
        cc.nu[i+1] <- ifelse(paste0("x",i) %in% names(c.nu) ,c.nu[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma, "nu" = cc.nu))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    nu.names <- c("nu.intercept", paste0("nu.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    names(c[["nu"]]) <- nu.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    coefs.nu <- rbind.data.frame(c[["nu"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- rownames(coefs.nu) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    colnames(coefs.nu) <- nu.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,c(2,4,6,7)] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,c(3,5,6)] != 0)
    resi$tp.nu <- rowSums(coefs.nu[,c(4,5,6)] != 0)
    
    resi$fp.mu <- rowSums(coefs.mu[,-c(1,2,4,6,7)] != 0)
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,3,5,6)] != 0)
    resi$fp.nu <- rowSums(coefs.nu[,-c(1,4,5,6)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2,4,6,7)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(3,5,6)] -  t(matrix(rep(sigmac, 1), nrow = 3)) ) ^2)
    resi$tp.nu.mse <- rowMeans((coefs.nu[,c(4,5,6)] - t(matrix(rep(nuc, 1), nrow = 3)) ) ^2)
    
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1,2,4,6,7)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,3,5,6)]^2)
    resi$fp.nu.mse <- rowMeans(coefs.nu[,-c(1,4,5,6)]^2)
    
    m <- cbind(m, coefs.mu[,c(2,4,6,7)], coefs.sigma[,c(3,5,6)],coefs.nu[,c(4,5,6)], resi)
    m1 <- rbind(m1, m)
  }
  
  # SDR3 bw noncyclic updating and CF = TRUE 
  if(TRUE){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .),
      update(f, nu ~ .)
    )
    
    fam <- ZANBI()
    
    # SDR3 noncyclic updating, CF = TRUE and early stopping selection via bic
    ptm <- proc.time()
    maxit <- 1000
    batch_ids <- min(round(2*nobs/3), 10000) 
    b <- sdr(formula = f, family = fam, data = d, maxit = maxit, updating = "noncyclic", 
             CF = TRUE, batch_ids = batch_ids)
    bic <- AIC(b, K = log(batch_ids))
    itermax <- max(20,which.min(bic))
    c <- coef(b, mstop = itermax, refit = TRUE)
    fu <- newformula(b, mstop = itermax)
    b <- sdr(formula = fu, family = fam, data = d, maxit = maxit, 
             coef_start = c, batch_ids = batch_ids, updating = "noncyclic",
             init = FALSE, CF = FALSE)
    time <- c(proc.time() - ptm)[3]
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    nu <- mean((nd$eta.nu - p$nu)^2)
    
    # transformation from linear predictor scale to dist parameter scale
    par.trafo <- function(list){
      list$mu <- exp(list$mu)
      list$sigma <- exp(list$sigma)
      list$nu <- 1/(1+exp(-list$nu))
      return(list)
    }
    p <- par.trafo(p)
    
    k = max(350, nd$y)
    mat <- fam$d(0, p)
    for(i in 1:k){
      pp <- fam$d(i, p)
      mat <- rbind(mat, pp)
      cat("\r");cat("i =",i, "/", k)
    }
    mat <- t(mat)
    colnames(mat) <- 0:k
    mean1 <- rowSums(mat %*% diag(c(0:k)))
    
    mse1 <- mean((nd$y - mean1)^2)
    crps1 <- try(rps(nd$y, mat),TRUE)
    if(!is.numeric(crps1)) crps1 <- NaN
    
    ll <- fam[["d"]]
    l <- try(sum(ll(nd$y, p, log = T)), T)
    if(!is.numeric(l)) l <- NaN
    
    
    # combine results
    m <- rbind(
      "sdr3_bw" = c("mu" = mu, "sigma" = sigma, "nu" = nu, "mse" = mse1, "crps" = crps1, "logLiK" = l, "mstop" = maxit, time)
    )
    
    m <- as.data.frame(m)
    rownames(m) <- NULL
    
    m$method <- "sdr3_bw"
    m$nobs <- nobs
    m$nnoise <- nnoise
    m$rho <- rho
    m$rep <- seed
    
    
    # Extract coefs and calc number of true and false positives
    c <- coef(b)
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      c.nu <- c[["nu"]]
      cc.mu <- cc.sigma <- cc.nu <-  rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      cc.nu[1] <- c.nu["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
        cc.nu[i+1] <- ifelse(paste0("x",i) %in% names(c.nu) ,c.nu[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma, "nu" = cc.nu))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    nu.names <- c("nu.intercept", paste0("nu.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    names(c[["nu"]]) <- nu.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    coefs.nu <- rbind.data.frame(c[["nu"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- rownames(coefs.nu) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    colnames(coefs.nu) <- nu.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,c(2,4,6,7)] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,c(3,5,6)] != 0)
    resi$tp.nu <- rowSums(coefs.nu[,c(4,5,6)] != 0)
    
    resi$fp.mu <- rowSums(coefs.mu[,-c(1,2,4,6,7)] != 0)
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,3,5,6)] != 0)
    resi$fp.nu <- rowSums(coefs.nu[,-c(1,4,5,6)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2,4,6,7)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(3,5,6)] -  t(matrix(rep(sigmac, 1), nrow = 3)) ) ^2)
    resi$tp.nu.mse <- rowMeans((coefs.nu[,c(4,5,6)] - t(matrix(rep(nuc, 1), nrow = 3)) ) ^2)
    
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1,2,4,6,7)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,3,5,6)]^2)
    resi$fp.nu.mse <- rowMeans(coefs.nu[,-c(1,4,5,6)]^2)
    
    m <- cbind(m, coefs.mu[,c(2,4,6,7)], coefs.sigma[,c(3,5,6)],coefs.nu[,c(4,5,6)], resi)
    m1 <- rbind(m1, m)
  }
  
  # SDR4 bw best subset updating and CF = TRUE
  if(TRUE){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .),
      update(f, nu ~ .)
    )
    
    fam <- ZANBI()
    
    # SDR2 best subset updating, CF =TRUE and early stopping selection via bic
    ptm <- proc.time()
    maxit <- 1000
    batch_ids <- min(round(2*nobs/3), 10000) 
    b <- sdr(formula = f, family = fam, data = d, maxit = maxit, updating = "bs", 
             CF = TRUE, batch_ids = batch_ids)
    bic <- AIC(b, K = log(batch_ids))
    itermax <- max(20,which.min(bic))
    c <- coef(b, mstop = itermax, refit = TRUE)
    fu <- newformula(b, mstop = itermax)
    b <- sdr(formula = fu, family = fam, data = d, maxit = maxit, 
             coef_start = c, batch_ids = batch_ids, updating = "bs",
             init = FALSE, CF = FALSE)
    time <- c(proc.time() - ptm)[3]
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    nu <- mean((nd$eta.nu - p$nu)^2)
    
    # transformation from linear predictor scale to dist parameter scale
    par.trafo <- function(list){
      list$mu <- exp(list$mu)
      list$sigma <- exp(list$sigma)
      list$nu <- 1/(1+exp(-list$nu))
      return(list)
    }
    p <- par.trafo(p)
    
    k = max(350, nd$y)
    mat <- fam$d(0, p)
    for(i in 1:k){
      pp <- fam$d(i, p)
      mat <- rbind(mat, pp)
      cat("\r");cat("i =",i, "/", k)
    }
    mat <- t(mat)
    colnames(mat) <- 0:k
    mean1 <- rowSums(mat %*% diag(c(0:k)))
    
    mse1 <- mean((nd$y - mean1)^2)
    crps1 <- try(rps(nd$y, mat),TRUE)
    if(!is.numeric(crps1)) crps1 <- NaN
    
    ll <- fam[["d"]]
    l <- try(sum(ll(nd$y, p, log = T)), T)
    if(!is.numeric(l)) l <- NaN
    
    
    # combine results
    m <- rbind(
      "sdr4_bw" = c("mu" = mu, "sigma" = sigma, "nu" = nu, "mse" = mse1, "crps" = crps1, "logLiK" = l, "mstop" = maxit, time)
    )
    
    m <- as.data.frame(m)
    rownames(m) <- NULL
    
    m$method <- "sdr4_bw"
    m$nobs <- nobs
    m$nnoise <- nnoise
    m$rho <- rho
    m$rep <- seed
    
    
    # Extract coefs and calc number of true and false positives
    c <- coef(b)
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      c.nu <- c[["nu"]]
      cc.mu <- cc.sigma <- cc.nu <-  rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      cc.nu[1] <- c.nu["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
        cc.nu[i+1] <- ifelse(paste0("x",i) %in% names(c.nu) ,c.nu[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma, "nu" = cc.nu))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    nu.names <- c("nu.intercept", paste0("nu.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    names(c[["nu"]]) <- nu.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    coefs.nu <- rbind.data.frame(c[["nu"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- rownames(coefs.nu) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    colnames(coefs.nu) <- nu.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,c(2,4,6,7)] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,c(3,5,6)] != 0)
    resi$tp.nu <- rowSums(coefs.nu[,c(4,5,6)] != 0)
    
    resi$fp.mu <- rowSums(coefs.mu[,-c(1,2,4,6,7)] != 0)
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,3,5,6)] != 0)
    resi$fp.nu <- rowSums(coefs.nu[,-c(1,4,5,6)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2,4,6,7)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(3,5,6)] -  t(matrix(rep(sigmac, 1), nrow = 3)) ) ^2)
    resi$tp.nu.mse <- rowMeans((coefs.nu[,c(4,5,6)] - t(matrix(rep(nuc, 1), nrow = 3)) ) ^2)
    
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1,2,4,6,7)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,3,5,6)]^2)
    resi$fp.nu.mse <- rowMeans(coefs.nu[,-c(1,4,5,6)]^2)
    
    m <- cbind(m, coefs.mu[,c(2,4,6,7)], coefs.sigma[,c(3,5,6)],coefs.nu[,c(4,5,6)], resi)
    m1 <- rbind(m1, m)
  }
  
  return(m1)
}



