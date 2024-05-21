library("bamlss")
library("gamlss.dist")
library("gamboostLSS")
library("ASL")
library("ff")
library("stagewise")
source("~/SCRATCH/stagewise/prg_R/Var_deselection.R")

### Define data generating process
## Parameter for linear predictors in dgp
cm1 = 1; cm2 = 2; cm3 = 0.5; cm4 = -1
cs3 = 0.5; cs4 = 0.75; cs5 = -0.3; cs6 = -0.5
muc <- c(cm1,cm2, cm3, cm4)
sigmac <- c(cs3,cs4, cs5, cs6)


dgp <- function(nobs = 6000, nnoise = 0, rho = 0, seed = rnorm(1)) {
  
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
  d$eta.mu <-   cm1*d$x1 + cm2 * d$x2 + cm3 * d$x3 + cm4*d$x4 
  d$eta.sigma <-  cs3 * d$x3 + cs4 * d$x4 + cs5 * d$x5 + cs6 * d$x6 
  d$y <- rGA(nobs, mu = exp(d$eta.mu), sigma = exp(d$eta.sigma)) + 10^-8
  # + 10^-8 to avoid errors related to rounding down to zero, 
  # for example 10^-50 gets rounded to 0 and gives error as 0 can not come from GA
  # mu is mean of GA
  # sigma*mu is standard deviation of GA
  
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
    
    fam <- GA
    
    ## (1) True model
    ftrue <- list(y ~ x1 + x2 + x3 + x4,
                  sigma ~ x3 + x4 + x5 + x6)
    
    ptm <- proc.time()
    b <- bamlss(ftrue, data = d, family = fam, sampler = FALSE)
    time <- c(proc.time() - ptm)[3]
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    
    p$mu <- exp(p$mu)
    p$sigma <- exp(p$sigma)
    
    crps_gamma <- function (y, mu, sigma) 
    {
      sd <- sigma * mu # sd of gamma dist
      scale <- sd^2 / mu
      shape <- mu^2 / sd^2
      
      p1 <- pgamma(y, shape, scale = scale)
      p2 <- pgamma(y, shape + 1, scale = scale)
      c <- y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1/beta(0.5, shape))
      return(mean(c))
    } 
    crps <- crps_gamma(nd$y, mu = p$mu, sigma = p$sigma)
    l <- sum(dGA(nd$y, mu = p$mu, sigma = p$sigma, log = T))
    mse <- mean((nd$y - p$mu)^2)
    # combine results
    m <- rbind(
      "likelihood" = c("mu" = mu, "sigma" = sigma, "mse" = mse, "crps" = crps, "logLiK" = l, "mstop" = maxit, time)
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
    c.mu <- c[nmu]
    c.sigma <- c[nsigma]
    nmu <- gsub("mu.p.","",nmu)
    nsigma <- gsub("sigma.p.","",nsigma)
    colnames(c.mu) <- nmu
    colnames(c.sigma) <- nsigma
    colnames(c.mu)[1] <- "(Intercept)"
    colnames(c.sigma)[1] <- "(Intercept)"
    c <- list("mu" = unlist(c.mu), "sigma" = unlist(c.sigma))
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      cc.mu <- cc.sigma <- rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,2:5] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,4:7] != 0)
    resi$fp.mu <- rowSums(coefs.mu[,-c(1:5)] != 0)
    
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,4:7)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2:5)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(4:7)] -  t(matrix(rep(sigmac, 1), nrow = 4)) ) ^2)
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1:5)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,4:7)]^2)
    m <- cbind(m, coefs.mu[,2:5], coefs.sigma[,4:7], resi)
    m1 <- m
  }
  
  # Thresdesc full batch
  if(nobs < 100000){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .)
    )
    
    fam <- GA
    
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
    
    p$mu <- exp(p$mu)
    p$sigma <- exp(p$sigma)
    
    crps_gamma <- function (y, mu, sigma) 
    {
      sd <- sigma * mu # sd of gamma dist
      scale <- sd^2 / mu
      shape <- mu^2 / sd^2
      
      p1 <- pgamma(y, shape, scale = scale)
      p2 <- pgamma(y, shape + 1, scale = scale)
      c <- y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1/beta(0.5, shape))
      return(mean(c))
    } 
    crps <- crps_gamma(nd$y, mu = p$mu, sigma = p$sigma)
    l <- sum(dGA(nd$y, mu = p$mu, sigma = p$sigma, log = T))
    mse <- mean((nd$y - p$mu)^2)
    
    # combine results
    m <- rbind(
      "thresdesc" = c("mu" = mu, "sigma" = sigma, "mse" = mse, "crps" = crps, "logLiK" = l, "mstop" = maxit_refit, time)
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
      cc.mu <- cc.sigma <- rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,2:5] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,4:7] != 0)
    resi$fp.mu <- rowSums(coefs.mu[,-c(1:5)] != 0)
    
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,4:7)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2:5)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(4:7)] -  t(matrix(rep(sigmac, 1), nrow = 4)) ) ^2)
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1:5)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,4:7)]^2)
    
    m <- cbind(m, coefs.mu[,2:5], coefs.sigma[,4:7], resi)
    m1 <- rbind(m1, m)
  }
  
  # Thresdesc batchwise
  if(TRUE){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .)
    )
    
    fam <- GA
    
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
    
    p$mu <- exp(p$mu)
    p$sigma <- exp(p$sigma)
    
    crps_gamma <- function (y, mu, sigma) 
    {
      sd <- sigma * mu # sd of gamma dist
      scale <- sd^2 / mu
      shape <- mu^2 / sd^2
      
      p1 <- pgamma(y, shape, scale = scale)
      p2 <- pgamma(y, shape + 1, scale = scale)
      c <- y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1/beta(0.5, shape))
      return(mean(c))
    } 
    crps <- crps_gamma(nd$y, mu = p$mu, sigma = p$sigma)
    l <- sum(dGA(nd$y, mu = p$mu, sigma = p$sigma, log = T))
    mse <- mean((nd$y - p$mu)^2)
    
    # combine results
    m <- rbind(
      "thresdesc_bw" = c("mu" = mu, "sigma" = sigma, "mse" = mse, "crps" = crps, "logLiK" = l, "mstop" = maxit_refit, time)
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
      cc.mu <- cc.sigma <- rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,2:5] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,4:7] != 0)
    resi$fp.mu <- rowSums(coefs.mu[,-c(1:5)] != 0)
    
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,4:7)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2:5)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(4:7)] -  t(matrix(rep(sigmac, 1), nrow = 4)) ) ^2)
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1:5)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,4:7)]^2)
    
    m <- cbind(m, coefs.mu[,2:5], coefs.sigma[,4:7], resi)
    m1 <- rbind(m1,m)
  }
  
  # Gradient boosting with 10 fold cross validation (gamboostLSS)
  if(nobs < 100000){
    fam <- GA
    
    # Gradient boosting with 10 fold cross validation
    f <- grep("x", names(d), value = TRUE)
    f <- paste("y ~", paste(f, collapse = "+"))
    f <- as.formula(f)
    
    ptm <- proc.time()
    maxitlss <- ifelse(nobs >= 5000, 3000, 2000)
    b <- glmboostLSS(f,  families = as.families(GA),
                     data = d, method = "noncyclic",
                     control = boost_control(mstop = maxitlss, nu = 0.1, trace = TRUE))
    
    set.seed(seed)
    folds <- mboost::cv(rep(1, nrow(d)), B = 10, type = "kfold")
    cvr <- cvrisk(b, folds = folds, papply = lapply)
    mstop(b) <- mstop(cvr)
    # b is input for variable deselection
    
    mstop <- mstop(b)[1]
    time <- c(proc.time() - ptm)[3]
    
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    
    p$mu <- exp(p$mu)
    p$sigma <- exp(p$sigma)
    
    crps_gamma <- function (y, mu, sigma) 
    {
      sd <- sigma * mu # sd of gamma dist
      scale <- sd^2 / mu
      shape <- mu^2 / sd^2
      
      p1 <- pgamma(y, shape, scale = scale)
      p2 <- pgamma(y, shape + 1, scale = scale)
      c <- y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1/beta(0.5, shape))
      return(mean(c))
    } 
    crps <- crps_gamma(nd$y, mu = p$mu, sigma = p$sigma)
    l <- sum(dGA(nd$y, mu = p$mu, sigma = p$sigma, log = T))
    mse <- mean((nd$y - p$mu)^2)
    
    # combine results
    m <- rbind(
      "boosting" = c("mu" = mu, "sigma" = sigma, "mse" = mse, "crps" = crps, "logLiK" = l, "mstop" = mstop, time)
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
      cc.mu <- cc.sigma <- rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,2:5] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,4:7] != 0)
    resi$fp.mu <- rowSums(coefs.mu[,-c(1:5)] != 0)
    
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,4:7)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2:5)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(4:7)] -  t(matrix(rep(sigmac, 1), nrow = 4)) ) ^2)
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1:5)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,4:7)]^2)
    
    m <- cbind(m, coefs.mu[,2:5], coefs.sigma[,4:7], resi)
    m1 <- rbind(m1, m)
  }
  # Var deselection, gradient boosting model is input
  # variables which contribute less than 1% of the total risk reduction are deselected, afterwards refitting until mstop
  if(nobs < 100000){
    fam <- GA
    
    f <- grep("x", names(d), value = TRUE)
    f <- paste("y ~", paste(f, collapse = "+"))
    f <- as.formula(f)
    
    timegradboost <- time
    ptm <- proc.time()
    b <- DeselectBoost(b, fam = as.families(GA), 
                       tau = 0.01)  
    
    time <- c(proc.time() - ptm)[3]
    time <- time + timegradboost
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- list()
    p$mu <- predict(b$mu,newdata=nd)
    p$sigma <- predict(b$sigma,newdata=nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    
    p$mu <- exp(p$mu)
    p$sigma <- exp(p$sigma)
    
    crps_gamma <- function (y, mu, sigma) 
    {
      sd <- sigma * mu # sd of gamma dist
      scale <- sd^2 / mu
      shape <- mu^2 / sd^2
      
      p1 <- pgamma(y, shape, scale = scale)
      p2 <- pgamma(y, shape + 1, scale = scale)
      c <- y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1/beta(0.5, shape))
      return(mean(c))
    } 
    crps <- crps_gamma(nd$y, mu = p$mu, sigma = p$sigma)
    l <- sum(dGA(nd$y, mu = p$mu, sigma = p$sigma, log = T))
    mse <- mean((nd$y - p$mu)^2)
    
    # combine results
    m <- rbind(
      "var_des" = c("mu" = mu, "sigma" = sigma, "mse" = mse, "crps" = crps, "logLiK" = l, "mstop" = NA, time)
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
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      cc.mu <- cc.sigma <- rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,2:5] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,4:7] != 0)
    resi$fp.mu <- rowSums(coefs.mu[,-c(1:5)] != 0)
    
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,4:7)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2:5)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(4:7)] -  t(matrix(rep(sigmac, 1), nrow = 4)) ) ^2)
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1:5)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,4:7)]^2)
    
    m <- cbind(m, coefs.mu[,2:5], coefs.sigma[,4:7], resi)
    m1 <- rbind(m1, m)
  }
  
  # bamlss stability selection
  # Estimate selection frequencies of variables on B = 100 resampled data sets
  if(nobs < 100000){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .)
    )
    
    fam <- GA
    
    # bamlss stability selection
    q = max(round(sqrt(6+nnoise)*2), 8)
    ptm <- proc.time()
    sel <- bamlss::stabsel(formula = f, data = d, q = q, B = 100, 
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
    nf <- stabf(sel, thr = 0.8, B = 100, par = c("mu", "sigma"))
    b <- bamlss(nf, data = d, family = fam)
    time <- c(proc.time() - ptm)[3]
    
    ######### Comparison #####################
    
    ## Predict linear predictors and calculate different metrics on validation data
    p <- predict(b, newdata = nd)
    
    # metrics
    mu <- mean((nd$eta.mu - p$mu)^2)
    sigma <- mean((nd$eta.sigma - p$sigma)^2)
    
    p$mu <- exp(p$mu)
    p$sigma <- exp(p$sigma)
    
    crps_gamma <- function (y, mu, sigma) 
    {
      sd <- sigma * mu # sd of gamma dist
      scale <- sd^2 / mu
      shape <- mu^2 / sd^2
      
      p1 <- pgamma(y, shape, scale = scale)
      p2 <- pgamma(y, shape + 1, scale = scale)
      c <- y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1/beta(0.5, shape))
      return(mean(c))
    } 
    crps <- crps_gamma(nd$y, mu = p$mu, sigma = p$sigma)
    l <- sum(dGA(nd$y, mu = p$mu, sigma = p$sigma, log = T))
    mse <- mean((nd$y - p$mu)^2)
    
    # combine results
    m <- rbind(
      "bamlss_stab" = c("mu" = mu, "sigma" = sigma, "mse" = mse, "crps" = crps, "logLiK" = l, "mstop" = NA, time)
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
    c.mu <- c[nmu]
    c.sigma <- c[nsigma]
    nmu <- gsub("mu.p.","",nmu)
    nsigma <- gsub("sigma.p.","",nsigma)
    names(c.mu) <- nmu
    names(c.sigma) <- nsigma
    names(c.mu)[1] <- "(Intercept)"
    names(c.sigma)[1] <- "(Intercept)"
    c.mu <- c.mu[-length(c.mu)] 
    c.sigma <- c.sigma[-length(c.sigma)] 
    c <- list("mu" = c.mu,"sigma" = c.sigma)
    
    coef.fun <- function(c){
      c.mu <- c[["mu"]]
      c.sigma <- c[["sigma"]]
      cc.mu <- cc.sigma <- rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,2:5] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,4:7] != 0)
    resi$fp.mu <- rowSums(coefs.mu[,-c(1:5)] != 0)
    
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,4:7)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2:5)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(4:7)] -  t(matrix(rep(sigmac, 1), nrow = 4)) ) ^2)
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1:5)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,4:7)]^2)
    
    m <- cbind(m, coefs.mu[,2:5], coefs.sigma[,4:7], resi)
    m1 <- rbind(m1, m)
  }
  
  
  
  # Full batch variant
  # SDR1 noncyclic updating, cap = 0 and early stopping selection via bic
  if(nobs < 100000){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .)
    )
    
    fam <- GA
    
    
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
    
    p$mu <- exp(p$mu)
    p$sigma <- exp(p$sigma)
    
    crps_gamma <- function (y, mu, sigma) 
    {
      sd <- sigma * mu # sd of gamma dist
      scale <- sd^2 / mu
      shape <- mu^2 / sd^2
      
      p1 <- pgamma(y, shape, scale = scale)
      p2 <- pgamma(y, shape + 1, scale = scale)
      c <- y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1/beta(0.5, shape))
      return(mean(c))
    } 
    crps <- crps_gamma(nd$y, mu = p$mu, sigma = p$sigma)
    l <- sum(dGA(nd$y, mu = p$mu, sigma = p$sigma, log = T))
    mse <- mean((nd$y - p$mu)^2)
    
    # combine results
    m <- rbind(
      "sdr1" = c("mu" = mu, "sigma" = sigma, "mse" = mse, "crps" = crps, "logLiK" = l, "mstop" = maxit, time)
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
      cc.mu <- cc.sigma <- rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,2:5] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,4:7] != 0)
    resi$fp.mu <- rowSums(coefs.mu[,-c(1:5)] != 0)
    
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,4:7)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2:5)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(4:7)] -  t(matrix(rep(sigmac, 1), nrow = 4)) ) ^2)
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1:5)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,4:7)]^2)
    
    m <- cbind(m, coefs.mu[,2:5], coefs.sigma[,4:7], resi)
    m1 <- rbind(m1, m)
  }
  
  # SDR2 best subset updating, CF = FALSE and early stopping selection via bic
  if(nobs < 100000){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .)
    )
    
    fam <- GA
    
    # SDR2 bs updating, cap = 0 and early stopping selection via bic
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
    
    p$mu <- exp(p$mu)
    p$sigma <- exp(p$sigma)
    
    crps_gamma <- function (y, mu, sigma) 
    {
      sd <- sigma * mu # sd of gamma dist
      scale <- sd^2 / mu
      shape <- mu^2 / sd^2
      
      p1 <- pgamma(y, shape, scale = scale)
      p2 <- pgamma(y, shape + 1, scale = scale)
      c <- y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1/beta(0.5, shape))
      return(mean(c))
    } 
    crps <- crps_gamma(nd$y, mu = p$mu, sigma = p$sigma)
    l <- sum(dGA(nd$y, mu = p$mu, sigma = p$sigma, log = T))
    mse <- mean((nd$y - p$mu)^2)
    
    # combine results
    m <- rbind(
      "sdr2" = c("mu" = mu, "sigma" = sigma, "mse" = mse, "crps" = crps, "logLiK" = l, "mstop" = maxit, time)
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
      cc.mu <- cc.sigma <- rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,2:5] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,4:7] != 0)
    resi$fp.mu <- rowSums(coefs.mu[,-c(1:5)] != 0)
    
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,4:7)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2:5)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(4:7)] -  t(matrix(rep(sigmac, 1), nrow = 4)) ) ^2)
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1:5)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,4:7)]^2)
    
    m <- cbind(m, coefs.mu[,2:5], coefs.sigma[,4:7], resi)
    m1 <- rbind(m1, m)
  }
  
  # SDR3 noncyclic updating and CF = TRUE 
  if(nobs < 100000){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .)
    )
    
    fam <- GA
    
    # SDR3 noncyclic updating and CF = TRUE
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
    
    p$mu <- exp(p$mu)
    p$sigma <- exp(p$sigma)
    
    crps_gamma <- function (y, mu, sigma) 
    {
      sd <- sigma * mu # sd of gamma dist
      scale <- sd^2 / mu
      shape <- mu^2 / sd^2
      
      p1 <- pgamma(y, shape, scale = scale)
      p2 <- pgamma(y, shape + 1, scale = scale)
      c <- y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1/beta(0.5, shape))
      return(mean(c))
    } 
    crps <- crps_gamma(nd$y, mu = p$mu, sigma = p$sigma)
    l <- sum(dGA(nd$y, mu = p$mu, sigma = p$sigma, log = T))
    mse <- mean((nd$y - p$mu)^2)
    
    # combine results
    m <- rbind(
      "sdr3" = c("mu" = mu, "sigma" = sigma, "mse" = mse, "crps" = crps, "logLiK" = l, "mstop" = maxit, time)
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
      cc.mu <- cc.sigma <- rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,2:5] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,4:7] != 0)
    resi$fp.mu <- rowSums(coefs.mu[,-c(1:5)] != 0)
    
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,4:7)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2:5)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(4:7)] -  t(matrix(rep(sigmac, 1), nrow = 4)) ) ^2)
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1:5)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,4:7)]^2)
    
    m <- cbind(m, coefs.mu[,2:5], coefs.sigma[,4:7], resi)
    m1 <- rbind(m1, m)
  }
  
  # SDR4 best subset updating and CF = TRUE
  if(nobs < 100000){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .)
    )
    
    fam <- GA
    
    # SDR4 bs updating and CF = TRUE
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
    
    p$mu <- exp(p$mu)
    p$sigma <- exp(p$sigma)
    
    crps_gamma <- function (y, mu, sigma) 
    {
      sd <- sigma * mu # sd of gamma dist
      scale <- sd^2 / mu
      shape <- mu^2 / sd^2
      
      p1 <- pgamma(y, shape, scale = scale)
      p2 <- pgamma(y, shape + 1, scale = scale)
      c <- y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1/beta(0.5, shape))
      return(mean(c))
    } 
    crps <- crps_gamma(nd$y, mu = p$mu, sigma = p$sigma)
    l <- sum(dGA(nd$y, mu = p$mu, sigma = p$sigma, log = T))
    mse <- mean((nd$y - p$mu)^2)
    
    # combine results
    m <- rbind(
      "sdr4" = c("mu" = mu, "sigma" = sigma, "mse" = mse, "crps" = crps, "logLiK" = l, "mstop" = maxit, time)
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
      cc.mu <- cc.sigma <- rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,2:5] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,4:7] != 0)
    resi$fp.mu <- rowSums(coefs.mu[,-c(1:5)] != 0)
    
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,4:7)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2:5)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(4:7)] -  t(matrix(rep(sigmac, 1), nrow = 4)) ) ^2)
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1:5)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,4:7)]^2)
    
    m <- cbind(m, coefs.mu[,2:5], coefs.sigma[,4:7], resi)
    m1 <- rbind(m1, m)
  }
  
  # Batchwise variant
  # SDR1 bw noncyclic updating, CF = FALSE and early stopping selection via bic
  if(TRUE){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .)
    )
    
    fam <- GA
    
    #if(nobs > 10000) d <- as.ffdf(d)
    
    # SDR1 bw noncyclic updating, cap = 0 and early stopping selection via bic
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
    
    p$mu <- exp(p$mu)
    p$sigma <- exp(p$sigma)
    
    crps_gamma <- function (y, mu, sigma) 
    {
      sd <- sigma * mu # sd of gamma dist
      scale <- sd^2 / mu
      shape <- mu^2 / sd^2
      
      p1 <- pgamma(y, shape, scale = scale)
      p2 <- pgamma(y, shape + 1, scale = scale)
      c <- y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1/beta(0.5, shape))
      return(mean(c))
    } 
    crps <- crps_gamma(nd$y, mu = p$mu, sigma = p$sigma)
    l <- sum(dGA(nd$y, mu = p$mu, sigma = p$sigma, log = T))
    mse <- mean((nd$y - p$mu)^2)
    
    # combine results
    m <- rbind(
      "sdr1_bw" = c("mu" = mu, "sigma" = sigma, "mse" = mse, "crps" = crps, "logLiK" = l, "mstop" = maxit, time)
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
      cc.mu <- cc.sigma <- rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,2:5] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,4:7] != 0)
    resi$fp.mu <- rowSums(coefs.mu[,-c(1:5)] != 0)
    
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,4:7)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2:5)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(4:7)] -  t(matrix(rep(sigmac, 1), nrow = 4)) ) ^2)
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1:5)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,4:7)]^2)
    
    m <- cbind(m, coefs.mu[,2:5], coefs.sigma[,4:7], resi)
    m1 <- rbind(m1, m)
  }
  
  # SDR2 bw best subset updating, CF = FALSE and early stopping selection via bic
  if(TRUE){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .)
    )
    
    fam <- GA
    
    #if(nobs > 10000) d <- as.ffdf(d)
    
    # SDR2 bw best subset updating, cap = 0 and early stopping selection via bic
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
    
    p$mu <- exp(p$mu)
    p$sigma <- exp(p$sigma)
    
    crps_gamma <- function (y, mu, sigma) 
    {
      sd <- sigma * mu # sd of gamma dist
      scale <- sd^2 / mu
      shape <- mu^2 / sd^2
      
      p1 <- pgamma(y, shape, scale = scale)
      p2 <- pgamma(y, shape + 1, scale = scale)
      c <- y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1/beta(0.5, shape))
      return(mean(c))
    } 
    crps <- crps_gamma(nd$y, mu = p$mu, sigma = p$sigma)
    l <- sum(dGA(nd$y, mu = p$mu, sigma = p$sigma, log = T))
    mse <- mean((nd$y - p$mu)^2)
    
    # combine results
    m <- rbind(
      "sdr2_bw" = c("mu" = mu, "sigma" = sigma, "mse" = mse, "crps" = crps, "logLiK" = l, "mstop" = maxit, time)
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
      cc.mu <- cc.sigma <- rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,2:5] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,4:7] != 0)
    resi$fp.mu <- rowSums(coefs.mu[,-c(1:5)] != 0)
    
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,4:7)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2:5)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(4:7)] -  t(matrix(rep(sigmac, 1), nrow = 4)) ) ^2)
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1:5)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,4:7)]^2)
    
    m <- cbind(m, coefs.mu[,2:5], coefs.sigma[,4:7], resi)
    m1 <- rbind(m1, m)
  }
  
  # SDR3 bw noncyclic updating and CF = TRUE 
  if(TRUE){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .)
    )
    
    fam <- GA
    
    #if(nobs > 10000) d <- as.ffdf(d)
    
    # SDR3 bw noncyclic updating and CF = TRUE
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
    
    p$mu <- exp(p$mu)
    p$sigma <- exp(p$sigma)
    
    crps_gamma <- function (y, mu, sigma) 
    {
      sd <- sigma * mu # sd of gamma dist
      scale <- sd^2 / mu
      shape <- mu^2 / sd^2
      
      p1 <- pgamma(y, shape, scale = scale)
      p2 <- pgamma(y, shape + 1, scale = scale)
      c <- y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1/beta(0.5, shape))
      return(mean(c))
    } 
    crps <- crps_gamma(nd$y, mu = p$mu, sigma = p$sigma)
    l <- sum(dGA(nd$y, mu = p$mu, sigma = p$sigma, log = T))
    mse <- mean((nd$y - p$mu)^2)
    
    # combine results
    m <- rbind(
      "sdr3_bw" = c("mu" = mu, "sigma" = sigma, "mse" = mse, "crps" = crps, "logLiK" = l, "mstop" = maxit, time)
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
      cc.mu <- cc.sigma <- rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,2:5] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,4:7] != 0)
    resi$fp.mu <- rowSums(coefs.mu[,-c(1:5)] != 0)
    
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,4:7)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2:5)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(4:7)] -  t(matrix(rep(sigmac, 1), nrow = 4)) ) ^2)
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1:5)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,4:7)]^2)
    
    m <- cbind(m, coefs.mu[,2:5], coefs.sigma[,4:7], resi)
    m1 <- rbind(m1, m)
  }
  
  # SDR4 bw best subset updating and CF = TRUE
  if(TRUE){
    maxit = NULL
    
    f <- as.formula(paste("y ~ ", paste0("x", 1:(6+nnoise), collapse = "+")))
    f <- list(
      f,
      update(f, sigma ~ .)
    )
    
    fam <- GA
    
    #if(nobs > 10000) d <- as.ffdf(d)
    
    # SDR4 bw best subset updating and CF = TRUE
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
    
    p$mu <- exp(p$mu)
    p$sigma <- exp(p$sigma)
    
    crps_gamma <- function (y, mu, sigma) 
    {
      sd <- sigma * mu # sd of gamma dist
      scale <- sd^2 / mu
      shape <- mu^2 / sd^2
      
      p1 <- pgamma(y, shape, scale = scale)
      p2 <- pgamma(y, shape + 1, scale = scale)
      c <- y * (2 * p1 - 1) - scale * (shape * (2 * p2 - 1) + 1/beta(0.5, shape))
      return(mean(c))
    } 
    crps <- crps_gamma(nd$y, mu = p$mu, sigma = p$sigma)
    l <- sum(dGA(nd$y, mu = p$mu, sigma = p$sigma, log = T))
    mse <- mean((nd$y - p$mu)^2)
    
    # combine results
    m <- rbind(
      "sdr4_bw" = c("mu" = mu, "sigma" = sigma, "mse" = mse, "crps" = crps, "logLiK" = l, "mstop" = maxit, time)
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
      cc.mu <- cc.sigma <- rep(0, 1+6+nnoise)
      cc.mu[1] <- c.mu["(Intercept)"]
      cc.sigma[1] <- c.sigma["(Intercept)"]
      
      for(i in 1:(length(cc.mu)-1)) {
        cc.mu[i+1] <- ifelse(paste0("x",i) %in% names(c.mu) ,c.mu[paste0("x",i)], 0 )
        cc.sigma[i+1] <- ifelse(paste0("x",i) %in% names(c.sigma) ,c.sigma[paste0("x",i)], 0 )
      }
      return(list("mu" = cc.mu, "sigma" = cc.sigma))
    }
    
    c <- coef.fun(c)
    
    mu.names <- c("mu.intercept", paste0("mu.x",1:(6+nnoise)))
    sigma.names <- c("sigma.intercept", paste0("sigma.x",1:(6+nnoise)))
    names(c[["mu"]]) <- mu.names
    names(c[["sigma"]]) <- sigma.names
    coefs.mu <- rbind.data.frame(c[["mu"]])
    coefs.sigma <- rbind.data.frame(c[["sigma"]])
    rownames(coefs.mu) <- rownames(coefs.sigma) <- NULL
    colnames(coefs.mu) <- mu.names
    colnames(coefs.sigma) <- sigma.names
    
    resi <- data.frame("tp.mu" = rowSums(coefs.mu[,2:5] != 0))
    resi$tp.sigma <- rowSums(coefs.sigma[,4:7] != 0)
    resi$fp.mu <- rowSums(coefs.mu[,-c(1:5)] != 0)
    
    resi$fp.sigma <- rowSums(coefs.sigma[,-c(1,4:7)] != 0)
    
    resi$tp.mu.mse <- rowMeans((coefs.mu[,c(2:5)] - t(matrix(rep(muc, 1), nrow = 4)) ) ^2)
    resi$tp.sigma.mse <- rowMeans((coefs.sigma[,c(4:7)] -  t(matrix(rep(sigmac, 1), nrow = 4)) ) ^2)
    resi$fp.mu.mse <- rowMeans(coefs.mu[,-c(1:5)]^2)
    resi$fp.sigma.mse <- rowMeans(coefs.sigma[,-c(1,4:7)]^2)
    
    m <- cbind(m, coefs.mu[,2:5], coefs.sigma[,4:7], resi)
    m1 <- rbind(m1, m)
  }
  
  return(m1)
}



