# to do
# build package
# R CMD REMOVE stagewise
# R CMD build stagewise
# R CMD INSTALL stagewise_ .tar.gz``


powerset <- function(set) {
  n <- length(set)
  masks <- 2^(1:n - 1)
  lapply(1:2^n - 1, function(u) set[bitwAnd(u, masks) != 0])
}

## Function to transform gamlss.family objects.
tF <- function(x, ...)
{
  if(is.function(x)) x <- x()
  if(!inherits(x, "gamlss.family")) stop('only "gamlss.family" objects can be transformed!')

  args <- list(...)
  bd <- if(is.null(args$bd)) 1 else args$bd
  args$bd <- NULL
  pr <- args$range
  check_range <- function(par) { return(par) }
  if(!is.null(pr)) {
    if(is.list(pr) | is.data.frame(pr)) {
      check_range <- function(par) {
        for(j in names(par)) {
          if(!is.null(pr[[j]])) {
            if(is.numeric(pr[[j]])) {
              par[[j]][par[[j]] < min(pr[[j]])] <- min(pr[[j]])
              par[[j]][par[[j]] > max(pr[[j]])] <- max(pr[[j]])
            }
          }
        }
        par
      }
    }
  }
  nx <- names(x$parameters)
  score <- hess <- initialize <- list()

  make_call <- function(fun) {
    fn <- deparse(substitute(fun), backtick = TRUE, width.cutoff = 500)
    nf <- names(formals(fun))
    if(length(nf) < 1) {
      call <- paste(fn, "()", sep = "")
    } else {
      call <- paste(fn, "(", if("y" %in% nf) "y," else "", sep = "")
      np <- nx[nx %in% nf]
      call <- paste(call, paste(np, '=', 'par$', np, sep = '', collapse = ','), sep = "")
      if("bd" %in% nf) {
        call <- paste(call, ",bd=", bd, sep = "")
      }
    }
    call <- parse(text = paste(call, ")", sep = ""))
    return(call)
  }

  if("mu" %in% nx) {
    mu.link <- make.link2(x$mu.link)
    mu.cs <- make_call(x$dldm)
    mu.hs <- make_call(x$d2ldm2)
    score$mu  <- function(y, par, ...) {
      par <- check_range(par)
      res <- eval(mu.cs) * mu.link$mu.eta(mu.link$linkfun(par$mu))
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    hess$mu <- function(y, par, ...) {
      par <- check_range(par)
      score <- eval(mu.cs)
      hess <- -1 * eval(mu.hs)
      eta <- mu.link$linkfun(par$mu)
      res <- drop(score * mu.link$mu.eta2(eta) + hess * mu.link$mu.eta(eta)^2)
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    if(!is.null(x$mu.initial)) {
      initialize$mu <- function(y, ...) {
        if(!is.null(attr(y, "contrasts"))) {
          if(!is.null(dim(y)))
            y <- y[, ncol(y)]
        }
        if(!is.null(bd))
          bd <- rep(bd, length.out = if(!is.null(dim(y))) nrow(y) else length(y))
        res <- eval(x$mu.initial)
        if(!is.null(dim(res)))
          res <- res[, 1]
        res
      }
    }
  }

  if("sigma" %in% nx) {
    sigma.link <- make.link2(x$sigma.link)
    sigma.cs <- make_call(x$dldd)
    sigma.hs <- make_call(x$d2ldd2)
    score$sigma  <- function(y, par, ...) {
      par <- check_range(par)
      res <- eval(sigma.cs) * sigma.link$mu.eta(sigma.link$linkfun(par$sigma))
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    hess$sigma <- function(y, par, ...) {
      par <- check_range(par)
      score <- eval(sigma.cs)
      hess <- -1 * eval(sigma.hs)
      eta <- sigma.link$linkfun(par$sigma)
      res <- drop(score * sigma.link$mu.eta2(eta) + hess * sigma.link$mu.eta(eta)^2)
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    if(!is.null(x$sigma.initial)) {
      initialize$sigma <- function(y, ...) {
        if(!is.null(bd))
          bd <- rep(bd, length.out = if(!is.null(dim(y))) nrow(y) else length(y))
        res <- eval(x$sigma.initial)
        if(!is.null(dim(res)))
          res <- res[, 1]
        res
      }
    }
  }

  if("nu" %in% nx) {
    nu.link <- make.link2(x$nu.link)
    nu.cs <- make_call(x$dldv)
    nu.hs <- make_call(x$d2ldv2)
    score$nu  <- function(y, par, ...) {
      par <- check_range(par)
      res <- eval(nu.cs) * nu.link$mu.eta(nu.link$linkfun(par$nu))
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    hess$nu <- function(y, par, ...) {
      par <- check_range(par)
      score <- eval(nu.cs)
      hess <- -1 * eval(nu.hs)
      eta <- nu.link$linkfun(par$nu)
      res <- drop(score * nu.link$mu.eta2(eta) + hess * nu.link$mu.eta(eta)^2)
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    if(!is.null(x$nu.initial)) {
      initialize$nu <- function(y, ...) {
        if(!is.null(bd))
          bd <- rep(bd, length.out = if(!is.null(dim(y))) nrow(y) else length(y))
        res <- eval(x$nu.initial)
        if(!is.null(dim(res)))
          res <- res[, 1]
        res
      }
    }
  }

  if("tau" %in% nx) {
    tau.link <- make.link2(x$tau.link)
    tau.cs <- make_call(x$dldt)
    tau.hs <- make_call(x$d2ldt2)
    score$tau  <- function(y, par, ...) {
      par <- check_range(par)
      res <- eval(tau.cs) * tau.link$mu.eta(tau.link$linkfun(par$tau))
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    hess$tau <- function(y, par, ...) {
      par <- check_range(par)
      score <- eval(tau.cs)
      hess <- -1 * eval(tau.hs)
      eta <- tau.link$linkfun(par$tau)
      res <- drop(score * tau.link$mu.eta2(eta) + hess * tau.link$mu.eta(eta)^2)
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    if(!is.null(x$tau.initial)) {
      initialize$tau <- function(y, ...) {
        if(!is.null(bd))
          bd <- rep(bd, length.out = if(!is.null(dim(y))) nrow(y) else length(y))
        res <- eval(x$tau.initial)
        if(!is.null(dim(res)))
          res <- res[, 1]
        res
      }
    }
  }

  dfun <- get(paste("d", x$family[1], sep = ""))
  pfun <- try(get(paste("p", x$family[1], sep = "")), silent = TRUE)
  qfun <- try(get(paste("q", x$family[1], sep = "")), silent = TRUE)
  rfun <- try(get(paste("r", x$family[1], sep = "")), silent = TRUE)

  nf <- names(formals(dfun))
  bdc <- "bd" %in% nf

  dc <- parse(text = paste('dfun(y,', paste(paste(nx, 'par$', sep = "="),
    nx, sep = '', collapse = ','), ',log=log,...',
    if(bdc) paste0(",bd=", bd) else NULL, ")", sep = ""))
  pc <- parse(text = paste('pfun(q,', paste(paste(nx, 'par$', sep = "="),
    nx, sep = '', collapse = ','), ',log=log,...',
    if(bdc) paste0(",bd=", bd) else NULL, ")", sep = ""))
  qc <- parse(text = paste('qfun(p,', paste(paste(nx, 'par$', sep = "="),
    nx, sep = '', collapse = ','), ',log=log,...',
    if(bdc) paste0(",bd=", bd) else NULL, ")", sep = ""))
  rc <- parse(text = paste('rfun(n,', paste(paste(nx, 'par$', sep = "="),
    nx, sep = '', collapse = ','), ',...',
    if(bdc) paste0(",bd=", bd) else NULL, ")", sep = ""))

  rval <- list(
    "family" = x$family[1],
    "names" = nx,
    "links" = unlist(x[paste(nx, "link", sep = ".")]),
    "score" = score,
    "hess" = hess,
    "d" = function(y, par, log = FALSE, ...) {
       par <- check_range(par)
       d <- try(eval(dc), silent = TRUE)
       if(inherits(d, "try-error"))
         d <- rep(NA, length(par[[1L]]))
       return(d)
    },
    "p" = if(!inherits(pfun, "try-error")) function(q, par, log = FALSE, ...) {
      par <- check_range(par)
      eval(pc)
    } else NULL,
    "q" = if(!inherits(qfun, "try-error")) function(p, par, log = FALSE, ...) {
      par <- check_range(par)
      eval(qc)
    } else NULL,
    "r" = if(!inherits(rfun, "try-error")) function(n, par, ...) {
      par <- check_range(par)
      eval(rc)
    } else NULL
  )
  names(rval$links) <- nx
  rval$valid.response <- x$y.valid
  rval$initialize <- initialize
  rval$type <- tolower(x$type)

  if(!is.null(x$mean)) {
    meanc <- make_call(x$mean)
    rval$mean  <- function(par, ...) {
      par <- check_range(par)
      res <- eval(meanc)
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
  } else {
    rval$mean <- function(par, ...) { par[[1L]] }
  }

  if(!is.null(x$variance)) {
    varc <- make_call(x$variance)
    rval$variance  <- function(par, ...) {
      par <- check_range(par)
      res <- eval(varc)
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
  } else {
    rval$variance <- function(par, ...) { par[[2L]] }
  }

  class(rval) <- "family.bamlss"
  rval
}

## Family completion.
complete.bamlss.family <- function(family)
{
  if(is.null(names(family$links)))
    names(family$links) <- family$names

  linkinv <- linkfun <- list()
  for(j in family$names) {
    link <- make.link2(family$links[j])
    linkinv[[j]] <- link$linkinv
    linkfun[[j]] <- link$linkfun
  }

  if(is.null(family$map2par)) {
    family$map2par <- function(eta) {
      for(j in family$names) {
        eta[[j]] <- linkinv[[j]](eta[[j]])
        eta[[j]][is.na(eta[[j]])] <- 0
        if(any(jj <- eta[[j]] == Inf))
          eta[[j]][jj] <- 10
        if(any(jj <- eta[[j]] == -Inf))
          eta[[j]][jj] <- -10
      }
      return(eta)
    }
  }

  if(is.null(family$mu)) {
    family$mu <- function(par) { make.link2(family$links[1])$linkinv(par[[1]]) }
  }

  if(is.null(family$loglik)) {
    if(!is.null(family$d)) {
      family$loglik <- function(y, par, ...) {
        logdens <- family$d(y, par, log = TRUE)
        if(any(i <- !is.finite(logdens))) {
          logdens[i] <- -100
        }
        return(sum(logdens, na.rm = TRUE))
      }
    }
  }

  err01 <- .Machine$double.eps^(1/3)
  err02 <- err01 * 2
  err11 <- .Machine$double.eps^(1/4)
  err12 <- err11 * 2

  if(is.null(family$score) & !is.null(family$d))
    family$score <- list()
  for(i in family$names) {
    if(is.null(family$score[[i]]) & !is.null(family$d)) {
      fun <- c(
        "function(y, par, ...) {",
        paste("  eta <- linkfun[['", i, "']](par[['", i, "']]);", sep = ""),
        paste("  par[['", i, "']] <- linkinv[['", i, "']](eta + err01);", sep = ""),
        "  d1 <- family$d(y, par, log = TRUE);",
        paste("  par[['", i, "']] <- linkinv[['", i, "']](eta - err01);", sep = ""),
        "  d2 <- family$d(y, par, log = TRUE);",
        "  return((d1 - d2) / err02)",
        "}"
      )
      family$score[[i]] <- eval(parse(text = paste(fun, collapse = "")))
      attr(family$score[[i]], "dnum") <- TRUE
    }
  }

  if(is.null(family$hess) & !is.null(family$d))
    family$hess <- list()
  for(i in family$names) {
    if(is.null(family$hess[[i]]) & !is.null(family$d)) {
      fun <- if(!is.null(attr(family$score[[i]], "dnum"))) {
        c(
          "function(y, par, ...) {",
          paste("  eta <- linkfun[['", i, "']](par[['", i, "']]);", sep = ""),
          paste("  par[['", i, "']] <- linkinv[['", i, "']](eta + err11);", sep = ""),
          paste("  d1 <- family$score[['", i, "']](y, par, ...);", sep = ""),
          paste("  par[['", i, "']] <- linkinv[['", i, "']](eta - err11);", sep = ""),
          paste("  d2 <- family$score[['", i, "']](y, par, ...);", sep = ""),
          "  return(-1 * (d1 - d2) / err12)",
          "}"
        )
      } else {
        c(
          "function(y, par, ...) {",
          paste("  eta <- linkfun[['", i, "']](par[['", i, "']]);", sep = ""),
          paste("  par[['", i, "']] <- linkinv[['", i, "']](eta + err01);", sep = ""),
          paste("  d1 <- family$score[['", i, "']](y, par, ...);", sep = ""),
          paste("  par[['", i, "']] <- linkinv[['", i, "']](eta - err01);", sep = ""),
          paste("  d2 <- family$score[['", i, "']](y, par, ...);", sep = ""),
          "  return(-1 * (d1 - d2) / err02)",
          "}"
        )
      }
      family$hess[[i]] <- eval(parse(text = paste(fun, collapse = "")))
    }
  }

  return(family)
}

## Second make.link function.
make.link2 <- function(link)
{
  if(is.null(link))
    link <- "identity"
  link0 <- link
  if(link0 == "tanhalf"){
    rval <- list(
      "linkfun" = function (mu) {
        tan(mu/2)},
      "linkinv" = function(eta) {
        2 * atan(eta)},
      "mu.eta" = function(eta) {
        2 / (eta^2 + 1)},
      "mu.eta2" = function(eta) {
        (-4 * eta ) / (eta^2 + 1)^2},
      "valideta" = function(eta) TRUE,
      "name" = "tanhalf"
      )
  } else {
    mu.eta2 <- function(x) {
      if(link0 == "identity") {
        x$mu.eta2 <- function(eta) rep.int(0, length(eta))
        return(x)
      }
      if(link0 == "log") {
        x$mu.eta2 <- function(eta) exp(eta)
        return(x)
      }
      if(link0 == "logit") {
        x$mu.eta2 <- function(eta) {
          eta <- exp(eta)
          return(-eta * (eta - 1) / (eta + 1)^3)
        }
        return(x)
      }
      if(link0 == "probit") {
        x$mu.eta2 <- function(eta) {
          -eta * dnorm(eta, mean = 0, sd = 1)
        }
        return(x)
      }
      if(link0 == "inverse") {
        x$mu.eta2 <- function(eta) {
          2 / (eta^3)
        }
        return(x)
      }
      if(link0 == "1/mu^2") {
        x$mu.eta2 <- function(eta) {
          0.75 / eta^(2.5)
        }
        return(x)
      }
      if(link0 == "sqrt") {
        x$mu.eta2 <- function(eta) { rep(2, length = length(eta)) }
        return(x)
      }
      x$mu.eta2 <- function(eta) rep.int(0, length(eta))
      ## warning(paste('higher derivatives of link "', link, '" not available!', sep = ''))
      return(x)
    }

    if(link %in% c("logit", "probit", "cauchit", "cloglog", "identity",
                   "log", "sqrt", "1/mu^2", "inverse")) {
      rval <- make.link(link)
    } else {
      rval <- switch(link,
        "rhogit" = list(
          "linkfun" = function(mu) { mu / sqrt(1 - mu^2) },
          "linkinv" = function(eta) {
              rval <- eta / sqrt(1 + eta^2)
              rval <- (abs(rval) - .Machine$double.eps) * sign(rval)
              rval
          },
          "mu.eta" = function(eta) { 1 / (1 + eta^2)^1.5 }
        ),
        "cloglog2" = list(
          "linkfun" = function(mu) { log(-log(mu)) },
          "linkinv" = function(eta) {
            pmax(pmin(1 - expm1(-exp(eta)), .Machine$double.eps), .Machine$double.eps)
          },
          "mu.eta" = function(eta) {
            eta <- pmin(eta, 700)
            pmax(-exp(eta) * exp(-exp(eta)), .Machine$double.eps)
          }
        ),
        "sigmoid" = list(
          "linkfun" = function(mu) {
            i <- mu <= -1
            if(any(i))
              mu[i] <- mu[i] <- -0.9999
            i <- mu >= 1
            if(any(i))
              mu[i] <- mu[i] <- 0.9999 
            -log(2/(mu + 1) - 1)
          },
          "linkinv" = function(eta) {
            tanh(eta/2)
          },
          "mu.eta" = function(eta) {
            0.5 / cosh(eta * 0.5)^2
          },
          "mu.eta2" = function(eta) {
            eta2 <- eta * 0.5
            -(0.5 * (2 * (sinh(eta2) * 0.5 * cosh(eta2)))/(cosh(eta2)^2)^2)
          }
        )
      )
    }

    rval <- mu.eta2(rval)
  }
  rval$name <- link
  rval
}


alpha2cap2 <- function(alpha = 0.01, nnobs = NULL, nvars = NULL, mean = 0) {
  t <- (1 - alpha)^(1/nvars)
  p <- function(cap) {
    pnorm(cap, mean = mean, sd = 1/sqrt(nnobs)) - pnorm(-cap, mean = mean, sd = 1/sqrt(nnobs))
  }
  vec <- 0:40000/1e+05
  s <- sapply(vec, FUN = function(i) abs(p(i) - t))
  cap <- vec[which.min(s)]
  return(cap)
}


## Stagewise distributional regression (SDR).
sdr <- function(formula, family = NULL, data = NULL, weights = NULL, batch_ids = NULL,
  cyclic = FALSE, subset = NULL, offset = NULL, contrasts = NULL, model = TRUE, CF = FALSE, cap = NULL,
  scalex = TRUE, boost = TRUE, eps = 0.01, nu = 0.1, ...) {
  
  
  stopifnot(requireNamespace("Formula"))
  stopifnot(requireNamespace("Matrix"))

  if (is.null(family))
    family <- gamlss.dist::NO
  if (is.function(family))
    family <- family()
  if (inherits(family, "gamlss.family"))
    family <- tF(family)
  family <- complete.bamlss.family(family)
  if (!is.list(formula))
    formula <- list(formula)
  if (length(formula) < length(family$names)) {
    for (j in (length(formula) + 1):length(family$names)) formula[[j]] <- ~1
  }
  formula <- rep(formula, length = length(family$names))
  names(formula) <- family$names

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
    names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)
  # print(mf)
  mfd <- list(x = list(), model = list())
  k <- 1
  for (i in names(formula)) {
    kk <- 1
    mfd$x[[i]] <- list()

    fi <- Formula::as.Formula(formula[[i]])

    fj <- formula(fi, lhs = 0, rhs = 1, drop = TRUE)
    vars <- all.vars(fj)
    if (length(vars) < 1) {
      ff <- ~1
    } else {
      ff <- as.formula(paste("~", paste(all.vars(fj), collapse = "+")))
    }

    mf[[2]] <- as.call(ff)

    mfd$model[[k]] <- eval(mf, parent.frame())

    # print(eval(mf, parent.frame()))
    mfd$x[[i]]$formula <- fj
    mfd$x[[i]]$terms <- terms(fj)
    mfd$x[[i]]$Formula <- fi

    if (length(rhs <- attr(mfd$x[[i]]$Formula, "rhs")) > 1) {
      mfd$x[[i]]$gfe <- list()
      for (j in 2:length(rhs)) {
        fj <- formula(mfd$x[[i]]$Formula, lhs = 0, rhs = j, drop = TRUE)
        ff <- as.formula(paste("~", paste(all.vars(fj), collapse = "+")))
        mf[[2]] <- as.call(ff)
        k <- k + 1
        mfd$model[[k]] <- eval(mf, parent.frame())
        mfd$x[[i]]$gfe[[j - 1]] <- list()
        mfd$x[[i]]$gfe[[j - 1]]$formula <- fj
        mfd$x[[i]]$gfe[[j - 1]]$fake.formula <- ff
        mfd$x[[i]]$gfe[[j - 1]]$terms <- terms(fj)
      }
    }

    k <- k + 1
  }


  fi <- Formula::as.Formula(formula[[1L]])
  fi <- formula(fi, lhs = 1, rhs = 0, drop = FALSE)
  mf[[2]] <- as.call(fi)

  mfd$y <- eval(mf, parent.frame())
  # return(mf) }

  mfd$model <- do.call("cbind", mfd$model)
  mfd$model <- mfd$model[, unique(names(mfd$model)), drop = FALSE]
  mfd$model <- na.omit(mfd$model)

  mfd$X <- list()
  if (scalex) {
    mfd$scaled <- list()
    mfd$Xscaled <- mfd$X
  }
  nvars <- NULL
  for (i in names(formula)) {
    tll <- attr(mfd$x[[i]]$terms, "term.labels")
    nvars <- c(nvars, length(tll))
    tls <- unlist(sapply(attr(mfd$x[[i]]$terms, "specials"), identity))
    if (!is.null(tls)) {
      tll <- tll[-tls]
      tls <- attr(mfd$x[[i]]$terms, "term.labels")[tls]
    }

    if (length(tll) < 1)
      tll <- "1"
    fl <- as.formula(paste("~", paste(tll, collapse = "+")))
    mfd$X[[i]] <- model.matrix(fl, data = mfd$model, contrasts)
    if (scalex)
      mfd$Xscaled[[i]] <- mfd$X[[i]]


    ## 0/1 scaling. changed to mittelwert abziehen und durch sd dividieren
    if ((ncol(mfd$X[[i]]) > 1) & scalex) {
      mfd$scaled[[i]] <- list()
      mfd$Xscaled[[i]] <- mfd$X[[i]]
      for (j in 2:ncol(mfd$X[[i]])) {
        mean.v <- mean(mfd$Xscaled[[i]][, j])
        sd.v <- sd(mfd$Xscaled[[i]][, j])

        mfd$scaled[[i]][[j - 1]] <- list(mean = mean.v, sd = sd.v)
        mfd$Xscaled[[i]][, j] <- (mfd$Xscaled[[i]][, j] - mfd$scaled[[i]][[j -
          1]]$mean)/mfd$scaled[[i]][[j - 1]]$sd
      }
    }
  }

  # nobs = batchsize 
  if (is.null(batch_ids)) nobs <- length(mfd$y[[1]])
  else if (is.numeric(batch_ids)) nobs <- batch_ids
  if (is.list(batch_ids)) nobs <- length(batch_ids[[1]])
  
  # Correlation Filtering
  if(CF){
    if(is.null(cap)){
      nvars <- max(nvars)
      
      if(length(mfd$y[[1]]) > nobs){
        # batchwise
        if(nobs <= 500) alpha = 0.1
        if(nobs <= 1000 & nobs > 500) alpha = 0.075
        if(nobs > 1000) alpha = 0.01
        
        cap <- alpha2cap2(alpha = alpha, nnobs = nobs, nvars = nvars, 
                          mean = 0.05)
      } else {
        # full sample
        if(nobs <= 500) alpha = 0.2
        if(nobs <= 1000 & nobs > 500) alpha = 0.1
        if(nobs > 1000) alpha = 0.01
        
        cap <- alpha2cap2(alpha = alpha, nnobs = nobs, nvars = nvars, 
                           mean = 0)
      }
      
    }
  } else cap <- 0
   
  if (scalex) {
    if (boost) {
      if (cyclic)
        mfd <- c(mfd, sdr.gradboostfit2(mfd$Xscaled, mfd$y[[1]], family = family, cap = cap, batch_ids = batch_ids,
                                        eps = eps, nu = nu, ...))
      if (!cyclic)
        mfd <- c(mfd, sdr.gradboostfit(mfd$Xscaled, mfd$y[[1]], family = family, cap = cap, batch_ids = batch_ids,
                                       eps = eps, nu = nu, ...))
      } else {
      mfd <- c(mfd, sdr.fit(mfd$Xscaled, mfd$y[[1]], family = family, cap = cap,batch_ids = batch_ids,
                            eps = eps, nu = nu,  ...))
    }
    nnn <- names(formula)
    mfd$unscaled.coef <- mfd[["coefficients"]]
    for (i in names(formula)) {
      nn <- names(mfd[["coefficients"]][[i]][1, ])
      if (length(nn) > 1) {
        if ("(Intercept)" %in% nn) {
          range <- 2:length(nn)
        } else {
          range <- 1:length(nn)
        }
        # rescaling
        mfd[["coefficients"]][[i]][, range] <- t(t(mfd[["coefficients"]][[i]][,
          range]) * (1/matrix(unlist(mfd[["scaled"]][[i]]), nrow = 2)[2,
          ]))
      }
    }
  } else {
    if (boost) {
      if (cyclic)
        mfd <- c(mfd, sdr.gradboostfit2(mfd$X, mfd$y[[1]], family = family, cap = cap, batch_ids = batch_ids,
                                        eps = eps, nu = nu, 
          ...))
      if (!cyclic)
        mfd <- c(mfd, sdr.gradboostfit(mfd$X, mfd$y[[1]], family = family, cap = cap,batch_ids = batch_ids,
                                       eps = eps,  nu = nu, 
          ...))
    } else {
      mfd <- c(mfd, sdr.fit(mfd$X, mfd$y[[1]], family = family, cap = cap, batch_ids = batch_ids,
                            eps = eps, nu = nu, ...))
    }
  }
  mfd$family <- family
  mfd$call <- cl
  mfd$formula <- formula
  mfd$nobs <- nobs
  

  if (!model)
    mfd["model"] <- NULL

  class(mfd) <- "stagewise"

  return(mfd)
}

## Main fitting functions.
# backfitting, no correlation filtering support 
sdr.fit <- function(X, y, family, maxit = 1000, start = NULL, 
                    eps = 0.01, nu = 0.1, 
                    eps_int = exp(seq(log(0.1), log(0.001), length = maxit)), 
                    nu_int = 0.05,
                    aic = FALSE, K = 0, full_grad = FALSE, cores = 1,
  batch_ids = NULL, verbose = TRUE, initialize = TRUE, int_ind = TRUE, plot = FALSE,
  steps = c(-1, 1), ts = F, coef_start = NULL, replace = F, comp = F,
  cap = 0, line.search = exp(seq(from = log(0.1), to = log(10), length.out = 10)),
  ...) {
  ia <- interactive()

  cap <- rep(cap, length.out = maxit)
  eps <- rep(eps, length.out = maxit)
  eps_int <- rep(eps_int, length.out = maxit)

  N <- nrow(X[[1L]])
  nx <- names(X)
  if (is.null(batch_ids)) {
    ind <- 1:N
    b.size <- N
    batch_ids <- lapply(1:maxit, function(...) ind)
  } else {
    if (is.numeric(batch_ids)) {
      ind <- 1:N
      b.size <- batch_ids
      batch_ids <- lapply(1:maxit, function(...) sample(ind, size = batch_ids,
        replace = replace))
    }
  }
  if (!is.list(batch_ids))
    stop("Argument batch_ids must be a list of indices!")
  if (length(batch_ids) != maxit)
    warning("Length of batch_ids != maxit, using batch_ids for setting maxit!")
  maxit <- length(batch_ids)
  tw <- length(strsplit(as.character(maxit), "")[[1]]) + 1L
  if (is.null(K))
    K <- log(length(y))

  ic.oos.list <- ic0.list <- ll.list <- ll0.list <- NULL
  beta <- list()
  for (j in nx) {
    beta[[j]] <- matrix(0, nrow = maxit, ncol = ncol(X[[j]]))
    colnames(beta[[j]]) <- colnames(X[[j]])
  }
  if (!is.null(coef_start)) {
    for (j in nx) {
      beta[[j]][1, ] <- coef_start[[j]]
    }
  }

  beta.grad <- beta

  if (initialize) {
    if (!is.null(family$initialize)) {
      betai <- list()
      nxx <- nx
      if (ts) {
        nxx <- "sigma"
      }
      for (j in nxx) {
        if (!is.null(family$initialize[[j]])) {
          linkfun <- make.link2(family$links[j])$linkfun
          beta[[j]][1L, "(Intercept)"] <- mean(linkfun(family$initialize[[j]](y,
          ...)), na.rm = TRUE)

        }
      }
    }
  }
  # print(beta[[1]][1L, '(Intercept)'])
  df <- sum(data.frame(beta)[1L, ] != 0)

  err01 <- .Machine$double.eps^(1/2)
  bs <- nrow(data.frame(batch_ids[1]))
  err02 <- err01 * 2 * bs
  # err02 is the denominator for central numeric differentiation.  b.size makes
  # gradient magnitute for different sample sizes comparable this translates
  # maximum log likelihood to maximum average loglikelihood
  # https://stats.stackexchange.com/questions/267847/motivation-for-average-log-likelihood

  ma <- function(x, n = 20) {
    ma1 <- filter(x, rep(1/n, n), sides = 1)
    ma2 <- rev(filter(rev(x), rep(1/n, n), sides = 1))
    ma3 <- ifelse(is.na(ma1), ma2, ma1)
    ma3
  }


  for (i in 2:maxit) {

    eta.oos <- eta <- val <- absval <- sgn <- list()

    ## Extract response.
    yi <- if (is.matrix(y)) {
      y[batch_ids[[i]], , drop = FALSE]
    } else {
      y[batch_ids[[i]]]
    }

    ## out of sample
    if (comp) {
      batch.oos <- setdiff(ind, batch_ids[[i]])
    } else {
      if (i != maxit) {
        i.oos <- i + 1
      } else {
        i.oos <- 1
      }
      batch.oos <- batch_ids[[i.oos]]
    }
    y.oos <- if (is.matrix(y)) {
      y[batch.oos, , drop = FALSE]
    } else {
      y[batch.oos]
    }


    for (j in nx) {
      ## Save last iteration.
      beta[[j]][i, ] <- beta[[j]][i - 1L, ]

      ## Setup linear predictor.
      eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% beta[[j]][i,
        ])

      ## out of sample
      eta.oos[[j]] <- drop(X[[j]][batch.oos, , drop = FALSE] %*% beta[[j]][i,
        ])

    }


    # out of sample
    ll.oos <- family$loglik(y.oos, family$map2par(eta.oos))

    ## Compute log-likelihood.
    df <- sum(data.frame(beta)[i, ] != 0)
    ll0 <- family$loglik(yi, family$map2par(eta))
    ic0 <- -2 * ll0 + K * df
    ic.oos <- -2 * ll.oos + K * df

    ll0.list <- c(ll0.list, ll0)
    ic0.list <- c(ic0.list, ic0)
    ic.oos.list <- c(ic.oos.list, ic.oos)

    if (verbose) {
      if (ia)
        cat("\r")
      cat("iter = ", formatC(i - 1L, width = tw, flag = " "), ", logLik = ",
        formatC(round(ll0, 4L), width = tw, flag = " "), ", df = ", formatC(df,
          width = tw, flag = " ") ,if (!ia)
          "\n" else NULL, sep = "")
    }
    bb = 20
    if (plot & (i%%10 == 0)) {
      if (i > bb + 4) {
        bic.min <- which.min(ma(ic0.list, order = bb))
        if (bic.min < bb + 1)
          bic.min <- bic.min + bb
      }
      par(mfrow = n2mfrow(length(nx) + 2))
      plot(y = ic.oos.list, x = 2:i, xlab = paste0(round(sign.list), 5), ylab = "BIC")
      if (i > bb + 4) {
        abline(v = bic.min)
        abline(v = bic.min - bb)
      }
      if (i > 5) {
        fit2 <- lowess(y = ic.oos.list, x = 2:i)
        lines(fit2)
      }
      plot(y = ll0.list, x = 2:i, xlab = "Iteration", ylab = "logLik")
      if (i > bb + 4) {
        abline(v = bic.min)
        abline(v = bic.min - bb)
      }
      if (i > 5) {
        fit2 <- lowess(y = ll0.list, x = 2:i)
        lines(fit2)
      }
      for (j in nx) {
        matplot(beta[[j]][1:i, ], type = "l", lty = 1, main = j, xlab = "Iteration",
          ylab = "Coefficients")
        if (i > bb + 4) {
          abline(v = bic.min)
          abline(v = bic.min - bb)
        }
      }
    }


    for (j in nx) {
      eta0 <- eta
      ## Get coefficients and setup.
      tbeta <- beta[[j]][i, ]
      ## Positive.
      tbeta[1] <- tbeta[1] + err01
      eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
      ll1 <- family$loglik(yi, family$map2par(eta))

      ## Negative
      tbeta[1] <- tbeta[1] - 2 * err01
      eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
      ll2 <- family$loglik(yi, family$map2par(eta))

      grad <- (ll1 - ll2)/err02
      grad <- eps_int[i] * grad
      if (i < 0.8 * maxit & abs(grad) <= nu_int * eps_int[i])
        grad <- nu_int * sign(grad) * eps_int[i]
      # eta[[j]] <- eta0[[j]]

      # intercept
      beta[[j]][i, 1] <- beta[[j]][i, 1] + grad
      tbeta <- beta[[j]][i, ]
      eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
      ll2 <- family$loglik(yi, family$map2par(eta))
      if (ll2 > ll0) {
        ll0 <- ll2
      } else {
        beta[[j]][i, 1] <- beta[[j]][i, 1] - grad
      }

    }

    eta0 <- eta
    sign.list <- pos.list <- rep(-Inf, length(nx))
    names(sign.list) <- names(pos.list) <- nx

    for (j in nx) {
      ## Get coefficients and setup.
      tbeta <- beta[[j]][i, ]
      nc <- length(tbeta)

      val[[j]] <- sgn[[j]] <- rep(-Inf, length = nc - 1)

      if (nc > 1) {
        # mclapply here?
        loop <- 2:nc

        for (jj in loop) {

          ## Positive.
          tbeta[jj] <- tbeta[jj] + err01
          eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
          ll1 <- family$loglik(yi, family$map2par(eta))

          ## Negative
          tbeta[jj] <- tbeta[jj] - 2 * err01
          eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
          ll2 <- family$loglik(yi, family$map2par(eta))

          grad <- (ll1 - ll2)/err02
          if (abs(grad) > eps[i])
          grad <- eps[i] * sign(grad)
          grad <- ifelse(i < 0.8 * maxit & abs(grad) < nu * eps[i], sign(grad) *
          nu * eps[i], grad)

          tbeta[jj] <- tbeta[jj] + err01 + grad

        }

        beta[[j]][i, ] <- tbeta
        eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% beta[[j]][i,
          ])

      }
    }
  }
  ## Compute log-likelihood.
  ll <- family$loglik(yi, family$map2par(eta))

  ## Extract 'out of sample' response.
  yi <- if (is.matrix(y)) {
    y[batch_ids[[1]], , drop = FALSE]
  } else {
    y[batch_ids[[1]]]
  }

  for (j in nx) {
    ## Setup linear predictor.
    eta[[j]] <- drop(X[[j]][batch_ids[[1]], , drop = FALSE] %*% beta[[j]][maxit,
      ])
  }
  ## Compute 'out of sample' log-likelihood.
  ll0 <- family$loglik(yi, family$map2par(eta))
  ll0.list <- c(ll0.list, ll0)

  if (verbose) {
    if (ia)
      cat("\r")
    cat("iter = ", formatC(i, width = tw, flag = " "), ", logLik = ", formatC(round(ll,
      4L), width = tw, flag = " "), "\n", sep = "")
  }

  rval <- list(coefficients = beta, logLik = ll0.list, maxit = maxit)

  return(rval)
}



# best subset gradboosting with correlation filtering
sdr.gradboostfit <- function(X, y, family, maxit = 1000, start = NULL, eps = 0.01, 
                             nu = 0.1, eps_int = exp(seq(log(0.1), log(0.001), length = maxit)), 
                             nu_int = 0.05,
                             aic = FALSE, K = 0, full_grad = FALSE, cores = 1,
  batch_ids = NULL, verbose = TRUE, initialize = TRUE, qstab = NULL, int_ind = TRUE,
  plot = FALSE, steps = c(-1, 1), ts = F, coef_start = NULL, replace = F,
  comp = F, cap = 0, length_ps = NULL, line.search = exp(seq(from = log(0.1), to = log(10),
    length.out = 10)), ...) {
  ia <- interactive()

  cap <- rep(cap, length.out = maxit)
  eps <- rep(eps, length.out = maxit)
  eps_int <- rep(eps_int, length.out = maxit)

  N <- nrow(X[[1L]])
  nx <- names(X)
  if (is.null(batch_ids)) {
    ind <- 1:N
    b.size <- N
    batch_ids <- lapply(1:maxit, function(...) ind)
  } else {
    if (is.numeric(batch_ids)) {
      ind <- 1:N
      b.size <- batch_ids
      batch_ids <- lapply(1:maxit, function(...) sample(ind, size = batch_ids,
        replace = replace))
    }
  }
  if (!is.list(batch_ids))
    stop("Argument batch_ids must be a list of indices!")
  if (length(batch_ids) != maxit)
    warning("Length of batch_ids != maxit, using batch_ids for setting maxit!")
  maxit <- length(batch_ids)
  tw <- length(strsplit(as.character(maxit), "")[[1]]) + 1L
  if (is.null(K))
    K <- log(length(y))

  ic.oos.list <- ic0.list <- ll.list <- ll0.list <- NULL
  beta <- list()
  for (j in nx) {
    beta[[j]] <- matrix(0, nrow = maxit, ncol = ncol(X[[j]]))
    colnames(beta[[j]]) <- colnames(X[[j]])
  }
  if (!is.null(coef_start)) {
    for (j in nx) {
      beta[[j]][1, ] <- coef_start[[j]]
    }
  }

  beta.grad <- beta
  powerset.list <- powerset(nx)[-1]
  if (!is.null(length_ps))
    powerset.list <- powerset.list[sapply(powerset.list, length) == length_ps]

  if (initialize) {
    if (!is.null(family$initialize)) {
      betai <- list()
      nxx <- nx
      if (ts) {
        nxx <- "sigma"
      }
      for (j in nxx) {
        if (!is.null(family$initialize[[j]])) {
          linkfun <- make.link2(family$links[j])$linkfun
          beta[[j]][1L, "(Intercept)"] <- mean(linkfun(family$initialize[[j]](y,
          ...)), na.rm = TRUE)

        }
      }
    }
  }
  # print(beta[[1]][1L, '(Intercept)'])
  df <- sum(data.frame(beta)[1L, ] != 0)

  err01 <- .Machine$double.eps^(1/2)
  err02 <- err01 * 2 * nrow(data.frame(batch_ids[1]))
  # err02 is the denominator for central numeric differentiation.  b.size makes
  # gradient magnitute for different sample sizes comparable this translates
  # maximum log likelihood to maximum average loglikelihood
  # https://stats.stackexchange.com/questions/267847/motivation-for-average-log-likelihood

  ma <- function(x, order = 20) {
    ma1 <- filter(x, rep(1/order, order), sides = 1)
    ma2 <- rev(filter(rev(x), rep(1/order, order), sides = 1))
    ma3 <- ifelse(is.na(ma1), ma2, ma1)
    ma3
  }


  for (i in 2:maxit) {

    eta.oos <- eta <- val <- absval <- sgn <- list()

    ## Extract response.
    yi <- if (is.matrix(y)) {
      y[batch_ids[[i]], , drop = FALSE]
    } else {
      y[batch_ids[[i]]]
    }

    ## out of sample
    if (comp) {
      batch.oos <- setdiff(ind, batch_ids[[i]])
    } else {
      if (i != maxit) {
        i.oos <- i + 1
      } else {
        i.oos <- 1
      }
      batch.oos <- batch_ids[[i.oos]]
    }
    y.oos <- if (is.matrix(y)) {
      y[batch.oos, , drop = FALSE]
    } else {
      y[batch.oos]
    }


    for (j in nx) {
      ## Save last iteration.
      beta[[j]][i, ] <- beta[[j]][i - 1L, ]

      ## Setup linear predictor.
      eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% beta[[j]][i,
        ])

      ## out of sample
      eta.oos[[j]] <- drop(X[[j]][batch.oos, , drop = FALSE] %*% beta[[j]][i,
        ])

    }


    # out of sample
    ll.oos <- family$loglik(y.oos, family$map2par(eta.oos))

    ## Compute log-likelihood.
    df <- sum(data.frame(beta)[i, ] != 0)
    ll0 <- family$loglik(yi, family$map2par(eta))
    ic0 <- -2 * ll0 + K * df
    ic.oos <- -2 * ll.oos + K * df

    ll0.list <- c(ll0.list, ll0)
    ic0.list <- c(ic0.list, ic0)
    ic.oos.list <- c(ic.oos.list, ic.oos)

    
    bb = 20
    if (plot & (i%%10 == 0)) {
      if (i > bb + 4) {
        bic.min <- which.min(ma(ic0.list, order = bb))
        if (bic.min < bb + 1)
          bic.min <- bic.min + bb
      }
      par(mfrow = n2mfrow(length(nx) + 2))
      plot(y = ic.oos.list, x = 2:i, xlab = "Iteration", ylab = "BIC")
      if (i > bb + 4) {
        abline(v = bic.min)
        abline(v = bic.min - bb)
      }
      if (i > 5) {
        fit2 <- lowess(y = ic.oos.list, x = 2:i)
        lines(fit2)
      }
      plot(y = ll0.list, x = 2:i, xlab = "Iteration", ylab = "logLik")
      if (i > bb + 4) {
        abline(v = bic.min)
        abline(v = bic.min - bb)
      }
      if (i > 5) {
        fit2 <- lowess(y = ll0.list, x = 2:i)
        lines(fit2)
      }
      for (j in nx) {
        matplot(beta[[j]][1:i, ], type = "l", lty = 1, main = j, xlab = "Iteration",
          ylab = "Coefficients")
        if (i > bb + 4) {
          abline(v = bic.min)
          abline(v = bic.min - bb)
        }
      }
    }


    for (j in nx) {
      eta0 <- eta
      ## Get coefficients and setup.
      tbeta <- beta[[j]][i, ]
      ## Positive.
      tbeta[1] <- tbeta[1] + err01
      eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
      ll1 <- family$loglik(yi, family$map2par(eta))

      ## Negative
      tbeta[1] <- tbeta[1] - 2 * err01
      eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
      ll2 <- family$loglik(yi, family$map2par(eta))

      grad <- (ll1 - ll2)/err02
      # grad <- ifelse(abs(grad) > eps_int[i], sign(grad)*eps_int[i],
      # ifelse(abs(grad) <= 0.1*eps_int[i], 0.1*sign(grad)*eps_int[i], grad))
      grad <- eps_int[i] * grad
      if (i < 0.8 * maxit & abs(grad) <= nu_int *  eps_int[i]) grad <- nu_int * sign(grad) * eps_int[i]
      # eta[[j]] <- eta0[[j]]

      # intercept
      beta[[j]][i, 1] <- beta[[j]][i, 1] + grad
      tbeta <- beta[[j]][i, ]
      eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
      ll2 <- family$loglik(yi, family$map2par(eta))
      if (ll2 > ll0) {
        ll0 <- ll2
      } else {
        beta[[j]][i, 1] <- beta[[j]][i, 1] - grad
      }

    }

    eta0 <- eta
    sign.list <- pos.list <- rep(-Inf, length(nx))
    names(sign.list) <- names(pos.list) <- nx

    for (j in nx) {

      ## Get coefficients and setup.
      tbeta <- beta[[j]][i, ]
      nc <- length(tbeta) - 1



      if (nc > 0) {
        eta[[j]] <- eta[[j]] + err01
        ll1 <- family$d(yi, family$map2par(eta), log = TRUE)

        ## Negative
        eta[[j]] <- eta[[j]] - 2 * err01
        ll2 <- family$d(yi, family$map2par(eta), log = TRUE)
        # if(T ) print(ll2)

        grad <- (ll1 - ll2)/(2 * err01)
        # print(mean(grad))
        eta[[j]] <- eta[[j]] + err01

        cc <- try(cor(grad, X[[j]][batch_ids[[i]], -1, drop = FALSE]), F)  #[,-1]))
        cc[is.na(cc)] <- 0
        if (is.numeric(cc)) {
          ## Select update
          jj <- which.max(abs(cc)) + 1
        } else {
          jj <- 1
        }  # interecept if cor gives error

        eta0 <- eta
        ## Get coefficients and setup.
        tbeta <- beta[[j]][i, ]

        if (max(abs(cc)) < cap[i]) {
          grad <- 0
        } else {
          ## Positive.
          tbeta[jj] <- tbeta[jj] + err01
          eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
          ll1 <- family$loglik(yi, family$map2par(eta))

          ## Negative
          tbeta[jj] <- tbeta[jj] - 2 * err01
          eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
          ll2 <- family$loglik(yi, family$map2par(eta))

          grad <- (ll1 - ll2)/err02
        }

        sign.list[[j]] <- grad
        pos.list[[j]] <- jj

      }

    }

    ps.final <- "no par"
    pset <- powerset.list
    for (j in nx) {
      if (sign.list[[j]] == 0 | sign.list[[j]] == -Inf) {
        ok <- !grepl(j, pset)
        pset <- pset[ok]
      }
    }

    if (length(pset) == 0) {
      
      for (j in nx) beta[[j]][i, ] <- beta[[j]][i - 1, ]
    } else {
      beta.final <- list()
      for (j in nx) beta.final[[j]] <- beta[[j]][i - 1, ]
      ic.oos.old <- ic.oos

      beta.wip <- list()
      for (j in nx) beta.wip[[j]] <- beta[[j]][i, ]
      for (l in 1:length(pset)) {
        # ps <- !is.na(pset[l,])
        sign.list2 <- sign.list[pset[[l]]]
        gnorm <- sqrt(sum(sign.list2^2))

        if (gnorm > eps[i]) {
          sign.list2 <- eps[i] * sign.list2/gnorm
        }
        for (ij in pset[[l]]) {
          sign.list2[ij] <- ifelse(i < 0.8 * maxit & abs(sign.list2[ij]) <
          nu * eps[i], sign(sign.list2[ij]) * nu * eps[i], sign.list2[ij])
        }

        for (j in pset[[l]]) {

          jj <- pos.list[[j]]
          grad <- sign.list2[[j]]
          
          beta[[j]][i, jj] <- beta[[j]][i, jj] + grad
          
          eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% beta[[j]][i,
          ])
        }


        df <- sum(data.frame(beta)[i, ] != 0)

        ## keep update only if oos information crit improves
        for (j in nx) eta.oos[[j]] <- drop(X[[j]][batch.oos, , drop = FALSE] %*%
          beta[[j]][i, ])
        ll.oos <- family$loglik(y.oos, family$map2par(eta.oos))
        ic.oos.new <- -2 * ll.oos + K * df
        if (ic.oos.new < ic.oos.old) {
          # keep current try
          for (j in nx) beta.final[[j]] <- beta[[j]][i, ]
          ps.final <- pset[[l]]
          ic.oos.old <- ic.oos.new
        }
        for (j in pset[[l]]) beta[[j]][i, ] <- beta.wip[[j]]

      }  # this bracket is for pset

      # save best beta
      for (j in nx) beta[[j]][i, ] <- beta.final[[j]]
      
    }

    if (verbose) {
      if (ia)
        cat("\r")
      cat("iter = ", formatC(i - 1L, width = tw, flag = " "), ", logLik = ",
          formatC(round(ll0, 4L), width = tw, flag = " "), ", df = ", formatC(df,
                                                                              width = tw, flag = " "), ", ", paste(ps.final, collapse = ", ") , "               " , if (!ia)
                                                                                "\n" else NULL, sep = "")
    }
  }

  ## Compute log-likelihood.
  ll <- family$loglik(yi, family$map2par(eta))

  ## Extract 'out of sample' response.
  yi <- if (is.matrix(y)) {
    y[batch_ids[[1]], , drop = FALSE]
  } else {
    y[batch_ids[[1]]]
  }

  for (j in nx) {
    ## Setup linear predictor.
    eta[[j]] <- drop(X[[j]][batch_ids[[1]], , drop = FALSE] %*% beta[[j]][maxit,
      ])
  }
  ## Compute 'out of sample' log-likelihood.
  ll0 <- family$loglik(yi, family$map2par(eta))
  ll0.list <- c(ll0.list, ll0)

  if (verbose) {
    if (ia)
      cat("\r")
    cat("iter = ", formatC(i, width = tw, flag = " "), ", logLik = ", formatC(round(ll,
      4L), width = tw, flag = " "), "\n", sep = "")
  }

  rval <- list(coefficients = beta, logLik = ll0.list, maxit = maxit)

  return(rval)
}


ma <- function(x, order = 20) {
  ma1 <- filter(x, rep(1/order, order), sides = 1)
  ma2 <- rev(filter(rev(x), rep(1/order, order), sides = 1))
  ma3 <- ifelse(is.na(ma1), ma2, ma1)
  ma3
}

# cyclical gradboosting with correlation filtering
sdr.gradboostfit2 <- function(X, y, family, maxit = 1000, start = NULL, eps = 0.01, nu = 0.1,
                              aic = FALSE, K = 0, full_grad = FALSE, cores = 1,
  batch_ids = NULL, verbose = TRUE, initialize = TRUE, qstab = NULL, int_ind = TRUE,
  plot = FALSE, steps = c(-1, 1), ts = F, coef_start = NULL, nu_int = 0.05, 
  eps_int = exp(seq(log(0.1), log(0.001), length = maxit)), replace = F,
  comp = F, cap = 0, line.search = exp(seq(from = log(0.1), to = log(10), length.out = 10)),
  ...) {
  ia <- interactive()

  cap <- rep(cap, length.out = maxit)
  eps <- rep(eps, length.out = maxit)
  eps_int <- rep(eps_int, length.out = maxit)

  N <- nrow(X[[1L]])
  nx <- names(X)
  if (is.null(batch_ids)) {
    ind <- 1:N
    b.size <- N
    batch_ids <- lapply(1:maxit, function(...) ind)
  } else {
    if (is.numeric(batch_ids)) {
      ind <- 1:N
      b.size <- batch_ids
      batch_ids <- lapply(1:maxit, function(...) sample(ind, size = batch_ids,
        replace = replace))
    }
  }
  if (!is.list(batch_ids))
    stop("Argument batch_ids must be a list of indices!")
  if (length(batch_ids) != maxit)
    warning("Length of batch_ids != maxit, using batch_ids for setting maxit!")
  maxit <- length(batch_ids)
  tw <- length(strsplit(as.character(maxit), "")[[1]]) + 1L
  if (is.null(K))
    K <- log(length(y))

  ic.oos.list <- ic0.list <- ll.list <- ll0.list <- NULL
  beta <- list()
  for (j in nx) {
    beta[[j]] <- matrix(0, nrow = maxit, ncol = ncol(X[[j]]))
    colnames(beta[[j]]) <- colnames(X[[j]])
  }
  if (!is.null(coef_start)) {
    for (j in nx) {
      beta[[j]][1, ] <- coef_start[[j]]
    }
  }

  beta.grad <- beta

  if (initialize) {
    if (!is.null(family$initialize)) {
      betai <- list()
      nxx <- nx
      if (ts) {
        nxx <- "sigma"
      }
      for (j in nxx) {
        if (!is.null(family$initialize[[j]])) {
          linkfun <- make.link2(family$links[j])$linkfun
          beta[[j]][1L, "(Intercept)"] <- mean(linkfun(family$initialize[[j]](y,
          ...)), na.rm = TRUE)

        }
      }
    }
  }
  # print(beta[[1]][1L, '(Intercept)'])
  df <- sum(data.frame(beta)[1L, ] != 0)

  err01 <- .Machine$double.eps^(1/2)
  err02 <- err01 * 2 * nrow(data.frame(batch_ids[1]))
  # err02 is the denominator for central numeric differentiation.  b.size makes
  # gradient magnitute for different sample sizes comparable this translates
  # maximum log likelihood to maximum average loglikelihood
  # https://stats.stackexchange.com/questions/267847/motivation-for-average-log-likelihood

  for (i in 2:maxit) {

    eta.oos <- eta <- val <- absval <- sgn <- list()

    ## Extract response.
    yi <- if (is.matrix(y)) {
      y[batch_ids[[i]], , drop = FALSE]
    } else {
      y[batch_ids[[i]]]
    }

    ## out of sample
    if (comp) {
      batch.oos <- setdiff(ind, batch_ids[[i]])
    } else {
      if (i != maxit) {
        i.oos <- i + 1
      } else {
        i.oos <- 1
      }
      batch.oos <- batch_ids[[i.oos]]
    }
    y.oos <- if (is.matrix(y)) {
      y[batch.oos, , drop = FALSE]
    } else {
      y[batch.oos]
    }


    for (j in nx) {
      ## Save last iteration.
      beta[[j]][i, ] <- beta[[j]][i - 1L, ]

      ## Setup linear predictor.
      eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% beta[[j]][i,
        ])

      ## out of sample
      eta.oos[[j]] <- drop(X[[j]][batch.oos, , drop = FALSE] %*% beta[[j]][i,
        ])

    }


    # out of sample
    ll.oos <- family$loglik(y.oos, family$map2par(eta.oos))

    ## Compute log-likelihood.
    df <- sum(data.frame(beta)[i, ] != 0)
    ll0 <- family$loglik(yi, family$map2par(eta))
    ic0 <- -2 * ll0 + K * df
    ic.oos <- -2 * ll.oos + K * df

    ll0.list <- c(ll0.list, ll0)
    ic0.list <- c(ic0.list, ic0)
    ic.oos.list <- c(ic.oos.list, ic.oos)

    if (verbose) {
      if (ia)
        cat("\r")
      cat("iter = ", formatC(i - 1L, width = tw, flag = " "), ", logLik = ",
        formatC(round(ll0, 4L), width = tw, flag = " "), ", df = ", formatC(df,
          width = tw, flag = " "), if (!ia)
          "\n" else NULL, sep = "")
    }
    bb = 20
    if (plot & (i%%10 == 0)) {
      if (i > bb + 4) {
        bic.min <- which.min(ma(ic0.list, order = bb))
        if (bic.min < bb + 1)
          bic.min <- bic.min + bb
      }
      par(mfrow = n2mfrow(length(nx) + 2))
      plot(y = ic.oos.list, x = 2:i, xlab = paste0(round(sign.list), 5), ylab = "BIC")
      if (i > bb + 4) {
        abline(v = bic.min)
        abline(v = bic.min - bb)
      }
      if (i > 5) {
        fit2 <- lowess(y = ic.oos.list, x = 2:i)
        lines(fit2)
      }
      plot(y = ll0.list, x = 2:i, xlab = "Iteration", ylab = "logLik")
      if (i > bb + 4) {
        abline(v = bic.min)
        abline(v = bic.min - bb)
      }
      if (i > 5) {
        fit2 <- lowess(y = ll0.list, x = 2:i)
        lines(fit2)
      }
      for (j in nx) {
        matplot(beta[[j]][1:i, ], type = "l", lty = 1, main = j, xlab = "Iteration",
          ylab = "Coefficients")
        if (i > bb + 4) {
          abline(v = bic.min)
          abline(v = bic.min - bb)
        }
      }
    }


    for (j in nx) {
      eta0 <- eta
      ## Get coefficients and setup.
      tbeta <- beta[[j]][i, ]
      ## Positive.
      tbeta[1] <- tbeta[1] + err01
      eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
      ll1 <- family$loglik(yi, family$map2par(eta))

      ## Negative
      tbeta[1] <- tbeta[1] - 2 * err01
      eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
      ll2 <- family$loglik(yi, family$map2par(eta))

      grad <- (ll1 - ll2)/err02
      # grad <- ifelse(abs(grad) > eps_int[i], sign(grad)*eps_int[i],
      # ifelse(abs(grad) <= 0.1*eps_int[i], 0.1*sign(grad)*eps_int[i], grad))
      grad <- eps_int[i] * grad
      if (i < 0.8 * maxit & abs(grad) <= nu_int * eps_int[i])
        grad <- nu_int * sign(grad) * eps_int[i]
      # eta[[j]] <- eta0[[j]]

      # intercept
      beta[[j]][i, 1] <- beta[[j]][i, 1] + grad
      tbeta <- beta[[j]][i, ]
      eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
      ll2 <- family$loglik(yi, family$map2par(eta))
      if (ll2 > ll0) {
        ll0 <- ll2
      } else {
        beta[[j]][i, 1] <- beta[[j]][i, 1] - grad
      }

    }

    eta0 <- eta
    sign.list <- pos.list <- rep(0, length(nx))
    names(sign.list) <- names(pos.list) <- nx

    for (j in nx) {

      ## Get coefficients and setup.
      tbeta <- beta[[j]][i, ]
      nc <- length(tbeta) - 1



      if (nc > 0) {
        eta[[j]] <- eta[[j]] + err01
        ll1 <- family$d(yi, family$map2par(eta), log = TRUE)

        ## Negative
        eta[[j]] <- eta[[j]] - 2 * err01
        ll2 <- family$d(yi, family$map2par(eta), log = TRUE)
        # if(T ) print(ll2)

        grad <- (ll1 - ll2)/(2 * err01)
        # print(mean(grad))
        eta[[j]] <- eta[[j]] + err01

        cc <- try(cor(grad, X[[j]][batch_ids[[i]], -1, drop = FALSE]), F)  #[,-1]))
        cc[is.na(cc)] <- 0
        if (is.numeric(cc)) {
          ## Select update
          jj <- which.max(abs(cc)) + 1
        } else {
          jj <- 1
        }  # interecept if cor gives error

        eta0 <- eta
        ## Get coefficients and setup.
        tbeta <- beta[[j]][i, ]

        if (max(abs(cc)) < cap[i]) {
          grad <- 0
        } else {
          ## Positive.
          tbeta[jj] <- tbeta[jj] + err01
          eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
          ll1 <- family$loglik(yi, family$map2par(eta))

          ## Negative
          tbeta[jj] <- tbeta[jj] - 2 * err01
          eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
          ll2 <- family$loglik(yi, family$map2par(eta))

          grad <- (ll1 - ll2)/err02
        }

        sign.list[[j]] <- grad
        pos.list[[j]] <- jj

      }

    }


    gnorm <- sqrt(sum(sign.list^2))

    if (gnorm > eps[i]) {
      sign.list <- eps[i] * sign.list/gnorm
    }
    for (ij in nx) {
      sign.list[ij] <- ifelse(i < 0.8 * maxit & sign.list[ij] != 0 & abs(sign.list[ij]) <
        nu * eps[i], sign(sign.list[ij]) * nu * eps[i], sign.list[ij])
    }
    for (j in nx) {
      if (sign.list[[j]] != 0) {
        jj <- pos.list[[j]]
        grad <- sign.list[[j]]
        beta[[j]][i, jj] <- beta[[j]][i, jj] + grad

        eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% beta[[j]][i,
          ])
      }
    }


    df <- sum(data.frame(beta)[i, ] != 0)

    ## keep update only if oos information crit improves
    for (j in nx) eta.oos[[j]] <- drop(X[[j]][batch.oos, , drop = FALSE] %*%
      beta[[j]][i, ])
    ll.oos <- family$loglik(y.oos, family$map2par(eta.oos))
    if (K >= 0)
      ic.oos.new <- -2 * ll.oos + K * df else ic.oos.new <- -Inf
    if (!ic.oos.new < ic.oos) {
      for (j in nx) {
        beta[[j]][i, ] <- beta[[j]][i - 1, ]
        # print('no update')
      }
    }



  }


  ## Compute log-likelihood.
  ll <- family$loglik(yi, family$map2par(eta))

  ## Extract 'out of sample' response.
  yi <- if (is.matrix(y)) {
    y[batch_ids[[1]], , drop = FALSE]
  } else {
    y[batch_ids[[1]]]
  }

  for (j in nx) {
    ## Setup linear predictor.
    eta[[j]] <- drop(X[[j]][batch_ids[[1]], , drop = FALSE] %*% beta[[j]][maxit,
      ])
  }
  ## Compute 'out of sample' log-likelihood.
  ll0 <- family$loglik(yi, family$map2par(eta))
  ll0.list <- c(ll0.list, ll0)

  if (verbose) {
    if (ia)
      cat("\r")
    cat("iter = ", formatC(i, width = tw, flag = " "), ", logLik = ", formatC(round(ll,
      4L), width = tw, flag = " "), "\n", sep = "")
  }

  rval <- list(coefficients = beta, logLik = ll0.list, maxit = maxit)

  return(rval)
}




## Model frame.
model.frame.stagewise <- function(formula, data = NULL, ...) {
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
  if (length(nargs) || is.null(formula$model) || !is.null(data)) {
    fcall <- formula$call
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
      names(fcall), 0L)
    fcall <- fcall[c(1L, m)]
    fcall$drop.unused.levels <- TRUE
    if (!is.null(data))
      fcall$data <- as.name(quote("data"))
    fcall[[1L]] <- quote(model.frame)
    mf <- list()
    k <- 1
    for (i in names(formula$formula)) {
      fi <- Formula::as.Formula(formula$formula[[i]])
      fj <- formula(fi, lhs = 0, rhs = 1, drop = TRUE)
      vars <- all.vars(fj)
      if (length(vars) < 1) {
        ff <- ~1
      } else {
        ff <- as.formula(paste("~", paste(all.vars(fj), collapse = "+")))
      }
      fcall[[2]] <- as.call(ff)
      mf[[k]] <- eval(fcall, if (is.null(data))
        parent.frame() else NULL)
      k <- k + 1
    }
    mf <- do.call("cbind", mf)
    mf <- mf[, unique(names(mf)), drop = FALSE]
    return(mf)
  } else {
    return(formula$model)
  }
}

## Model matrix.
model.matrix.stagewise <- function(object, ...) {
  mf <- model.frame(object, ...)

  X <- list()
  for (i in names(object$formula)) {
    tll <- attr(object$x[[i]]$terms, "term.labels")
    tls <- unlist(sapply(attr(object$x[[i]]$terms, "specials"), identity))
    if (!is.null(tls)) {
      tll <- tll[-tls]
      tls <- attr(object$x[[i]]$terms, "term.labels")[tls]
    }

    if (length(tll) < 1)
      tll <- "1"
    fl <- as.formula(paste("~", paste(tll, collapse = "+")))
    X[[i]] <- model.matrix(fl, data = mf)

    # ## 0/1 scaling. fix to mittelwert abziehen und durch sd dividieren
    # if((ncol(X[[i]]) > 1) & !is.null(object$scaled)) { for(j in
    # 2:ncol(X[[i]])) { #xr <- object$scaled[[i]][[j - 1]]$max -
    # object$scaled[[i]][[j - 1]]$min X[[i]][, j] <- (X[[i]][, j] -
    # mean(X[[i]][, j]) ) / sd(X[[i]][, j]) # (X[[i]][, j] -
    # object$scaled[[i]][[j - 1]]$min) / xr } }
  }
  return(X)
}

## Predict method.
predict.stagewise <- function(object, newdata = NULL, type = c("link", "parameter"), drop = TRUE, mstart = NULL, mstop = NULL, model = NULL,
   ...) {
  # X1 <- model.matrix(object, data = newdata) X <- object[['X']] for(i in
  # names(object$formula)){ if(ncol(X[[i]]) == 1){ # only intercept (or if no
  # int spec, then only one variable in model matrix) X[[i]] <-
  # X[[i]][1:nrow(X1[[i]])]*0 X[[i]] <- X[[i]]+X1[[i]] } else { # normal case,
  # multiple columns in model matrix X[[i]] <- X[[i]][1:nrow(X1[[i]]),]*0 for(j
  # in colnames(X1[[i]])){ if(!is.null(ncol(X[[i]]))) { X[[i]][,j] <-
  # X[[i]][,j]+X1[[i]][,j] } else { X[[i]][j] <- X[[i]][j]+X1[[i]][j] } } } }
  X <- model.matrix(object, data = newdata, ...)
  fun <- "mean" # matrix alternative not working properly
  if (!is.null(mstop)) {
    if (length(mstop) > 2L)
      stop("Argument mstop must be a single number!")
  }
  maxit <- nrow(object$coefficients[[1]])
  if (is.null(mstop))
    mstop <- maxit
  if (is.null(mstart))
    mstart <- min(round(maxit/2), 300)
  if (is.null(model))
    model <- names(object$formula)
  if (is.character(model)) {
    model <- match.arg(model, names(object$formula), several.ok = TRUE)
  }
  if (is.null(type))
    type <- "link"
  if (length(type) > 1)
    type <- type[1]

  eta <- list()
  if (type != "terms") {

    for (j in model) {
      eta[[j]] <- if (fun == "mean") {
        0
      }
      if (fun == "matrix") {
        matrix(0, nrow = nrow(X[[j]]), ncol = mstop - mstart + 1)
      }
      iterations <- mstart:mstop
      for (i in iterations) {
        if (fun == "mean") {
          eta[[j]] <- eta[[j]] + drop(X[[j]] %*% object$coefficients[[j]][i,
          ])
        }
        if (fun == "matrix") {
          eta[[j]][, i] <- drop(X[[j]] %*% object$coefficients[[j]][i, ])
        }
      }
      if (fun == "mean")
        eta[[j]] <- eta[[j]]/(mstop - mstart + 1)
      if (type != "link") {
        linkinv <- make.link2(object$family$links[j])$linkinv
        eta[[j]] <- linkinv(eta[[j]])
      }
    }
    if ((length(eta) < 2L) & drop)
      eta <- eta[[1L]]
  }
  
  return(eta)
}

# Extract coefficients.
coef.stagewise <- function(object, model = NULL, refit = FALSE,
                           mstop = NULL, mstart = mstop, ...) {
  coef_fun <- function(...){
    if (!is.null(mstop)) {
      if (length(mstop) > 2L)
        stop("Argument mstop must be a single number!")
    }
    keep = NULL 
    if (is.null(model))
      model <- names(object$formula)
    if (is.character(model))
      model <- match.arg(model, names(object$formula), several.ok = TRUE) else model <- names(object$formula)[model]
      maxit <- nrow(object$coefficients[[1]])
      if (is.null(keep))
        keep <- maxit else keep <- floor(maxit * (1 - keep))
      coef <- list()
      for (j in model) coef[[j]] <- apply(object$coefficients[[j]][if (is.null(mstop))
        keep:maxit else mstart:mstop, , drop = FALSE], 2, mean)
      if (length(coef) < 2L)
        coef <- coef[[1L]]
      return(coef)
  }
  
  if(refit){
    if (is.null(mstop)){ 
      mstop <- nrow(object$coefficients[[1]])
      cat("Last iteration is used as 'mstop' is not provided.")
    }
    if(is.null(mstart)) mstart <- mstop
    if(!mstart == mstop) cat("'mstart' is ignored. Coefs are based only on 'mstop'")
    
    if (!is.null(object[["unscaled.coef"]])) {
      # exist unscaled object?
      coef <- object[["unscaled.coef"]][names(object[["unscaled.coef"]])]
      nc <- names(coef)
      coef <- lapply(nc, FUN = function(i) {
        coef[[i]] <- coef[[i]][mstop, ]
      })
      
    } else {
      coef <- coef_fun(mstart = mstop, mstop = mstop, ...)
      nc <- names(coef)
    }
    names(coef) <- nc
    
    coef <- lapply(nc, FUN = function(i) {
      ind <- coef[[i]] != 0
      ind[1] <- TRUE  # Intercept is always true
      coef[[i]] <- coef[[i]][ind]
    })
    names(coef) <- nc
  } else {
    coef <- coef_fun(...)
  }
  
  return(coef)
}


# coef.stagewise <- function(object, ...) {
# 
# }

### update function for refitting
# newformula
newformula <- function(object, mstop = NULL, name = NULL) {
  if (is.null(mstop)){
    mstop <- nrow(object$coefficients[[1]])
    cat("Last iteration is used as 'mstop' is not provided.")
  }
    
  coef <- coef.stagewise(object = object, mstop = mstop, refit = TRUE)
  nc <- names(coef)
  coef <- lapply(nc, FUN = function(i) {
    xx <- setdiff(names(coef[[i]]), "(Intercept)")
    if (length(xx) == 0)
      xx <- 1
    if (i == nc[1]) {
      if (!is.null(name)) {
        f <- as.formula(paste(name, " ~ ", paste(xx, collapse = "+")))
      } else {
        f <- as.formula(paste("y ~ ", paste(xx, collapse = "+")))
      }
      
    } else {
      f <- as.formula(paste(" ~ ", paste(xx, collapse = "+")))
    }
    
    return(f)
  })
  names(coef) <- nc 
  
  return(coef)
}












## Summary method.
summary.stagewise <- function(object,
                              digits = max(3, getOption("digits") - 3), mstart = round(0.5*object$maxit),
                              mstop = object$maxit, ...) {
  x <- object
  # for computing parameter summary
  parsum <- function(d, vec = mstart:mstop) {
    dd <- t(d)[, 1:4]
    if (ncol(d) == 1) {
      dd[1] <- mean(d[vec, ])
      dd[2:4] <- quantile(d[vec, ], c(0.025, 0.5, 0.975))
      dd <- matrix(dd, ncol = 4)
      rownames(dd) <- "(Intercept)"
    } else {
      if(length(vec) == 1){
        dd[, 1] <- t(apply(d[c(vec,vec), ], 2, mean))
        dd[, 2:4] <- t(apply(d[c(vec,vec), ], 2, quantile, c(0.025, 0.5, 0.975)))
        # dd[,5] <- t(d[nrow(d),])
      } else {
        dd[, 1] <- t(apply(d[vec, ], 2, mean))
        dd[, 2:4] <- t(apply(d[vec, ], 2, quantile, c(0.025, 0.5, 0.975)))
        # dd[,5] <- t(d[nrow(d),])  
      }
      
    }
    colnames(dd) <- c("mean", "2.5%", "50%", "97.5%")
    return(dd)
  }
  
  
  cat("\nCall:\n")
  print(x$call)
  cat("---\n")
  print(x$family, full = FALSE)
  cat("*---\n")
  
  cat("\nmstart:\n")
  x$mstart <- mstart
  print(x$mstart)
  cat("---\n")
  
  cat("\nmstop:\n")
  x$mstop <- mstop
  print(x$mstop)
  cat("---\n")
  
  for (i in names(x$formula)) {
    if (length(names(x$formula)) > 1) {
      cat("\nFormula", i, ":\n")
      cat("---\n")
      print(unname(x$formula[i]))
      cat("---\n")
    }
    if (length(names(x$formula)) == 1) {
      cat("\nFormula mu:\n")
      cat("---\n")
      print(x$formula)
      cat("---\n")
    }
    if (!is.null(x$coefficients[[i]])) {
      cat("-\n")
      cat("Parametric coefficients:\n")
      x$parsum[[i]] <- round(parsum(x$coefficients[[i]]), digits)
      print(x$parsum[[i]])
      # printCoefmat(parsum(x$coefficients[[i]], digits = digits)
      cat("---\n")
    }
  }
  
  cat("\nlogLik-mstart-mstop:\n")
  print(x$logLik[c(x$mstart, x$mstop)])
  cat("---\n")
  return(invisible(x[6:12]))
  
}


## From bamlss.
print_bamlss_formula <- function(x, ...) {
  if(!inherits(x, "list") & !inherits(x, "bamlss.formula")) {
    print(x)
  } else {
    nx <- names(x)
    if(is.null(nx))
      nx <- as.character(1:length(x))
    for(i in seq_along(x)) {
      cat("Formula ", nx[i], ":\n---\n", sep = "")
      if(inherits(x[[i]], "list") & "h1" %in% names(x[[i]])) {
        for(j in seq_along(x[[i]])) {
          cat("h", j, ": ", sep = "")
          attr(x[[i]][[j]], "name") <- NULL
          attr(x[[i]][[j]]$formula, ".Environment") <- NULL
          if(is.character(x[[i]][[j]]$formula)) {
            cat(x[[i]][[j]]$formula, "\n")
          } else print(x[[i]][[j]]$formula, showEnv = FALSE)
        }
      } else {
        attr(x[[i]], "name") <- NULL
        attr(x[[i]]$formula, "name") <- NULL
        attr(x[[i]]$formula, ".Environment") <- NULL
        if("formula" %in% names(x[[i]])) {
          if(is.character(x[[i]]$formula))
            cat(x[[i]]$formula, "\n")
          else
            print(x[[i]]$formula, showEnv = FALSE)
        } else print(x[[i]])
      }
      if(i < length(x))
      cat("\n")
    }
  }
  invisible(NULL)
}

##### Print needs rework Print method for summary.
print.summary.stagewise <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("---\n")
  print(x$family, full = FALSE)
  cat("*---\n")

  cat("\nIterations:\n")
  print(x$maxit)
  cat("---\n")

  cat("\nburnin:\n")
  print(x$burnin)
  cat("---\n")


  for (i in names(x$formula)) {
    print_bamlss_formula(x$formula[i])
    if (!is.null(x$parsum[[i]])) {
      cat("-\n")
      cat("Parametric coefficients:\n")
      print(x$parsum[[i]])
      # printCoefmat(parsum(x$coefficients[[i]], digits = digits)
      cat("---\n")
    }
  }
  cat("\nlogLik:\n")
  print(x$logLik)
  cat("---\n")

}





logLik.stagewise <- function(object, mstart = 1, mstop = object$maxit,
  all = TRUE, ...) {
  ll <- object$logLik[mstart:mstop]

  ### degrees of freedom
  df <- rep(0, mstop - mstart + 1)
  beta <- object$coefficients
  for (i in 1:length(beta)) {
    if (length(data.frame(beta[i])) != 1)
      df <- df + rowSums(ifelse(data.frame(beta[i])[mstart:mstop, ] ==
        0, 0, 1))
  }
  if (!all) {
    ll <- mean(ll)
    df <- mean(df)
  }
  attr(ll, "df") <- df
  class(ll) <- c("logLik", "logLik.stagewise")
  return(ll)
}

AIC.stagewise <- function(object, K = 2, ...) {
  ll <- logLik(object, ...)
  ic <- as.numeric(-2 * ll + K * attr(ll, "df"))
  
  class(ic) <- "stagewise_AIC"
  return(ic)
}

BIC.stagewise <- function(object, ...) {
  AIC.stagewise(object, K = log(object$nobs), ...)
}

## Plotting method.  terms plot needs more work
plot.stagewise <- function(x, which = c("all", "coefficients", "AIC"), K = 2,
                           bw = 0, spar = TRUE, ...) {
  if(length(which) > 1) which <- which[1]
  if (which == "all") {
    which <- match.arg(which)
    ll <- logLik(x, ...)
    ic <- as.numeric(-2 * ll + K * attr(ll, "df"))
    bic.min <- ifelse(bw > 0, which.min(ma(ic, order = bw)), which.min(ic))
    nx <- names(x$formula)
    if (spar)
      par(mfrow = n2mfrow(length(nx) + 2))
    
    plot(y = ic, x = 1:length(ic), xlab = "Iteration", ylab = "AIC")
    abline(v = bic.min)
    abline(v = bic.min - bw)
    
    plot(y = ll, x = 1:length(ll), xlab = "Iteration", ylab = "logLik")
    abline(v = bic.min)
    abline(v = bic.min - bw)
    
    
    for (j in nx) {
      matplot(x$coefficients[[j]], type = "l", lty = 1, main = j, xlab = "Iteration",
              ylab = "Coefficients")
      abline(v = bic.min, lty = 2)
      abline(v = bic.min - bw, lty = 2)
    }
  }
  
  
  if (which == "coefficients") {
    nx = names(x$formula)
    dimpar = n2mfrow(length(nx))
    par(mfrow = dimpar, mar = c(6, 3, 2, 6))
    
    ll <- logLik(x, ...)
    ic <- as.numeric(-2 * ll + K * attr(ll, "df"))
    bic.min <- ifelse(bw > 0, which.min(ma(ic, order = bw)), which.min(ic))
    
    for (j in nx) {
      matplot(x$coefficients[[j]], type = "l", lty = 1, cex.main = 1.5,
              cex.lab = 1.5, cex.axis = 1.3, xlab = "Iteration", ylab = "", main = j)
      # for(i in colnames(x$coefficients[[j]][,-1]))
      abline(v = bic.min, lty = 2)
      abline(v = bic.min - bw, lty = 2)
      l <- nrow(x$coefficients[[j]])
      what <- x$coefficients[[j]][l, ] != 0
      nam <- colnames(x$coefficients[[j]])[what]
      plab <- x$coefficients[[j]][l, nam]
      rang <- abs(diff(par("usr")[3:4]))
      o <- order(plab, decreasing = TRUE)
      nam1 <- nam <- nam[o]
      plab <- plab[o]
      if(length(plab) > 1){
        for (i in 1:(length(plab)-1)) {
          dp <- abs(plab[i] - plab[i + 1])/rang
          if (dp <= 0.025) {
            nam[i + 1] <- paste(c(nam[i], nam[i + 1]), collapse = ",")
            nam[i] <- ""
          }
        }
      } 
      
      for (i in 1:length(nam)) {
        pos <- plab[i]
        if (pos != 0)
          text(x = par("usr")[2], labels = nam[i], y = pos, pos = 4, xpd = TRUE,
               cex = 1.1)
      }
    }
    
  }
  
  
  if (which == "AIC"){
    ll <- logLik(x, ...)
    ic <- as.numeric(-2 * ll + K * attr(ll, "df"))
    bic.min <- ifelse(bw > 0, which.min(ma(ic, order = bw)), which.min(ic))
    
    plot(y = ic, x = 1:length(ic), xlab = "Iteration", ylab = "AIC")
    abline(v = bic.min, lty = 2)
    abline(v = bic.min - bw, lty = 2)
    
    }
  
 
  return(invisible(NULL))
}


residuals.stagewise <- function(object, type = c("quantile", "response"), nsamps = NULL,
  ...) {
  family <- object$family
  if (!is.null(family$residuals)) {
    res <- family$residuals(object, type = type, nsamps = nsamps, ...)
    if (length(class(res)) < 2) {
      if (inherits(res, "numeric"))
        class(res) <- c("stagewise_residuals", class(res))
    }
  } else {
    type <- match.arg(type)
    y <- NULL
    if (!is.null(nsamps)) {
      y <- nsamps[, names(object[["y"]])]
      par <- predict(object, newdata = nsamps, drop = FALSE, ...)
    } else {
      y <- unlist(object$y)
      par <- predict(object, newdata = NULL, drop = FALSE, ...)
    }
    # if (!is.null(object$y)) { y <- if (is.data.frame(object$y)) { if
    # (ncol(object$y) < 2) { object$y[[1]] } else object$y } else { object$y }
    # } if (!is.null(nd <- list(...)$newdata)) { #is.null(nsamps)){ # rn <-
    # names(object[['y']]) y <- nd[[rn]] if (is.null(y)) stop(paste('the
    # response', rn, 'is not available in newdata!')) #rm(nd) n <- if
    # (is.null(dim(y))) length(y) else nrow(y) } if (is.null(y)) { rn <-
    # names(object[['y']]) y <- object$y[[rn]] } y <- y[-nas]if (is.null(y))
    # stop('response variable is missing, cannot compute residuals!')

    nas <- attr(par, "na.action")
    if (!is.null(nas)) {
      if (is.null(dim(y))) {
      } else {
        y <- y[-nas, ]
      }
    }
    nod <- is.null(dim(par[[1L]]))
    for (j in family$names) {
      if (!nod)
        par[[j]] <- as.matrix(par[[j]])
      par[[j]] <- make.link2(family$links[j])$linkinv(par[[j]])
    }
    if (type == "quantile") {
      if (is.null(family$p)) {
        type <- "response"
        warning(paste("no $p() function in family '", family$family, "', cannot compute quantile residuals, computing response resdiuals instead!",
          sep = ""))
      } else {
        discrete <- FALSE
        if (!is.null(family$type)) {
          if (tolower(family$type) == "discrete")
          discrete <- TRUE
        }
        if (family$family == "binomial")
          discrete <- TRUE
        if (discrete) {
          ymin <- min(y, na.rm = TRUE)
          a <- family$p(ifelse(y == ymin, y, y - 1), par)
          a <- ifelse(y == ymin, 0, a)
          b <- family$p(y, par)
          u <- runif(length(y), a, b)
          u <- ifelse(u > 0.999999, u - 1e-16, u)
          u <- ifelse(u < 1e-06, u + 1e-16, u)
          res <- qnorm(u)
        } else {
          prob <- family$p(y, par)
          res <- qnorm(prob)
          if (any(isnf <- !is.finite(res))) {
          warning("non finite quantiles from probabilities, set to NA!")
          res[isnf] <- NA
          }
        }
        attr(res, "type") <- "Quantile"
      }
    }
    if (type == "response") {
      mu <- if (is.null(family$mu)) {
        function(par, ...) {
          par[[1]]
        }
      } else family$mu
      res <- y - mu(par)
      attr(res, "type") <- "Response"
    }
    class(res) <- c("stagewise_residuals", class(res))
  }
  if (any(j <- !is.finite(res)))
    res[j] <- NA
  return(res)
}


c95 <- function (x) 
{
  qx <- quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  return(c(qx[1], Mean = mean(x, na.rm = TRUE), qx[2]))
}

plot.stagewise_residuals <- function(x, which = c("hist-resid", "qq-resid", "wp"),
  spar = TRUE, ...) {
  which.match <- c("hist-resid", "qq-resid", "wp")
  if (!is.character(which)) {
    if (any(which > 3L))
      which <- which[which <= 3L]
    which <- which.match[which]
  } else which <- which.match[pmatch(tolower(which), which.match)]
  if (length(which) > length(which.match) || !any(which %in% which.match))
    stop("argument which is specified wrong!")
  if (is.null(dim(x)))
    x <- matrix(x, ncol = 1)
  nc <- ncol(x)
  cn <- colnames(x)
  if (nc > 10) {
    nc <- 1
    cn <- NULL
  }
  add <- list(...)$add
  if (is.null(add))
    add <- FALSE
  if (add)
    spar <- FALSE
  if (spar) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow = n2mfrow(length(which) * nc))
  }
  type <- attr(x, "type")
  for (j in 1:nc) {
    for (w in which) {
      args <- list(...)
      if (w == "hist-resid") {
        if (ncol(x) > 1) {
          x2 <- rowMeans(x, na.rm = TRUE)
        } else {
          x2 <- x
        }
        rdens <- density(as.numeric(x2), na.rm = TRUE)
        rh <- hist(as.numeric(x2), plot = FALSE)
        args$ylim <- c(0, max(c(rh$density, rdens$y)))
        args$freq <- FALSE
        args$x <- as.numeric(x2)
        args <- delete.args("hist.default", args, package = "graphics", not = c("xlim",
          "ylim"))
        if (is.null(args$xlab))
          args$xlab <- if (is.null(type))
          "Residuals" else paste(type, "residuals")
        if (is.null(args$ylab))
          args$ylab <- "Density"
        if (is.null(args$main))
          args$main <- paste("Histogramm and density", if (!is.null(cn[j]))
          paste(":", cn[j]) else NULL)
        ok <- try(do.call("hist", args))
        if (!inherits(ok, "try-error"))
          lines(rdens)
        box()
      }
      if (w == "qq-resid") {
        if (ncol(x) > 1) {
          x2 <- t(apply(x, 1, c95))
          args$x <- NULL
          args$plot.it <- FALSE
          args <- delete.args("qqnorm.default", args, package = "stats",
          not = c("col", "pch", "cex"))
          if (is.null(args$main))
          args$main <- paste("Normal Q-Q plot", if (!is.null(cn[j]))
            paste(":", cn[j]) else NULL)
          args$y <- x2[, "Mean"]
          mean <- do.call(qqnorm, args)
          args$y <- x2[, "2.5%"]
          lower <- do.call(qqnorm, args)
          args$y <- x2[, "97.5%"]
          upper <- do.call(qqnorm, args)
          ylim <- range(c(as.numeric(mean$y), as.numeric(lower$y), as.numeric(upper$y)),
          na.rm = TRUE)
          args$plot.it <- TRUE
          if (is.null(args$ylim))
          args$ylim <- ylim
          args$y <- x2[, "Mean"]
          mean <- do.call(qqnorm, args)
          if (is.null(args$ci.col))
          args$ci.col <- 1
          if (is.null(args$ci.lty))
          args$ci.lty <- 2
          lines(lower$x[order(lower$x)], lower$y[order(lower$x)], lty = args$ci.lty,
          col = args$ci.col)
          lines(upper$x[order(upper$x)], upper$y[order(upper$x)], lty = args$ci.lty,
          col = args$ci.col)
          args$y <- x2[, "Mean"]
          qqline(args$y)
        } else {
          args$y <- as.numeric(x)
          args$x <- NULL
          args <- delete.args("qqnorm.default", args, package = "stats",
          not = c("col", "pch", "xlim", "ylim", "cex"))
          if (is.null(args$main))
          args$main <- paste("Normal Q-Q plot", if (!is.null(cn[j]))
            paste(":", cn[j]) else NULL)
          args$plot.it <- !add
          ok <- try(do.call(qqnorm, args))
          if (add) {
          args <- delete.args("points.default", list(...), package = "graphics",
            not = c("col", "pch", "cex"))
          points(ok$x, ok$y, pch = args$pch, col = args$col, cex = args$cex)
          } else {
          if (!inherits(ok, "try-error"))
            qqline(args$y)
          }
        }
      }
      if (w == "wp") {
        xlo <- xup <- NULL
        if (ncol(x) > 1) {
          x2 <- t(apply(x, 1, c95))
          xlo <- x2[, "2.5%"]
          xup <- x2[, "97.5%"]
          x2 <- x2[, "Mean"]
        } else {
          x2 <- x
        }
        d <- qqnorm(x2, plot = FALSE)
        probs <- c(0.25, 0.75)
        y3 <- quantile(x2, probs, type = 7, na.rm = TRUE)
        x3 <- qnorm(probs)
        slope <- diff(y3)/diff(x3)
        int <- y3[1L] - slope * x3[1L]
        d$y <- d$y - (int + slope * d$x)
        if (!is.null(xlo)) {
          d2 <- qqnorm(xlo, plot = FALSE)
          d$ylo <- d2$y - d2$x
          d$xlo <- d2$x
          d2 <- qqnorm(xup, plot = FALSE)
          d$yup <- d2$y - d2$x
          d$xup <- d2$x
        }
        level <- 0.95
        xlim <- max(abs(d$x), na.rm = TRUE)
        xlim <- c(-xlim, xlim)
        ylim <- max(abs(c(as.numeric(d$y), as.numeric(d$ylo), as.numeric(d$yup))),
          na.rm = TRUE)
        ylim <- c(-ylim, ylim)
        if (!is.null(args$ylim2))
          ylim <- args$ylim2
        if (!is.null(args$xlim2))
          xlim <- args$xlim2
        z <- seq(xlim[1] - 10, xlim[2] + 10, 0.25)
        p <- pnorm(z)
        se <- (1/dnorm(z)) * (sqrt(p * (1 - p)/length(d$y)))
        low <- qnorm((1 - level)/2) * se
        high <- -low
        args <- list(...)
        if (is.null(args$col))
          args$col <- 1
        if (is.null(args$pch))
          args$pch <- 1
        if (is.null(args$cex))
          args$cex <- 1
        if (is.null(args$ylab))
          args$ylab <- "Deviation"
        if (is.null(args$xlab))
          args$xlab <- "Unit normal quantile"
        if (add) {
          points(d$x, d$y, col = args$col, pch = args$pch, cex = args$cex)
        } else {
          if (is.null(args$main))
          args$main <- paste("Worm plot", if (!is.null(cn[j]))
            paste(":", cn[j]) else NULL)
          plot(d$x, d$y, ylab = args$ylab, xlab = args$xlab, main = args$main,
          xlim = xlim, ylim = ylim, col = NA, type = "n")
          grid(lty = "solid")
          abline(0, 0, lty = 2, col = "lightgray")
          abline(0, 1e+05, lty = 2, col = "lightgray")
          lines(z, low, lty = 2)
          lines(z, high, lty = 2)
          points(d$x, d$y, col = args$col, pch = args$pch, cex = args$cex)
        }
        if (!is.null(xlo)) {
          if (is.null(args$ci.col))
          args$ci.col <- 4
          if (is.null(args$ci.lty))
          args$ci.lty <- 2
          i <- order(d$xlo)
          lines(d$ylo[i] ~ d$xlo[i], lty = args$ci.lty, col = args$ci.col)
          i <- order(d$xup)
          lines(d$yup[i] ~ d$xup[i], lty = args$ci.lty, col = args$ci.col)
        }
      }
    }
  }
  return(invisible(NULL))
}



delete.args <- function(fun = NULL, args = NULL, not = NULL, package = NULL) {
  if (is.character(fun) & !is.null(package))
    fun <- eval(parse(text = paste(package, paste(rep(":", 3), collapse = ""),
      fun, sep = "")))
  nf <- names(formals(fun))
  na <- names(args)
  for (elmt in na) if (!elmt %in% nf) {
    if (!is.null(not)) {
      if (!elmt %in% not)
        args[elmt] <- NULL
    } else args[elmt] <- NULL
  }
  return(args)
}


rps <- function(obs, pred) {
  nr <- nrow(pred)
  test <- apply(pred, 1, sum)
  id <- is.finite(obs) & is.finite(test)  # & test <= 1 & test > 0.8
  obs <- obs[id]
  pred <- matrix(pred[id, ], nrow = nr)
  OBS <- matrix(0, nrow = length(obs), ncol = ncol(pred))
  k = nrow(OBS)
  k1 = ncol(OBS)
  k2 = ncol(pred)
  for (i in 1:nrow(OBS)) {
    OBS[i, 1 + obs[i]] <- 1
    # cat('\r');cat('i =',i, '-', k, '/// 0 /',k1 ,'/// 0 /',k2 )
  }
  OBS2 <- OBS
  for (i in 1:ncol(OBS)) {
    OBS2[, i] <- apply(matrix(OBS[, 1:i], nrow = nr), 1, sum)
    cat("\r")
    cat(k, "-", k, "/// i =", i, "-", k1, "/// 0 /", k2)
  }
  PRED <- OBS
  for (i in 1:ncol(pred)) {
    PRED[, i] <- apply(matrix(pred[, 1:i], nrow = nr), 1, sum)
    cat("\r")
    cat(k, "-", k, "///", k1, " / ", k1, "///", i, "-", k2)
  }
  RPS <- mean(apply((PRED - OBS2)^2, 1, sum))/(ncol(pred) - 1)

  return(RPS)
}


ma <- function(x, order = 20) {
  ma1 <- filter(x, rep(1/order, order), sides = 1)
  ma2 <- rev(filter(rev(x), rep(1/order, order), sides = 1))
  ma3 <- ifelse(is.na(ma1), ma2, ma1)
  ma3
}


# formula updating after variable selection step
f.update <- function(model, mstop, bb = 0, names = NULL, max.char = 24) {
  l <- length(model[["coefficients"]])
  f <- list()
  for (i in 1:l) {
    dd <- 2:ncol(model[["coefficients"]][[i]])
    if (mstop > bb)
      int <- (mstop - bb):mstop
    if (mstop <= bb)
      int <- mstop:(mstop + bb)

    coef <- t(colMeans(model[["coefficients"]][[i]][int, dd]))

    x <- colnames(coef)[abs(coef) > 0]
    nn <- c(1:max.char, "inf_cin", "cin", "cloud", "no_cloud", paste0(0, 1:max.char))
    if (!is.null(names))
      for (j in names) {
        if (any(paste0(j, nn) %in% x))
          x <- c(setdiff(x, paste0(j, nn)), j)
      }

    if (length(x) == 0)
      x = "1"

    f[i] <- ifelse(i == 1, paste(names(model[["y"]]), " ~ ", paste(x, collapse = "+")),
      paste(" ~ ", paste(x, collapse = "+")))

  }

  f1 <- lapply(f, function(i) as.formula(i))
  return(f1)
}


f2.update <- function(model, mstop, bb = 0) {
  dd <- 2:ncol(model[["coefficients"]][["mu"]])
  mu.coef <- t(colMeans(model[["coefficients"]][["mu"]][mstop - bb:mstop, dd]))
  sigma.coef <- t(colMeans(model[["coefficients"]][["sigma"]][mstop - bb:mstop,
    dd]))
  # nu.coef <- t(colMeans(model[['coefficients']][['nu']][mstop-bb:mstop,dd]))

  mu.x <- colnames(mu.coef)[abs(mu.coef) > 0]
  sigma.x <- colnames(sigma.coef)[abs(sigma.coef) > 0]
  # nu.x <- colnames(mu.coef)[abs(nu.coef) > 0]

  if (length(mu.x) == 0)
    mu.x = "1"
  if (length(sigma.x) == 0)
    sigma.x = "1"
  # if(length(nu.x) == 0) nu.x = '1'

  f.mu <- as.formula(paste("y ~ ", paste(mu.x, collapse = "+")))
  f.sigma <- as.formula(paste(" ~ ", paste(sigma.x, collapse = "+")))
  # f.nu <- as.formula(paste(' ~ ', paste(nu.x, collapse = '+')))

  f <- list(f.mu, f.sigma)
  return(f)
}


f3.update <- function(model, mstop, bb = 0) {
  dd <- 2:ncol(model[["coefficients"]][["mu"]])
  mu.coef <- t(colMeans(model[["coefficients"]][["mu"]][mstop - bb:mstop, dd]))
  sigma.coef <- t(colMeans(model[["coefficients"]][["sigma"]][mstop - bb:mstop,
    dd]))
  nu.coef <- t(colMeans(model[["coefficients"]][["nu"]][mstop - bb:mstop, dd]))

  mu.x <- colnames(mu.coef)[abs(mu.coef) > 0]
  sigma.x <- colnames(sigma.coef)[abs(sigma.coef) > 0]
  nu.x <- colnames(nu.coef)[abs(nu.coef) > 0]

  if (length(mu.x) == 0)
    mu.x = "1"
  if (length(sigma.x) == 0)
    sigma.x = "1"
  if (length(nu.x) == 0)
    nu.x = "1"

  f.mu <- as.formula(paste("y ~ ", paste(mu.x, collapse = "+")))
  f.sigma <- as.formula(paste(" ~ ", paste(sigma.x, collapse = "+")))
  f.nu <- as.formula(paste(" ~ ", paste(nu.x, collapse = "+")))

  f <- list(f.mu, f.sigma, f.nu)
  return(f)
}

f4.update <- function(model, mstop, bb = 0) {
  dd <- 2:ncol(model[["coefficients"]][["mu"]])
  mu.coef <- t(colMeans(model[["coefficients"]][["mu"]][mstop - bb:mstop, dd]))
  sigma.coef <- t(colMeans(model[["coefficients"]][["sigma"]][mstop - bb:mstop,
    dd]))
  nu.coef <- t(colMeans(model[["coefficients"]][["nu"]][mstop - bb:mstop, dd]))
  tau.coef <- t(colMeans(model[["coefficients"]][["tau"]][mstop - bb:mstop, dd]))

  mu.x <- colnames(mu.coef)[abs(mu.coef) > 0]
  sigma.x <- colnames(sigma.coef)[abs(sigma.coef) > 0]
  nu.x <- colnames(nu.coef)[abs(nu.coef) > 0]
  tau.x <- colnames(tau.coef)[abs(tau.coef) > 0]

  if (length(mu.x) == 0)
    mu.x = "1"
  if (length(sigma.x) == 0)
    sigma.x = "1"
  if (length(nu.x) == 0)
    nu.x = "1"
  if (length(tau.x) == 0)
    tau.x = "1"

  f.mu <- as.formula(paste("y ~ ", paste(mu.x, collapse = "+")))
  f.sigma <- as.formula(paste(" ~ ", paste(sigma.x, collapse = "+")))
  f.nu <- as.formula(paste(" ~ ", paste(nu.x, collapse = "+")))
  f.tau <- as.formula(paste(" ~ ", paste(tau.x, collapse = "+")))

  f <- list(f.mu, f.sigma, f.nu, f.tau)
  return(f)
}


rootogram <- function(model, newdata = NULL, counts = NULL, maxit, bb = 20, max.k = NULL,
  max.ylim = NULL, main = NULL, score = F) {

  range <- c(maxit - bb, maxit)
  print(range)

  if (is.null(newdata)) {
    newdata <- model[["model"]]
    newdata <- data.frame(newdata, model$y)
  }

  pred <- predict(object = model, newdata = newdata, type = "parameter", mstart = range[1],
    mstop = range[2])


  mat <- model[["family"]]$d(0, pred)
  k = ifelse(is.null(max.k), max(30, newdata$y), max.k)
  for (i in 1:k) {
    p <- model[["family"]]$d(i, pred)
    mat <- rbind(mat, p)
    cat("\r")
    cat("i =", i, "/", k)
  }
  if (score) {
    rankprob <- rps(counts, t(mat))
    print(paste0("rps = ", rankprob))
  }

  mat <- rowSums(mat)
  mat <- sqrt(mat)
  y <- rep(0, k + 1)

  if (is.null(counts))
    for (i in 1:(k + 1)) y[i] <- sqrt(sum(newdata$y == i - 1))
  if (!is.null(counts))
    for (i in 1:(k + 1)) y[i] <- sqrt(sum(counts == i - 1))

  diff <- mat - y
  names(diff) <- names(mat) <- NULL

  ylim = c(min(diff, 0), max(y, mat))
  if (!is.null(max.ylim))
    ylim[2] <- max.ylim
  barplot(diff ~ c(0:k), ylab = expression(sqrt(frequency)), ylim = ylim, main = main,
    xlab = "flash counts", space = 0)
  x <- c(0:k) + 0.5
  lines(y ~ x, col = "red", type = "l")
  lines(mat ~ x, col = "blue", type = "l")
  legend("topright", horiz = F, bty = "n", legend = c(expression(sqrt(pred)), expression(sqrt(obs)),
    expression(sqrt(pred) - sqrt(obs))), cex = 0.9, col = c("blue", "red", "grey"),
    title = NULL, lwd = 1)


}





prob.stop <- function(model, cyclic = FALSE, nnoise) {
  iter.first.p <- NULL
  for (i in names(model[["coefficients"]])) {
    ifp <- which.max(rowSums(model[["coefficients"]][[i]][, paste0("p", 1:(6 +
      nnoise))] != 0) != 0)
    iter.first.p <- c(iter.first.p, ifp)
  }
  r <- iter.first.p - 1
  if (!cyclic)
    r <- min(r)

  return(r)
}


prob.sample <- function(data, name = NULL) {
  nam <- colnames(data)
  if (!is.null(name)) {
    nam <- paste0(name, 1:length(nam))
  }
  data <- t(sample(data.frame(t(data))))
  rownames(data) <- NULL
  colnames(data) <- nam
  return(data)
}


