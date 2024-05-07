powerset <- function(set) {
  n <- length(set)
  masks <- 2 ^ (1:n - 1)
  lapply(1:2 ^ n - 1, function(u) set[bitwAnd(u, masks) != 0])
}

## Function to transform gamlss.family objects.
tF <- function(x, ...) {
  if (is.function(x))
    x <- x()
  if (!inherits(x, "gamlss.family"))
    stop('only "gamlss.family" objects can be transformed!')
  
  args    <- list(...)
  bd      <- if (is.null(args$bd)) 1 else args$bd
  args$bd <- NULL
  pr      <- args$range
  check_range <- function(par) {
    return(par)
  }
  if (!is.null(pr)) {
    if (is.list(pr) | is.data.frame(pr)) {
      check_range <- function(par) {
        for (j in names(par)) {
          if (!is.null(pr[[j]])) {
            if (is.numeric(pr[[j]])) {
              par[[j]][par[[j]] < min(pr[[j]])] <- min(pr[[j]])
              par[[j]][par[[j]] > max(pr[[j]])] <- max(pr[[j]])
            }
          }
        }
        return(par)
      }
    }
  }
  nx    <- names(x$parameters)
  score <- hess <- initialize <- list()
  
  make_call <- function(fun) {
    fn <- deparse(substitute(fun),
                  backtick = TRUE,
                  width.cutoff = 500)
    nf <- names(formals(fun))
    if (length(nf) < 1) {
      call <- paste(fn, "()", sep = "")
    } else {
      call <- paste(fn, "(", if ("y" %in% nf) "y," else "", sep = "")
      np <- nx[nx %in% nf]
      call <- paste(call,
                    paste(np, '=', 'par$', np, sep = '', collapse = ','),
                    sep = "")
      if ("bd" %in% nf) {
        call <- paste(call, ",bd=", bd, sep = "")
      }
    }
    call <- parse(text = paste(call, ")", sep = ""))
    return(call)
  }
  
  if ("mu" %in% nx) {
    mu.link <- make.link2(x$mu.link)
    mu.cs   <- make_call(x$dldm)
    mu.hs   <- make_call(x$d2ldm2)

    score$mu  <- function(y, par, ...) {
      par <- check_range(par)
      res <- eval(mu.cs) * mu.link$mu.eta(mu.link$linkfun(par$mu))
      if (!is.null(dim(res)))
        res <- res[, 1]
      res
    }

    hess$mu <- function(y, par, ...) {
      par <- check_range(par)
      score <- eval(mu.cs)
      hess <- -1 * eval(mu.hs)
      eta <- mu.link$linkfun(par$mu)
      res <-
        drop(score * mu.link$mu.eta2(eta) + hess * mu.link$mu.eta(eta) ^ 2)
      if (!is.null(dim(res)))
        res <- res[, 1]
      res
    }

    if (!is.null(x$mu.initial)) {
      initialize$mu <- function(y, ...) {
        if (!is.null(attr(y, "contrasts"))) {
          if (!is.null(dim(y))) y <- y[, ncol(y)]
        }
        if (!is.null(bd)) {
          bd <- rep(bd, length.out = if (!is.null(dim(y))) nrow(y) else length(y))
        }
        res <- eval(x$mu.initial)
        if (!is.null(dim(res))) res <- res[, 1]
        return(res)
      }
    }
  }
  
  if ("sigma" %in% nx) {
    sigma.link <- make.link2(x$sigma.link)
    sigma.cs  <- make_call(x$dldd)
    sigma.hs  <- make_call(x$d2ldd2)

    score$sigma  <- function(y, par, ...) {
      par <- check_range(par)
      res <-
        eval(sigma.cs) * sigma.link$mu.eta(sigma.link$linkfun(par$sigma))
      if (!is.null(dim(res)))
        res <- res[, 1]
      return(res)
    }

    hess$sigma <- function(y, par, ...) {
      par <- check_range(par)
      score <- eval(sigma.cs)
      hess <- -1 * eval(sigma.hs)
      eta <- sigma.link$linkfun(par$sigma)
      res <- drop(score * sigma.link$mu.eta2(eta) + hess * sigma.link$mu.eta(eta)^2)
      if (!is.null(dim(res))) res <- res[, 1]
      return(res)
    }

    if (!is.null(x$sigma.initial)) {
      initialize$sigma <- function(y, ...) {
        if (!is.null(bd)) {
          bd <- rep(bd, length.out = if (!is.null(dim(y))) nrow(y) else length(y))
        }
        res <- eval(x$sigma.initial)
        if (!is.null(dim(res))) res <- res[, 1]
        return(res)
      }
    }
  }
  
  if ("nu" %in% nx) {
    nu.link <- make.link2(x$nu.link)
    nu.cs   <- make_call(x$dldv)
    nu.hs   <- make_call(x$d2ldv2)

    score$nu <- function(y, par, ...) {
      par <- check_range(par)
      res <- eval(nu.cs) * nu.link$mu.eta(nu.link$linkfun(par$nu))
      if (!is.null(dim(res))) res <- res[, 1]
      return(res)
    }

    hess$nu <- function(y, par, ...) {
      par   <- check_range(par)
      score <- eval(nu.cs)
      hess  <- -1 * eval(nu.hs)
      eta   <- nu.link$linkfun(par$nu)
      res   <- drop(score * nu.link$mu.eta2(eta) + hess * nu.link$mu.eta(eta) ^ 2)
      if (!is.null(dim(res))) res <- res[, 1]
      return(res)
    }

    if (!is.null(x$nu.initial)) {
      initialize$nu <- function(y, ...) {
        if (!is.null(bd)) {
          bd <- rep(bd, length.out = if (!is.null(dim(y))) nrow(y) else length(y))
        }
        res <- eval(x$nu.initial)
        if (!is.null(dim(res))) res <- res[, 1]
        return(res)
      }
    }
  }
  
  if ("tau" %in% nx) {
    tau.link <- make.link2(x$tau.link)
    tau.cs   <- make_call(x$dldt)
    tau.hs   <- make_call(x$d2ldt2)

    score$tau <- function(y, par, ...) {
      par <- check_range(par)
      res <- eval(tau.cs) * tau.link$mu.eta(tau.link$linkfun(par$tau))
      if (!is.null(dim(res))) res <- res[, 1]
      return(res)
    }

    hess$tau <- function(y, par, ...) {
      par   <- check_range(par)
      score <- eval(tau.cs)
      hess  <- -1 * eval(tau.hs)
      eta   <- tau.link$linkfun(par$tau)
      res   <- drop(score * tau.link$mu.eta2(eta) + hess * tau.link$mu.eta(eta) ^ 2)
      if (!is.null(dim(res))) res <- res[, 1]
      return(res)
    }

    if (!is.null(x$tau.initial)) {
      initialize$tau <- function(y, ...) {
        if (!is.null(bd)) {
          bd <- rep(bd, length.out = if (!is.null(dim(y))) nrow(y) else length(y))
        }
        res <- eval(x$tau.initial)
        if (!is.null(dim(res))) res <- res[, 1]
        return(res)
      }
    }
  }
  
  dfun <- get(paste("d", x$family[1], sep = ""))
  pfun <- try(get(paste("p", x$family[1], sep = "")), silent = TRUE)
  qfun <- try(get(paste("q", x$family[1], sep = "")), silent = TRUE)
  rfun <- try(get(paste("r", x$family[1], sep = "")), silent = TRUE)
  
  nf   <- names(formals(dfun))
  bdc  <- "bd" %in% nf
  
  dc <- parse(text = paste("dfun(y, ",
              paste(paste(nx, "par$", sep = "="), nx, sep = "", collapse = ","),
              ", log = log, ...",
              if (bdc) paste0(", bd = ", bd) else NULL, ")", sep = ""))

  pc <- parse(text = paste("pfun(q, ",
              paste(paste(nx, "par$", sep = "="), nx, sep = "", collapse = ","),
              ", log = log, ...",
              if (bdc) paste0(",bd = ", bd) else NULL, ")", sep = ""))

  qc <- parse(text = paste("qfun(p, ",
              paste(paste(nx, "par$", sep = "="), nx, sep = "", collapse = ","),
              ", log = log, ...",
              if (bdc) paste0(",bd = ", bd) else NULL, ")", sep = ""))

  rc <- parse(text = paste("rfun(n, ",
              paste(paste(nx, 'par$', sep = "="), nx, sep = "", collapse = ","),
              ", ...",
              if (bdc) paste0(",bd=", bd) else NULL, ")", sep = ""))
  
  rval <- list(
    "family" = x$family[1],
    "names"  = nx,
    "links"  = unlist(x[paste(nx, "link", sep = ".")]),
    "score"  = score,
    "hess"   = hess,
    "d" = function(y, par, log = FALSE, ...) {
      par <- check_range(par)
      d <- try(eval(dc), silent = TRUE)
      if (inherits(d, "try-error"))
        d <- rep(NA, length(par[[1L]]))
      return(d)
    },
    "p" = if (!inherits(pfun, "try-error"))
      function(q, par, log = FALSE, ...) {
        par <- check_range(par)
        eval(pc)
      } else
        NULL,
    "q" = if (!inherits(qfun, "try-error"))
      function(p, par, log = FALSE, ...) {
        par <- check_range(par)
        eval(qc)
      } else
        NULL,
    "r" = if (!inherits(rfun, "try-error"))
      function(n, par, ...) {
        par <- check_range(par)
        eval(rc)
      } else
        NULL
  )

  names(rval$links)   <- nx
  rval$valid.response <- x$y.valid
  rval$initialize     <- initialize
  rval$type           <- tolower(x$type)
  
  if (!is.null(x$mean)) {
    meanc <- make_call(x$mean)
    rval$mean  <- function(par, ...) {
      par <- check_range(par)
      res <- eval(meanc)
      if (!is.null(dim(res))) res <- res[, 1]
      return(res)
    }
  } else {
    rval$mean <- function(par, ...) {
      return(par[[1L]])
    }
  }
  
  if (!is.null(x$variance)) {
    varc <- make_call(x$variance)
    rval$variance  <- function(par, ...) {
      par <- check_range(par)
      res <- eval(varc)
      if (!is.null(dim(res))) res <- res[, 1]
      return(res)
    }
  } else {
    rval$variance <- function(par, ...) {
      return(par[[2L]])
    }
  }
  
  class(rval) <- "family.bamlss"
  return(rval)
}

## Family completion.
complete.bamlss.family <- function(family) {

  if (is.null(names(family$links)))
    names(family$links) <- family$names
  
  linkinv <- linkfun <- list()
  for (j in family$names) {
    link <- make.link2(family$links[j])
    linkinv[[j]] <- link$linkinv
    linkfun[[j]] <- link$linkfun
  }
  
  if (is.null(family$map2par)) {
    family$map2par <- function(eta) {
      for (j in family$names) {
        eta[[j]] <- linkinv[[j]](eta[[j]])
        eta[[j]][is.na(eta[[j]])] <- 0
        if (any(jj <- eta[[j]] == Inf))  eta[[j]][jj] <- 10
        if (any(jj <- eta[[j]] == -Inf)) eta[[j]][jj] <- -10
      }
      return(eta)
    }
  }
  
  if (is.null(family$mu)) {
    family$mu <- function(par) {
        make.link2(family$links[1])$linkinv(par[[1]])
    }
  }
  
  if (is.null(family$loglik)) {
    if (!is.null(family$d)) {
      family$loglik <- function(y, par, ...) {
        logdens <- family$d(y, par, log = TRUE)
        if (any(i <- !is.finite(logdens))) {
          logdens[i] <- -100
        }
        return(sum(logdens, na.rm = TRUE))
      }
    }
  }
  
  err01 <- .Machine$double.eps^(1 / 3)
  err02 <- err01 * 2
  err11 <- .Machine$double.eps^(1 / 4)
  err12 <- err11 * 2
  
  if (is.null(family$score) & !is.null(family$d)) family$score <- list()

  for (i in family$names) {
    if (is.null(family$score[[i]]) & !is.null(family$d)) {
      fun <- c("function(y, par, ...) {",
               paste("  eta <- linkfun[['", i, "']](par[['", i, "']]);", sep = ""),
               paste("  par[['", i, "']] <- linkinv[['", i, "']](eta + err01);", sep = ""),
               "  d1 <- family$d(y, par, log = TRUE);",
               paste("  par[['", i, "']] <- linkinv[['", i, "']](eta - err01);", sep = ""),
               "  d2 <- family$d(y, par, log = TRUE);",
               "  return((d1 - d2) / err02)",
               "}")
      family$score[[i]] <- eval(parse(text = paste(fun, collapse = "")))
      attr(family$score[[i]], "dnum") <- TRUE
    }
  }
  
  if (is.null(family$hess) & !is.null(family$d)) family$hess <- list()

  for (i in family$names) {
    if (is.null(family$hess[[i]]) & !is.null(family$d)) {
      fun <- if (!is.null(attr(family$score[[i]], "dnum"))) {
        c("function(y, par, ...) {",
          paste0("  eta <- linkfun[['", i, "']](par[['", i, "']]);"),
          paste0("  par[['", i, "']] <- linkinv[['", i, "']](eta + err11);"),
          paste0("  d1 <- family$score[['", i, "']](y, par, ...);", sep = ""),
          paste0("  par[['", i, "']] <- linkinv[['", i, "']](eta - err11);"),
          paste0("  d2 <- family$score[['", i, "']](y, par, ...);"),
          "  return(-1 * (d1 - d2) / err12)",
          "}")
      } else {
        c("function(y, par, ...) {",
          paste0("  eta <- linkfun[['", i, "']](par[['", i, "']]);"),
          paste0("  par[['", i, "']] <- linkinv[['", i, "']](eta + err01);"),
          paste0("  d1 <- family$score[['", i, "']](y, par, ...);"),
          paste0("  par[['", i, "']] <- linkinv[['", i, "']](eta - err01);"),
          paste0("  d2 <- family$score[['", i, "']](y, par, ...);"),
          "  return(-1 * (d1 - d2) / err02)",
          "}")
      }
      family$hess[[i]] <- eval(parse(text = paste(fun, collapse = "")))
    }
  }
  
  return(family)
}

## Second make.link function.
make.link2 <- function(link) {
  if (is.null(link)) link <- "identity"

  link0 <- link
  if (link0 == "tanhalf") {
    rval <- list(
      "linkfun" = function (mu) {
        tan(mu / 2)
      },
      "linkinv" = function(eta) {
        2 * atan(eta)
      },
      "mu.eta" = function(eta) {
        2 / (eta^2 + 1)
      },
      "mu.eta2" = function(eta) {
        (-4 * eta) / (eta^2 + 1)^2
      },
      "valideta" = function(eta) {
        TRUE
      },
      "name" = "tanhalf"
    )
  } else {
    mu.eta2 <- function(x) {
      if (link0 == "identity") {
        x$mu.eta2 <- function(eta)
          rep.int(0, length(eta))
        return(x)
      }
      if (link0 == "log") {
        x$mu.eta2 <- function(eta)
          exp(eta)
        return(x)
      }
      if (link0 == "logit") {
        x$mu.eta2 <- function(eta) {
          eta <- exp(eta)
          return(-eta * (eta - 1) / (eta + 1)^3)
        }
        return(x)
      }
      if (link0 == "probit") {
        x$mu.eta2 <- function(eta) {
          -eta * dnorm(eta, mean = 0, sd = 1)
        }
        return(x)
      }
      if (link0 == "inverse") {
        x$mu.eta2 <- function(eta) {
          2 / (eta^3)
        }
        return(x)
      }
      if (link0 == "1/mu^2") {
        x$mu.eta2 <- function(eta) {
          0.75 / eta^(2.5)
        }
        return(x)
      }
      if (link0 == "sqrt") {
        x$mu.eta2 <- function(eta) {
          rep(2, length = length(eta))
        }
        return(x)
      }
      x$mu.eta2 <- function(eta)
        rep.int(0, length(eta))
      ## warning(paste('higher derivatives of link "', link, '" not available!', sep = ''))
      return(x)
    }
    
    if (link %in% c("logit", "probit", "cauchit", "cloglog", "identity", "log", "sqrt", "1/mu^2", "inverse")) {
      rval <- make.link(link)
    } else {
      rval <- switch(
        link,
        "rhogit" = list(
          "linkfun" = function(mu) {
            mu / sqrt(1 - mu^2)
          },
          "linkinv" = function(eta) {
            rval <- eta / sqrt(1 + eta^2)
            rval <- (abs(rval) - .Machine$double.eps) * sign(rval)
            rval
          },
          "mu.eta" = function(eta) {
            1 / (1 + eta^2)^1.5
          }
        ),
        "cloglog2" = list(
          "linkfun" = function(mu) {
            log(-log(mu))
          },
          "linkinv" = function(eta) {
            pmax(pmin(1 - expm1(-exp(eta)), .Machine$double.eps),
                 .Machine$double.eps)
          },
          "mu.eta" = function(eta) {
            eta <- pmin(eta, 700)
            pmax(-exp(eta) * exp(-exp(eta)), .Machine$double.eps)
          }
        ),
        "sigmoid" = list(
          "linkfun" = function(mu) {
            i <- mu <= -1
            if (any(i)) mu[i] <- mu[i] <- -0.9999
            i <- mu >= 1
            if (any(i)) mu[i] <- mu[i] <- +0.9999
            - log(2 / (mu + 1) - 1)
          },
          "linkinv" = function(eta) {
            tanh(eta / 2)
          },
          "mu.eta" = function(eta) {
            0.5 / cosh(eta * 0.5)^2
          },
          "mu.eta2" = function(eta) {
            eta2 <- eta * 0.5 - (0.5 * (2 * (sinh(eta2) * 0.5 * cosh(eta2))) / (cosh(eta2)^2)^2)
          }
        )
      )
    }
    
    rval <- mu.eta2(rval)
  }
  rval$name <- link

  return(rval)
}


alpha2cap2 <- function(alpha = 0.01, nnobs = NULL, nvars = NULL, mean = 0) {
    t <- (1 - alpha) ^ (1 / nvars)
    if(mean != 0){
    p <- function(cap) {
      pnorm(cap, mean = mean, sd = sqrt(nnobs)/(nnobs-1)) - pnorm(-cap, mean = mean, sd = sqrt(nnobs)/(nnobs-1))
    }
    vec <- 0:40000 / 1e+05
    s <- sapply(vec, FUN = function(i) abs(p(i) - t))
    cap <- vec[which.min(s)]
    } else {
      cap <- qnorm((1 + t) / 2) * sqrt(nnobs) / (nnobs - 1)
    }
    return(cap)
}


## Stagewise distributional regression (SDR).
sdr <- function(formula,
                family = NULL,
                data = NULL,
                batch_ids = NULL,
                updating = c("thresdesc", "bs", "cyclic", "noncyclic"),
                light = FALSE, # TRUE means no data matrix gets returned
                CF = TRUE,
                cap = NULL,
                caps = NULL,
                scalex = TRUE,
                refitting = TRUE,
                eps = 0.01,
                nu = 0.1,
                ncaps = 10,
                quick_ffdf = FALSE,
                ...) {

    stopifnot(requireNamespace("Formula"))
    stopifnot(requireNamespace("Matrix"))

    # Sanity check for input args

    # Checking and manipulating formula (if needed)
    # * If character of length 1 we assume it should be a formula and convert it.
    # * If list but elements contain character length 1, coerce to formula.
    if (is.character(formula) & length(formula) == 1) formula <- as.formula(formula)
    if (is.list(formula)) formula <- lapply(formula, function(x) if (is.character(x) & length(x) == 1) as.formula(x) else x)
    stopifnot("argument `formula` must be a formula or list of formulae objects" =
              inherits(formula, "formula") || (is.list(formula) && all(sapply(formula, function(x) inherits(x, "formula")))))

    stopifnot("argument `batch_ids` incorrect" =
              is.null(batch_ids) || is.numeric(batch_ids) || (is.list(batch_ids) && all(sapply(batch_ids, is.numeric))))
    updating <- match.arg(updating)
    stopifnot("argument `light` must be logical TRUE or FALSE"     = isTRUE(light) || isFALSE(light))
    stopifnot("argument `CF` must be logical TRUE or FALSE"        = isTRUE(CF) || isFALSE(CF))
    stopifnot("argument `cap` must be NULL or single numeric"      = is.null(cap) || (is.numeric(cap) && length(cap) == 1))
    stopifnot("argument `caps` must be NULL or numeric"            = is.null(caps) || is.numeric(cap))
    stopifnot("argument `scalex` must be logical TRUE or FALSE"    = isTRUE(scalex) || isFALSE(scalex))
    stopifnot("argument `refitting` must be logical TRUE or FALSE" = isTRUE(refitting) || isFALSE(refitting))

    stopifnot("argument `eps` must be positive numeric of length 1" = is.numeric(eps) && length(eps) == 1 && eps > 0)
    stopifnot("argument `nu` must be positive numeric of length 1"  = is.numeric(nu) && length(nu) == 1 && nu > 0)
    stopifnot("argument `ncaps` must be positive numeric of length 1"  = is.numeric(ncaps) && length(ncaps) == 1 && ncaps > 0)
    stopifnot("argument `quick_ffdf` must be logical TRUE or FALSE"    = isTRUE(quick_ffdf) || isFALSE(quick_ffdf))

    if (quick_ffdf)                         scalex <- FALSE

    # Preparing family object
    if (is.null(family))                    family <- gamlss.dist::NO
    if (is.function(family))                family <- family()
    if (inherits(family, "gamlss.family"))  family <- tF(family)
    family <- complete.bamlss.family(family)

    # Ensure formula is a list containing a formula for each parameter of the family
    if (!is.list(formula)) formula <- list(formula)
    if (length(formula) < length(family$names)) {
      for (j in (length(formula) + 1):length(family$names))
        formula[[j]] <- ~ 1
    }
    formula <- rep(formula, length = length(family$names))
    names(formula) <- family$names
    
    # Starting to prepare the objects for estimation
    mfd <- list(varnames = list(), formula = list())
    terms <- NULL
    # list of model formula
    # Intercept only
    for (i in names(formula)) {
      fi <- as.Formula(formula[[i]])
      mfd$varnames[[i]] <-  unname(attr(terms(fi), "term.labels"))
      if(!length(mfd$varnames[[i]]) > 0){
        mfd$formula[[i]] <- ~ 1
      } else {
        mfd$formula[[i]] <- as.formula(paste(" ~ ", paste0(mfd$varnames[[i]], collapse = "+")))
      } 
    }

    mfd$y <- y <- all.vars(formula(as.Formula(formula[[1]]),
                                   lhs = 1, rhs = 0, drop = TRUE))

    # All variables that are needed
    vars <- unique(unlist(mfd$varnames))
    if(!length(vars) > 0){
      formula_all <- ~ 1
    } else {
      formula_all <- as.formula(paste(" ~ ", paste0(vars, collapse = "+")))
    }
    
    if(inherits(data, "big.matrix")){
      quick_ffdf <- TRUE
      light <- TRUE
      scalex <- FALSE
    }

    # delete bm.csv on exit
    # quick_ffdf means data is in format of modelframe
    if(!quick_ffdf){
      if (inherits(data, "ffdf") ){
        stopifnot(requireNamespace("ff"))
        outfile <- tempfile(fileext = ".csv")
        
        prepare_bigmem <- function(infile, outfile, formula, BATCHSIZE = 10000, overwrite = FALSE) {
          if (!overwrite & file.exists(outfile)) {
            stop("output file \"", outfile, "\" already exists (overwrite = FALSE)")
          } else if (file.exists(outfile)) {
            file.remove(outfile)
          }
          counter <- 0
          fun <- function(x, f, file) {
            x <- cbind(x[,y], model.matrix(object = f, data = x))
            colnames(x)[1:length(y)] <- y
            if (counter == 0) {
              write.table(x, file = file, col.names = TRUE, row.names = FALSE, sep = ",")
            } else {
              write.table(x, file = file, col.names = FALSE, row.names = FALSE, sep = ",", append = TRUE)
            }
            counter <<- counter + nrow(x)
            cat(counter, "\r")
            
          }
          ff::ffrowapply(fun(infile[i1:i2, ], f = formula, file = outfile), X = infile, BATCHSIZE = BATCHSIZE)
          
          return(ff::read.csv.ffdf(file = outfile))
        }
        
        data <- prepare_bigmem(infile = data, outfile = outfile, f = formula_all, overwrite = TRUE)
        #TODO: read.csv.ffdf(file = outfile)#read.big.matrix(outfile, header = TRUE, type = "double")#prepare_bigmem(infile = data, outfile = outfile, f = formula_all, overwrite = TRUE)
      } else {
        if(is.vector(data)) {
          data <- as.matrix(data)
          colnames(data) <- y
        }
        data <- cbind(data[,y], model.matrix(data = data.frame(data), object = formula_all))
        colnames(data)[1:length(y)] <- y
        data <- as.matrix(data)
      }
    }
    
    if(!light & !inherits(data, "ffdf")){
      mfd$y <- as.matrix(data[,y])
      colnames(mfd$y) <- y
      mfd$X <- data[,vars]
    }

    x_mu <- x_sd <- NULL
    if(length(vars) > 0 & scalex){
      # Scaling
      if (is.matrix(data)) {
        if(length(vars) > 1){
          x_mu <- apply(data[,vars], 2, mean)
          x_sd <- apply(data[,vars], 2, sd)
        } else {
          x_mu <- mean(data[,vars])
          x_sd <- sd(data[,vars])
        }  
      } else if (inherits(data, "ffdf")) {
        # Now we need to 'batch' trough the big.matrix and calculate the mean
        bm_mu <- function(x, BATCHSIZE = 20000) {
          i1   <- seq(1, nrow(x), by = BATCHSIZE) # Starting index
          rval <- lapply(i1, function(i) {
            ii <- seq.int(i, pmin(i + BATCHSIZE - 1, nrow(x))) # Indices
            if(length(vars > 1)) {
              ret <- apply(x[ii, vars], 2, sum)
            } else {
              ret <- sum(x[ii, vars])
            }
            return(ret)
          })
          return(colSums(do.call(rbind, rval)) / nrow(x))
        }
        
        bm_sd <- function(x, mu, BATCHSIZE = 20000) {
          i1 <- seq(1, nrow(x), by = BATCHSIZE) # Starting index
          fn <- function(i) {
            ii <- seq.int(i, pmin(i + BATCHSIZE - 1, nrow(x))) # Indices
            z  <- x[ii, vars]
            return(sapply(colnames(z), function(n) sum((z[, n] - mu[[n]])^2)))
          }
          rval <- lapply(i1, fn)
          return(sqrt(colSums(do.call(rbind, rval)) / nrow(x)))
        }
        
        x_mu <- bm_mu(data)
        x_sd <- bm_sd(data, mu = x_mu)
      }
    } 
    
    if (is.matrix(data) | inherits(data, "ffdf") | inherits(data, "big.matrix")) {
      ndata <- nrow(data)
    } else {
      ndata <- 1
    }

    # nobs = batchsize
    if (is.null(batch_ids)){
      nobs <- ndata
    } else if (is.numeric(batch_ids)){
      nobs <- batch_ids
    } else if (is.list(batch_ids)){
      nobs <- length(batch_ids[[1]])
    }
    
    # max nvars in dist parameters
    nvars <- max(sapply(names(formula), FUN = function(i) length(mfd$varnames[[i]])))
    
    # Correlation Filtering
    # todo: adapt for different nobs in parameter
    if (CF & nvars != 0) {
      if (is.null(cap)) {
        cap <- alpha2cap2(alpha = 0.05, nnobs = nobs, nvars = nvars, mean = 0)
        cap <- min(max(cap, 0.075), 0.175) # Limiting cap to [0.075 - 0.175]
      }
    } else {
      cap <- 0
    }
    
    if (updating == "thresdesc") {
      if (is.null(caps)) {
        if(nobs < 1000 )                      caps <- seq(max(0.15, min(cap + 0.15, 0.4)), min(0.1,max(cap - 0.15, 0.03)), length.out = ncaps)
        else if(nobs >= 1000 & nobs <= 10000) caps <- seq(max(0.15, min(cap + 0.08, 0.4)), min(0.1,max(cap - 0.12, 0.03)), length.out = ncaps)
        else                                  caps <- seq(max(0.15, min(cap + 0.10, 0.4)), min(0.1,max(cap - 0.10, 0.03)), length.out = ncaps)
        cat("\n", "Proposed thresholds: ", round(caps, 4), "\n")
      }
      
      mfd <- c(mfd,
               sdr.thresdesc(
                 data = data,
                 y = y,
                 vardist = mfd$varnames,
                 family = family,
                 caps = caps,
                 batch_ids = batch_ids,
                 eps = eps,
                 nu = nu,
                 refitting = refitting,
                 x_mu = x_mu,
                 x_sd = x_sd,
                 nvars = nvars,
                 scalex = scalex,
                 vars = vars,
                 ...
               )
             ) # end mfd <- ...
    } else if (updating == "bs") {
      mfd$cap <- cap
      mfd <- c(mfd,
               sdr.gradboostfit(
                 data = data,
                 y = y,
                 vardist = mfd$varnames,
                 family = family,
                 cap = cap,
                 batch_ids = batch_ids,
                 eps = eps,
                 nu = nu,
                 refitting = refitting,
                 x_mu = x_mu,
                 x_sd = x_sd,
                 nvars = nvars,
                 scalex = scalex,
                 vars = vars,
                 ...
               )
            ) # end mfd <- ...
    } else if (updating == "noncyclic") {
      mfd$cap <- cap
      mfd <- c(mfd,
               sdr.gradboostfit(
                 data = data,
                 y = y,
                 vardist = mfd$varnames,
                 family = family,
                 cap = cap,
                 batch_ids = batch_ids,
                 eps = eps,
                 nu = nu,
                 refitting = refitting,
                 x_mu = x_mu,
                 x_sd = x_sd,
                 nvars = nvars,
                 scalex = scalex,
                 vars = vars,
                 length_ps = 1,
                 ...
               )
             ) # end mdf <- 
    } else if (updating == "cyclic") {
      mfd$cap <- cap
      mfd <- c(mfd,
               sdr.gradboostfit(
                 data = data,
                 y = y,
                 vardist = mfd$varnames,
                 family = family,
                 cap = cap,
                 batch_ids = batch_ids,
                 eps = eps,
                 nu = nu,
                 refitting = refitting,
                 x_mu = x_mu,
                 x_sd = x_sd,
                 nvars = nvars,
                 scalex = scalex,
                 vars = vars,
                 ...
               )
             ) # end mfd <- ...
    } else {
        stop("Whoops, we should never end up here.")
    }
    
    if (scalex) {
      for (i in names(formula)) {
        nn <- names(mfd[["coefficients"]][[i]][1, ])
        if (length(nn) > 1) {
          # rescaling
          range <- mfd$varnames[[i]]
          mfd[["coefficients"]][[i]][, range] <- t(t(mfd[["coefficients"]][[i]][, range]) * (1 / x_sd[range]))
          
          for(j in range){
            mfd[["coefficients"]][[i]][, "(Intercept)"] <- mfd[["coefficients"]][[i]][, "(Intercept)"] - mfd[["coefficients"]][[i]][, j] * x_mu[j]
          }  
        }
      }
    } # end scalex
    
    
    mfd$family <- family
    mfd$nobs <- nobs
    
    # If outfile was created, delete it
    if (exists("outfile") && file.exists(outfile)) unlink(outfile)

  # Delete formula environment or else function environment gets copied 
  # into fomula at end of function
  for(x in 1:length(mfd$formula)){
    environment(mfd$formula[[x]]) <- emptyenv()
  }
  
  # Adding class to mfd and return
  class(mfd) <- "sdr"
  mfd$call <- match.call()                      
  
  return(mfd)
}



# best subset gradboosting with correlation filtering
# noncyclic is a special case with length_ps = 1
sdr.thresdesc <- function(data,
                          family,
                          maxit = 200,
                          refitting = TRUE,
                          maxit_refit = 2*maxit,
                          vardist = NULL,
                          eps = 0.01,
                          nu = 0.1,
                          #eps_int = exp(seq(log(0.1), log(0.001), length = ifelse(refitting, maxit+maxit_refit, maxit))),
                          eps_int = c(exp(seq(log(0.1), log(0.001), length = maxit)), rep(0.001, maxit_refit)),
                          nu_int = 0.05,
                          oos_batch = "next",
                          aic = FALSE,
                          K = 0,
                          batch_ids = NULL,
                          verbose = TRUE,
                          initialize = TRUE,
                          plot = FALSE,
                          coef_start = NULL,
                          length_ps = NULL,
                          cap = 0.1,
                          x_mu = NULL,
                          x_sd = NULL,
                          nvars = NULL,
                          scalex = TRUE,
                          y = NULL,
                          vars = NULL,
                          caps = seq(0.2,0.05, length.out = 10),
                          ...) {

    ia <- interactive()
    if (!refitting) maxit_refit <- 0
    
    # scaling is done on each respective batch separately
    myscale <- function(x, x_mu, x_sd) {
      for (n in colnames(x)) x[, n] <- (x[, n] - x_mu[[n]]) / x_sd[[n]]
      return(x)
    }
    
    maxit1  <- maxit + maxit_refit
    eps     <- rep(eps, length.out = maxit1)
    eps_int <- rep(eps_int, length.out = maxit1)
    
    N       <- nrow(data)
    nx      <- names(vardist)
    ind     <- 1:N
    
    if (is.null(batch_ids)) {
      ind       <- 1:N
      b.size    <- N
      batch_ids <- lapply(1:maxit1, function(...) ind)
    } else if (is.numeric(batch_ids)) {
      ind       <- 1:N
      b.size    <- batch_ids
      batch_ids <- lapply(1:maxit1, function(...) sample(ind, size = batch_ids, replace = FALSE))
    } # else batch_ids is a list of indices
    
    if (!is.list(batch_ids))
      stop("Argument batch_ids must be a list of indices!")
    if (length(batch_ids) != maxit1)
      warning("Length of batch_ids != maxit+maxit_refit, using batch_ids for setting maxit+maxit_refit!")
    maxit1 <- length(batch_ids)

    # number of plots
    tw <- length(strsplit(as.character(maxit1), "")[[1]]) + 1L
    if(!exists("b.size")) b.size <- length(batch_ids[[1]])
    # K is penalty for IC
    if (is.null(K)) K <- log(b.size)
    
    beta <- list()
    intvar <- list()
    for (j in nx) {
      if(length(vardist[[j]]) == 0){
        intvar[[j]] <- "(Intercept)"
        ncolu <- 1
      } else {
        intvar[[j]] <- c("(Intercept)",vardist[[j]])
        ncolu <- length(intvar[[j]])
      }
      beta[[j]] <- matrix(0, nrow = maxit1, ncol = ncolu)
      colnames(beta[[j]]) <- intvar[[j]]
    }
    if (!is.null(coef_start)) {
      initialize <- FALSE
      for (j in nx) {
        beta[[j]][1, ] <- coef_start[[j]]
      }
    }
    
    #beta.grad <- beta
    powerset.list <- powerset(nx)[-1]
    if (!is.null(length_ps))
      powerset.list <- powerset.list[sapply(powerset.list, length) == length_ps]

    # used for formatting cat() output if verbose = TRUE
    ll <- which.max(sapply(powerset.list, length))
    pset_fmt <- paste0("%-", max(sapply(powerset.list[ll], function(x) nchar(paste(x, collapse = ", ")))), "s")                                      
    
    if (initialize) {
      if (!is.null(family$initialize)) {
        betai <- list()
        for (j in nx) {
          if (!is.null(family$initialize[[j]])) {
            linkfun <- make.link2(family$links[j])$linkfun
            beta[[j]][1L, "(Intercept)"] <- mean(linkfun(family$initialize[[j]](data[batch_ids[[1]], y], )), na.rm = TRUE)
            
          }
        }
      }
    }
    df    <- sum(sapply(beta, function(b) sum(b[1, ] != 0)))
    err01 <- .Machine$double.eps ^ (1 / 2)
    err02 <- err01 * 2 * b.size
    # err02 is the denominator for central numeric differentiation.  b.size makes
    # gradient magnitute for different sample sizes comparable this translates
    # maximum log likelihood to maximum average loglikelihood
    # https://stats.stackexchange.com/questions/267847/motivation-for-average-log-likelihood
    
    ma <- function(x, order = 20) {
      ma1 <- filter(x, rep(1 / order, order), sides = 1)
      ma2 <- rev(filter(rev(x), rep(1 / order, order), sides = 1))
      return(ifelse(is.na(ma1), ma2, ma1))
    }
    
    bs_fun <- function(coef_list, yi, X, yoos, Xoos, nu = 0.1, nu_int = 0.05, eps = 0.1, eps_int = 0.1,
                       intvar, vardist, b.size = length(X[[1]]), cap = 0, i = 1) {
      
      coef_list_new <- coef_list
      eta <- etaoos <- list()
      
      for (j in nx) {
        ## Setup linear predictor.
        etaoos[[j]] <- drop(Xoos[[j]] %*% coef_list_new[[j]])
        eta[[j]]    <- drop(X[[j]] %*% coef_list_new[[j]])
      }
      
      lloos <- family$loglik(yoos, family$map2par(etaoos))
      
      # best subset intercept updating
      sign.list <- setNames(rep(-Inf, length(nx)), nx)
      
      for (j in nx) {
       
        ## Get coefficients and setup.
        tbeta    <- coef_list_new[[j]]
        ## Positive.
        tbeta[1] <- tbeta[1] + err01
        eta[[j]] <- drop(X[[j]] %*% tbeta)
        ll1 <- family$loglik(yi, family$map2par(eta))
        
        ## Negative
        tbeta[1] <- tbeta[1] - 2 * err01
        eta[[j]] <- drop(X[[j]] %*% tbeta)
        ll2 <- family$loglik(yi, family$map2par(eta))
        
        grad <- (ll1 - ll2) / err02
        sign.list[[j]] <- grad
        
        tbeta <- coef_list_new[[j]]
      }
      
      pset <- powerset.list
      for (j in nx) {
        if (sign.list[[j]] == -Inf) {
          ok <- !grepl(j, pset)
          pset <- pset[ok]
        }
      }
      
      if (length(pset) != 0) {
        beta.wip <- beta.final <- list()
        for (j in nx){
          beta.wip[[j]] <- beta.final[[j]] <- coef_list_new[[j]]
        }  
        # ic.old is the old information criterion insample
        # without penalty for complexity, i.e. -2*logLik. K = 0 is default.
        ic.old <- -2 * lloos 
        
        for (l in 1:length(pset)) {
          # ps <- !is.na(pset[l,])
          sign.list2 <- sign.list[pset[[l]]]
          gnorm <- sqrt(sum(sign.list2 ^ 2))
          
          if (gnorm > eps_int) {
            sign.list2 <- eps_int * sign.list2 / gnorm
          }

          for (ij in pset[[l]]) {
            sign.list2[ij] <- ifelse(i < 0.8 * maxit1 & abs(sign.list2[ij]) < nu_int * eps_int,
                                     sign(sign.list2[ij]) * nu_int * eps_int,
                                     sign.list2[ij])
          }
          
          for (j in pset[[l]]) {
            grad <- sign.list2[[j]]
            coef_list_new[[j]][1] <- coef_list_new[[j]][1] + grad
          }
          
          
          ## keep update only if oos information crit improves
          for (j in nx)
            etaoos[[j]] <- drop(Xoos[[j]] %*% coef_list_new[[j]])

          ll <- family$loglik(yoos, family$map2par(etaoos))
          ic.new <- -2 * ll 
          if (ic.new < ic.old) {
            # keep current try
            for (j in nx)
              beta.final[[j]] <- coef_list_new[[j]]
            ps.final <- pset[[l]]
            ic.old   <- ic.new
          }
          for (j in pset[[l]])
            coef_list_new[[j]] <- beta.wip[[j]]
          
        } # this bracket is for pset
        
        # save best beta
        for (j in nx){
          coef_list_new[[j]] <- beta.final[[j]]
          eta[[j]]           <- drop(X[[j]] %*% coef_list_new[[j]])
        }
      }
      
      
      sign.list <- pos.list <- setNames(rep(-Inf, length(nx)), nx)
      
      for (j in nx) {
        ## Get coefficients and setup.
        nc <- length(intvar[[j]]) - 1
        
        if (nc > 0) {
          
          eta[[j]] <- eta[[j]] + err01
          ll1      <- family$d(yi, family$map2par(eta), log = TRUE)
          
          ## Negative
          eta[[j]] <- eta[[j]] - 2 * err01
          ll2      <- family$d(yi, family$map2par(eta), log = TRUE)
          
          grad     <- (ll1 - ll2) / (2 * err01)
          eta[[j]] <- eta[[j]] + err01
          
          cc <- try(cor(grad, X[[j]][, vardist[[j]], drop = FALSE]), F)
          #[,-1]))
          cc[is.na(cc)] <- 0
          if (is.numeric(cc)) {
            ## Select update
            jj <- which.max(abs(cc)) + 1
          } else {
            jj <- 1
          } # interecept if cor gives error
          
          
          if (max(abs(cc)) < cap) {
            grad <- 0
          } else {
            # average partial derivative with respect to best variable
            grad <- t(grad) %*% X[[j]][, intvar[[j]][jj], drop = FALSE] / b.size 
          }
          
          sign.list[[j]] <- grad
          pos.list[[j]]  <- jj # position in intvar
        }
      } # end j in nx loop
                                    
      ps.final <- "no par"
      pset     <- powerset.list
      for (j in nx) {
        if (sign.list[[j]] == 0 | sign.list[[j]] == -Inf) {
          ok   <- !grepl(j, pset)
          pset <- pset[ok]
        }
      }
      
      if (length(pset) != 0) {
        ## beta.final <- list()
        for (j in nx)
          beta.final[[j]] <- coef_list_new[[j]]
        #ic.old <- ic
        
        beta.wip <- list()
        for (j in nx)
          beta.wip[[j]] <- coef_list_new[[j]]
        for (l in 1:length(pset)) {
          # ps <- !is.na(pset[l,])
          sign.list2 <- sign.list[pset[[l]]]
          gnorm      <- sqrt(sum(sign.list2 ^ 2))
          
          if (gnorm > eps) {
            sign.list2 <- eps * sign.list2 / gnorm
          }
          for (ij in pset[[l]]) {
            sign.list2[ij] <- ifelse(i < 0.8 * maxit1 & abs(sign.list2[ij]) < nu * eps,
                                     sign(sign.list2[ij]) * nu * eps,
                                     sign.list2[ij])
          }
          
          for (j in pset[[l]]) {
            jj   <- pos.list[[j]]
            grad <- sign.list2[[j]]
            
            coef_list_new[[j]][intvar[[j]][jj]] <- coef_list_new[[j]][intvar[[j]][jj]] + grad
          }
          
          ## keep update only if oos information crit improves
          for (j in nx)
            etaoos[[j]] <- drop(Xoos[[j]] %*% coef_list_new[[j]])
          ll <- family$loglik(yoos, family$map2par(etaoos))
          ic.new <- -2 * ll 
          if (ic.new < ic.old) {
            # keep current try
            for (j in nx)
              beta.final[[j]] <- coef_list_new[[j]]
            ps.final <- pset[[l]]
            ic.old   <- ic.new
          }
          for (j in pset[[l]])
            coef_list_new[[j]] <- beta.wip[[j]]
          
        } # this bracket is for pset
        
        # save best coefs
        for (j in nx)
          coef_list_new[[j]] <- beta.final[[j]]
        
      }
      return(c(coef_list_new, "ic.old" = ic.old, "ps.final" = ps.final))
    }
    
    #caps <- seq(0.2, 0.05, length.out = 5) # c(0.2, 0.15, 0.1, 0.5)
    beta.list <- setNames(lapply(caps, FUN = function(x) beta), caps)
    intvar1   <- intvar
    vardist1  <- vardist
    ll.list   <- list()
    df.list   <- cs.list <- NULL
    #seed <- rnorm(1)
    #vardist = mfd$varnames

    for (c in 1:length(caps)) {
      ic.oos.list <- ic0.list   <- NULL
      ll0.list    <- lloos.list <- 0
      beta    <- beta.list[[paste(caps[c])]]
      cap     <- caps[c]
      intvar  <- intvar1
      vardist <- vardist1
      cap     <- c(rep(cap, length.out = maxit), rep(0, length.out = maxit_refit))
      
      # # reordering, helpfull for small data to avoid bad starting batches
      # batch_ids[1:maxit] <- batch_ids[sample(1:maxit)]
      # batch_ids[(maxit+1):maxit1] <- batch_ids[sample((maxit+1):maxit1)]
      
      for (i in 2:maxit1) {
        # deselect non selected variables for updating if maxit_refit > 0
        if (i == maxit +1) {
          for (j in nx) {
            cond <- beta[[j]][maxit,] != 0
            cond["(Intercept)"] <- TRUE
            intvar[[j]]  <- intvar[[j]][cond]
            vardist[[j]] <- vardist[[j]][cond[-1]]
            # beta[[j]][maxit+1,] <- beta[[j]][maxit,] <- beta.list[[paste(caps[1])]][[j]][1,]
          }
        }
        
        ## out of sample
        if (oos_batch == "next") {
          i.oos <- if (i != maxit1) i + 1 else 1
        } else if (oos_batch == "same") {
          i.oos <- i
        } else if (oos_batch == "random") {
          i.oos <- sample(maxit1, 1)
        }
        
        batch.oos <- batch_ids[[i.oos]]
        
        ## Extract response.
        yi    <- as.matrix(data[batch_ids[[i]],y])
        y.oos <- as.matrix(data[batch.oos,y])
        
        # draw batchwise from big.matrix and scale
        XX    <- as.matrix(data[batch_ids[[i]],vars])
        XXoos <- as.matrix(data[batch.oos,vars])
        if (scalex) {
          XX    <- myscale(XX, x_mu = x_mu, x_sd = x_sd)
          XXoos <- myscale(XXoos, x_mu = x_mu, x_sd = x_sd)
        }
        X <- Xoos <- list()
        for (j in nx) {
          ## Save last iteration.
          beta[[j]][i, ]      <- beta[[j]][i - 1L, ]
          X[[j]]              <- cbind( "(Intercept)" = rep(1,b.size), XX[,setdiff(colnames(beta[[j]]), "(Intercept)" )])
          Xoos[[j]]           <- cbind( "(Intercept)" = rep(1,b.size),XXoos[,setdiff(colnames(beta[[j]]), "(Intercept)" )])
          colnames(Xoos[[j]]) <- colnames(X[[j]]) <- c("(Intercept)", setdiff(colnames(beta[[j]]), "(Intercept)" ))
        }
       
        coef_list <- list()
        for (j in nx) {
          coef_list[[j]] <- beta[[j]][i-1, ]
        }
        
        coef_list_new <- bs_fun(coef_list = coef_list, yi = yi, X = X, yoos = y.oos, Xoos = Xoos,
                                nu = nu, nu_int = nu_int, eps = eps[i], eps_int = eps_int[i],
                                intvar = intvar, vardist = vardist, b.size = nrow(X[[1]]), cap = cap[i])
        
        for (j in nx) {
          beta[[j]][i, ] <- coef_list_new[[j]]
        }

        # Selected best subset
        ps.final <- unlist(coef_list_new[grep("ps.final", names(coef_list_new))])
        
        # #Nonzero coefs
        df <- sum(sapply(beta, function(b) sum(b[i, ] != 0)))
        # perhaps change to:
        # ic.old <- coef_list_new[["ic.old"]] + K * df
        
        # BIC = -2*loglik_batch + log(batch.size) * df              ####* (nrow(X[[1]])/N)
        ic.old     <- coef_list_new[["ic.old"]] + K * df #* (nrow(X[[1]])/N)
        lloos.list <- c(lloos.list,  ic.old)
        ll0.list   <- c(lloos.list,  -coef_list_new[["ic.old"]]/2)
        
        if (verbose) {
          if (ia) cat("\r")
          cat("thres = ",  formatC(cap[1], width = tw, flag = " "),
              ", iter = ", formatC(i, width = tw, flag = " "),
              ", AIC = ",  formatC(round(ic.old, 4L), width = tw, flag = " "),
              ", df = ",   formatC(df, width = tw, flag = " "),
              ", ",        sprintf(pset_fmt, paste(ps.final, collapse = ", ")),
              if (!ia) "\n" else NULL,
              sep = "")
        }
        ## plot AIC and coef-paths
        if (plot & (i %% 10 == 0)) {
          
          par(mfrow = n2mfrow(length(nx) + 1))
          plot(y = lloos.list, x = 1:i, xlab = "Iteration", ylab = "AIC")
          
          if (i > 5) {
            fit2 <- lowess(y = lloos.list, x = 1:i)
            lines(fit2)
          }
          for (j in nx) {
            matplot(beta[[j]][1:i, ], type = "l", lty = 1, main = j,
                    xlab = "Iteration", ylab = "Coefficients")
          }

        }
        
      }
      
      
      ll.list[[c]] <- ll0.list
      
      if (verbose) {
        cs <- mean(lloos.list[(maxit1-30):maxit1])
        
        if (ia) cat("\r")
        cat(paste0("Thres = ",round(caps[c],4), ", df = ",df, ", AIC = ", round(cs,3)), "\n", sep = "")
      }
      
      # save betas
      beta.list[[paste(cap[1])]] <- beta
      
      # Continue the selection steps with smaller thresholds 
      if (cap[1] != caps[length(caps)]) {
        for (j in nx) {
          beta.list[[paste(caps[c + 1])]][[j]][1,] <- beta[[j]][maxit, ]
        }
      }
      
      df.list <- c(df.list, df)
      cs.list <- c(cs.list, cs)
      #if(cmax[c] < 0) break
    }
    
    
    final <- which.min(cs.list)
    print(paste("selected thres = ", round(caps[final], 3)))
    for (j in nx) beta[[j]] <- beta.list[[final]][[j]][(maxit + 1):maxit1, ]
    wmax <- maxit1

    return(list(coefficients = beta,
                logLik = ll0.list[[final]],
                maxit = list(var_selection = maxit, refitting = maxit_refit),
                iter = maxit_refit, cap = caps[final]))
} # end of function: sdr.thresdesc

# Helper Function to Draw Batch IDs (row ids)
#
# Intended for package internal use only; thus no sanity checks.
# x is the user object (NULL, single numeric, or a list
# of prepared indices for the iterations). "i" is the
# batch we would like to draw, N the total sample size.
get_batch_ids <- function(x, i, N) {
    if (is.null(x))         seq_len(N) # simply 1:N
    else if (is.list(x))    x[[i]]     # i'th user batch
    else if (is.numeric(x)) sample(seq_len(N), x, replace = FALSE) # random
    else stop("problems drawing batch IDs")
}


# best subset gradboosting with correlation filtering
# noncyclic is a special case with length_ps = 1
sdr.gradboostfit <- function(data,
                             family,
                             maxit = 1000,
                             refitting = TRUE,
                             maxit_refit = 0,
                             vardist = NULL,
                             eps = 0.01,
                             nu = 0.1,
                             eps_int = exp(seq(log(0.1), log(0.001), length = ifelse(refitting, maxit+maxit_refit, maxit))),
                             nu_int = 0.05,
                             oos_batch = "next",
                             aic = FALSE,
                             K = 0,
                             batch_ids = NULL,
                             verbose = TRUE,
                             initialize = TRUE,
                             plot = FALSE,
                             coef_start = NULL,
                             length_ps = NULL,
                             cap = 0.1,
                             x_mu = NULL,
                             x_sd = NULL,
                             nvars = NULL,
                             scalex = TRUE,
                             y = NULL,
                             vars = NULL,
                             ...) {
  
  ia <- interactive()
  if(!refitting) maxit_refit <- 0
  
  # scaling is done on each respective batch separately
  myscale <- function(x, x_mu, x_sd) {
    for (n in colnames(x)) x[, n] <- (x[, n] - x_mu[[n]]) / x_sd[[n]]
    return(x)
  }
  
  maxit1 <- maxit + maxit_refit
  cap <- c(rep(cap, length.out = maxit), rep(0,length.out = maxit_refit))
  eps <- rep(eps, length.out = maxit1)
  eps_int <- rep(eps_int, length.out = maxit1)
  
  N  <- nrow(data)
  nx <- names(vardist)
  
  # If the user provided a list of ids for the batches: check
  # if that matches the number of requested iterations + refitting.
  # Throw warning if not.
  
  if (is.list(batch_ids) && length(batch_ids) != maxit1) {
    warning("Length of batch_ids != maxit + maxit_refit; using length of batch_ids provided for setting maxit + maxit_refit!")
    maxit1 <- length(batch_ids)
  }
  
  # Draw initial set of batch IDs; `bids` (batch ids) will
  # be extended during iteration to keep track of ids.
  bids <- list(initial = get_batch_ids(batch_ids, i = 1, N = N))
  
  # TODO(R): Why "plots" and what do we actually do here?
  # number of plots
  ##ORIG#   tw <- length(strsplit(as.character(maxit1), "")[[1]]) + 1L
  tw <- nchar(as.character(maxit1)) + 1
  if (!exists("b.size")) b.size <- length(bids$initial)
  
  # K is penalty for IC
  if (is.null(K)) K <- log(b.size)
  
  ic.oos.list <- ic0.list <- ll.list <- ll0.list <- NULL
  ll.oos.afterint.list <- 0
  beta <- list()
  intvar <- list()
  for (j in nx) {
    if(length(vardist[[j]]) == 0){
      intvar[[j]] <- "(Intercept)"
      ncolu <- 1
    } else {
      intvar[[j]] <- c("(Intercept)",vardist[[j]])
      ncolu <- length(intvar[[j]])
    }
    beta[[j]] <- matrix(0, nrow = maxit1, ncol = ncolu)
    colnames(beta[[j]]) <- intvar[[j]]
  }
  if (!is.null(coef_start)) {
    initialize <- FALSE
    for (j in nx) {
      beta[[j]][1, ] <- coef_start[[j]]
    }
  }
  
  #beta.grad <- beta
  powerset.list <- powerset(nx)[-1]
  if (!is.null(length_ps))
    powerset.list <-
    powerset.list[sapply(powerset.list, length) == length_ps]
  
  if (initialize) {
    if (!is.null(family$initialize)) {
      betai <- list()
      for (j in nx) {
        if (!is.null(family$initialize[[j]])) {
          linkfun <- make.link2(family$links[j])$linkfun
          beta[[j]][1L, "(Intercept)"] <-
            mean(linkfun(family$initialize[[j]](data[bids$initial, y])), na.rm = TRUE)
        }
      }
    }
  }
  
  # Calculating degrees of freedom ...
  df <- sum(sapply(beta, function(b) { sum(b[1, ] != 0) }))
  
  err01 <- .Machine$double.eps^(1 / 2)
  err02 <- err01 * 2 * length(bids$initial)
  # err02 is the denominator for central numeric differentiation.  b.size makes
  # gradient magnitute for different sample sizes comparable this translates
  # maximum log likelihood to maximum average loglikelihood
  # https://stats.stackexchange.com/questions/267847/motivation-for-average-log-likelihood
  
  ma <- function(x, order = 20) {
    ma1 <- filter(x, rep(1 / order, order), sides = 1)
    ma2 <- rev(filter(rev(x), rep(1 / order, order), sides = 1))
    return(ifelse(is.na(ma1), ma2, ma1))
  }
  
  # If oos_batch is 'next' we would like to use the next batch for
  # validation. Thus we need to keep two sets of batch IDs and store
  # the batch ids of the current ($current) as well as the next ($next)
  # iteration into our `bids` list.
  #
  # After each iteration $current is overwritten with $next and a new
  # the new batch IDs for the $next iteration are added.
  # When reaching maxit1, $next will be set to $initial (first batch).
  bids$current = bids$initial # Set current to batch 1,
  bids$`next`  = get_batch_ids(batch_ids, i = 2, N = N) # next to batch 2
  
  # Iterating; looping 2:maxit1
  for (i in 2:maxit1) {
    
    # Update batch ids; write next -> current and draw
    # a new next (initial if i == maxit, else i + 1).
    bids$current <- bids$`next`
    bids$`next`  <- if (i == maxit1) bids$initial else get_batch_ids(batch_ids, i = i + 1, N = N)
    
    # deselect non selected variables for updating if maxit_refit > 0
    if (i == maxit +1) {
      for(j in nx){
        cond <- beta[[j]][maxit,] != 0
        cond["(Intercept)"] <- TRUE
        intvar[[j]]  <- intvar[[j]][cond]
        vardist[[j]] <- vardist[[j]][cond[-1]]
      }
    }
    
    eta.oos <- eta <- val <- absval <- sgn <- list()
    
    ## Extract response.
    yi <- as.matrix(data[bids$current, y])
    
    # Draw batch IDs for out of sample batch
    # Overwrite $oos in the bids list.
    bids$oos <- if (oos_batch == "next") {
      bids$`next`
    } else if (oos_batch == "random") {
      get_batch_ids(batch_ids, i = sample(maxit1, 1), N = N)
    } else if (oos_batch == "same") {
      bids$current
    } else stop("Whoops, unknown oos_batch case; we should never end up here")
    
    y.oos <- as.matrix(data[bids$oos, y])
    # draw batchwise from big.matrix and scale
    XX    <- as.matrix(data[bids$current, vars])
    XXoos <- as.matrix(data[bids$oos, vars])
    if (scalex){
      XX <- myscale(XX, x_mu = x_mu, x_sd = x_sd)
      XXoos <- myscale(XXoos, x_mu = x_mu, x_sd = x_sd)
    }
    X <- Xoos <- list()
    
    for (j in nx) {
      ## Save last iteration.
      beta[[j]][i, ] <- beta[[j]][i - 1L, ]
      if(ncol(XX) == 1){
        if(ncol(beta[[j]]) == 1 ){
          X[[j]] <- cbind( "(Intercept)" = rep(1,b.size))
          Xoos[[j]] <- cbind( "(Intercept)" = rep(1,b.size))
        } else {
          X[[j]] <- cbind( "(Intercept)" = rep(1,b.size), XX)
          Xoos[[j]] <- cbind( "(Intercept)" = rep(1,b.size),XXoos)
        }
        colnames(Xoos[[j]]) <- colnames(X[[j]]) <- colnames(beta[[j]]) 
      } else {
        X[[j]] <- cbind( "(Intercept)" = rep(1,b.size), XX[,setdiff(colnames(beta[[j]]), "(Intercept)" )])
        Xoos[[j]] <- cbind( "(Intercept)" = rep(1,b.size),XXoos[,setdiff(colnames(beta[[j]]), "(Intercept)" )])
        colnames(Xoos[[j]]) <- colnames(X[[j]]) <- c("(Intercept)", setdiff(colnames(beta[[j]]), "(Intercept)" ))
      }
      ## Setup linear predictor.
      eta[[j]] <-
        drop(X[[j]] %*% beta[[j]][i, ])
      
      ## out of sample
      eta.oos[[j]] <-
        drop(Xoos[[j]] %*% beta[[j]][i, ])
    }
    
    
    # out of sample
    ll.oos <- family$loglik(y.oos, family$map2par(eta.oos))
    ll.oos.afterint <- ll.oos
    
    ## Compute log-likelihood.
    df <- sum(sapply(beta, function(b) {
      sum(b[i, ] != 0)
    }))
    ll0 <- family$loglik(yi, family$map2par(eta))
    ic0 <- -2 * ll0 + K * df
    ic.oos <- -2 * ll.oos + K * df
    
    ll0.list <- c(ll0.list, ll0)
    ic0.list <- c(ic0.list, ic0)
    ic.oos.list <- c(ic.oos.list, ic.oos)
    
    
    bb = 20
    if (plot & (i %% 10 == 0)) {
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
      plot( y = ll0.list, x = 2:i, xlab = "Iteration", ylab = "logLik")
      
      if (i > bb + 4) {
        abline(v = bic.min)
        abline(v = bic.min - bb)
      }
      if (i > 5) {
        fit2 <- lowess(y = ll0.list, x = 2:i)
        lines(fit2)
      }
      for (j in nx) {
        matplot(beta[[j]][1:i, ], type = "l", lty = 1, main = j, xlab = "Iteration", ylab = "Coefficients")
        
        if (i > bb + 4) {
          abline(v = bic.min)
          abline(v = bic.min - bb)
        }
      }
    }
    
    
    
    sign.list <- rep(-Inf, length(nx))
    names(sign.list) <- nx
    
    for (j in nx) {
      eta0 <- eta
      ## Get coefficients and setup.
      tbeta <- beta[[j]][i, ]
      ## Positive.
      tbeta[1] <- tbeta[1] + err01
      eta[[j]] <- drop(X[[j]] %*% tbeta)
      ll1      <- family$loglik(yi, family$map2par(eta))
      
      ## Negative
      tbeta[1] <- tbeta[1] - 2 * err01
      eta[[j]] <- drop(X[[j]] %*% tbeta)
      ll2      <- family$loglik(yi, family$map2par(eta))
      
      grad <- (ll1 - ll2) / err02
      sign.list[[j]] <- grad
    }
    
    pset <- powerset.list
    for (j in nx) {
      if (sign.list[[j]] == -Inf) {
        ok <- !grepl(j, pset)
        pset <- pset[ok]
      }
    }
    
    if (length(pset) != 0) {
      beta.final <- list()
      for (j in nx)
        beta.final[[j]] <- beta[[j]][i, ]
      ic.oos.old <- ic.oos
      
      beta.wip <- list()
      for (j in nx)
        beta.wip[[j]] <- beta[[j]][i, ]
      
      for (l in 1:length(pset)) {
        # ps <- !is.na(pset[l,])
        sign.list2 <- sign.list[pset[[l]]]
        gnorm <- sqrt(sum(sign.list2 ^ 2))
        
        if (gnorm > eps_int[i]) {
          sign.list2 <- eps_int[i] * sign.list2 / gnorm
        }
        for (ij in pset[[l]]) {
          sign.list2[ij] <- ifelse(
            i < 0.8 * maxit1 & abs(sign.list2[ij]) <
              nu_int * eps_int[i],
            sign(sign.list2[ij]) * nu_int * eps_int[i],
            sign.list2[ij]
          )
        }
        
        for (j in pset[[l]]) {
          grad <- sign.list2[[j]]
          beta[[j]][i, 1] <- beta[[j]][i, 1] + grad
          
          eta[[j]] <- drop(X[[j]] %*% beta[[j]][i, ])
        }
        
        df <- sum(sapply(beta, function(b) {
          sum(b[i, ] != 0)
        }))
        ##df <- sum(data.frame(beta)[i, ] != 0)
        
        ## keep update only if oos information crit improves
        for (j in nx)
          eta.oos[[j]] <- drop(Xoos[[j]] %*% beta[[j]][i, ])
        
        ll.oos <- family$loglik(y.oos, family$map2par(eta.oos))
        ic.oos.new <- -2 * ll.oos + K * df
        if (ic.oos.new < ic.oos.old) {
          # keep current try
          for (j in nx)
            beta.final[[j]] <- beta[[j]][i, ]
          ps.final <- pset[[l]]
          ic.oos.old <- ic.oos.new
          ll.oos.afterint <- ll.oos
        }
        for (j in pset[[l]])
          beta[[j]][i, ] <- beta.wip[[j]]
        
      }  # this bracket is for pset
      
      # save best beta
      for (j in nx)
        beta[[j]][i, ] <- beta.final[[j]]
      
    }
    
    eta0 <- eta
    sign.list <- pos.list <- rep(-Inf, length(nx))
    names(sign.list) <- names(pos.list) <- nx
    
    for (j in nx) {
      ## Get coefficients and setup.
      tbeta <- beta[[j]][i, ]
      nc <- length(intvar[[j]]) - 1
      
      if (nc > 0) {
        eta[[j]] <- eta[[j]] + err01
        ll1 <- family$d(yi, family$map2par(eta), log = TRUE)
        
        ## Negative
        eta[[j]] <- eta[[j]] - 2 * err01
        ll2 <- family$d(yi, family$map2par(eta), log = TRUE)
        # if(T ) print(ll2)
        
        grad <- (ll1 - ll2) / (2 * err01)
        # print(mean(grad))
        eta[[j]] <- eta[[j]] + err01
        
        cc <- try(cor(grad, X[[j]][, vardist[[j]], drop = FALSE]), FALSE)
        
        #[,-1]))
        cc[is.na(cc)] <- 0
        if (is.numeric(cc)) {
          ## Select update
          jj <- which.max(abs(cc)) + 1
        } else {
          jj <- 1
        }  # interecept if cor gives error
        
        eta0  <- eta
        ## Get coefficients and setup.
        tbeta <- beta[[j]][i, ]
        
        if (max(abs(cc)) < cap[i]) {
          grad <- 0
        } else {
          # average partial derivative with respect to best variable
          grad <- t(grad) %*% X[[j]][,intvar[[j]][jj] , drop = FALSE]/ b.size 
          # ## Positive.
          # tbeta[jj] <- tbeta[jj] + err01
          # eta[[j]] <-
          #   drop(X[batch_ids[[i]],intvar[[j]] , drop = FALSE] %*% tbeta)
          # ll1 <- family$loglik(yi, family$map2par(eta))
          # 
          # ## Negative
          # tbeta[jj] <- tbeta[jj] - 2 * err01
          # eta[[j]] <-
          #   drop(X[batch_ids[[i]],intvar[[j]] , drop = FALSE] %*% tbeta)
          # ll2 <- family$loglik(yi, family$map2par(eta))
          # 
          # grad <- (ll1 - ll2) / err02
        }
        
        sign.list[[j]] <- grad
        pos.list[[j]]  <- jj # position in intvar
        
      }
      
    }
    
    # Max number of characters in case all parameters are set;
    # used for formatting cat() output if verbose = TRUE
    pset_fmt <- paste0("%-", max(sapply(pset, function(x) nchar(paste(x, collapse = ", ")))), "s")
    
    ps.final <- "no par"
    for (j in nx) {
      if (sign.list[[j]] == 0 | sign.list[[j]] == -Inf) {
        ok   <- !grepl(j, pset)
        pset <- pset[ok]
      }
    }
    
    if (length(pset) == 0) {
      # for (j in nx) beta[[j]][i, ] <- beta[[j]][i - 1, ]
    } else {
      ## beta.final <- list()
      for (j in nx)
        beta.final[[j]] <- beta[[j]][i, ]
      #ic.oos.old <- ic.oos
      
      beta.wip <- list()
      for (j in nx)
        beta.wip[[j]] <- beta[[j]][i, ]
      for (l in 1:length(pset)) {
        # ps <- !is.na(pset[l,])
        sign.list2 <- sign.list[pset[[l]]]
        gnorm      <- sqrt(sum(sign.list2 ^ 2))
        
        if (gnorm > eps[i]) {
          sign.list2 <- eps[i] * sign.list2 / gnorm
        }
        for (ij in pset[[l]]) {
          sign.list2[ij] <- ifelse(i < 0.8 * maxit1 & abs(sign.list2[ij]) < nu * eps[i],
                                   sign(sign.list2[ij]) * nu * eps[i],
                                   sign.list2[ij])
        }
        
        for (j in pset[[l]]) {
          jj   <- pos.list[[j]]
          grad <- sign.list2[[j]]
          
          beta[[j]][i, intvar[[j]][jj]] <- beta[[j]][i, intvar[[j]][jj]] + grad
          
          eta[[j]] <- drop(X[[j]] %*% beta[[j]][i, ])
        }
        
        ## df <- sum(data.frame(beta)[i, ] != 0)
        df <- sum(sapply(beta, function(b) {
          sum(b[i, ] != 0)
        }))
        
        ## keep update only if oos information crit improves
        for (j in nx)
          eta.oos[[j]] <- drop(Xoos[[j]] %*% beta[[j]][i, ])
        
        ll.oos <- family$loglik(y.oos, family$map2par(eta.oos))
        ic.oos.new <- -2 * ll.oos + K * df
        if (ic.oos.new < ic.oos.old) {
          # keep current try
          for (j in nx)
            beta.final[[j]] <- beta[[j]][i, ]
          
          ps.final <- pset[[l]]
          ic.oos.old <- ic.oos.new
          ll.oos.aftercoef <- ll.oos
        }
        for (j in pset[[l]])
          beta[[j]][i, ] <- beta.wip[[j]]
        
      }  # this bracket is for pset
      
      # save best beta
      for (j in nx)
        beta[[j]][i, ] <- beta.final[[j]]
      
    }
    ll.oos.afterint.list <- c(ll.oos.afterint.list, ll.oos.afterint)
    if (verbose) {
      if (ia) cat("\r")
      cat(
        "iter = ",     formatC(i, width = tw, flag = " "),
        ", logLik = ", formatC(round(ll0, 4L), width = tw, flag = " "),
        ", df = ",     formatC(df, width = tw, flag = " "),
        ", ",          sprintf(pset_fmt, paste(ps.final, collapse = ", ")),
        if (!ia) "\n" else NULL,
        sep = ""
      )
    }
    
  } # End of loop over 2:maxit1
  
  
  ## Compute log-likelihood.
  ll <- family$loglik(yi, family$map2par(eta))
  
  ## Extract 'out of sample' response.
  ## bids$oos is the last out-of-sample batch index vector.
  ## If oos_batch = 'next' is used this one should be the
  ## first one again (identical to bids$initial).
  
  y.oos <- as.matrix(data[bids$oos, y])
  XX    <- as.matrix(data[bids$oos, vars])
  if (scalex) XX <- myscale(XX, x_mu = x_mu, x_sd = x_sd)
  
  X <- list()
  for (j in nx) {
    if (ncol(XX) == 1){
      if (ncol(beta[[j]]) == 1 ) {
        Xoos[[j]] <- cbind("(Intercept)" = rep(1, b.size))
      } else {
        Xoos[[j]] <- cbind("(Intercept)" = rep(1, b.size), XXoos)
      }
      colnames(Xoos[[j]]) <- colnames(beta[[j]]) 
    } else {
      Xoos[[j]]   <- cbind("(Intercept)" = rep(1, b.size),
                           XXoos[, setdiff(colnames(beta[[j]]), "(Intercept)" )])
      colnames(Xoos[[j]]) <- c("(Intercept)", setdiff(colnames(beta[[j]]), "(Intercept)" ))
    }
  } 
  for (j in nx) {
    ## Setup linear predictor.
    eta[[j]] <- drop(Xoos[[j]] %*% beta[[j]][maxit1, ])
  }
  
  ## Compute 'out of sample' log-likelihood.
  ll0      <- family$loglik(y.oos, family$map2par(eta))
  ll0.list <- c(ll0.list, ll0)
  
  if (verbose) {
    if (ia) cat("\r")
    cat("iter = ", formatC(i, width = tw, flag = " "), ", ",
        "logLik = ", formatC(round(ll, 4L), width = tw, flag = " "), "\n", sep = "")
  }
  
  return(list(coefficients = beta,
              logLik       = ll0.list,
              logLik.afterint = ll.oos.afterint.list,
              maxit        = list(var_selection = maxit, refitting = maxit_refit)))
} # end of function: sdr.gradboostfit

# cyclic gradboosting with correlation filtering
sdr.gradboostfit2 <- function(data,
                              family,
                              maxit = 1000,
                              refitting = TRUE,
                              maxit_refit = 0,
                              vardist = NULL,
                              eps = 0.01,
                              nu = 0.1,
                              eps_int = exp(seq(log(0.1), log(0.001), length = ifelse(refitting, maxit + maxit_refit, maxit))),
                              nu_int = 0.05,
                              oos_batch = "next",
                              aic = FALSE,
                              K = 0,
                              batch_ids = NULL,
                              verbose = TRUE,
                              initialize = TRUE,
                              plot = FALSE,
                              coef_start = NULL,
                              cap = 0.1,
                              x_mu = NULL,
                              x_sd = NULL,
                              nvars = NULL,
                              scalex = TRUE,
                              y = NULL,
                              vars = NULL,
                              ...) {
    ia <- interactive()
    if(!refitting) maxit_refit <- 0
    
    # scaling is done on each respective batch separately
    myscale <- function(x, x_mu, x_sd) {
      for (n in colnames(x)) x[, n] <- (x[, n] - x_mu[[n]]) / x_sd[[n]]
      return(x)
    }
    
    maxit1  <- maxit + maxit_refit
    cap     <- c(rep(cap, length.out = maxit), rep(0,length.out = maxit_refit))
    eps     <- rep(eps, length.out = maxit1)
    eps_int <- rep(eps_int, length.out = maxit1)
    
    N <- nrow(data)
    nx <- names(vardist)
    if (is.null(batch_ids)) {
      ind <- 1:N
      b.size <- N
      batch_ids <- lapply(1:maxit1, function(...) ind)
    } else {
      if (is.numeric(batch_ids)) {
        ind <- 1:N
        b.size <- batch_ids
        batch_ids <- lapply(1:maxit1, function(...) sample(ind, size = batch_ids, replace = FALSE))
      }
    }

    if (!is.list(batch_ids))
      stop("Argument batch_ids must be a list of indices!")

    if (length(batch_ids) != maxit1)
      warning("Length of batch_ids != maxit+maxit_refit, using batch_ids for setting maxit+maxit_refit!")

    maxit1 <- length(batch_ids)
    # number of plots
    tw <- length(strsplit(as.character(maxit1), "")[[1]]) + 1L
    if(!exists("b.size")) b.size <- length(batch_ids[[1]])

    # K is penalty for IC
    if (is.null(K)) K <- log(b.size)
    
    ic.oos.list <- ic0.list <- ll.list <- ll0.list <- NULL
    beta   <- list()
    intvar <- list()
    for (j in nx) {
      if(length(vardist[[j]]) == 0){
        intvar[[j]] <- "(Intercept)"
        ncolu <- 1
      } else {
        intvar[[j]] <- c("(Intercept)",vardist[[j]])
        ncolu <- length(intvar[[j]])
      }
      beta[[j]] <- matrix(0, nrow = maxit1, ncol = ncolu)
      colnames(beta[[j]]) <- intvar[[j]]
    }
    if (!is.null(coef_start)) {
      initialize <- FALSE
      for (j in nx) {
        beta[[j]][1, ] <- coef_start[[j]]
      }
    }
    
    if (initialize) {
      if (!is.null(family$initialize)) {
        betai <- list()
        for (j in nx) {
          if (!is.null(family$initialize[[j]])) {
            linkfun <- make.link2(family$links[j])$linkfun
            beta[[j]][1L, "(Intercept)"] <-
              mean(linkfun(family$initialize[[j]](data[batch_ids[[1]],y],...
              )), na.rm = TRUE)
            
          }
        }
      }
    }
    df <- sum(sapply(beta, function(b) {
      sum(b[1, ] != 0)
    }))
    err01 <- .Machine$double.eps^(1 / 2)
    err02 <- err01 * 2 * nrow(data.frame(batch_ids[1]))
    # err02 is the denominator for central numeric differentiation.  b.size makes
    # gradient magnitute for different sample sizes comparable this translates
    # maximum log likelihood to maximum average loglikelihood
    # https://stats.stackexchange.com/questions/267847/motivation-for-average-log-likelihood
    
    # TODO(R): Die funktion gibt es igendwie 5x an verschiedensten stellen!
    ma <- function(x, order = 20) {
      ma1 <- filter(x, rep(1 / order, order), sides = 1)
      ma2 <- rev(filter(rev(x), rep(1 / order, order), sides = 1))
      return(ifelse(is.na(ma1), ma2, ma1))
    }
    
    for (i in 2:maxit1) {
      # deselect non selected variables for updating
      if(i == maxit +1){
        for(j in nx){
          cond <- beta[[j]][maxit,] != 0
          cond["(Intercept)"] <- TRUE
          intvar[[j]]  <- intvar[[j]][cond]
          vardist[[j]] <- vardist[[j]][cond[-1]]
        }
      }
      
      eta.oos <- eta <- val <- absval <- sgn <- list()
      
      ## Extract response.
      yi <- as.matrix(data[batch_ids[[i]],y])
      
      ## out of sample
      if (oos_batch == "next") {
        if (i != maxit1) {
          i.oos <- i + 1
        } else {
          i.oos <- 1
        }
      }
      if (oos_batch == "same") {
        i.oos <- i
      }
      if (oos_batch == "random") {
        i.oos <- sample(maxit1, 1)
      }
      
      batch.oos <- batch_ids[[i.oos]]
      y.oos     <- as.matrix(data[batch.oos,y])
      # draw batchwise from big.matrix and scale
      XX        <- as.matrix(data[batch_ids[[i]], vars])
      XXoos     <- as.matrix(data[batch.oos,vars])
      if(scalex){
        XX <- myscale(XX, x_mu = x_mu, x_sd = x_sd)
        XXoos <- myscale(XXoos, x_mu = x_mu, x_sd = x_sd)
      }

      X <- Xoos <- list()
      for (j in nx) {
        ## Save last iteration.
        beta[[j]][i, ] <- beta[[j]][i - 1L, ]
        X[[j]]         <- cbind( "(Intercept)" = rep(1,b.size), XX[,setdiff(colnames(beta[[j]]), "(Intercept)" )])
        Xoos[[j]]      <- cbind( "(Intercept)" = rep(1,b.size),XXoos[,setdiff(colnames(beta[[j]]), "(Intercept)" )])
        colnames(Xoos[[j]]) <- colnames(X[[j]]) <- c("(Intercept)", setdiff(colnames(beta[[j]]), "(Intercept)" ))

        ## Setup linear predictor.
        eta[[j]]     <- drop(X[[j]] %*% beta[[j]][i, ])
        ## out of sample
        eta.oos[[j]] <-  drop(Xoos[[j]] %*% beta[[j]][i, ])
      }
      
      # out of sample
      ll.oos <- family$loglik(y.oos, family$map2par(eta.oos))
      
      ## Compute log-likelihood.
      df <- sum(sapply(beta, function(b) {
        sum(b[i, ] != 0)
      }))
      ll0    <- family$loglik(yi, family$map2par(eta))
      ic0    <- -2 * ll0 + K * df
      ic.oos <- -2 * ll.oos + K * df
      
      ll0.list    <- c(ll0.list, ll0)
      ic0.list    <- c(ic0.list, ic0)
      ic.oos.list <- c(ic.oos.list, ic.oos)
      
      
      bb = 20
      if (plot & (i %% 10 == 0)) {
        if (i > bb + 4) {
          bic.min <- which.min(ma(ic0.list, order = bb))
          if (bic.min < bb + 1) bic.min <- bic.min + bb
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
          matplot(beta[[j]][1:i, ], type = "l", lty = 1, main = j, xlab = "Iteration", ylab = "Coefficients")

          if (i > bb + 4) {
            abline(v = bic.min)
            abline(v = bic.min - bb)
          }
        }
      }
      
      
      
      # Cyclic updating of intercepts
      for (j in nx) {
        eta0 <- eta
        ## Get coefficients and setup.
        tbeta <- beta[[j]][i, ]
        ## Positive.
        tbeta[1] <- tbeta[1] + err01
        eta[[j]] <- drop(X[[j]] %*% tbeta)
        ll1 <- family$loglik(yi, family$map2par(eta))
        
        ## Negative
        tbeta[1] <- tbeta[1] - 2 * err01
        eta[[j]] <- drop(X[[j]] %*% tbeta)
        ll2 <- family$loglik(yi, family$map2par(eta))
        
        grad <- (ll1 - ll2) / err02
        
        if(abs(grad) > eps_int[i]){
          grad <- sign(grad) * eps_int[i]
        } else if (i < 0.8 * maxit & abs(grad) <= nu_int * eps_int[i]) {
          grad <- nu_int * sign(grad) * eps_int[i]
        }
        # eta[[j]] <- eta0[[j]]
        
        # intercept
        beta[[j]][i, 1] <- beta[[j]][i, 1] + grad
        tbeta <- beta[[j]][i, ]
        ## out of sample
        eta.oos[[j]] <- drop(Xoos[[j]] %*% tbeta)
        
        ll2 <- family$loglik(y.oos, family$map2par(eta.oos))
        if (ll2 > ll0) {
          ll0 <- ll2
        } else {
          beta[[j]][i, 1] <- beta[[j]][i, 1] - grad
        }
        eta[[j]] <- drop(X[[j]] %*% beta[[j]][i, ])
      }
      
      
      names(sign.list) <- names(pos.list) <- nx
      
      for (j in nx) {
        ## Get coefficients and setup.
        tbeta <- beta[[j]][i, ]
        nc    <- length(intvar[[j]]) - 1
        
        if (nc > 0) {
          eta[[j]] <- eta[[j]] + err01
          ll1      <- family$d(yi, family$map2par(eta), log = TRUE)
          
          ## Negative
          eta[[j]] <- eta[[j]] - 2 * err01
          ll2      <- family$d(yi, family$map2par(eta), log = TRUE)
          # if(T ) print(ll2)
          
          grad     <- (ll1 - ll2) / (2 * err01)
          # print(mean(grad))
          eta[[j]] <- eta[[j]] + err01
          
          cc <- try(cor(grad, X[[j]][, vardist[[j]], drop = FALSE]), FALSE)
          #[,-1]))
          cc[is.na(cc)] <- 0
          if (is.numeric(cc)) {
            ## Select update
            jj <- which.max(abs(cc)) + 1
          } else {
            jj <- 1
          }  # interecept if cor gives error
          
          eta0  <- eta
          ## Get coefficients and setup.
          tbeta <- beta[[j]][i, ]
          
          if (max(abs(cc)) < cap[i]) {
            grad <- 0
          } else {
            # average partial derivative with respect to best variable
            grad <- t(grad) %*% X[[j]][,intvar[[j]][jj] , drop = FALSE]/ b.size 
            
            if(abs(grad) > eps[i]){
              grad <- sign(grad) * eps[i]
            } else if (i < 0.8 * maxit & abs(grad) <= nu * eps[i]) {
              grad <- nu * sign(grad) * eps[i]
            }
            # var update
            beta[[j]][i, intvar[[j]][jj]] <- beta[[j]][i, intvar[[j]][jj]] + grad
            tbeta <- beta[[j]][i, ]
            ## out of sample
            eta.oos[[j]] <- drop(Xoos[[j]] %*% tbeta)
            
            ll2 <- family$loglik(y.oos, family$map2par(eta.oos))
            if (ll2 > ll0) {
              ll0 <- ll2
            } else {
              beta[[j]][i, intvar[[j]][jj]] <- beta[[j]][i, intvar[[j]][jj]] - grad
            }
            eta[[j]] <- drop(X[[j]] %*% beta[[j]][i, ])
          }
          
        }
        
      }
      
      if (verbose) {
        if (ia) cat("\r")
        cat("iter = ", formatC(i, width = tw, flag = " "),
            ", logLik = ", formatC(round(ll0, 4L), width = tw, flag = " "),
            ", df = ", formatC(df, width = tw, flag = " "),
            ", ", sprintf(pset_fmt, paste(ps.final, collapse = ", ")),
        if (!ia) "\n" else NULL, sep = "")
      }
    }
    
    
    ## Compute log-likelihood.
    ll <- family$loglik(yi, family$map2par(eta))
    
    ## Extract 'out of sample' response.
    yoos <- as.matrix(data[batch_ids[[1]],y])
    
    XX   <- as.matrix(data[batch_ids[[1]],vars])
    XX   <- myscale(XX, x_mu = x_mu, x_sd = x_sd)
    
    X    <- list()

    for (j in nx) {
      Xoos[[j]] <- cbind( "(Intercept)" = rep(1,b.size), XX[,setdiff(colnames(beta[[j]]), "(Intercept)" )])
    } 

    for (j in nx) {
      ## Setup linear predictor.
      eta[[j]] <-
        drop(Xoos[[j]] %*% beta[[j]][maxit1, ])
    }

    ## Compute 'out of sample' log-likelihood.
    ll0 <- family$loglik(yoos, family$map2par(eta))
    ll0.list <- c(ll0.list, ll0)
    
    if (verbose) {
      if (ia) cat("\r")
      cat("iter = ", formatC(i, width = tw, flag = " "),
          ", logLik = ", formatC(round(ll, 4L), width = tw, flag = " "), "\n", sep = "")
    }
    
    rval <- list(coefficients = beta,
                 logLik = ll0.list,
                 maxit = list(var_selection = maxit, refitting = maxit_refit))
    return(rval)
} # end of function: sdr.gradboostfit2


ma <- function(x, order = 20) {
  ma1 <- filter(x, rep(1 / order, order), sides = 1)
  ma2 <- rev(filter(rev(x), rep(1 / order, order), sides = 1))
  return(ifelse(is.na(ma1), ma2, ma1))
}

# cyclical gradboosting with correlation filtering
sdr.gradboostfit2_old <- function(X,
                                  y,
                                  family,
                                  maxit = 1000,
                                  start = NULL,
                                  eps = 0.01,
                                  nu = 0.1,
                                  aic = FALSE,
                                  K = 0,
                                  full_grad = FALSE,
                                  cores = 1,
                                  batch_ids = NULL,
                                  verbose = TRUE,
                                  initialize = TRUE,
                                  qstab = NULL,
                                  int_ind = TRUE,
                                  plot = FALSE,
                                  steps = c(-1, 1),
                                  ts = F,
                                  coef_start = NULL,
                                  nu_int = 0.05,
                                  eps_int = exp(seq(log(0.1), log(0.001), length = maxit)),
                                  replace = F,
                                  comp = F,
                                  cap = 0,
                                  line.search = exp(seq(
                                    from = log(0.1),
                                    to = log(10),
                                    length.out = 10
                                  )),
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
        batch_ids <- lapply(1:maxit, function(...) sample(ind, size = batch_ids, replace = replace))
      } 
    }

    if (!is.list(batch_ids))
      stop("Argument batch_ids must be a list of indices!")
    if (length(batch_ids) != maxit)
      warning("Length of batch_ids != maxit, using batch_ids for setting maxit!")

    maxit <- length(batch_ids)
    tw    <- length(strsplit(as.character(maxit), "")[[1]]) + 1L
    if (is.null(K)) K <- log(length(y))
    
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
        if (ts) nxx <- "sigma"

        for (j in nxx) {
          if (!is.null(family$initialize[[j]])) {
            linkfun <- make.link2(family$links[j])$linkfun
            beta[[j]][1L, "(Intercept)"] <- mean(linkfun(family$initialize[[j]](y, ...)), na.rm = TRUE)
          }
        }
      }
    }
    # print(beta[[1]][1L, '(Intercept)'])
    df <- sum(data.frame(beta)[1L, ] != 0)
    
    err01 <- .Machine$double.eps^(1 / 2)
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
        eta[[j]]     <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% beta[[j]][i, ])
        
        ## out of sample
        eta.oos[[j]] <- drop(X[[j]][batch.oos, , drop = FALSE] %*% beta[[j]][i, ])
      }
      
      
      # out of sample
      ll.oos <- family$loglik(y.oos, family$map2par(eta.oos))
      
      ## Compute log-likelihood.
      df     <- sum(data.frame(beta)[i, ] != 0)
      ll0    <- family$loglik(yi, family$map2par(eta))
      ic0    <- -2 * ll0 + K * df
      ic.oos <- -2 * ll.oos + K * df
      
      ll0.list    <- c(ll0.list, ll0)
      ic0.list    <- c(ic0.list, ic0)
      ic.oos.list <- c(ic.oos.list, ic.oos)
      
      if (verbose) {
        if (ia) cat("\r")
        cat("iter = ", formatC(i - 1L, width = tw, flag = " "),
            ", logLik = ", formatC(round(ll0, 4L), width = tw, flag = " "),
            ", df = ", formatC(df, width = tw, flag = " "),
        if (!ia) "\n" else NULL, sep = "")
      }
      bb = 20
      if (plot & (i %% 10 == 0)) {
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
          matplot(beta[[j]][1:i, ], type = "l", lty = 1, main = j, xlab = "Iteration", ylab = "Coefficients")

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
        ll1      <- family$loglik(yi, family$map2par(eta))
        
        ## Negative
        tbeta[1] <- tbeta[1] - 2 * err01
        eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
        ll2      <- family$loglik(yi, family$map2par(eta))
        
        grad <- (ll1 - ll2) / err02
        # grad <- ifelse(abs(grad) > eps_int[i], sign(grad)*eps_int[i],
        # ifelse(abs(grad) <= 0.1*eps_int[i], 0.1*sign(grad)*eps_int[i], grad))
        grad <- eps_int[i] * grad
        if (i < 0.8 * maxit & abs(grad) <= nu_int * eps_int[i])
          grad <- nu_int * sign(grad) * eps_int[i]
        # eta[[j]] <- eta0[[j]]
        
        # intercept
        beta[[j]][i, 1] <- beta[[j]][i, 1] + grad
        tbeta    <- beta[[j]][i, ]
        eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
        ll2 <- family$loglik(yi, family$map2par(eta))
        if (ll2 > ll0) {
          ll0 <- ll2
        } else {
          beta[[j]][i, 1] <- beta[[j]][i, 1] - grad
        }
      }
      
      eta0             <- eta
      sign.list        <- pos.list <- rep(0, length(nx))
      names(sign.list) <- names(pos.list) <- nx
      
      for (j in nx) {
        ## Get coefficients and setup.
        tbeta <- beta[[j]][i, ]
        nc <- length(tbeta) - 1
        
        if (nc > 0) {
          eta[[j]] <- eta[[j]] + err01
          ll1      <- family$d(yi, family$map2par(eta), log = TRUE)
          
          ## Negative
          eta[[j]] <- eta[[j]] - 2 * err01
          ll2      <- family$d(yi, family$map2par(eta), log = TRUE)
          # if(T ) print(ll2)
          
          grad     <- (ll1 - ll2) / (2 * err01)
          # print(mean(grad))
          eta[[j]] <- eta[[j]] + err01
          
          cc <- try(cor(grad, X[[j]][batch_ids[[i]], -1, drop = FALSE]), F)
          #[,-1]))
          cc[is.na(cc)] <- 0
          if (is.numeric(cc)) {
            ## Select update
            jj <- which.max(abs(cc)) + 1
          } else {
            jj <- 1
          } # interecept if cor gives error
          
          eta0  <- eta
          ## Get coefficients and setup.
          tbeta <- beta[[j]][i, ]
          
          if (max(abs(cc)) < cap[i]) {
            grad <- 0
          } else {
            ## Positive.
            tbeta[jj] <- tbeta[jj] + err01
            eta[[j]]  <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
            ll1       <- family$loglik(yi, family$map2par(eta))
            
            ## Negative
            tbeta[jj] <- tbeta[jj] - 2 * err01
            eta[[j]]  <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% tbeta)
            ll2       <- family$loglik(yi, family$map2par(eta))
            
            grad <- (ll1 - ll2) / err02
          }
          
          sign.list[[j]] <- grad
          pos.list[[j]]  <- jj
          
        }
        
      } # end for j in nx
      
      
      gnorm <- sqrt(sum(sign.list ^ 2))
      
      if (gnorm > eps[i]) {
        sign.list <- eps[i] * sign.list / gnorm
      }

      for (ij in nx) {
        sign.list[ij] <- ifelse(i < 0.8 * maxit & sign.list[ij] != 0 & abs(sign.list[ij]) < nu * eps[i],
                                sign(sign.list[ij]) * nu * eps[i],
                                sign.list[ij])
      }
      for (j in nx) {
        if (sign.list[[j]] != 0) {
          jj               <- pos.list[[j]]
          grad             <- sign.list[[j]]
          beta[[j]][i, jj] <- beta[[j]][i, jj] + grad
          
          eta[[j]] <- drop(X[[j]][batch_ids[[i]], , drop = FALSE] %*% beta[[j]][i, ])
        }
      }
      
      
      df <- sum(data.frame(beta)[i, ] != 0)
      
      ## keep update only if oos information crit improves
      for (j in nx)
        eta.oos[[j]] <- drop(X[[j]][batch.oos, , drop = FALSE] %*% beta[[j]][i, ])

      ll.oos <- family$loglik(y.oos, family$map2par(eta.oos))
      if (K >= 0)
        ic.oos.new <- -2 * ll.oos + K * df
      else
        ic.oos.new <- -Inf

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
      eta[[j]] <- drop(X[[j]][batch_ids[[1]], , drop = FALSE] %*% beta[[j]][maxit, ])
    }
    ## Compute 'out of sample' log-likelihood.
    ll0      <- family$loglik(yi, family$map2par(eta))
    ll0.list <- c(ll0.list, ll0)
    
    if (verbose) {
      if (ia) cat("\r")
      cat("iter = ", formatC(i, width = tw, flag = " "),
          ", logLik = ", formatC(round(ll, 4L), width = tw, flag = " "), "\n", sep = "")
    }
    
    return(list(coefficients = beta,
                logLik = ll0.list,
                maxit = maxit))
} # end of function: sdr.gradboostfit2_old




## Model frame.
model.frame.sdr <- function(formula, data = NULL, ...) {
  dots  <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
  if (length(nargs) || is.null(formula$model) || !is.null(data)) {
    fcall <- formula$call
    m     <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(fcall), 0L)
    fcall <- fcall[c(1L, m)]
    fcall$drop.unused.levels <- TRUE
    if (!is.null(data))
      fcall$data <- as.name(quote("data"))
    fcall[[1L]] <- quote(model.frame)

    mf <- list()
    k  <- 1
    for (i in names(formula$formula)) {
      fi <- Formula::as.Formula(formula$formula[[i]])
      fj <- formula(fi, lhs = 0, rhs = 1, drop = TRUE)
      vars <- all.vars(fj)
      if (length(vars) < 1) {
        ff <- ~ 1
      } else {
        ff <- as.formula(paste("~", paste(all.vars(fj), collapse = "+")))
      }

      fcall[[2]] <- as.call(ff)
      mf[[k]] <- eval(fcall, if (is.null(data)) parent.frame() else NULL)
      k <- k + 1
    }
    mf <- do.call("cbind", mf)
    mf <- mf[, unique(names(mf)), drop = FALSE]
    return(mf)
  } else {
    return(formula$model)
  }
} # end of function: model.frame.sdr


## Model matrix.
model.matrix.sdr <- function(object, ...) {
  mf <- model.frame(object, ...)
  
  X <- list()
  for (i in names(object$formula)) {
    tll <- attr(object$x[[i]]$terms, "term.labels")
    tls <- unlist(sapply(attr(object$x[[i]]$terms, "specials"), identity))
    if (!is.null(tls)) {
      tll <- tll[-tls]
      tls <- attr(object$x[[i]]$terms, "term.labels")[tls]
    }
    
    if (length(tll) < 1) tll <- "1"
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
} # end of function: model.matrix.sdr


## Predict method.
predict.sdr <-
  function(object,
           newdata = NULL,
           type = c("link", "parameter"),
           drop = TRUE,
           mstart = NULL,
           mstop = NULL,
           iterations = NULL,
           model = NULL,
           what = c("mean", "matrix"),
           ...) {
    
    # Needed for sdr.treshdesc variant
    if(is.null(mstop) & "iter" %in% names(object) ) mstop = object$iter
    if(is.null(mstart) & "iter" %in% names(object) ) mstart = object$iter
    
    if(is.null(newdata)) newdata <- object$X
    # All variables that are needed
    vars <- unique(unlist(object$varnames))
    if(!length(vars) > 0){
      X <- model.matrix(object = ~ 1, data = data.frame(newdata))
    } else {
      formula_all <- as.formula(paste(" ~ ", paste0(vars, collapse = "+")))
      X <- model.matrix(object = formula_all, data = data.frame(newdata))
    }
    
    if (!is.null(mstop)) {
      if (length(mstop) > 2L)
        stop("Argument mstop must be a single number!")
    }
    maxit <- nrow(object$coefficients[[1]])
    if (is.null(mstop))
      mstop <- maxit
    if (is.null(mstart))
      mstart <- mstop
    if (is.null(model))
      model <- names(object$formula)
    if (is.character(model)) {
      model <- match.arg(model, names(object$formula), several.ok = TRUE)
    }
    if (is.null(type))
      type <- "link"
    if (length(type) > 1)
      type <- type[1]
    if (is.null(what))
      what <- "mean"
    if (length(what) > 1)
      what <- what[1]
    if(is.null(iterations))
      iterations <- mstart:mstop
    
    eta <- list()
    
    if(what == "mean"){
      for (j in model) {
        eta[[j]] <- 0
        
        for (i in iterations) {
          eta[[j]] <- eta[[j]] + drop(as.matrix(X[,colnames(object$coefficients[[j]])]) %*% object$coefficients[[j]][i, ])
        }
        eta[[j]] <- eta[[j]] / (mstop - mstart + 1)
        if (type != "link") {
          linkinv <- make.link2(object$family$links[j])$linkinv
          eta[[j]] <- linkinv(eta[[j]])
        }
      }
      if ((length(eta) < 2L) & drop)
        eta <- eta[[1L]]
    } else if(what == "matrix") {
      
      for (j in model) {
        eta[[j]] <- matrix(0, nrow = nrow(X), ncol = length(iterations))
        
        for (i in 1:length(iterations)) {
          eta[[j]][,i] <- drop(as.matrix(X[,colnames(object$coefficients[[j]])]) %*% object$coefficients[[j]][iterations[i], ])
          #str(eta[[j]][,i])
          if (type != "link") {
            linkinv <- make.link2(object$family$links[j])$linkinv
            eta[[j]][,i] <- linkinv(eta[[j]][,i])
          }
        }
        
      }
      # if ((length(eta) < 2L) & drop)
      #   eta <- eta[[1L]]
    } else {
      stop("what must be either mean or matrix!")
    }
    
    return(eta)
  }

# Extract coefficients.
coef.sdr <- function(object,
                           model = NULL,
                           refit = FALSE,
                           mstop = NULL,
                           mstart = mstop,
                           ...) {
  
  # Needed for sdr.treshdesc variant
  if (is.null(mstop) & "iter" %in% names(object)) mstop <- object$iter
  
  coef_fun <- function(...) {
    if (!is.null(mstop)) {
      if (length(mstop) > 2L)
        stop("Argument mstop must be a single number!")
    }
    keep = NULL
    if (is.null(model))
      model <- names(object$formula)
    if (is.character(model))
      model <- match.arg(model, names(object$formula), several.ok = TRUE)
    else
      model <- names(object$formula)[model]
    maxit <- nrow(object$coefficients[[1]])
    if (is.null(keep))
      keep <- maxit
    else
      keep <- floor(maxit * (1 - keep))
    coef <- list()
    for (j in model)
      coef[[j]] <- apply(object$coefficients[[j]][if (is.null(mstop)) keep:maxit
                                                  else mstart:mstop, , drop = FALSE], 2, mean)
    #if (length(coef) < 2L)
    #  coef <- coef[[1L]]

    return(coef)
  }
  
  if (refit) {
    if (is.null(mstop)) {
      mstop <- nrow(object$coefficients[[1]])
      cat("Last iteration is used as 'mstop' is not provided.")
    }
    if (is.null(mstart))
      mstart <- mstop
    if (!mstart == mstop)
      cat("'mstart' is ignored. Coefs are based only on 'mstop'")
    
    coef <- coef_fun(mstart = mstop, mstop = mstop, ...)
    nc   <- names(coef)
    
    names(coef) <- nc
    
    coef <- lapply(nc,
      FUN = function(i) {
        ind <- coef[[i]] != 0
        ind[1] <- TRUE  # Intercept is always true
        coef[[i]] <- coef[[i]][ind]
      }
    )
    names(coef) <- nc
  } else {
    coef <- coef_fun(...)
  }
  
  return(coef)
}


# coef.sdr <- function(object, ...) {
#
# }

### update function for refitting
# newformula
newformula <- function(object,
                       mstop = NULL,
                       name = NULL) {
  if (is.null(mstop)) {
    mstop <- nrow(object$coefficients[[1]])
    cat("Last iteration is used as 'mstop' is not provided.")
  }
  
  coef <- coef.sdr(object = object, mstop = mstop, refit = TRUE)
  nc   <- names(coef)
  coef <- lapply(nc,
    FUN = function(i) {
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
    }
  )
  names(coef) <- nc
  
  return(coef)
}




## Summary method.
# Umschreiben so dass iter automatisch uebernommen wird. Parameterverteilung
# interessiert dann keinen mehr
summary.sdr <- function(object,
                              digits = max(3, getOption("digits") - 3),
                              mstart = round(0.5 * length(object$logLik)),
                              mstop  = length(object$logLik),
                              ...) {
  # for computing parameter summary
  parsum <- function(d, vec = mstart:mstop) {
    # TODO(R): Fails if you only have less than 4 iterations!
    dd <- t(d)[, 1:4]
    if (ncol(d) == 1) {
      dd[1]   <- mean(d[vec, ])
      dd[2:4] <- quantile(d[vec, ], c(0.025, 0.5, 0.975))
      dd      <- matrix(dd, ncol = 4)
      rownames(dd) <- "(Intercept)"
    } else {
      if (length(vec) == 1) {
        dd[, 1]   <- t(apply(d[c(vec, vec), ], 2, mean))
        dd[, 2:4] <- t(apply(d[c(vec, vec), ], 2, quantile, c(0.025, 0.5, 0.975)))
        # dd[,5] <- t(d[nrow(d),])
      } else {
        dd[, 1]   <- t(apply(d[vec, ], 2, mean))
        dd[, 2:4] <- t(apply(d[vec, ], 2, quantile, c(0.025, 0.5, 0.975)))
        # dd[,5] <- t(d[nrow(d),])
      }
      
    }
    colnames(dd) <- c("mean", "2.5%", "50%", "97.5%")
    return(dd)
  }
  
  
  cat("\nCall:\n")
  print(object$call)
  cat("---\n")
  print(object$family, full = FALSE)
  cat("*---\n")
  
  cat("\nmstart:\n")
  object$mstart <- mstart
  print(object$mstart)
  cat("---\n")
  
  cat("\nmstop:\n")
  object$mstop <- mstop
  print(object$mstop)
  cat("---\n")
  
  for (i in names(object$formula)) {
    if (length(names(object$formula)) > 1) {
      cat("\nFormula", i, ":\n")
      cat("---\n")
      print(unname(object$formula[i]))
      cat("---\n")
    }
    if (length(names(object$formula)) == 1) {
      cat("\nFormula mu:\n")
      cat("---\n")
      print(object$formula)
      cat("---\n")
    }
    if (!is.null(object$coefficients[[i]])) {
      cat("Parametric coefficients:\n")
      object$parsum[[i]] <- round(parsum(object$coefficients[[i]]), digits)
      print(object$parsum[[i]])
      # printCoefmat(parsum(object$coefficients[[i]], digits = digits)
      cat("---\n")
    }
  }
  
  cat("\nlogLik-mstart-mstop:\n")
  print(object$logLik[c(object$mstart, object$mstop)])
  cat("---\n")
  return(invisible(object[6:12]))
}


## From bamlss.
print_bamlss_formula <- function(x, ...) {
  if (!inherits(x, "list") & !inherits(x, "bamlss.formula")) {
    print(x)
  } else {
    nx <- names(x)
    if (is.null(nx))
      nx <- as.character(1:length(x))
    for (i in seq_along(x)) {
      cat("Formula ", nx[i], ":\n---\n", sep = "")
      if (inherits(x[[i]], "list") & "h1" %in% names(x[[i]])) {
        for (j in seq_along(x[[i]])) {
          cat("h", j, ": ", sep = "")
          attr(x[[i]][[j]], "name") <- NULL
          attr(x[[i]][[j]]$formula, ".Environment") <- NULL
          if (is.character(x[[i]][[j]]$formula)) {
            cat(x[[i]][[j]]$formula, "\n")
          } else {
            print(x[[i]][[j]]$formula, showEnv = FALSE)
          }
        }
      } else {
        attr(x[[i]], "name") <- NULL
        attr(x[[i]]$formula, "name") <- NULL
        attr(x[[i]]$formula, ".Environment") <- NULL
        if ("formula" %in% names(x[[i]])) {
          if (is.character(x[[i]]$formula)) {
            cat(x[[i]]$formula, "\n")
          } else {
            print(x[[i]]$formula, showEnv = FALSE)
          }
        } else {
          print(x[[i]])
        }
      }
      if (i < length(x)) cat("\n")
    }
  }
  invisible(NULL)
}

##### Print needs rework Print method for summary.
print.summary.sdr <- function(x, ...) {
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





logLik.sdr <- function(object,
                             mstart = 1,
                             mstop = length(object$logLik),
                             all = TRUE,
                             ...) {
  ll <- object$logLik[mstart:mstop]
  
  ### degrees of freedom
  df <- rep(0, mstop - mstart + 1)
  beta <- object$coefficients
  for (i in 1:length(beta)) {
    if (length(data.frame(beta[i])) != 1)
      df <- df + rowSums(ifelse(data.frame(beta[i])[mstart:mstop, ] == 0, 0, 1))
  }
  if (!all) {
    ll <- mean(ll)
    df <- mean(df)
  }
  attr(ll, "df") <- df
  class(ll) <- c("logLik", "logLik.sdr")
  return(ll)
}

AIC.sdr <- function(object, K = 2, ...) {
  ll <- logLik(object, ...)
  ic <- as.numeric(-2 * ll + K * attr(ll, "df"))
  
  class(ic) <- "sdr_AIC"
  return(ic)
}

BIC.sdr <- function(object, ...) {
  structure(AIC.sdr(object, K = log(object$nobs), ...),
            class = "sdr_BIC")
}

plot.sdr <- function(x,
                           which = c("all", "coefficients", "AIC"),
                           K = 2,
                           bw = 0,
                           spar = TRUE,
                           ...) {
  if (length(which) > 1)
    which <- which[1]
  if (which == "all") {
    ll <- logLik(x, ...)
    ic <- as.numeric(-2 * ll + K * attr(ll, "df"))
    bic.min <-
      ifelse(bw > 0, which.min(ma(ic, order = bw)), which.min(ic))
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
      matplot(x$coefficients[[j]], type = "l", lty = 1, main = j, xlab = "Iteration", ylab = "Coefficients")
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
      matplot(x$coefficients[[j]], type = "l", lty = 1, cex.main = 1.5, cex.lab = 1.5,
              cex.axis = 1.3, xlab = "Iteration", ylab = "", main = j)
      # for(i in colnames(x$coefficients[[j]][,-1]))
      abline(v = bic.min, lty = 2)
      abline(v = bic.min - bw, lty = 2)
      l    <- nrow(x$coefficients[[j]])
      what <- x$coefficients[[j]][l, ] != 0
      nam  <- colnames(x$coefficients[[j]])[what]
      plab <- x$coefficients[[j]][l, nam]
      rang <- abs(diff(par("usr")[3:4]))
      o    <- order(plab, decreasing = TRUE)
      nam1 <- nam <- nam[o]
      plab <- plab[o]

      if (length(plab) > 1) {
        for (i in 1:(length(plab) - 1)) {
          dp <- abs(plab[i] - plab[i + 1]) / rang
          if (dp <= 0.025) {
            nam[i + 1] <- paste(c(nam[i], nam[i + 1]), collapse = ",")
            nam[i] <- ""
          }
        }
      }
      
      for (i in 1:length(nam)) {
        pos <- plab[i]
        if (pos != 0)
          text(x = par("usr")[2], labels = nam[i], y = pos, pos = 4, xpd = TRUE, cex = 1.1)
      }
    }
    
  }
  
  if (which == "AIC") {
    ll <- logLik(x, ...)
    ic <- as.numeric(-2 * ll + K * attr(ll, "df"))
    bic.min <- ifelse(bw > 0, which.min(ma(ic, order = bw)), which.min(ic))
    
    plot(y = ic, x = 1:length(ic), xlab = "Iteration", ylab = "AIC")
    abline(v = bic.min, lty = 2)
    abline(v = bic.min - bw, lty = 2)
  }
  
  return(invisible(NULL))
}


residuals.sdr <- function(object,
                                type = c("quantile", "response"),
                                nsamps = NULL,
                                ...) {
  family <- object$family
  ynam   <- colnames(object[["y"]])
  if(is.null(ynam)) ynam <- object[["y"]]
  
  if (!is.null(family$residuals)) {
    res <- family$residuals(object, type = type, nsamps = nsamps, ...)
    if (length(class(res)) < 2) {
      if (inherits(res, "numeric"))
        class(res) <- c("sdr_residuals", class(res))
    }
  } else {
    type <- match.arg(type)
    y    <- NULL
    if (!is.null(nsamps)) {
      y   <- nsamps[, ynam]
      par <- predict(object, newdata = nsamps, drop = FALSE, ...)
    } else {
      y   <- unlist(object$y)
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
        # TODO: Better if (!is.null(dim(y))) and remove this empty condition 
      } else {
        y <- y[-nas, ]
      }
    }
    nod <- is.null(dim(par[[1L]]))
    for (j in family$names) {
      if (!nod) par[[j]] <- as.matrix(par[[j]])
      par[[j]] <- make.link2(family$links[j])$linkinv(par[[j]])
    }
    if (type == "quantile") {
      if (is.null(family$p)) {
        type <- "response"
        warning(paste("no $p() function in family '", family$family,
                "', cannot compute quantile residuals, computing response resdiuals instead!",
                sep = ""))
      } else {
        discrete <- FALSE
        if (!is.null(family$type)) {
          if (tolower(family$type) == "discrete") discrete <- TRUE
        }
        if (family$family == "binomial") discrete <- TRUE

        if (discrete) {
          ymin <- min(y, na.rm = TRUE)
          a    <- family$p(ifelse(y == ymin, y, y - 1), par)
          a    <- ifelse(y == ymin, 0, a)
          b    <- family$p(y, par)
          u    <- runif(length(y), a, b)
          u    <- ifelse(u > 0.999999, u - 1e-16, u)
          u    <- ifelse(u < 1e-06, u + 1e-16, u)
          res  <- qnorm(u)
        } else {
          prob <- family$p(y, par)
          res  <- qnorm(prob)
          if (any(isnf <- !is.finite(res))) {
            warning("non finite quantiles from probabilities, set to NA!")
            res[isnf] <- NA
          }
        }
        attr(res, "type") <- "Quantile"
      }
    }

    if (type == "response") {
      mu  <- if (is.null(family$mu)) function(par, ...) par[[1]] else family$mu
      res <- y - mu(par)
      attr(res, "type") <- "Response"
    }

    class(res) <- c("sdr_residuals", class(res))
  }

  if (any(j <- !is.finite(res))) res[j] <- NA
  return(res)
}


c95 <- function (x) {
  qx <- quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  return(c(qx[1], Mean = mean(x, na.rm = TRUE), qx[2]))
}

plot.sdr_residuals <- function(x,
                                     which = c("hist-resid", "qq-resid", "wp"),
                                     spar = TRUE,
                                     ...) {
  which.match <- c("hist-resid", "qq-resid", "wp")
  if (!is.character(which)) {
    if (any(which > 3L)) which <- which[which <= 3L]
    which <- which.match[which]
  } else
    which <- which.match[pmatch(tolower(which), which.match)]
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
        rdens     <- density(as.numeric(x2), na.rm = TRUE)
        rh        <- hist(as.numeric(x2), plot = FALSE)
        args$ylim <- c(0, max(c(rh$density, rdens$y)))
        args$freq <- FALSE
        args$x    <- as.numeric(x2)
        args      <- delete.args("hist.default",
                                 args,
                                 package = "graphics",
                                 not = c("xlim", "ylim"))
        if (is.null(args$xlab))
          args$xlab <- if (is.null(type)) "Residuals"
        else
          paste(type, "residuals")

        if (is.null(args$ylab))
          args$ylab <- "Density"

        if (is.null(args$main)) {
          args$main <- paste("Histogramm and density", if (!is.null(cn[j])) paste(":", cn[j]) else NULL)
        }
        ok <- try(do.call("hist", args))
        if (!inherits(ok, "try-error")) {
            lines(rdens)
        }
        box()
      }
      if (w == "qq-resid") {
        if (ncol(x) > 1) {
          x2 <- t(apply(x, 1, c95))
          args$x <- NULL
          args$plot.it <- FALSE
          args <- delete.args("qqnorm.default",
                              args,
                              package = "stats",
                              not = c("col", "pch", "cex"))
          if (is.null(args$main)) {
            args$main <- paste("Normal Q-Q plot", if (!is.null(cn[j])) paste(":", cn[j]) else NULL)
          }
          args$y <- x2[, "Mean"]
          mean   <- do.call(qqnorm, args)
          args$y <- x2[, "2.5%"]
          lower  <- do.call(qqnorm, args)
          args$y <- x2[, "97.5%"]
          upper  <- do.call(qqnorm, args)
          ylim   <- range(c(as.numeric(mean$y),
                            as.numeric(lower$y),
                            as.numeric(upper$y)), na.rm = TRUE)

          args$plot.it <- TRUE
          if (is.null(args$ylim)) args$ylim <- ylim
          args$y <- x2[, "Mean"]
          mean   <- do.call(qqnorm, args)
          if (is.null(args$ci.col)) args$ci.col <- 1
          if (is.null(args$ci.lty)) args$ci.lty <- 2
          lines(lower$x[order(lower$x)], lower$y[order(lower$x)],
                lty = args$ci.lty, col = args$ci.col)
          lines(upper$x[order(upper$x)], upper$y[order(upper$x)],
                lty = args$ci.lty, col = args$ci.col)
          args$y <- x2[, "Mean"]
          qqline(args$y)
        } else {
          args$y <- as.numeric(x)
          args$x <- NULL
          args <- delete.args("qqnorm.default",
                              args,
                              package = "stats",
                              not = c("col", "pch", "xlim", "ylim", "cex"))
          if (is.null(args$main)) {
            args$main <- paste("Normal Q-Q plot", if (!is.null(cn[j])) paste(":", cn[j]) else NULL)
          }
          args$plot.it <- !add
          ok <- try(do.call(qqnorm, args))
          if (add) {
            args <- delete.args("points.default",
                                list(...),
                                package = "graphics",
                                not = c("col", "pch", "cex"))
            points(ok$x, ok$y, pch = args$pch, col = args$col, cex = args$cex)
          } else {
            if (!inherits(ok, "try-error")) qqline(args$y)
          }
        }
      }

      if (w == "wp") {
        xlo <- xup <- NULL
        if (ncol(x) > 1) {
          x2  <- t(apply(x, 1, c95))
          xlo <- x2[, "2.5%"]
          xup <- x2[, "97.5%"]
          x2  <- x2[, "Mean"]
        } else {
          x2 <- x
        }

        d     <- qqnorm(x2, plot = FALSE)
        probs <- c(0.25, 0.75)
        y3    <- quantile(x2, probs, type = 7, na.rm = TRUE)
        x3    <- qnorm(probs)
        slope <- diff(y3) / diff(x3)
        int   <- y3[1L] - slope * x3[1L]
        d$y   <- d$y - (int + slope * d$x)
        if (!is.null(xlo)) {
          d2    <- qqnorm(xlo, plot = FALSE)
          d$ylo <- d2$y - d2$x
          d$xlo <- d2$x
          d2    <- qqnorm(xup, plot = FALSE)
          d$yup <- d2$y - d2$x
          d$xup <- d2$x
        }
        level <- 0.95
        xlim  <- max(abs(d$x), na.rm = TRUE)
        xlim  <- c(-xlim, xlim)
        ylim  <- max(abs(c(as.numeric(d$y),
                           as.numeric(d$ylo),
                           as.numeric(d$yup))), na.rm = TRUE)

        ylim <- c(-ylim, ylim)
        if (!is.null(args$ylim2)) ylim <- args$ylim2
        if (!is.null(args$xlim2)) xlim <- args$xlim2
        z    <- seq(xlim[1] - 10, xlim[2] + 10, 0.25)
        p    <- pnorm(z)
        se   <- (1 / dnorm(z)) * (sqrt(p * (1 - p) / length(d$y)))
        low  <- qnorm((1 - level) / 2) * se
        high <- -low

        args <- list(...)
        if (is.null(args$col))  args$col <- 1
        if (is.null(args$pch))  args$pch <- 1
        if (is.null(args$cex))  args$cex <- 1
        if (is.null(args$ylab)) args$ylab <- "Deviation"
        if (is.null(args$xlab)) args$xlab <- "Unit normal quantile"
        if (add) {
          points(d$x, d$y, col = args$col, pch = args$pch, cex = args$cex)
        } else {
          if (is.null(args$main)) {
            args$main <- paste("Worm plot", if (!is.null(cn[j])) paste(":", cn[j]) else NULL)
          }
          plot(d$x, d$y, ylab = args$ylab, xlab = args$xlab, main = args$main, xlim = xlim, ylim = ylim, col = NA, type = "n")
          grid(lty = "solid")
          abline(0, 0, lty = 2, col = "lightgray")
          abline(0, 1e+05, lty = 2, col = "lightgray")
          lines(z, low, lty = 2)
          lines(z, high, lty = 2)
          points(d$x, d$y, col = args$col, pch = args$pch, cex = args$cex)
        }

        if (!is.null(xlo)) {
          if (is.null(args$ci.col)) args$ci.col <- 4
          if (is.null(args$ci.lty)) args$ci.lty <- 2
          i <- order(d$xlo)
          lines(d$ylo[i] ~ d$xlo[i], lty = args$ci.lty, col = args$ci.col)
          i <- order(d$xup)
          lines(d$yup[i] ~ d$xup[i], lty = args$ci.lty, col = args$ci.col)
        }
      }
    }
  }

  return(invisible(NULL))
} # end of function: plot.sdr_residuals



delete.args <- function(fun = NULL,
                        args = NULL,
                        not = NULL,
                        package = NULL) {
  if (is.character(fun) & !is.null(package))
    fun <- eval(parse(text = paste(
                package, paste(rep(":", 3), collapse = ""),
                fun, sep = "")))

  nf <- names(formals(fun))
  na <- names(args)
  for (elmt in na) {
    if (!elmt %in% nf) {
      if (!is.null(not)) {
        if (!elmt %in% not) args[elmt] <- NULL
      } else {
        args[elmt] <- NULL
      }
    }
  }

  return(args)
}


rps <- function(obs, pred) {
  nr   <- nrow(pred)
  test <- apply(pred, 1, sum)
  id   <- is.finite(obs) & is.finite(test)  # & test <= 1 & test > 0.8
  obs  <- obs[id]
  pred <- matrix(pred[id, ], nrow = nr)
  OBS  <- matrix(0, nrow = length(obs), ncol = ncol(pred))
  k    <- nrow(OBS)
  k1   <- ncol(OBS)
  k2   <- ncol(pred)
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
  RPS <- sum(apply((PRED - OBS2) ^ 2, 1, sum)) / (ncol(pred) - 1)
  
  return(RPS)
}


ma <- function(x, order = 20) {
  ma1 <- filter(x, rep(1 / order, order), sides = 1)
  ma2 <- rev(filter(rev(x), rep(1 / order, order), sides = 1))
  return(ifelse(is.na(ma1), ma2, ma1))
}


# formula updating after variable selection step
f.update <- function(model,
                     mstop,
                     bb = 0,
                     names = NULL,
                     max.char = 24) {
    l <- length(model[["coefficients"]])
    f <- list()
    for (i in 1:l) {
      dd <- 2:ncol(model[["coefficients"]][[i]])
      #if (mstop > bb)  int <- (mstop - bb):mstop
      #if (mstop <= bb) int <- mstop:(mstop + bb)
      int <- if (mstop > bb) (mstop - bb):mstop else mstop:(mstop + bb)
      
      coef <- t(colMeans(model[["coefficients"]][[i]][int, dd]))
      
      x  <- colnames(coef)[abs(coef) > 0]
      nn <- c(1:max.char, "inf_cin", "cin", "cloud", "no_cloud", paste0(0, 1:max.char))
      if (!is.null(names)) {
        for (j in names) {
          if (any(paste0(j, nn) %in% x)) x <- c(setdiff(x, paste0(j, nn)), j)
        }
      }
      
      if (length(x) == 0) x = "1"
      
      f[i] <- ifelse(i == 1, paste(names(model[["y"]]), " ~ ", paste(x, collapse = "+")),
                     paste(" ~ ", paste(x, collapse = "+")))
      
    }
    
    f1 <- lapply(f, function(i) as.formula(i))
    return(f1)
  }


f2.update <- function(model, mstop, bb = 0) {
  dd         <- 2:ncol(model[["coefficients"]][["mu"]])
  mu.coef    <- t(colMeans(model[["coefficients"]][["mu"]][mstop - bb:mstop, dd]))
  sigma.coef <- t(colMeans(model[["coefficients"]][["sigma"]][mstop - bb:mstop,
                                                  dd]))
  # nu.coef <- t(colMeans(model[['coefficients']][['nu']][mstop-bb:mstop,dd]))
  
  mu.x    <- colnames(mu.coef)[abs(mu.coef) > 0]
  sigma.x <- colnames(sigma.coef)[abs(sigma.coef) > 0]
  # nu.x <- colnames(mu.coef)[abs(nu.coef) > 0]
  
  if (length(mu.x) == 0)    mu.x = "1"
  if (length(sigma.x) == 0) sigma.x = "1"
  # if(length(nu.x) == 0) nu.x = '1'
  
  f.mu    <- as.formula(paste("y ~ ", paste(mu.x, collapse = "+")))
  f.sigma <- as.formula(paste(" ~ ",  paste(sigma.x, collapse = "+")))
  # f.nu <- as.formula(paste(' ~ ', paste(nu.x, collapse = '+')))
  
  f <- list(f.mu, f.sigma)
  return(f)
}


f3.update <- function(model, mstop, bb = 0) {
  dd         <- 2:ncol(model[["coefficients"]][["mu"]])
  mu.coef    <- t(colMeans(model[["coefficients"]][["mu"]][mstop - bb:mstop, dd]))
  sigma.coef <- t(colMeans(model[["coefficients"]][["sigma"]][mstop - bb:mstop,
                                                  dd]))
  nu.coef <- t(colMeans(model[["coefficients"]][["nu"]][mstop - bb:mstop, dd]))
  
  mu.x    <- colnames(mu.coef)[abs(mu.coef) > 0]
  sigma.x <- colnames(sigma.coef)[abs(sigma.coef) > 0]
  nu.x    <- colnames(nu.coef)[abs(nu.coef) > 0]
  
  if (length(mu.x) == 0)     mu.x    <- "1"
  if (length(sigma.x) == 0)  sigma.x <- "1"
  if (length(nu.x) == 0)     nu.x    <- "1"
  
  f.mu    <- as.formula(paste("y ~ ", paste(mu.x, collapse = "+")))
  f.sigma <- as.formula(paste(" ~ ",  paste(sigma.x, collapse = "+")))
  f.nu    <- as.formula(paste(" ~ ",  paste(nu.x, collapse = "+")))
  
  f <- list(f.mu, f.sigma, f.nu)
  return(f)
}

f4.update <- function(model, mstop, bb = 0) {
  dd         <- 2:ncol(model[["coefficients"]][["mu"]])
  mu.coef    <- t(colMeans(model[["coefficients"]][["mu"]][mstop - bb:mstop, dd]))
  sigma.coef <- t(colMeans(model[["coefficients"]][["sigma"]][mstop - bb:mstop,
                                                  dd]))
  nu.coef  <- t(colMeans(model[["coefficients"]][["nu"]][mstop - bb:mstop, dd]))
  tau.coef <- t(colMeans(model[["coefficients"]][["tau"]][mstop - bb:mstop, dd]))
  
  mu.x    <- colnames(mu.coef)[abs(mu.coef) > 0]
  sigma.x <- colnames(sigma.coef)[abs(sigma.coef) > 0]
  nu.x    <- colnames(nu.coef)[abs(nu.coef) > 0]
  tau.x   <- colnames(tau.coef)[abs(tau.coef) > 0]
  
  if (length(mu.x) == 0)     mu.x    <- "1"
  if (length(sigma.x) == 0)  sigma.x <- "1"
  if (length(nu.x) == 0)     nu.x    <- "1"
  if (length(tau.x) == 0)    tau.x   <- "1"
  
  f.mu    <- as.formula(paste("y ~ ", paste(mu.x, collapse = "+")))
  f.sigma <- as.formula(paste(" ~ ",  paste(sigma.x, collapse = "+")))
  f.nu    <- as.formula(paste(" ~ ",  paste(nu.x, collapse = "+")))
  f.tau   <- as.formula(paste(" ~ ",  paste(tau.x, collapse = "+")))
  
  return(list(f.mu, f.sigma, f.nu, f.tau))
}


rootogram <- function(model,
                      newdata = NULL,
                      counts = NULL,
                      maxit,
                      bb = 20,
                      max.k = NULL,
                      max.ylim = NULL,
                      main = NULL,
                      score = F) {
    range <- c(maxit - bb, maxit)
    #print(range)
    
    if (is.null(newdata)) {
      newdata <- model[["model"]]
      newdata <- data.frame(newdata, model$y)
    }
    
    pred <- predict(object = model, newdata = newdata, type = "parameter",
                    mstart = range[1], mstop = range[2])
    
    mat <- model[["family"]]$d(0, pred)
    k   <- ifelse(is.null(max.k), max(30, newdata$y), max.k)
    for (i in 1:k) {
      p   <- model[["family"]]$d(i, pred)
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
    y   <- rep(0, k + 1)
    
    if (is.null(counts)) {
      for (i in 1:(k + 1))
        y[i] <- sqrt(sum(newdata$y == i - 1))
    }
    if (!is.null(counts)) {
      for (i in 1:(k + 1))
        y[i] <- sqrt(sum(counts == i - 1))
    }
    
    diff <- mat - y
    names(diff) <- names(mat) <- NULL
    
    ylim = c(min(diff, 0), max(y, mat))
    if (!is.null(max.ylim))
      ylim[2] <- max.ylim
    barplot(diff ~ c(0:k), ylab = expression(sqrt(frequency)), ylim = ylim,
            main = main, xlab = "flash counts", space = 0)
    x <- c(0:k) + 0.5
    lines(y ~ x, col = "red", type = "l")
    lines(mat ~ x, col = "blue", type = "l")
    legend("topright", horiz = F, bty = "n",
           legend = c( expression(sqrt(pred)), expression(sqrt(obs)), expression(sqrt(pred) - sqrt(obs))),
           cex = 0.9, col = c("blue", "red", "grey"), title = NULL, lwd = 1)
    
} # end of function: rootogram



prob.stop <- function(model, cyclic = FALSE, nnoise) {
  iter.first.p <- NULL
  for (i in names(model[["coefficients"]])) {
    ifp <- which.max(rowSums(model[["coefficients"]][[i]][, paste0("p", 1:(6 + nnoise))] != 0) != 0)
    iter.first.p <- c(iter.first.p, ifp)
  }
  r <- iter.first.p - 1
  if (!cyclic) r <- min(r)
  
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
