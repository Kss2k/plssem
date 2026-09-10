PRIOR_OP <- ":~"

mcmc_pls <- function(syntax,
                     data,
                     ...,
                     chains = 1L,
                     iter = 2000,
                     warmup = floor(iter / 2),
                     parallel = "no",
                     ncores = chains,
                     iseed = runif(1, min = 100000, max = 999999),
                     sampler = c("Metropolis-Hastings", "Gibbs"),
                     verbose = interactive(),

                     warm.start = TRUE,
                     noise.correction = TRUE,
                     noise.correction.R = min(max(floor(warmup/4L), 200), 1000),

                     # capture
                     bootstrap = NULL,
                     consistent = FALSE,
                     mcpls = NULL,
                     probit = FALSE,

                     # Advanced stuff
                     rng.R = 20,
                     rng.s.start = 0.20,
                     rng.tune.pct = 0.1,
                     rng.acceptance.rate = 0.44,
                     N = 20000,
                     acceptance.rate = \(d) 0.234 + 0.21 / d,
                     Q = list(
                       r = \(x, s) as.vector(mvtnorm::rmvnorm(n = 1, mean = x, sigma = s)),
                       d = \(x, y, s, log = TRUE) mvtnorm::dmvnorm(matrix(y, nrow = 1), mean = x, sigma = s, log = log)
                     )) {

  sampler <- match.arg(tolower(sampler), c("metropolis-hastings", "gibbs"))

  pls_stopif(
    length(rng.R) != 1L || !is.finite(rng.R) || rng.R < 2,
    "`rng.R` must be at least 2!"
  )
  pls_stopif(
    length(rng.s.start) != 1L || !is.finite(rng.s.start) || rng.s.start <= 0,
    "`rng.s.start` must be positive!"
  )
  pls_stopif(
    length(rng.tune.pct) != 1L || !is.finite(rng.tune.pct) ||
      rng.tune.pct < 0 || rng.tune.pct > 1,
    "`rng.tune.pct` must be between 0 and 1!"
  )
  pls_stopif(
    length(rng.acceptance.rate) != 1L || !is.finite(rng.acceptance.rate) ||
      rng.acceptance.rate <= 0 || rng.acceptance.rate >= 1,
    "`rng.acceptance.rate` must be between 0 and 1!"
  )

  # Parse priors specified in the model syntax.
  input <- modsem::modsemify(syntax, parentheses.as.string = TRUE)

  inputModel   <- input[input$op != PRIOR_OP, , drop = FALSE]
  inputModel[isFunc(inputModel$mod), "mod"] <- ""

  inputPriors0 <- input[isFunc(input$mod), , drop = FALSE]
  inputPriors1 <- input[input$op == PRIOR_OP, , drop = FALSE]
  inputPriors0 <- addReverseCovariancesToParTable(inputPriors0)


  fit.mc <- pls(
    syntax = parTableToSyntax(inputModel),
    data   = data,
    bootstrap = TRUE, # neccessary
    consistent = consistent,
    mcpls = warm.start,
    probit = probit,
    ...
  )

  if (warm.start) {
    fit0 <- combinedModel(fit.mc)@status$fit0
  } else {
    fit0 <- fit.mc
  }

  ordered <- combinedModel(fit0)@info$ordered
  thresholdStruct0 <- combinedModel(fit0)@thresholdStruct
  data <- modelData(fit0)
  vars <- colnames(data)

  pls_stopif(!is.null(attr(data, "cluster")),
    "Bayesian estimation of Multilevel/Mixed-Effects",
    "models is not supported (yet)!"
  )

  parTableAll <- getParTableEstimates(fit0)
  parTable <- getFreeParamsTable(fit0)

  parTable$par <- paste0(parTable$lhs, parTable$op, parTable$rhs)
  parTableAll$par <- paste0(parTableAll$lhs, parTableAll$op, parTableAll$rhs)

  pars <- parTable[parTable$is.free, "par"]
  thr.pars <- parTableAll[parTableAll$op == "|", "par"]

  boot.probs <- fit.mc@boot$boot.probs
  coef0 <- coef(fit0, use.labels = FALSE)[pars]
 
  if (warm.start) {
    vcov.mc <- vcov(fit.mc, use.labels = FALSE)
    vcov0 <- stats::cov(fit.mc@boot$boot, use = "complete.obs")

    vcov.mc <- vcov.mc[pars, pars, drop = FALSE]
    vcov0 <- vcov0[pars, pars, drop = FALSE]

    coef.mc <- coef(fit.mc, use.labels = FALSE)[pars]

  } else {
    vcov0 <- vcov(fit0, use.labels = FALSE)
    vcov0 <- vcov0[pars, pars, drop = FALSE]
  }

  labs <- stats::setNames(
    names(coef(fit0, use.labels = TRUE)),
    names(coef(fit0, use.labels = FALSE))
  )

  priorsIndirect <- getPriorFunctions(
    nm = inputPriors1$lhs, exprs = inputPriors1$rhs
  )

  priorsDirect <- getPriorFunctions(
    nm = getParNamesFromParTable(inputPriors0),
    exprs = inputPriors0$mod
  )

  priors <- emptyNamedList(pars)

  for (par in pars) {
    lab <- tryCatch(labs[[par]], error = \(.e) NULL)

    if (par %in% names(priorsDirect))
      priors[[par]] <- priorsDirect[[par]]
    else if (!is.null(lab) && lab %in% names(priorsIndirect))
      priors[[par]] <- priorsIndirect[[lab]]
    else
      priors[[par]] <- Uniform.pInf.nInf # flat prior
  }

  L <- function(x, rng = autoCorrelatedRNG(R = rng.R), W = 0) {
    fit.sim <- fit0

    parTablex <- parTable[c("lhs", "op", "rhs", "est", "is.free")]
    parTablex[parTablex$is.free, "est"] <- x

    N1 <- min(max(1, floor(rng$prop * N)), N)
    N0 <- N - N1

    # Both simulations use the same size so their prefixes remain stable as
    # rng$prop changes. The PLS estimator is run once on the combined sample.
    sim0 <- simulateDataParTable(parTable = parTablex, N = N, seed = rng$rng0)
    sim1 <- simulateDataParTable(parTable = parTablex, N = N, seed = rng$rng1)

    sim.ov <- rbind(
      sim0$ov[seq_len(N0), , drop = FALSE],
      sim1$ov[seq_len(N1), , drop = FALSE]
    )

    if (length(ordered)) {
      thresholdStruct <- thresholdStruct0

      idx0 <- withSeed(sample(NROW(boot.probs), 1L), seed = rng$rng0)
      idx1 <- withSeed(sample(NROW(boot.probs), 1L), seed = rng$rng1)

      probs0 <- boot.probs[idx0, , drop = TRUE]
      probs1 <- boot.probs[idx1, , drop = TRUE]

      probs <- (1 - rng$prop) * probs0 + rng$prop * probs1

      thresholdStruct@proportions <- probs

      sim.ov <- ordinalizeDataFrame(
        df = sim.ov, thresholdStruct = thresholdStruct,
        return.thr = TRUE
      )

      fit.sim@thresholdStruct <- attr(sim.ov, "thresholdStruct")
    }

    Y <- Rfast::standardise(as.matrix(sim.ov[vars]))
    S <- Rfast::cova(Y)

    # Update observed-data (lowest-order) model input
    modelData(fit.sim)  <- Y
    indCorrMatrix(fit.sim) <- S

    fit.y <- estimatePLS_Inner(fit.sim)

    y <- coef(fit.y, use.labels = FALSE)
    l <- mvtnorm::dmvnorm(
      matrix(y[pars], nrow = 1), mean = coef0, sigma = vcov0 - W, log = TRUE
    )

    attr(l, "y") <- y[pars]
    if (length(ordered))
      attr(l, "thresholds") <- y[thr.pars]

    if (N1 >= N0) {
      attr(l, "lower") <- sim1$lower
      attr(l, "upper") <- sim1$upper
    } else {
      attr(l, "lower") <- sim0$lower
      attr(l, "upper") <- sim0$upper
    }

    l
  }

  P <- function(x) {
    p <- numeric(length(x))
    for (i in seq_along(x)) {
      par <- pars[[i]]
      p[[i]] <- priors[[par]](x[[i]])
    }
    sum(log(p))
  }

  if (warm.start) {
    pls_msg_note("Using naive MAP estimates as warm start...")

    objective <- function(x) {
      
      - P(x) - mvtnorm::dmvnorm(
        matrix(x, nrow = 1),
        mean = coef.mc,
        sigma = vcov0,
        log = TRUE
      )
    }

    start <- tryCatch(stats::nlminb(
        start = coef.mc,
        objective = objective,
        control = list(iter.max = 1000, eval.max = 2000)
      )$par,
      error = \(e) coef0
    )

  } else {
    start <- coef0
  }

  if (noise.correction) {
    Y <- matrix(NA_real_, noise.correction.R, length(pars))

    pls_msg_note("Correcting for sampling error...")
    for (b in seq_len(noise.correction.R)) {
      rng.b <- autoCorrelatedRNG(R = rng.R)
      lb <- L(coef0, rng = rng.b)
      Y[b,] <- attr(lb, "y")
    }

    W <- stats::cov(Y)
    alpha.W <- safeSubtractCovariance(V = vcov0, W = W)$fraction

    pls_warnif(alpha.W <= 0.01,
      "Noise correction seems to be unstable..."
    )

  } else {
    W <- diag(length(pars))
    alpha.W <- 0

  }

  acceptProposal <- function(log.ratio)
    !is.na(log.ratio) && log(stats::runif(1L)) <= min(0, log.ratio)

  pls_stopif(warmup >= iter, "warmup must be less than iter!")
  x  <- coef0
  S0 <- vcov0
  S  <- S0

  if (sampler == "metropolis-hastings") {
    blocks <- list(seq_along(x))

  } else {
    split <- as.data.frame(splitParameterNames(pars))
    oblocks <- lapply(ordered, function(ord) {
      which(split$lhs == ord & split$op == "|")
    })

    mblock <- list(which(split$op == "=~"))
    pblock <- list(which(split$op == "~"))
    cblock <- list(setdiff( # covariances
      seq_along(x),
      c(unlist(oblocks), unlist(mblock), unlist(pblock))
    ))

    blocks <- c(oblocks, mblock, pblock, cblock)
  }


  # One log proposal-SD multiplier per block
  # If S is already a reasonable estimate of the target covariance,
  # 2.38 / sqrt(d) is a good starting point.
  log.scale <- vapply(
    X = blocks,
    FUN.VALUE = numeric(1),
    FUN = function(idx) {
      d <- length(idx)
      if (d > 0) log(2.38 / sqrt(d)) else NA_real_
    }
  )
  target.acceptance <- vapply(
    blocks, \(idx) acceptance.rate(length(idx)), numeric(1L)
  )
  target.rng.acceptance <- rng.acceptance.rate

  runChain <- function(chain, p = printf) {
    # Number of adaptation steps for each parameter block and the RNG state.
    adapt.n <- integer(length(blocks))
    adapt.n.rng <- 0L
    acceptances.rng <- logical(0L)
    acceptances.par <- logical(0L)
    log.rng.s <- log(rng.s.start)

    samples <- matrix(
      NA, nrow = iter, ncol = length(x) + length(thr.pars),
      dimnames = list(NULL, c(pars, thr.pars))
    )

    rng <- autoCorrelatedRNG(R = rng.R)
    Lx <- L(x, rng = rng, W = alpha.W * W)
    Px <- P(x)
    thresholds.x <- attr(Lx, "thresholds")

    for (i in seq_len(iter)) {
      mode <- if (i > warmup) "sampling" else "warmup"
      tune.rng <- i <= warmup &&
        floor(i * rng.tune.pct) > floor((i - 1L) * rng.tune.pct)

      if (verbose && (i == 1 || i %% max(1, floor(iter / 10)) == 0)) {
        w <- nchar(as.character(iter))
        c <- nchar(as.character(chains))

        fstring <- paste0(
          "Chain: %", c, "d, ",
          "iter: %", w, "d/%", w, "d, ",
          "mode: %8s, ",
          "ar(par): %.3f, ",
          "ar(rng): %.3f...\n"
        )
      
        acc.rate.par <- recentAcceptance(acceptances.par, n = 500)
        acc.rate.rng <- recentAcceptance(acceptances.rng, n = 250)

        p(sprintf(
          fstring,
          chain, i, iter, mode, acc.rate.par, acc.rate.rng
        ))
      }

      if (i <= warmup && i %% 10 == 0) {
        # Update S
        wpct <- warmup / iter
        n    <- iter - warmup
        n1   <- floor(i * wpct)
        n0   <- max(0, n - n1)
        sub  <- tail(samples[seq_len(i), pars, drop = FALSE], n = n1)

        if (NROW(sub) > 10) {
          S1 <- stats::cov(sub, use = "complete.obs")
          S <- ((n0 - 1) * S0 + (n1 * 1) * S1) / (n0 + n1 - 2)
        }
      }

      if (tune.rng) {
        # Symmetric: Q(rng.star|rng)=Q(rng|rng.star)
        rng.star <- updateAutoCorrelatedRNG(
          rng, rho = (rng$rho + rnorm(1L, sd = exp(log.rng.s))) %% rng$R
        )

        Lx.star <- L(x, rng = rng.star, W = alpha.W * W)
        accept <- acceptProposal(Lx.star - Lx)
        acceptances.rng <- c(acceptances.rng, accept)
        adapt.n.rng <- adapt.n.rng + 1L

        # Robbins-Monro learning rate
        eta <- 0.5 / (10 + adapt.n.rng)^0.6
        log.rng.s <- log.rng.s +
          eta * (as.numeric(accept) - target.rng.acceptance)
        log.rng.s <- min(log(0.5), max(log(1 / N), log.rng.s))

        if (accept) {
          rng <- rng.star
          Lx <- Lx.star
        }

      } else {
        for (block in seq_along(blocks)) {
          x.star <- x
          idx <- blocks[[block]]
          d <- length(idx)

          if (d <= 0)
            next

          scale <- exp(log.scale[[block]])
          S.star.b <- scale^2 * S[idx, idx, drop = FALSE]
          x.star[idx] <- Q$r(x = x[idx], s = S.star.b)

          # The wrapped random-walk proposal for rho is symmetric.
          rng.star <- updateAutoCorrelatedRNG(
            rng,
            rho = (rng$rho + rnorm(1L, sd = exp(log.rng.s))) %% rng$R
          )

          Lx.star <- L(x.star, rng = rng.star, W = alpha.W * W)

          # L() returns the bounds used to constrain the proposed parameters.
          x.star[idx] <- pmin(x.star[idx], attr(Lx.star, "upper")[idx])
          x.star[idx] <- pmax(x.star[idx], attr(Lx.star, "lower")[idx])

          q.star <- Q$d(
            x = x[idx], y = x.star[idx], s = S.star.b, log = TRUE
          )
          q.x <- Q$d(
            x = x.star[idx], y = x[idx], s = S.star.b, log = TRUE
          )
          Px.star <- P(x.star)

          log.ratio <-
            (q.x - q.star) + (Lx.star + Px.star) - (Lx + Px)
          accept <- acceptProposal(log.ratio)
          acceptances.par <- c(acceptances.par, accept)

          if (i <= warmup) {
            adapt.n[[block]] <- adapt.n[[block]] + 1L

            # Robbins-Monro learning rate
            eta <- 0.5 / (10 + adapt.n[[block]])^0.6
            log.scale[[block]] <- log.scale[[block]] +
              eta * (as.numeric(accept) - target.acceptance[[block]])
          }

          if (accept) {
            thresholds.x <- attr(Lx.star, "thresholds")
            x[idx] <- x.star[idx]
            Px <- Px.star
            Lx <- Lx.star
            rng <- rng.star
          }
        }
      }

      samples[i, ] <- c(x, thresholds.x)
    }

    attr(samples, "rng.s") <- exp(log.rng.s)
    attr(samples, "rng.acceptance.rate") <- recentAcceptance(
      acceptances.rng, n = length(acceptances.rng)
    )
    samples
  }

  workers <- if (parallel == "no") 1L else ncores
  if (workers <= 1L) {
    set.seed(iseed)
    results <- lapply(seq_len(chains), runChain)

  } else {
    oldPlan <- future::plan()
    on.exit(future::plan(oldPlan), add = TRUE)

    if (parallel == "multicore" && .Platform$OS.type == "windows") {
      pls_msg_warn(paste0(
        "The `boot.parallel = 'multicore'` option is not supported on Windows.\n",
        "Falling back to `boot.parallel = 'multisession'`."
      ))
      parallel <- "multisession"
    }

    if (parallel == "multicore") {
      future::plan(future::multicore, workers = workers)
    } else {
      future::plan(future::multisession, workers = workers)
    }

    if (verbose) {
      livePrint <- progressr::make_progression_handler(
        name = "mcmc",
        reporter = list(
          update = function(config, state, progression, ...) {
            if (length(state$message) && nzchar(state$message)) {
              cat(state$message)
              flush.console()
            }
          }
        )
      )

      results <- progressr::with_progress({
        report.every <- max(1L, floor(iter / 10))

        report.iter <- which(
          seq_len(iter) == 1L |
          seq_len(iter) %% report.every == 0L
        )

        p <- progressr::progressor(steps = chains * length(report.iter))

        future.apply::future_lapply(
          X = seq_len(chains),
          FUN = \(chain) runChain(chain, p = p),
          future.seed = iseed,
          future.packages = "plssem"
        )

      }, handlers = livePrint, enable = TRUE, delay_stdout = FALSE)

    } else {

      results <- future.apply::future_lapply(
        X = seq_len(chains),
        FUN = runChain,
        future.seed = iseed,
        future.packages = "plssem"
      )
    }
  }
  
  samples <- do.call(
    rbind, lapply(
      X = results,
      FUN = \(samples) samples[(warmup+1):(iter), , drop = FALSE]
    )
  )

  coef1 <- apply(samples, MARGIN = 2, FUN = mean, na.rm = TRUE)
  vcov1 <- stats::cov(samples, use = "complete.obs")

  parTablex <- parTable[c("lhs", "op", "rhs", "est", "is.free")]
  parTablex[parTablex$is.free, "est"] <- coef1[pars]

  fit0.combined <- combinedModel(fit0)

  fit.out <- updateModelFromFreeParTableMC(
    parTable        = parTablex,
    model           = fit0.combined,
    mc.reps         = N,
    thresholdStruct = thresholdStruct0,
    ordered         = ordered,
    seed            = NULL,
    clusterSizes    = NULL,
    clusterName     = NULL,
    full            = TRUE,
    retry           = TRUE
  )

  replace <- intersect(pars, fit.out@params$names)

  fit.out@boot$samples <- plssemMatrix(samples)
  fit.out@boot$boot <- plssemMatrix(samples)

  # add samples to bootstrap results
  fit.out@boot$chains.all <- results
  fit.out@boot$chains.sample <- lapply(
    X = results, FUN = \(X) X[(warmup+1):(iter), , drop = FALSE]
  )
  fit.out@boot$chains.warmup <- lapply(
    X = results, FUN = \(X) X[seq_len(warmup), , drop = FALSE]
  )

  fit.out@boot$vcov <- vcov1[replace, replace, drop = FALSE]

  names(fit.out@params$se) <- fit.out@params$names
  fit.out@params$se[] <- NA_real_
  fit.out@params$se[replace] <- sqrt(diag(vcov1))[replace]

  fit.out@status$fit0 <- fit0

  fit.out@status$iterations <- iter
  fit.out@status$warmups <- warmup

  fit.out@parTable <- getParTableEstimates(fit.out)

  fit.out@parTable <- addMCMC_DiagnosticsParTable(
    parTable = fit.out@parTable,
    chains = fit.out@boot$chains.sample
  )

  fit.out
}


integerSeeds <- function(n = 1) {
  floor(runif(n, min = 100000, max = 999999))
}


autoCorrelatedRNG <- function(seed = NULL, rho = 500.5, R = 1000) {
  rng.seeds <- withSeed(integerSeeds(R), seed = seed)

  updateAutoCorrelatedRNG(list(
    rng  = rng.seeds,
    rho  = NULL,
    idx0 = NULL,
    idx1 = NULL,
    prop = NULL,
    rng0 = NULL,
    rng1 = NULL,
    R    = R
  ), rho = rho)
}


updateAutoCorrelatedRNG <- function(rng, rho = rng$rho) {
  rho  <- rho %% rng$R
  base <- floor(rho)

  idx0 <- base + 1L
  idx1 <- idx0 %% rng$R + 1L
  prop <- rho - base

  rng$rho  <- rho
  rng$prop <- prop
  rng$idx0 <- idx0
  rng$idx1 <- idx1
  rng$rng0 <- rng$rng[[idx0]]
  rng$rng1 <- rng$rng[[idx1]]

  refreshInactiveRNG(rng)
}


refreshInactiveRNG <- function(rng) {
  inactive <- setdiff(seq_len(rng$R), c(rng$idx0, rng$idx1))
  rng$rng[inactive] <- integerSeeds(length(inactive))
  rng
}


withSeed <- function(expr, seed = NULL) {
  if (!is.null(seed) && exists(".Random.seed")) {
    .Random.seed.orig <- .Random.seed
    on.exit(.Random.seed <<- .Random.seed.orig)
  }

  if (!is.null(seed))
    set.seed(seed)

  expr
}


recentAcceptance <- function(x, n = 100L) {
  if (!length(x))
    return(NA_real_)

  mean(tail(x, min(n, length(x))))
}


safeSubtractCovariance <- function(V, W, tol = 1e-8) {
  alpha <- 1

  while (alpha > 1e-6) {
    corrected <- V - alpha * W

    if (!inherits(try(chol(corrected), silent = TRUE), "try-error"))
      return(list(vcov = corrected, fraction = alpha))

    alpha <- alpha / 2
  }

  list(vcov = V, fraction = 0)
}


isFunc <- function(exprs) {
  grepl("[A-z\\._0-9]\\(.*\\)", exprs)
}


getPriorFunctions <- function(nm, exprs) {
  if (!length(nm) || !length(exprs))
    return(NULL)

  stats::setNames(lapply(
    X = exprs,
    FUN = function(expr) {
      fexpr <- paste0("(function(x) x |> ", expr, ")")
      f <- eval(parse(text=fexpr))
      attr(f, "prior.label") <- expr
      f
    }
  ), nm = nm)
}


addMCMC_DiagnosticsParTable <- function(parTable, chains) {
  k <- length(chains)
  n <- NROW(chains[[1L]])

  cchains <- do.call(cbind, chains)

  parTable$rhat <- NA_real_
  parTable$ess.tail <- NA_real_
  parTable$ess.bulk <- NA_real_

  for (i in seq_len(NROW(parTable))) {
    row <- parTable[i,,drop=TRUE]
    par <- paste0(row$lhs, row$op, row$rhs)

    if (!par %in% colnames(chains[[1L]]))
      next

    chains.par <- do.call(cbind,
      lapply(chains, FUN = \(chain) chain[,par,drop=FALSE])
    )
   
    if (k > 1) rhat.par <- posterior::rhat(chains.par)
    else       rhat.par <- NA_real_

    cchain.par   <- cchains[,par,drop=TRUE]
    se.par       <- stats::sd(cchain.par, na.rm = TRUE)
    z.par        <- parTable[i, "est"] / se.par
    ci.lower.par <- quantile(cchain.par, probs = 0.025)
    ci.upper.par <- quantile(cchain.par, probs = 0.975)
    ess.bulk.par <- posterior::ess_bulk(chains.par)
    ess.tail.par <- posterior::ess_tail(chains.par)

    # non-symmetric P-value (Mplus note: https://www.statmodel.com/download/FAQ-Bootstrap%20-%20Pvalue.pdf)
    M.par        <- sum(cchain.par, na.rm = TRUE)
    B.par        <- sum(!is.na(cchain.par))
    p.value.par  <- 2 * min(M.par/B.par, 1 - M.par/B.par) # Two-sided (non-symmetric)

    parTable[i, "rhat"] <- rhat.par
    parTable[i, "ess.bulk"] <- ess.bulk.par
    parTable[i, "ess.tail"] <- ess.tail.par
    parTable[i, "se"] <- se.par
    parTable[i, "z"]  <- 
    parTable[i, "ci.lower"] <- ci.lower.par
    parTable[i, "ci.upper"] <- ci.upper.par
    parTable[i, "pvalue"]   <- p.value.par
  }

  pls_msg_warn("Z-stats are not computed correctly for MCMC models (yet)")


  parTable
}


Uniform.pInf.nInf <- function(...) {
  1 # flat prior across all
}
attr(Uniform.pInf.nInf, "prior.label") <- "Uniform[-Inf,+Inf]"
