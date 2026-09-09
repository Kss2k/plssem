PRIOR_OP <- ":~"
PRIOR_OP_ALIAS <- "==" # override this operator, as it's not used for anything
                       # this can be fixed properly in modsem later


mcmc_pls <- function(syntax,
                     data,
                     ...,
                     priors = list(),
                     chains = 1L,
                     iter = 2000,
                     warmup = floor(iter / 2),
                     parallel = "no",
                     ncores = chains,
                     iseed = runif(1, min = 100000, max = 999999),
                     sampler = c("Metropolis-Hastings", "Gibbs"),
                     verbose = interactive(),

                     # capture
                     bootstrap = NULL,
                     consistent = FALSE,
                     mcpls = FALSE,
                     probit = FALSE,

                     # Advanced stuff
                     rng.R = 20,
                     rng.s.start = 0.05,
                     rng.tune.pct   = 0.1, # spend 10% of the warmup tuning rng sampling
                     rng.tune.times = 10, # tune it 10 times
                     N = 20000,
                     acceptance.rate = \(d) 0.234 + 0.21 / d, # acceptance ratio by the number of dimensions in a block
                     Q = list(
                       r = \(x, s) as.vector(mvtnorm::rmvnorm(n = 1, mean = x, sigma = s)),
                       d = \(x, y, s, log = TRUE) mvtnorm::dmvnorm(matrix(y, nrow = 1), mean = x, sigma = s, log = log)
                     )) {

  sampler <- match.arg(
    tolower(sampler), c("metropolis-hastings", "gibbs")
  )

  # have priors been specified in the syntax?
  syntax <- stringr::str_replace_all(syntax, stringr::fixed(PRIOR_OP), PRIOR_OP_ALIAS)
  input <- modsem::modsemify(syntax, parentheses.as.string = TRUE)
  input[input$op == PRIOR_OP_ALIAS, "op"] <- PRIOR_OP

  inputModel  <- input[input$op != PRIOR_OP, , drop = FALSE]
  inputPriors <- input[input$op == PRIOR_OP, , drop = FALSE]

  if (NROW(inputPriors)) {
    priorsSyntax <- lapply(
      stats::setNames(
        paste0("(function(x) x |> ", inputPriors$rhs, ")"),
        nm = inputPriors$lhs
      ),
      FUN = \(expr) eval(parse(text=expr))
    )

  } else {
    priorsSyntax <- NULL

  }

  fit0 <- pls(
    syntax = parTableToSyntax(inputModel),
    data   = data,
    bootstrap = TRUE, # neccessary
    consistent = consistent,
    mcpls = mcpls,
    probit = probit,
    ...
  )

  ordered <- combinedModel(fit0)@info$ordered
  thresholdStruct0 <- combinedModel(fit0)@thresholdStruct
  data <- modelData(fit0)
  vars <- colnames(data)

  pls_stopif(!is.null(attr(data, "cluster")),
    "Bayesian estimation of Multilevel/Mixed-Effects",
    "models is not supported (yet)!"
  )

  parTableAll <- parameter_estimates(fit0)
  parTable <- getFreeParamsTable(fit0)

  parTable$par <- paste0(parTable$lhs, parTable$op, parTable$rhs)
  parTableAll$par <- paste0(parTableAll$lhs, parTableAll$op, parTableAll$rhs)

  pars <- parTable[parTable$is.free, "par"]
  thr.pars <- parTableAll[parTableAll$op == "|", "par"]

  boot.probs <- fit0@boot$boot.probs
  vcov <- vcov(fit0, use.labels = FALSE)[pars, pars, drop = FALSE]
  coef <- coef(fit0, use.labels = FALSE)[pars]

  labs <- stats::setNames(
    names(coef(fit0, use.labels = TRUE)),
    names(coef(fit0, use.labels = FALSE))
  )

  missingPriors <- setdiff(pars, names(priors))
  for (missing in missingPriors) {
    lab <- tryCatch(labs[[missing]], error = \(.e) NULL)

    if (!is.null(lab) && lab %in% names(priorsSyntax))
      priors[[missing]] <- priorsSyntax[[lab]]
    else
      priors[[missing]] <- \(...) 1 # flat prior
  }

  L <- function(x, rng = autoCorrelatedRNG(R = rng.R)) {
    fit.sim <- fit0

    parTablex <- parTable[c("lhs", "op", "rhs", "est", "is.free")]
    parTablex[parTablex$is.free, "est"] <- x

    N1 <- min(max(1, floor(rng$prop * N)), N)
    N0 <- N - N1

    if (N0 > 0) {
      sim0 <- simulateDataParTable(
        parTable = parTablex,
        N = N,
        seed = rng$rng0
      )

      sim.ov0 <- sim0$ov[seq_len(N0),,drop=FALSE]
    } else {
      sim.ov0 <- NULL
    }
    
    if (N1 > 0) {
      sim1 <- simulateDataParTable(
        parTable = parTablex,
        N = N,
        seed = rng$rng1
      )

      sim.ov1 <- sim1$ov[seq_len(N1),,drop=FALSE]
    } else {
      sim.ov1 <- NULL
    }

    sim.ov <- rbind(sim.ov0, sim.ov1)

    if (length(ordered)) {
      thresholdStruct <- thresholdStruct0

      idx0 <- withSeed(sample(NROW(boot.probs), 1), seed = rng$rng0)
      idx1 <- withSeed(sample(NROW(boot.probs), 1), seed = rng$rng1)

      probs0 <- boot.probs[idx0,,drop=TRUE]
      probs1 <- boot.probs[idx1,,drop=TRUE]

      # mix
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

    y <- coef(fit.y, use.labels=FALSE)
    l <- mvtnorm::dmvnorm(matrix(y[pars], nrow = 1), mean = coef, sigma = vcov, log = TRUE)

    attr(l, "y") <- y[pars]
    if (length(ordered))
      attr(l, "thresholds") <- y[thr.pars]

    if (N1 >= N0) {
      # inherited from sim with largest N
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

  # for (chain in chains) {
  pls_stopif(warmup >= iter, "warmup must be less than iter!")
  x  <- coef
  S0 <- vcov
  S  <- S0

  if (sampler == "metropolis-hastings") {
    blocks <- list(seq_along(x))

  } else {
    split <- as.data.frame(splitParameterNames(pars))
    op <- split$op

    oblocks <- lapply(ordered, function(ord) {
      which(split$lhs == ord & op == "|")
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

  runChain <- function(chain, p = printf) {
    # Tuning settings
    tune             <- "block"
    block.tune.times <- rng.tune.times
    block.tune.pct   <- 1 - rng.tune.pct
    block.tune.iters <- max(floor(block.tune.pct * warmup / block.tune.times), 1)
    rng.tune.iters   <- max(floor(rng.tune.pct * warmup / rng.tune.times), 1)
    tuned.times      <- 0

    # Number of adaptation steps for each block/rng
    adapt.n           <- integer(length(blocks))
    adapt.n.rng       <- 0
    acceptances.rng   <- numeric(0)
    acceptances.par   <- numeric(0)
    log.rng.s         <- log(rng.s.start)

    samples <- matrix(
      NA, nrow = iter, ncol = length(x) + length(thr.pars),
      dimnames = list(NULL, c(pars, thr.pars))
    )

    rng <- autoCorrelatedRNG(R = rng.R)

    Lx <- NULL
    Px <- NULL

    for (i in seq_len(iter)) {
      mode <- if (i > warmup) "sampling" else "warmup"
      tuned.times <- tuned.times + 1

      if (i > warmup) {
        tune <- "block"
      } else if (i <= warmup && tune == "block" && tuned.times > block.tune.iters) {
        tune <- "rng"
        tuned.times <- 0
      } else if (i <= warmup && tune == "rng" && tuned.times > rng.tune.iters) {
        tune <- "block"
        tuned.times <- 0
      }

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
      
        acc.rate.par <- recentAcceptance(acceptances.par)
        acc.rate.rng <- recentAcceptance(acceptances.rng)

        p(sprintf(
          fstring,
          chain, i, iter, mode, acc.rate.par, acc.rate.rng
        ))
      }

      if (i %% 10 == 0) {
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

      if (tune == "rng") {
        # Symmetric: Q(rng.star|rng)=Q(rng|rng.star)
        rng.star <- updateAutoCorrelatedRNG(
          rng, rho = (rng$rho + rnorm(1L, sd = exp(log.rng.s))) %% rng$R
        )

        if (is.null(Lx)) Lx <- L(x, rng = rng)
        Lx.star <- L(x, rng = rng.star)

        a <- min(log(1), Lx.star - Lx)
        k <- log(runif(1, min = 0, max = 1))

        accept <- k<=a && !is.na(k<=a)
        acceptances.rng <- c(acceptances.rng, accept)
      
        target <- acceptance.rate(1)
        adapt.n.rng <- adapt.n.rng + 1L

        # Robbins-Monro learning rate
        eta <- 0.5 / (10 + adapt.n.rng)^0.6
        log.rng.s <- log.rng.s + eta * (as.numeric(accept) - target)

        if (accept) {
          rng <- rng.star
          Lx <- Lx.star
        }

      } else for (block in seq_along(blocks)) {
        x.star <- x
        idx <- blocks[[block]]
        d <- length(idx) # dimension

        if (d <= 0)
          next # nothing to do

        # sample proposal
        scale <- exp(log.scale[[block]])
        S.star.b <- scale^2 * S[idx, idx, drop = FALSE]
        x.star[idx] <- Q$r(x = x[idx], s = S.star.b)

        # Symmetric: Q(rng.star|rng)=Q(rng|rng.star)
        rng.star <- updateAutoCorrelatedRNG(
          rng, rho = (rng$rho + rnorm(1L, sd = exp(log.rng.s))) %% rng$R
        )

        # Evaluate Lx(x*), Lx(x*) automatically truncates x* (if necessary)
        if (is.null(Lx)) Lx <- L(x, rng = rng)
        Lx.star <- L(x.star, rng = rng.star)

        # Lx(x*) automatically truncates x* (if necessary)
        x.star[idx] <- pmin(x.star[idx], attr(Lx.star, "upper")[idx])
        x.star[idx] <- pmax(x.star[idx], attr(Lx.star, "lower")[idx])

        q.star <- Q$d(x = x[idx], y = x.star[idx], s = S.star.b, log = TRUE)
        q.x    <- Q$d(x = x.star[idx], y = x[idx], s = S.star.b, log = TRUE)

        if (is.null(Px)) Px <- P(x)
        Px.star <- P(x.star)

        a <- min(log(1), (q.x - q.star) + (Lx.star + Px.star) - (Lx + Px)) # q.x/q.star = 1 for symmetric distributions
        k <- log(runif(1, min = 0, max = 1))

        accept <- k<=a && !is.na(k<=a)
        acceptances.par <- c(acceptances.par, accept)
      
        if (tune == "block") {
          target <- acceptance.rate(d)
          adapt.n[[block]] <- adapt.n[[block]] + 1L

          # Robbins-Monro learning rate
          eta <- 0.5 / (10 + adapt.n[[block]])^0.6
          log.scale[[block]] <- log.scale[[block]] + eta * (as.numeric(accept) - target)
        }

        if (accept) {
          thresholds.x <- attr(Lx.star, "thresholds")
          x[idx] <- x.star[idx]

          Px <- Px.star
          Lx <- Lx.star
          rng <- rng.star

        } else {
          thresholds.x <- attr(Lx, "thresholds")
        }
      }

      samples[i,] <- c(x, thresholds.x)
    }

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

  list(
    fit = fit0,
    results = results,
    samples = samples
  )
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
