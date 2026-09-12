PRIOR_OP <- ":~"

#' Bayesian Partial Least Squares Structural Equation Modeling
#'
#' Fits a Bayesian PLS-SEM model using a pseudo-marginal Markov chain Monte
#' Carlo sampler. The estimator is an extension of the MC-PLSc estimator.
#' Priors can be specified directly for parameters in the model syntax, or
#' specified on parameter labels with the `:~` operator.
#'
#' @param syntax Character string with \code{lavaan}-style model syntax describing
#'   both measurement (\code{=~}) and structural (\code{~}) relations. Random effects are
#'   specified with \code{(term | cluster)} statements.
#'
#' @param data A \code{data.frame} or coercible object containing the manifest
#'   indicators referenced in \code{syntax}. Ordered factors are automatically
#'   detected, but can also be supplied explicitly through \code{ordered}.
#'
#' @param ... Additional arguments passed to [pls()], such as `ordered` and
#'   `boot.R`.
#'
#' @param chains Positive integer giving the number of MCMC chains.
#'
#' @param iter Positive integer giving the total number of iterations per chain,
#'   including warmup.
#'
#' @param warmup Non-negative integer giving the number of initial iterations
#'   discarded from each chain.
#'
#' @param parallel The type of parallel operation to be used (if any). The
#'   default is \code{"no"}. \code{"multisession"} runs the chains using multiple
#'   background \code{R} sessions (works on all platforms), while \code{"multicore"}
#'   uses forked processes (not available on Windows). \code{"snow"} is kept for
#'   backwards compatibility and is treated as an alias for \code{"multisession"}.
#'   Internally this is implemented using the \code{future} package 
#'
#' @param ncores Positive integer giving the number of parallel workers.
#'
#' @param iseed Integer seed used for the chains
#'
#' @param sampler MCMC blocking scheme. `"Metropolis-Hastings"` proposes all
#'   parameters jointly; `"Gibbs"` uses groups of related parameters as blocks.
#'
#' @param verbose Should verbose output be printed?
#'
#' @param point.estimate Posterior point estimate used in the returned model.
#'   Either `"median"` or `"mean"`.
#'
#' @param warm.start Logical; use an MC-PLS fit to initialize the sampler and its
#'   proposal covariance? Defaults to \code{TRUE}.
#'
#' @param noise.correction Logical; estimate and correct for Monte Carlo noise in
#'   the synthetic likelihood?
#'
#' @param noise.correction.R Positive integer giving the number of simulations
#'   used to estimate Monte Carlo noise.
#'
#' @param bootstrap Reserved argument. Bootstrap estimation is performed
#'   internally because it is required to construct the synthetic likelihood.
#'
#' @param consistent Logical; request the PLSc consistency correction in the
#'   underlying PLS fits. Should in general be set to \code{FALSE}.
#'
#' @param mcpls Reserved argument; MC-PLS use is controlled by `warm.start`.
#'
#' @param probit Logical; use probit factor scores in the underlying PLS fits?
#'   Should in general be set to \code{FALSE}.
#'
#' @param rng.s.start Initial scale used when correlating proposed and current
#'   simulation innovations. For advanced users.
#'
#' @param rng.tune.pct Fraction determining how often the innovation scale is
#'   tuned during warmup. For advanced users.
#'
#' @param rng.acceptance.rate Target acceptance rate for innovation proposals.
#'   For advanced users.
#'
#' @param mc.reps Positive integer giving the simulated sample size used to approximate
#'   the likelihood.
#'
#' @param acceptance.rate Function mapping a proposal-block dimension to its
#'   target acceptance rate. For advanced users.
#'
#' @param Q List with functions `r` and `d` for drawing from and evaluating the
#'   parameter proposal distribution. For advanced users.
#'
#' @return A fitted `PlsModel` object. Posterior draws are available in the
#'   model's bootstrap results, including retained, warmup, and per-chain draws.
#'
#' @seealso [pls()]
#'
#' @examples
#' \dontrun{
#' m <- '
#'   X =~ load * x1 + load * x2 + load * x3
#'   Z =~ load * z1 + load * z2 + load * z3
#'   Y =~ load * y1 + load * y2 + load * y3
#' 
#'   Y ~ "dnorm(.4, .1)" * X +
#'      "dnorm(.35, .1)" * Z +
#'      "dnorm(.45, .1)" * X:Z +
#'      "dnorm(0, .005)" * X:X
#' 
#'   load :~ dnorm(.8, .5)
#' '
#'
#' set.seed(23942)
#' fit <- bpls(
#'   m, modsem::oneInt, boot.R = 500, warmup = 5000, iter = 10000,
#'   parallel = "multisession", chains = 2
#' )
#'
#' summary(fit)
#' }
#'
#' @export
bpls <- function(syntax,
                 data,
                 ...,
                 chains = 1L,
                 iter = 2000,
                 warmup = floor(iter / 2),
                 parallel = c("no", "multicore", "multisession", "snow"),
                 ncores = chains,
                 iseed = stats::runif(1, min = 100000, max = 999999),
                 sampler = c("Metropolis-Hastings", "Gibbs"),
                 verbose = interactive(),
                 point.estimate = c("median", "mean"),

                 warm.start = TRUE,
                 noise.correction = TRUE,
                 noise.correction.R = min(max(floor(warmup/4L), 200), 1000),

                 # capture
                 bootstrap = NULL,
                 consistent = FALSE,
                 mcpls = NULL,
                 probit = FALSE,

                 # Advanced stuff
                 rng.s.start = 0.20,
                 rng.tune.pct = 0.1,
                 rng.acceptance.rate = 0.44,
                 mc.reps = 20000,
                 acceptance.rate = \(d) 0.234 + 0.21 / d,
                 Q = list(
                   r = \(x, s) as.vector(mvnfast::rmvn(n = 1, mu = x, sigma = s)),
                   d = \(x, y, s, log = TRUE) mvnfast::dmvn(matrix(y, nrow = 1), mu = x, sigma = s, log = log)
                 )) {
  # Check arguments
  if (is.null(parallel)) parallel <- "no"
  parallel <- match.arg(parallel, c("no", "multicore", "multisession", "snow"))
  if (parallel == "snow") parallel <- "multisession"
  point.estimate <- match.arg(tolower(point.estimate), c("median", "mean"))
  sampler <- match.arg(tolower(sampler), c("metropolis-hastings", "gibbs"))

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

  is.hi.ord <- isTRUE(combinedModel(fit0)@info$is.high.ord)
  ordered <- combinedModel(fit0)@info$ordered
  thresholdStruct0 <- combinedModel(fit0)@thresholdStruct
  data <- modelData(fit0)
  vars <- colnames(data)

  if (isMLM(fit0)) {
    clusterSizes <- as.numeric(table(attr(data, "cluster")))
    clusterName  <- colnames(attr(data, "cluster"))
  } else {
    clusterSizes <- NULL
    clusterName  <- NULL
  }

  parTableAll <- getParTableEstimates(fit0)
  parTable <- getFreeParamsTable(fit0)

  parTable$par <- paste0(parTable$lhs, parTable$op, parTable$rhs)
  parTableAll$par <- paste0(parTableAll$lhs, parTableAll$op, parTableAll$rhs)

  pars <- parTable[parTable$is.free, "par"]
  thr.pars <- parTableAll[parTableAll$op == "|", "par"]

  empirical.vpars <- intersect(
    getParNamesFromParTable(parTableAll),
    getEmpiricalVarParsParTable(
      parTable, clusterSizes = clusterSizes, clusterName = clusterName
    )
  )

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
      priors[[par]] <- \(...) 1 # flat/no prior
  }

  compiled.info <- NULL

  L <- function(x, innovations = NULL, W = 0) {
    fit.sim <- fit0

    parTablex <- parTable[c("lhs", "op", "rhs", "est", "is.free")]
    parTablex[parTablex$is.free, "est"] <- x

    sim <- simulateDataParTable(
      parTable                = parTablex,
      N                       = mc.reps,
      check.hi.ord            = is.hi.ord,
      clusterSizes            = clusterSizes,
      clusterName             = clusterName,
      collect.empirical.vpars = TRUE,
      innovations             = innovations,
      return.innovations      = TRUE,
      compiled.info           = compiled.info
    )

    if (is.null(compiled.info))
      compiled.info <<- sim$compiled.info

    sim.ov <- sim$ov
    innovations <- sim$innovations

    if (length(ordered)) {
      thresholdStruct <- thresholdStruct0

      ordinal.key <- "ordinal-bootstrap"
      if (is.null(innovations$blocks[[ordinal.key]]))
        innovations$blocks[[ordinal.key]] <- stats::rnorm(1L)

      u <- stats::pnorm(innovations$blocks[[ordinal.key]][[1L]])
      idx <- min(NROW(boot.probs), floor(u * NROW(boot.probs)) + 1L)
      thresholdStruct@proportions <- boot.probs[idx, , drop = TRUE]

      sim.ov <- ordinalizeDataFrame(
        df = sim.ov, thresholdStruct = thresholdStruct,
        return.thr = TRUE
      )

      fit.sim@thresholdStruct <- attr(sim.ov, "thresholdStruct")
    }

    Y <- Rfast::standardise(as.matrix(sim.ov[vars]))
    S <- Rfast::cova(Y)

    if (!is.null(sim$cluster))
      attr(Y, "cluster") <- sim$cluster

    # Update observed-data (lowest-order) model input
    modelData(fit.sim)  <- Y
    indCorrMatrix(fit.sim) <- S

    fit.y <- estimatePLS_Inner(fit.sim)

    y <- coef(fit.y, use.labels = FALSE)
    l <- mvnfast::dmvn(
      matrix(y[pars], nrow = 1), mu = coef0, sigma = vcov0 - W, log = TRUE
    )

    attr(l, "y") <- y[pars]
    attr(l, "innovations") <- innovations
    if (length(ordered))
      attr(l, "thresholds") <- y[thr.pars]

    attr(l, "empirical.vpars") <- sim$empirical.vpars[empirical.vpars]
    attr(l, "lower") <- sim$lower[parTablex$is.free]
    attr(l, "upper") <- sim$upper[parTablex$is.free]

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
    objective <- function(x) {
      
      - P(x) - mvnfast::dmvn(
        matrix(x, nrow = 1),
        mu = coef.mc,
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
      lb <- L(coef0)
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
  x  <- start
  S0 <- if (warm.start) vcov.mc else vcov0
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
      NA, nrow = iter,
      ncol = length(x) + length(thr.pars) + length(empirical.vpars),
      dimnames = list(NULL, c(pars, thr.pars, empirical.vpars))
    )

    Lx <- L(x, W = alpha.W * W)
    innovations <- attr(Lx, "innovations")
    Px <- P(x)
    thresholds.x <- attr(Lx, "thresholds")
    empirical.vpars.x <- attr(Lx, "empirical.vpars")

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
        sub  <- utils::tail(samples[seq_len(i), pars, drop = FALSE], n = n1)

        if (NROW(sub) > 10) {
          S1 <- stats::cov(sub, use = "complete.obs")
          S <- ((n0 - 1) * S0 + (n1 * 1) * S1) / (n0 + n1 - 2)
        }
      }

      if (tune.rng) {
        innovations.star <- perturbInnovations(
          innovations, scale = exp(log.rng.s)
        )

        Lx.star <- L(x, innovations = innovations.star, W = alpha.W * W)
        accept <- acceptProposal(Lx.star - Lx)
        acceptances.rng <- c(acceptances.rng, accept)
        adapt.n.rng <- adapt.n.rng + 1L

        # Robbins-Monro learning rate
        eta <- 0.5 / (10 + adapt.n.rng)^0.6
        log.rng.s <- log.rng.s +
          eta * (as.numeric(accept) - target.rng.acceptance)
        log.rng.s <- min(log(0.5), max(log(1 / mc.reps), log.rng.s))

        if (accept) {
          innovations <- innovations.star
          Lx <- Lx.star
          thresholds.x <- attr(Lx.star, "thresholds")
          empirical.vpars.x <- attr(Lx.star, "empirical.vpars")
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

          innovations.star <- perturbInnovations(
            innovations, scale = exp(log.rng.s)
          )

          Lx.star <- L(
            x.star, innovations = innovations.star, W = alpha.W * W
          )

          # L() returns the bounds used to constrain the proposed parameters.
          # L() has already constrained the parameters when evaluation the
          # log likelihood. So we should constrain our parameters as well
          x.star[idx] <- pmin(x.star[idx], attr(Lx.star, "upper")[idx])
          x.star[idx] <- pmax(x.star[idx], attr(Lx.star, "lower")[idx])

          q.star <- Q$d(x = x[idx], y = x.star[idx], s = S.star.b, log = TRUE)
          q.x <- Q$d(x = x.star[idx], y = x[idx], s = S.star.b, log = TRUE)
          Px.star <- P(x.star)

          log.ratio <- (q.x - q.star) + (Lx.star + Px.star) - (Lx + Px)
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
            empirical.vpars.x <- attr(Lx.star, "empirical.vpars")
            x[idx] <- x.star[idx]
            Px <- Px.star
            Lx <- Lx.star
            innovations <- innovations.star
          }
        }
      }

      samples[i, ] <- c(x, thresholds.x, empirical.vpars.x)
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
              utils::flush.console()
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

  if (point.estimate == "median") .agg <- stats::median
  else                            .agg <- mean

  coef1 <- apply(samples, MARGIN = 2, FUN = .agg, na.rm = TRUE)
  vcov1 <- stats::cov(samples, use = "complete.obs")

  parTablex <- parTable[c("lhs", "op", "rhs", "est", "is.free")]
  parTablex[parTablex$is.free, "est"] <- coef1[pars]

  fit0.combined <- combinedModel(fit0)

  fit.out <- updateModelFromFreeParTableMC(
    parTable        = parTablex,
    model           = fit0.combined,
    mc.reps         = mc.reps,
    thresholdStruct = thresholdStruct0,
    ordered         = ordered,
    seed            = NULL,
    clusterSizes    = clusterSizes,
    clusterName     = clusterName,
    full            = TRUE,
    retry           = TRUE
  )

  replace <- intersect(pars, fit.out@params$names)

  fit.out@boot$samples <- plssemMatrix(samples)
  fit.out@boot$boot <- plssemMatrix(samples)
  fit.out@info$estimator <- paste0("B", fit.out@info$estimator)

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
    chains = fit.out@boot$chains.sample,
    priors = priors
  )

  fit.out
}


recentAcceptance <- function(x, n = 100L) {
  if (!length(x))
    return(NA_real_)

  mean(utils::tail(x, min(n, length(x))))
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


addMCMC_DiagnosticsParTable <- function(parTable, chains, priors = list()) {
  k <- length(chains)
  n <- NROW(chains[[1L]])

  cchains <- do.call(cbind, chains)

  parTable$rhat <- NA_real_
  parTable$prior <- NA_character_
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
    ci.lower.par <- stats::quantile(cchain.par, probs = 0.025)
    ci.upper.par <- stats::quantile(cchain.par, probs = 0.975)
    ess.bulk.par <- posterior::ess_bulk(chains.par)
    ess.tail.par <- posterior::ess_tail(chains.par)

    # non-symmetric P-value (Mplus note: https://www.statmodel.com/download/FAQ-Bootstrap%20-%20Pvalue.pdf)
    M.par        <- sum(cchain.par>0, na.rm = TRUE)
    B.par        <- sum(!is.na(cchain.par))
    p.value.par  <- 2 * min(M.par/B.par, 1 - M.par/B.par) # Two-sided (non-symmetric)

    if (par %in% names(priors)) {
      prior.par <- attr(priors[[par]], "prior.label")
      parTable[i, "prior"] <- ifelse(is.null(prior.par), NA_character_, prior.par)
    }

    parTable[i, "rhat"] <- rhat.par
    parTable[i, "ess.bulk"] <- ess.bulk.par
    parTable[i, "ess.tail"] <- ess.tail.par
    parTable[i, "se"] <- se.par
    parTable[i, "z"]  <- z.par
    parTable[i, "ci.lower"] <- ci.lower.par
    parTable[i, "ci.upper"] <- ci.upper.par
    parTable[i, "pvalue"]   <- p.value.par
  }

  parTable
}
