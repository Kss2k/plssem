mcpls <- function(
  fit0,
  p.start                     = fit0@info$mc.args$p.start,
  min.iter                    = fit0@info$mc.args$min.iter,
  max.iter                    = fit0@info$mc.args$max.iter,
  mc.reps                     = fit0@info$mc.args$mc.reps,
  rng.seed                    = fit0@info$mc.args$rng.seed,
  tol                         = fit0@info$mc.args$tol,
  fixed.seed                  = fit0@info$mc.args$fixed.seed,
  verbose                     = fit0@info$verbose,
  polyak.juditsky             = fit0@info$mc.args$polyak.juditsky,
  fn.args                     = fit0@info$mc.args$fn.args,
  pj.extrapolate              = fit0@info$mc.args$pj.extrapolate,
  delta.jacobian              = fit0@info$mc.args$delta.se && fit0@info$boot$bootstrap,
  delta.fixed.seed            = TRUE,
  delta.jacobian.k            = fit0@info$mc.args$delta.jacobian.k,
  diag.secant                 = fit0@info$mc.args$diag.secant,
  small.sample                = fit0@info$mc.args$small.sample,
  small.sample.max.k          = fit0@info$mc.args$small.sample.max.k,
  small.sample.point.estimate = fit0@info$mc.args$small.sample.point.estimate,
  parallel                    = fit0@info$boot$parallel,
  ncores                      = fit0@info$boot$ncores,
  ...
) {
  fit0.base <- fit0
  fit0.combined <- combinedModel(fit0.base)

  # Is the base fit admissible?
  pls_warnif(!is_admissible(fit0.combined),
    "Base fit is inadmissible!",
    "The MC-PLS algorithm might not converge to a proper solution!"
  )

  data      <- fit0.base@data
  n         <- NROW(data)
  vars      <- colnames(data)
  ordered   <- fit0@info$ordered
  is.probit <- fit0@info$is.probit
  is.hi.ord <- isTRUE(fit0.combined@info$is.high.ord)
  thresholdStruct0 <- fit0.combined@thresholdStruct
  estimator <- fit0.combined@info$path.estimator

  # Residual-covariance handling:
  #   reduced: Residual covariances are treated as constrained parameters
  #   full:    Residual covariances are treated as free parameters
  mc.rescov <- fit0.base@info$mc.args$rescov

  use.full.rescov <- switch(mc.rescov,
    full    = TRUE,
    reduced = FALSE,
    auto    = estimator == "gls",
    pls_msg_stop("Unrecognized value for `mc.rescov` argument:", mc.rescov)
  )
  
  if (small.sample) {
    pls_stopif(
      !is.numeric(small.sample.max.k) ||
      length(small.sample.max.k) != 1L ||
      !is.finite(small.sample.max.k),
      "`mc.small.sample.max.k` must be a single finite number."
    )

    max.k.int <- round(small.sample.max.k)

    pls_warnif(abs(small.sample.max.k - max.k.int) > 1e-12,
      "`mc.small.sample.max.k` should be an integer!",
      sprintf("Using `small.sample.max.k = %d`", max.k.int)
    )

    k <- max(max.k.int, 1)
    mc.reps <- min(k * n, max(mc.reps - mc.reps %% n, n))
    mc.reps.k <- max(floor(mc.reps / n), 1)
  } else {
    mc.reps.k <- 1L
  }

  par0 <- getFreeParamsTable(fit0.combined)

  if (use.full.rescov) {
    etas      <- par0$lhs[par0$op == "~"]
    is.rescov <- (par0$lhs %in% etas | par0$rhs %in% etas) &
                 par0$lhs != par0$rhs & par0$op == "~~"
    par0$is.free[is.rescov] <- TRUE

  } else {
    pls_warnif(estimator == "gls",
      "Residual covariances are not bias-corrected in `reduced` mode.",
      "Pass `mc.rescov = \"full\"` to estimate them as free parameters."
    )
  }

  par1 <- par0[c("lhs", "op", "rhs", "est", "is.free")]

  if (fixed.seed && is.null(rng.seed)) {
    rng.seed <- floor(stats::runif(1L, min = 0, max = 9999999))
    if (verbose) pls_msg_note(sprintf("Using fixed seed %i...", rng.seed))
  }

  .parTable <- function(p) {
    parx <- par1
    parx[parx$is.free, "est"] <- p
    parx
  }

  # `seed` defaults to `rng.seed`, but the Jacobian replicates (see
  # `calcMcJacobians()`) each use their own seed
  .simulate <- function(p, standardize = FALSE, seed = rng.seed) {
    simulateDataParTable(
      parTable     = .parTable(p),
      N            = mc.reps,
      seed         = seed,
      check.hi.ord = is.hi.ord,
      standardize  = standardize,
      full         = use.full.rescov
    )
  }

  .f <- function(p, thresholdStruct = thresholdStruct0, sim = NULL, seed = rng.seed) {

    if (is.null(sim)) {
      par1[par1$is.free, "est"] <- p
      sim <- simulateDataParTable(
        parTable     = par1,
        N            = mc.reps,
        seed         = seed,
        check.hi.ord = is.hi.ord,
        full         = use.full.rescov
      )
    }

    # sim.ov  <- ordinalizeDataFrame(
    #   df = sim$ov, thresholdStruct = thresholdStruct
    # )
    # sim.mat <- as.matrix(sim.ov[vars])

    fit.sim <- fit0.base
    modelStatusIsQuick(fit.sim) <- TRUE

    free <- par0$is.free
    nk <- max(floor(NROW(sim$ov) / mc.reps.k), 1)

    .estimate <- function(i) {
      offset <- (i - 1) * nk
      idx <- (offset+1):(offset+nk)

      sim.ov  <- ordinalizeDataFrame(
        df = sim$ov[idx,,drop=FALSE], thresholdStruct = thresholdStruct
      )

      X <- Rfast::standardise(as.matrix(sim.ov[vars]))

      if (is.probit) S <- getCorrMat(X, probit = TRUE, ordered = ordered)
      else           S <- Rfast::cova(X)

      # Update observed-data (lowest-order) model input
      modelData(fit.sim)  <- X
      indCorrMatrix(fit.sim) <- S

      # Thresholds are not part of the root equation. Avoid recomputing them on
      # every Robbins-Monro iteration.
      fit2 <- estimatePLS_Inner(fit.sim)
      par2 <- getFreeParamsTable(combinedModel(fit2))

      eps <- par2$est - par0$est
      eps[free]
    }

    out <- averageMcReplicates(
      k = mc.reps.k,
      point.estimate = small.sample.point.estimate,
      fun = .estimate
    )

    attr(out, "lower") <- sim$lower[free]
    attr(out, "upper") <- sim$upper[free]
    out
  }

  .g <- function(p, thresholdStruct = thresholdStruct0, sim = NULL, seed = rng.seed) {
    fit <- updateModelFromFreeParTableMC(
      parTable         = .parTable(p),
      model            = fit0.combined,
      mc.reps          = mc.reps,
      thresholdStruct  = thresholdStruct,
      ordered          = ordered,
      seed             = seed,
      sim              = sim,
      params.only      = TRUE,
      full             = use.full.rescov
    )

    fit@params$values
  }

  .fg <- function(p, thresholdStruct = thresholdStruct0, sim = NULL, seed = rng.seed) {
    if (is.null(sim))
      sim <- .simulate(p, standardize = TRUE, seed = seed)

    list(f = .f(p, thresholdStruct = thresholdStruct, sim = sim, seed = seed),
         g = .g(p, thresholdStruct = thresholdStruct, sim = sim, seed = seed))
  }

  # Starting parameters
  ok.start <- !is.null(p.start) && length(p.start) == sum(par1$is.free)
  p        <- if (ok.start) p.start else par1[par1$is.free, "est"]
  names(p) <- getParNamesFromParTable(par1)[par1$is.free]

  # If p.start was not supplied, check if any parameters were specified
  # in the model syntax/parameter table. We add the reversed covariances
  # to the partable, in case the user has specified starting values for the
  # covariances
  parTableInput <- addReverseCovariancesToParTable(
    fit0.combined@parTableInput
  )

  if (!ok.start && any(!is.na(parTableInput$start))) {
    idx <- which(!is.na(parTableInput$start))

    start <- parTableInput[idx, "start"]
    names(start) <- getParNamesFromParTable(parTableInput)[idx]

    update <- intersect(names(start), names(p))
    p[update] <- start[update]
  }

  lower <- getMcLowerBounds(par1)
  upper <- getMcUpperBounds(par1)

  mcfit <- solveMcRoot(
    p               = p,
    f               = .f,
    lower           = lower,
    upper           = upper,
    tol             = tol,
    min.iter        = min.iter,
    max.iter        = max.iter,
    verbose         = verbose,
    polyak.juditsky = polyak.juditsky,
    pj.extrapolate  = pj.extrapolate,
    fn.args         = fn.args,
    diag.secant     = diag.secant,
    ...
  )

  if (!mcfit$ok)
    modelStatus(fit0.combined)$is.admissible <- FALSE

  par1[par1$is.free, "est"] <- as.vector(mcfit$root)

  fit1.combined <- updateModelFromFreeParTableMC(
    parTable        = par1,
    model           = fit0.combined,
    mc.reps         = mc.reps,
    thresholdStruct = thresholdStruct0,
    ordered         = ordered,
    seed            = rng.seed,
    full            = use.full.rescov,
    retry           = TRUE
  )

  if (delta.jacobian) {

    if (is.null(delta.jacobian.k))
      delta.jacobian.k <- floor(fit0@info$boot$R / 50)

    pls_stopif(
      !length(delta.jacobian.k) || !is.finite(delta.jacobian.k[1L]) ||
      delta.jacobian.k[[1L]] <= 0,
      "`mc.delta.jacobian.k` must be a positive integer or `NULL`."
    )

    if (verbose) pls_msg_note("Calculating Jacobian...")

    nm <- paste0(par1$lhs, par1$op, par1$rhs)
    p0 <- stats::setNames(mcfit$root, nm[par1$is.free])
    p1 <- fit1.combined@params$values

    # Delta-method SEs assume `p0` is close to the root.
    history.f <- mcfit$history.f[, names(p0), drop = FALSE]
    history.f <- history.f[stats::complete.cases(history.f), , drop = FALSE]

    # Only use the tail (steady-state) half of the trajectory: the early,
    # far-from-root iterations have their own large, systematic swings on top
    # of MC noise, which would otherwise swamp the steady-state behaviour.
    n.hist   <- NROW(history.f)
    tail.idx <- ceiling(n.hist / 2):n.hist
    tail.f   <- history.f[tail.idx, , drop = FALSE]
    n.tail   <- NROW(tail.f)

    # Successive residuals are autocorrelated (`p` moves slowly) and under
    # `mc.fixed.seed = TRUE` consecutive iterations share the same MC error.
    # Here we use a batch-means estimator, where we split the tail into `n.batch`
    # blocks long enough to break the autocorrelation, and treat the block
    # means as approximately independent replicates. Batches are kept at least
    # `5` iterations long; with `mc.min.iter = 50` a typical run only leaves
    # ~25 steady-state iterations, and demanding longer batches would silently
    # disable the check. The `t` quantile below compensates for the resulting
    # small number of batches.
    n.batch <- max(min(floor(sqrt(n.tail)), floor(n.tail / 5L)), 0L)

    if (n.batch < 3L) {
      # Too few iterations to say anything about the steady state
      resid     <- stats::setNames(rep(NA_real_, length(p0)), names(p0))
      bad.resid <- rep(FALSE, length(p0))

    } else {
      b       <- floor(n.tail / n.batch)
      batch.f <- rowsum(
        tail.f[seq_len(n.batch * b), , drop = FALSE],
        group = rep(seq_len(n.batch), each = b)
      ) / b

      resid    <- colMeans(batch.f)
      resid.se <- apply(batch.f, MARGIN = 2L, FUN = stats::sd) / sqrt(n.batch)
      resid.se[!is.finite(resid.se)] <- Inf # not reliable

      # Bonferroni-adjusted t-score - `n.batch` is small, so the normal
      # quantile would be too tight
      p.criterion <- 0.001
      resid.tol   <- stats::qt(
        1 - 0.5 * p.criterion / length(p0), df = n.batch - 1L
      )

      # Require the residual to be both statistically and practically
      # non-zero. With `mc.reps` large the sampling error is tiny, so a
      # residual well inside `mc.tol` can be "significant" without mattering.
      bad.resid <- abs(resid) > resid.tol * resid.se & abs(resid) > tol
    }

    bad.pars <- names(p0)[bad.resid]
    max.res  <- if (any(bad.resid)) max(abs(resid[bad.resid])) else NA_real_

    pls_warnif(
      any(bad.resid),
      "The MC-PLS residuals did not settle around zero for:",
      paste0(bad.pars, collapse = ", "),
      sprintf("(largest |mean residual| = %.4g).", max.res),
      "Delta-method standard errors might be unreliable for these parameters.",
      "Consider decreasing `mc.tol`, increasing `mc.max.iter`, or using",
      "bootstrap standard errors instead (`mc.delta.se = FALSE`)."
    )

    delta.jacobian.k <- delta.jacobian.k[[1L]]

    if (delta.fixed.seed) {
      seeds <- as.list(floor(stats::runif(delta.jacobian.k, min = 0, max = 9999999)))
    } else {
      seeds <- rep(list(rng.seed), delta.jacobian.k)
    }

    # Seed for the parallel workers. The replicate seeds are drawn above, so a
    # given `parallel`/`ncores` setting is reproducible. Serial and parallel
    # runs are not bit-identical though: the workers use L'Ecuyer-CMRG streams,
    # so the simulated data sets differ (by Monte-Carlo noise only) -- the same
    # applies to `bootstrap()`.
    jac.iseed <- floor(stats::runif(1L, min = 0, max = 9999999))

    JAC <- calcMcJacobians(
      .fg             = \(p, seed) .fg(p, thresholdStruct = thresholdStruct0, seed = seed),
      .probs          = \(seed) calcMcThresholdJacobians(
        .f              = .f,
        .simulate       = .simulate,
        p0              = p0,
        p1              = p1,
        thresholdStruct = thresholdStruct0,
        seed            = seed
      ),
      probs.names     = names(thresholdStruct0@proportions),
      seeds           = seeds,
      p0              = p0,
      p1              = p1,
      lower           = lower,
      upper           = upper,
      verbose         = verbose,
      parallel        = parallel,
      ncores          = ncores,
      iseed           = jac.iseed
    )

    fit1.combined@params$Jacobian0 <- JAC$J0
    fit1.combined@params$Jacobian1 <- JAC$J1
    fit1.combined@params$JacobianProbs0 <- JAC$Jp
    fit1.combined@params$JacobianProbs1 <- JAC$Gp
  }

  fit1.combined@params$mcpls.history <- plssemMatrix(mcfit$history.p)
  fit1.combined@params$mcpls.history.f <- plssemMatrix(mcfit$history.f)
  fit1.combined@status$par0 <- par0
  fit1.combined@status$fit0 <- fit0.base

  fit1.combined@status$iterations    <- mcfit$iter
  fit1.combined@info$mc.args$p.start <- as.vector(mcfit$root)

  fit0.base@combinedModel        <- fit1.combined
  fit0.base@status$iterations    <- mcfit$iter
  fit0.base@info$mc.args$p.start <- as.vector(mcfit$root)
  fit0.base
}


averageMcReplicates <- function(k, fun, point.estimate = "mean", catch = FALSE) {
  results <- vector("list", k)
  error   <- NULL

  for (i in seq_len(k)) {
    results[[i]] <- if (!catch) fun(i) else tryCatch(fun(i), error = \(e) {
      if (is.null(error)) error <<- conditionMessage(e)
      NULL
    })
  }

  failed <- vapply(results, FUN.VALUE = logical(1L), FUN = is.null)
  pls_stopif(all(failed),
    "The estimation failed for all the simulated data sets!",
    "Message (first failure):", error
  )

  THETA <- matrix(NA_real_, nrow = k, ncol = length(results[[which(!failed)[1L]]]))
  for (i in which(!failed))
    THETA[i, ] <- results[[i]]

  switch(point.estimate,
    median = colMedians(THETA, na.rm = TRUE),
    colMeans(THETA, na.rm = TRUE) # mean (default)
  )
}


ordinalizeDataFrame <- function(df, thresholdStruct) {
  nm <- colnames(df)
  ordered <- thresholdStruct@ordered
  probs   <- thresholdStruct@proportions
  indices <- thresholdStruct@indices

  quickdf(stats::setNames(
    lapply(nm, FUN = function(v) {
      if (v %in% ordered) ordinalizeVectorCpp(df[[v]], probs = probs[indices[[v]]])
      else df[[v]]
    }),
    nm = nm
  ))
}


getFreeParamsTable <- function(model) {
  model <- combinedModel(model)
  parTable <- getParTableEstimates(
    model, rm.tmp.ov = FALSE, clean.tmp.ind = FALSE
  )

  lhs <- parTable$lhs
  op  <- parTable$op
  rhs <- parTable$rhs

  inds.b <- rhs[op == "<~"]
  etas <- lhs[op == "~"]
  is.rescov <- (
    (lhs %in% etas | rhs %in% etas) & lhs != rhs & op == "~~"
  )

  cond1 <- !(lhs == rhs & op == "~~" & !grepl("~", rhs))
  cond2 <- !((isIntTermVariable(lhs) | isIntTermVariable(rhs)) & op == "~~")
  cond3 <- !op %in% c("~1", "|", ":=")
  cond4 <- !(lhs %in% inds.b & op == "~~") & !(rhs %in% inds.b & op == "~~")
  cond  <- cond1 & cond2 & cond3 & cond4

  out <- parTable[cond, , drop = FALSE]
  attr(out, "cond") <- cond

  out$is.free <- (op != "<~" & !is.rescov)[cond]
  out
}


isIntTermVariable <- function(x) {
  # Check if x is a intTerm variable name. However, it can be a parameter label
  # with an interaction term (e.g., "Y~X:Z")
  grepl(":", x) & !grepl("~", x)
}


updateModelFromFreeParTableMC <- function(parTable,
                                          model,
                                          mc.reps,
                                          thresholdStruct,
                                          ordered,
                                          seed = NULL,
                                          sim = NULL,
                                          params.only = FALSE,
                                          full = FALSE,
                                          retry = FALSE,
                                          n.retry = 5) {
  if (is.null(sim)) {
    sim <- simulateDataParTable(
      parTable     = parTable,
      N            = mc.reps,
      seed         = seed,
      check.hi.ord = model@info$is.high.ord,
      standardize  = TRUE,
      full         = full
    )

    if (retry && !sim$is.admissible) {
      sim0 <- sim

      for (i in seq_len(n.retry)) {
        if (is.null(seed)) seed.i <- NULL
        else               seed.i <- seed + i

        sim.i <- simulateDataParTable(
          parTable     = parTable,
          N            = mc.reps,
          seed         = seed.i,
          check.hi.ord = model@info$is.high.ord,
          standardize  = TRUE,
          full         = full
        )

        if (sim.i$is.admissible) {
          sim <- sim.i
          break
        }
      }

      # if we failed, revert to the original
      if (!sim$is.admissible)
        sim <- sim0
    }
  }

  SC     <- Rfast::cova(as.matrix(sim$all))
  ovs    <- colnames(model@matrices$S)
  lvsc   <- colnames(model@matrices$C)
  mode.a <- model@info$mode.a
  mode.b <- model@info$mode.b

  if (!params.only) {
    sim.ord <- ordinalizeDataFrame(
      df = sim$ov, thresholdStruct = thresholdStruct
    )

    model@matrices$S.ord.expected <- cov2cor(Rfast::cova(as.matrix(sim.ord)[, ovs]))
    model@matrices$S.ord.observed <- cov2cor(Rfast::cova(model@data[, ovs]))
    model@matrices$sim.ov.cont    <- sim$ov
    model@matrices$sim.ov.ord     <- sim.ord
  }

  model@matrices$S  <- SC[ovs,  ovs,  drop = FALSE]
  model@matrices$C  <- SC[lvsc, lvsc, drop = FALSE]
  model@matrices$SC <- SC[c(ovs, lvsc), c(ovs, lvsc), drop = FALSE]

  fit            <- model@fit
  fitMeasurement <- fit$fitMeasurement
  fitStructural  <- fit$fitStructural
  fitCov         <- fit$fitCov
  fitTheta       <- fit$fitTheta

  select       <- model@matrices$select
  selectLambda <- select$lambda
  selectGamma  <- select$gamma
  selectCov    <- select$cov
  selectTheta  <- select$theta

  vlhs <- parTable$lhs
  vop  <- parTable$op
  vrhs <- parTable$rhs

  getpar <- function(lhs, op, rhs) {
    cond <- vlhs == lhs & vop == op & vrhs == rhs
    if (op == "~~")
      cond <- cond | (vlhs == rhs & vop == op & vrhs == lhs)
    par <- parTable[cond, "est"]
    if (!length(par)) NA_real_ else par[[1L]]
  }

  cn <- colnames(fitMeasurement)
  rn <- rownames(fitMeasurement)

  mode.a <- intersect(mode.a, cn)
  mode.b <- intersect(mode.b, cn)

  if (!length(mode.a) && !length(mode.b)) {
    pls_msg_warn(
      "mode.a and mode.b are missing! This is likely a bug!"
    )

    # fallback
    mode.a <- cn
  }

  # Mode A
  for (lv in mode.a) {
    inds.lv <- rn[selectLambda[,lv, drop = TRUE]]

    for (ov in inds.lv) {

      par <- getpar(lhs = lv, op = "=~", rhs = ov)
      if (is.na(par)) par <- tryCatchNA(SC[ov, lv])

      fitMeasurement[ov, lv] <- par

      if (selectTheta[ov, ov])
        fitTheta[ov, ov] <- max(0, 1 - par^2)
    }
  }

  # Mode B
  for (lv in mode.b) {
    inds.lv <- rn[selectLambda[,lv, drop = TRUE]]

    for (ov in inds.lv) {

      par <- getpar(lhs = lv, op = "<~", rhs = ov)
      if (is.na(par)) par <- tryCatchNA(SC[ov, lv]) # this is likely a bad fallback
                                                    # but currently it's ok, since
      fitMeasurement[ov, lv] <- par                 # mode b isn't supported for more than
                                                    # one indicator per composite
      if (selectTheta[ov, ov])
        fitTheta[ov, ov] <- 1
    }
  }

  for (dep in colnames(fitStructural)) for (indep in rownames(fitStructural)) {
    if (!selectGamma[indep, dep]) next
    par <- getpar(lhs = dep, op = "~", rhs = indep)
    fitStructural[indep, dep] <- par
  }

  k          <- NCOL(fitCov)
  C          <- SC[colnames(fitStructural), colnames(fitStructural), drop = FALSE]
  projCov.mc <- t(fitStructural) %*% C %*% fitStructural
  diag(C)    <- diag(C) - diag(projCov.mc)

  for (i in seq_len(k)) for (j in seq_len(i)) {
    lhs <- colnames(fitCov)[[i]]
    rhs <- rownames(fitCov)[[j]]

    if (!selectCov[lhs, rhs]) next

    par <- getpar(lhs = lhs, op = "~~", rhs = rhs)
    if (is.na(par)) par <- tryCatchNA(C[lhs, rhs])
    fitCov[i, j] <- fitCov[j, i] <- par
  }

  model@fit$fitMeasurement    <- fitMeasurement
  model@fit$fitStructural     <- fitStructural
  model@fit$fitCov            <- fitCov
  model@fit$fitTheta          <- fitTheta
  model@status$is.admissible  <- (
    model@status$is.admissible && sim$is.admissible
  )

  model@status$mcpls.update.args <- list(
    parTable        = parTable,
    model           = model,
    mc.reps         = mc.reps,
    thresholdStruct = thresholdStruct,
    ordered         = ordered,
    seed            = seed,
    params.only     = params.only,
    full            = full,
    retry           = retry
  )

  model@thresholdStruct <- updateThresholds(
    thr = thresholdStruct, sim.cont = sim$ov
  )

  refreshModelParams(model, update.names = TRUE)
}


resampleMCPLS_Fit <- function(model, ...) {
  args         <- model@status$mcpls.update.args
  new.args     <- list(...)
  args[names(new.args)] <- new.args
  do.call(updateModelFromFreeParTableMC, args)
}


thresholdJacobian <- function(thresholdStruct, sim.cont = NULL, eps = 1e-3,
                              zero.tol = .Machine$double.eps^0.5) {

  # if we have no simulated data, we assume the normal distribution
  if (is.null(sim.cont)) {
    thr <- thresholdStruct@thresholds
    probs <- thresholdStruct@proportions

    out <- diag(1 / stats::dnorm(thr), nrow = length(thr))
    dimnames(out) <- list(names(thr), names(probs))

    return(out)
  }

  # Get empirical finite difference jacobian. Each threshold only depends on
  # its own (cumulative) proportion, so all of the proportions can be perturbed
  # at once. The steps are bounded by the neighbouring proportions of the same
  # variable (and 0 and 1), such that the perturbed proportions are still
  # increasing (a requirement for computing the quantiles).
  probs <- thresholdStruct@proportions
  step  <- rep(eps, length(probs))

  for (ord in thresholdStruct@ordered) {
    idx  <- thresholdStruct@indices[[ord]]
    gaps <- diff(c(0, probs[idx], 1))

    step[idx] <- pmin(eps, 0.45 * gaps[-length(gaps)], 0.45 * gaps[-1L])
  }

  # no stable finite-difference step (e.g., an empty category)
  unstable <- step <= zero.tol
  step[unstable] <- 0

  p0 <- probs - step
  p1 <- probs + step

  # update thresholdStruct
  T0 <- T1 <- thresholdStruct
  T0@proportions <- p0
  T1@proportions <- p1

  t0 <- updateThresholds(T0, sim.cont = sim.cont)@thresholds
  t1 <- updateThresholds(T1, sim.cont = sim.cont)@thresholds

  d <- (t1 - t0) / (p1 - p0)
  d[unstable] <- 0

  out <- diag(d, nrow = length(t1))
  dimnames(out) <- list(names(t0), names(p0))

  out
}


# Estimates the Jacobians used for the (implicit) delta-method standard errors
#   J0 = df/dp, where f(p) are the (naive) statistics of the root equation
#   J1 = dg/dp, where g(p) are the values of all the reported parameters
#
# Optionally, `.probs(seed)` returns the Jacobians w.r.t. the threshold
# proportions
calcMcJacobians <- function(.fg,
                            seeds,
                            p0,
                            p1,
                            parallel,
                            ncores,
                            verbose,
                            iseed,
                            lower       = -Inf,
                            upper       = Inf,
                            eps         = 5e-3,
                            .probs      = NULL,
                            probs.names = NULL) {

  k <- length(seeds)

  J0 <- matrix(
    0,
    nrow = length(p0), ncol = length(p0),
    dimnames = list(names(p0), names(p0))
  )

  J1 <- matrix(
    0,
    nrow = length(p1), ncol = length(p0),
    dimnames = list(names(p1), names(p0))
  )

  Jp <- matrix(
    0,
    nrow = length(p0), ncol = length(probs.names),
    dimnames = list(names(p0), probs.names)
  )

  Gp <- matrix(
    0,
    nrow = length(p1), ncol = length(probs.names),
    dimnames = list(names(p1), probs.names)
  )

  tasks.k <- function(k) {
    par.tasks <- lapply(
      seq_along(p0),
      FUN = \(i) list(k = k, type = "par", index = i)
    )

    if (is.null(.probs) || !length(probs.names)) par.tasks
    else c(par.tasks, list(list(k = k, type = "probs")))
  }

  tasks <- unlist(
    lapply(X = seq_len(k), FUN = tasks.k),
    recursive = FALSE
  )

  do.task <- function(task) {
    seed <- seeds[[task$k]]

    if (task$type == "probs")
      return(.probs(seed))

    points <- boundedParameterFiniteDiffPoints(
      x = p0, i = task$index, eps = eps, lower = lower, upper = upper
    )

    fg.p <- .fg(points$plus,  seed = seed)
    fg.m <- .fg(points$minus, seed = seed)

    list(
      J0 = (fg.p$f - fg.m$f) / points$denominator,
      J1 = (fg.p$g - fg.m$g) / points$denominator
    )
  }

  results <- plapply(
    X        = tasks,
    FUN      = do.task,
    parallel = parallel,
    ncores   = ncores,
    verbose  = verbose,
    iseed    = iseed,
    label    = "Jacobian"
  )

  for (j in seq_along(tasks)) {
    task <- tasks[[j]]
    res  <- results[[j]]

    if (task$type == "probs") {
      Jp <- Jp + res$Jp / k
      Gp <- Gp + res$Gp / k

    } else {
      J0[, task$index] <- J0[, task$index] + res$J0 / k
      J1[, task$index] <- J1[, task$index] + res$J1 / k

    }
  }

  list(J0 = J0, J1 = J1, Jp = Jp, Gp = Gp)
}


calcMcThresholdJacobians <- function(.f, .simulate, p0, p1, thresholdStruct,
                                     eps = 5e-3, seed = NULL,
                                     sim = .simulate(p0, standardize = TRUE, seed = seed),
                                     sim.cont = sim$ov) {
  probs0 <- thresholdStruct@proportions

  Jp <- matrix(
    0,
    nrow = length(p0), ncol = length(probs0),
    dimnames = list(names(p0), names(probs0))
  )

  Gp <- matrix(
    0,
    nrow = length(p1), ncol = length(probs0),
    dimnames = list(names(p1), names(probs0))
  )

  for (i in seq_along(probs0)) {
    points <- boundedProbabilityFiniteDiffPoints(probs0, i = i, eps = eps)
    if (is.null(points)) {
      pls_msg_warn(
        "Skipping a threshold-probability derivative because no stable ",
        "finite-difference step is available for: ", names(probs0)[[i]]
      )
      next
    }

    # check bounds
    T0 <- T1 <- thresholdStruct
    T1@proportions <- points$plus
    T0@proportions <- points$minus

    Jp[,i] <- (
      .f(p0, thresholdStruct = T1, sim = sim) -
      .f(p0, thresholdStruct = T0, sim = sim)
    ) / points$denominator
  }

  # Jacobian probs->thresholds
  # use larger eps for better numerical stability
  T <- thresholdJacobian(thresholdStruct, sim.cont = sim.cont, eps = 2 * eps)

  thr.rows  <- intersect(rownames(Gp), rownames(T))
  prob.cols <- intersect(colnames(Gp), colnames(T))

  Gp[thr.rows, prob.cols] <- T[thr.rows, prob.cols, drop = FALSE]

  list(Jp = Jp, Gp = Gp)
}


boundedParameterFiniteDiffPoints <- function(x, i, eps, lower = -Inf, upper = Inf) {
  lower <- rep_len(lower, length(x))
  upper <- rep_len(upper, length(x))

  plus <- minus <- x
  plus[i]  <- min(x[i] + eps, upper[i])
  minus[i] <- max(x[i] - eps, lower[i])
  denominator <- plus[i] - minus[i]

  pls_stopif(
    !is.finite(denominator) || denominator <= .Machine$double.eps^0.5,
    "Unable to calculate a finite-difference derivative at a parameter bound."
  )

  list(plus = plus, minus = minus, denominator = denominator)
}


boundedProbabilityFiniteDiffPoints <- function(probs, i, eps = 1e-3,
                                               tol = .Machine$double.eps^0.5) {
  par <- names(probs)[[i]]
  var <- stringr::str_split_i(par, pattern = "\\|", i = 1L)

  idx <- which(grepl(paste0("^", var, "\\|P[0-9]+$"), names(probs)))
  probs.x <- probs[idx] # keep only probabilities for the relevant variable

  ix <- which(idx == i)
  bound <- \(j) if (j < 1) 0 else if (j > length(probs.x)) 1 else probs.x[j]

  lower <- bound(ix - 1)
  upper <- bound(ix + 1)
  step <- min(eps, 0.45 * (probs.x[ix] - lower), 0.45 * (upper - probs.x[ix]))

  if (!is.finite(step) || step <= tol)
    return(NULL)

  plus <- minus <- probs
  plus[i] <- probs[i] + step
  minus[i] <- probs[i] - step

  list(plus = plus, minus = minus, denominator = 2 * step)
}


getMcLowerBounds <- function(par, tol = 1e-3) {
  parf  <- par[par$is.free, , drop = FALSE]
  lower <- rep(-Inf, NROW(parf))

  lower[parf$op == "=~"]                        <- tol - 1
  lower[parf$op == "~~" & parf$lhs == parf$rhs] <- tol

  lower
}


getMcUpperBounds <- function(par, tol = 1e-3) {
  parf  <- par[par$is.free, , drop = FALSE]
  upper <- rep(Inf, NROW(parf))

  upper[parf$op == "=~"] <- 1 - tol

  upper
}
