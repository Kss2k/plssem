mcpls <- function(
  fit0,
  p.start          = fit0@info$mc.args$p.start,
  min.iter         = fit0@info$mc.args$min.iter,
  max.iter         = fit0@info$mc.args$max.iter,
  mc.reps          = fit0@info$mc.args$mc.reps,
  rng.seed         = fit0@info$mc.args$rng.seed,
  tol              = fit0@info$mc.args$tol,
  fixed.seed       = fit0@info$mc.args$fixed.seed,
  verbose          = fit0@info$verbose,
  polyak.juditsky  = fit0@info$mc.args$polyak.juditsky,
  fn.args          = fit0@info$mc.args$fn.args,
  pj.extrapolate   = fit0@info$mc.args$pj.extrapolate,
  delta.jacobian   = fit0@info$mc.args$delta.se && fit0@info$boot$bootstrap,
  delta.fixed.seed = TRUE,
  delta.jacobian.k = fit0@info$mc.args$delta.jacobian.k,
  ...
) {
  fit0.base <- fit0
  fit0.combined <- combinedModel(fit0.base)

  data      <- fit0.base@data
  vars      <- colnames(data)
  ordered   <- fit0@info$ordered
  is.hi.ord <- isTRUE(fit0.combined@info$is.high.ord)
  thresholdStruct0 <- fit0.combined@thresholdStruct

  par0 <- getFreeParamsTable(fit0.combined)
  par1 <- par0[c("lhs", "op", "rhs", "est", "is.free")]

  if (fixed.seed && is.null(rng.seed)) {
    rng.seed <- floor(stats::runif(1L, min = 0, max = 9999999))
    if (verbose) pls_msg_note(sprintf("Using fixed seed %i...", rng.seed))
  }

  if (isMLM(fit0)) {
    clusterSizes <- as.numeric(table(attr(data, "cluster")))
    clusterName  <- colnames(attr(data, "cluster"))

  } else {
    clusterSizes <- NULL
    clusterName  <- NULL

  }

  .parTable <- function(p) {
    parx <- par1
    parx[parx$is.free, "est"] <- p
    parx
  }

  .simulate <- function(p, standardize = FALSE) {
    simulateDataParTable(
      parTable     = .parTable(p),
      N            = mc.reps,
      seed         = rng.seed,
      check.hi.ord = is.hi.ord,
      clusterSizes = clusterSizes,
      clusterName  = clusterName,
      standardize  = standardize
    )
  }

  .f <- function(p, thresholdStruct = thresholdStruct0, sim = NULL) {

    if (is.null(sim)) {
      par1[par1$is.free, "est"] <- p
      sim <- simulateDataParTable(
        parTable     = par1,
        N            = mc.reps,
        seed         = rng.seed,
        check.hi.ord = is.hi.ord,
        clusterSizes = clusterSizes,
        clusterName  = clusterName
      )
    }

    sim.ov  <- ordinalizeDataFrame(
      df = sim$ov, thresholdStruct = thresholdStruct
    )

    fit.sim <- fit0.base
    X       <- Rfast::standardise(as.matrix(sim.ov[vars]))
    S       <- Rfast::cova(X)

    if (!is.null(sim$cluster))
      attr(X, "cluster") <- sim$cluster

    # Update observed-data (lowest-order) model input
    modelData(fit.sim)  <- X
    indCorrMatrix(fit.sim) <- S

    # Thresholds are not part of the root equation. Avoid recomputing them on
    # every Robbins-Monro iteration.
    fit2 <- estimatePLS_Inner(fit.sim)
    par2 <- getFreeParamsTable(combinedModel(fit2))

    eps <- par2$est - par0$est + sim$penalty
    eps[par0$is.free]
  }

  .g <- function(p, thresholdStruct = thresholdStruct0, sim = NULL) {
    fit <- updateModelFromFreeParTableMC(
      parTable         = .parTable(p),
      model            = fit0.combined,
      mc.reps          = mc.reps,
      thresholdStruct  = thresholdStruct,
      ordered          = ordered,
      seed             = rng.seed,
      clusterSizes     = clusterSizes,
      clusterName      = clusterName,
      sim              = sim,
      params.only      = TRUE
    )

    fit@params$values
  }

  .fg <- function(p, thresholdStruct = thresholdStruct0, sim = NULL) {
    if (is.null(sim))
      sim <- .simulate(p, standardize = TRUE)

    list(f = .f(p, thresholdStruct = thresholdStruct, sim = sim),
         g = .g(p, thresholdStruct = thresholdStruct, sim = sim))
  }

  ok.start <- !is.null(p.start) && length(p.start) == sum(par1$is.free)
  p        <- if (ok.start) p.start else par1[par1$is.free, "est"]

  lower <- getMcLowerBounds(par1)
  upper <- getMcUpperBounds(par1)

  if (polyak.juditsky && !pj.extrapolate) {
    # If we're not using a Nonlinear Regression to solve for the convergence
    # point, we will get a biased root with Polyak-Juditsky averaging, if
    # we have no warmup
    if (verbose) pls_msg_note("Warming up...")

    mcfit <- robbinsMonro1951(
      p               = p,
      f               = .f,
      tol             = 10 * tol,
      min.iter        = 5L,
      max.iter        = 20L,
      verbose         = verbose,
      polyak.juditsky = FALSE,
      fn.args         = fn.args,
      lower           = lower,
      upper           = upper,
      ...
    )

    p <- as.vector(mcfit$root)
  }

  mcfit <- robbinsMonro1951(
    p               = p,
    f               = .f,
    tol             = tol,
    min.iter        = min.iter,
    max.iter        = max.iter,
    verbose         = verbose,
    polyak.juditsky = polyak.juditsky,
    fn.args         = fn.args,
    pj.extrapolate  = pj.extrapolate,
    lower           = lower,
    upper           = upper,
    ...
  )

  iter <- mcfit$iter
  if (iter >= max.iter && !polyak.juditsky) {
    pls_msg_warn(
      "Maximum number of (initial) iterations reached!\n",
      sprintf("Attempting to use Polyak Juditsky averaging...")
    )

    mcfit <- robbinsMonro1951(
      p               = as.vector(mcfit$root),
      f               = .f,
      tol             = tol,
      min.iter        = min.iter,
      max.iter        = max.iter,
      verbose         = verbose,
      polyak.juditsky = TRUE,
      fn.args         = fn.args,
      pj.extrapolate  = pj.extrapolate,
      lower           = lower,
      upper           = upper,
      ...
    )

    iter <- iter + mcfit$iter
  }

  # Check status of (last) mcfit
  if (mcfit$iter >= max.iter) {
    pls_msg_warn(
      "Maximum number of iterations reached!\n",
      "Parameter estimates might be unreliable!"
    )

    modelStatus(fit0.combined)$is.admissible <- FALSE
  }

  par1[par1$is.free, "est"] <- as.vector(mcfit$root)

  fit1.combined <- updateModelFromFreeParTableMC(
    parTable        = par1,
    model           = fit0.combined,
    mc.reps         = mc.reps,
    thresholdStruct = thresholdStruct0,
    ordered         = ordered,
    seed            = rng.seed,
    clusterSizes    = clusterSizes,
    clusterName     = clusterName
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


    probs0 <- thresholdStruct0@proportions
    nm <- paste0(par1$lhs, par1$op, par1$rhs)
    p0 <- stats::setNames(mcfit$root, nm[par1$is.free])
    p1 <- fit1.combined@params$values

    if (verbose) {
      pb <- utils::txtProgressBar(
        min     = 0,
        max     = delta.jacobian.k * (length(p) + length(probs0)),
        initial = 0,
        style   = 3,
        file    = stderr()
      )

      on.exit(close(pb), add = TRUE)

    } else {
      pb <- NULL

    }

    delta.jacobian.k <- delta.jacobian.k[[1L]]
    J0 <- J1 <- Jp <- Gp <- 0

    for (i in seq_len(delta.jacobian.k)) {

      if (delta.fixed.seed)
        rng.seed <- floor(stats::runif(1L, min = 0, max = 9999999))

      JAC  <- calcMcJacobians(
        .fg             = .fg,
        .f              = .f,
        .simulate       = .simulate,
        p0              = p0,
        p1              = p1,
        thresholdStruct = thresholdStruct0,
        lower           = lower,
        upper           = upper,
        progressBar     = pb,
        k               = i
      )

      J0.i <- JAC$J0
      J1.i <- JAC$J1
      Jp.i <- JAC$Jp
      Gp.i <- JAC$Gp

      J1 <- J1 + J1.i / delta.jacobian.k
      J0 <- J0 + J0.i / delta.jacobian.k
      Jp <- Jp + Jp.i / delta.jacobian.k
      Gp <- Gp + Gp.i / delta.jacobian.k
    }

    fit1.combined@params$Jacobian0 <- J0
    fit1.combined@params$Jacobian1 <- J1
    fit1.combined@params$JacobianProbs0 <- Jp
    fit1.combined@params$JacobianProbs1 <- Gp
  }


  fit1.combined@status$iterations    <- mcfit$iter
  fit1.combined@info$mc.args$p.start <- as.vector(mcfit$root)

  fit0.base@combinedModel        <- fit1.combined
  fit0.base@status$iterations    <- mcfit$iter
  fit0.base@info$mc.args$p.start <- as.vector(mcfit$root)
  fit0.base
}


ordinalize <- function(x, probs) {
  probs  <- sort(probs[probs < 1])
  breaks <- collapse::fquantile(x, probs = probs)
  findInterval(x, vec = breaks)
}


ordinalizeDataFrame <- function(df, thresholdStruct) {
  nm <- colnames(df)
  ordered <- thresholdStruct@ordered
  probs   <- thresholdStruct@proportions
  indices <- thresholdStruct@indices

  quickdf(stats::setNames(
    lapply(nm, FUN = function(v) {
      if (v %in% ordered) ordinalize(df[[v]], probs = probs[indices[[v]]])
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

  inds.b <- rhs[parTable$op == "<~"]

  cond1 <- !(lhs == rhs & op == "~~" & !grepl("~", rhs))
  cond2 <- !((isIntTermVariable(lhs) | isIntTermVariable(rhs)) & op == "~~")
  cond3 <- op != "~1"
  cond4 <- !(lhs %in% inds.b & op == "~~") & !(rhs %in% inds.b & op == "~~")
  cond5 <- op != "|"
  cond  <- cond1 & cond2 & cond3 & cond4 & cond5

  out <- parTable[cond, , drop = FALSE]
  attr(out, "cond") <- cond

  out$is.free <- out$op != "<~"
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
                                          clusterSizes = NULL,
                                          clusterName = NULL,
                                          sim = NULL,
                                          params.only = FALSE) {
  if (is.null(sim)) {
    sim <- simulateDataParTable(
      parTable     = parTable,
      N            = mc.reps,
      seed         = seed,
      check.hi.ord = model@info$is.high.ord,
      clusterSizes = clusterSizes,
      clusterName  = clusterName,
      standardize  = TRUE
    )
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

  model@matrices$S              <- SC[ovs,  ovs,  drop = FALSE]
  model@matrices$C              <- SC[lvsc, lvsc, drop = FALSE]
  model@matrices$SC             <- SC[c(ovs, lvsc), c(ovs, lvsc), drop = FALSE]

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
  C          <- SC[colnames(fitStructural), colnames(fitStructural)]
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

  if (isMLM(model)) {
    params <- parTableToParams(parTable)
    modelFitLmer(model)$values <- params$values
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
    clusterSizes    = clusterSizes,
    clusterName     = clusterName,
    params.only     = params.only
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

  # Get empirical finite difference jacobian
  probs <- thresholdStruct@proportions

  p0 <- probs - eps
  p1 <- probs + eps

  # check bounds
  p0[p0 < 0] <- probs[p0 < 0]
  p1[p1 > 1] <- probs[p1 > 1]

  # update thresholdStruct
  T0 <- T1 <- thresholdStruct
  T0@proportions <- p0
  T1@proportions <- p1

  t0 <- updateThresholds(T0, sim.cont = sim.cont)@thresholds
  t1 <- updateThresholds(T1, sim.cont = sim.cont)@thresholds

  # get step sizes (potentially affected by clamping above) from p1 and p0
  out <- diag((t1 - t0) / (p1 - p0), nrow = length(t1))
  dimnames(out) <- list(names(t0), names(p0))

  out
}


calcMcJacobians <- function(.fg, .f, .simulate, p0, p1,
                            thresholdStruct, lower = -Inf, upper = Inf,
                            eps = 5e-3, progressBar = NULL, k = 1) {
  probs0 <- thresholdStruct@proportions

  J0 <- matrix(
    NA_real_,
    nrow = length(p0), ncol = length(p0),
    dimnames = list(names(p0), names(p0))
  )

  J1 <- matrix(
    NA_real_,
    nrow = length(p1), ncol = length(p0),
    dimnames = list(names(p1), names(p0))
  )

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

  for (i in seq_along(p0)) {
    points <- boundedParameterFiniteDiffPoints(
      x = p0, i = i, eps = eps, lower = lower, upper = upper
    )

    fg.p <- .fg(points$plus, thresholdStruct = thresholdStruct)
    fg.m <- .fg(points$minus, thresholdStruct = thresholdStruct)
    J0[,i] <- (fg.p$f - fg.m$f) / points$denominator
    J1[,i] <- (fg.p$g - fg.m$g) / points$denominator

    if (!is.null(progressBar)) {
      utils::setTxtProgressBar(
        progressBar, (k - 1) * (length(p0) + length(probs0)) + i
      )
    }
  }

  offset <- length(p0)
  sim0 <- if (length(probs0)) .simulate(p0, standardize = TRUE) else NULL

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
      .f(p0, thresholdStruct = T1, sim = sim0) -
      .f(p0, thresholdStruct = T0, sim = sim0)
    ) / points$denominator

    if (!is.null(progressBar)) {
      utils::setTxtProgressBar(
        progressBar, (k - 1) * (length(p0) + length(probs0)) + offset + i
      )
    }
  }

  # Jacobian probs->thresholds
  if (length(probs0)) {
    # use larger eps for better numerical stability
    T <- thresholdJacobian(thresholdStruct, sim.cont = sim0$ov, eps = 2 * eps)

    thr.rows <- intersect(rownames(Gp), rownames(T))
    prob.cols <- intersect(colnames(Gp), colnames(T))

    Gp[thr.rows, prob.cols] <- T[thr.rows, prob.cols, drop = FALSE]
  }

  list(J0 = J0, J1 = J1, Jp = Jp, Gp = Gp)
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
