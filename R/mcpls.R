mcpls <- function(
  fit0,
  p.start         = fit0@info$mc.args$p.start,
  min.iter        = fit0@info$mc.args$min.iter,
  max.iter        = fit0@info$mc.args$max.iter,
  mc.reps         = fit0@info$mc.args$mc.reps,
  rng.seed        = fit0@info$mc.args$rng.seed,
  tol             = fit0@info$mc.args$tol,
  fixed.seed      = fit0@info$mc.args$fixed.seed,
  verbose         = fit0@info$verbose,
  polyak.juditsky = fit0@info$mc.args$polyak.juditsky,
  fn.args         = fit0@info$mc.args$fn.args,
  ...
) {
  fit0.base <- fit0
  fit0.combined <- combinedModel(fit0.base)

  data      <- fit0.base@data
  vars      <- colnames(data)
  ordered   <- fit0@info$ordered
  is.hi.ord <- isTRUE(fit0.combined@info$is.high.ord)

  PROBS <- getPROBS(data = data, ordered = ordered)

  par0 <- getFreeParamsTable(fit0.combined)
  par1 <- par0[c("lhs", "op", "rhs", "est", "is.free")]

  if (fixed.seed && is.null(rng.seed)) {
    rng.seed <- floor(stats::runif(1L, min = 0, max = 9999999))
    if (verbose) message(sprintf("Using fixed seed %i...", rng.seed))
  }

  if (isMLM(fit0)) {
    clusterSizes <- as.numeric(table(attr(data, "cluster")))
    clusterName  <- colnames(attr(data, "cluster"))

  } else {
    clusterSizes <- NULL
    clusterName  <- NULL

  }

  .f <- function(p) {
    par1[par1$is.free, "est"] <- p

    sim <- simulateDataParTable(
      parTable     = par1,
      N            = mc.reps,
      seed         = rng.seed,
      check.hi.ord = is.hi.ord,
      clusterSizes = clusterSizes,
      clusterName  = clusterName
    )

    sim.ov  <- ordinalizeDataFrame(
      df      = sim$ov,
      PROBS   = PROBS,
      ordered = ordered
    )

    fit.sim <- fit0.base
    X       <- Rfast::standardise(as.matrix(sim.ov[vars]))
    S       <- Rfast::cova(X)

    if (!is.null(sim$cluster))
      attr(X, "cluster") <- sim$cluster

    # Update observed-data (lowest-order) model input
    modelData(fit.sim)  <- X
    indCorrMatrix(fit.sim) <- S

    fit2 <- estimatePLS_Inner(fit.sim)
    par2 <- getFreeParamsTable(combinedModel(fit2))

    eps <- par2$est - par0$est + sim$penalty
    eps[par0$is.free]
  }

  ok.start <- !is.null(p.start) && length(p.start) == sum(par1$is.free)
  p        <- if (ok.start) p.start else par1[par1$is.free, "est"]

  mcfit <- robbinsMonro1951(
    p               = p,
    f               = .f,
    tol             = tol,
    min.iter        = min.iter,
    max.iter        = max.iter,
    verbose         = verbose,
    polyak.juditsky = polyak.juditsky,
    fn.args         = fn.args,
    ...
  )

  iter <- mcfit$iter
  if (iter >= max.iter && !fixed.seed) {
    rng.seed <- floor(stats::runif(1L, min = 0, max = 9999999))
    warning2("Maximum number of (initial) iterations reached!\n",
             sprintf("Attempting to use fixed seed %i...", rng.seed))

    mcfit <- robbinsMonro1951(
      p               = as.vector(mcfit$root),
      f               = .f,
      tol             = tol,
      min.iter        = min.iter,
      max.iter        = max.iter,
      verbose         = verbose,
      polyak.juditsky = polyak.juditsky,
      fn.args         = fn.args,
      ...
    )

    iter <- iter + mcfit$iter
  }

  # Check status of (last) mcfit
  if (mcfit$iter >= max.iter) {
    warning2(
      "Maximum number of iterations reached!\n",
      "Parameter estimates might be unreliable!"
    )

    modelStatus(fit0.combined)$is.admissible <- FALSE
  }

  par1[par1$is.free, "est"] <- as.vector(mcfit$root)

  fit1.combined <- updateModelFromFreeParTableMC(
    parTable     = par1,
    model        = fit0.combined,
    mc.reps      = mc.reps,
    PROBS        = PROBS,
    ordered      = ordered,
    seed         = rng.seed,
    clusterSizes = clusterSizes,
    clusterName  = clusterName
  )

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


ordinalizeDataFrame <- function(df, PROBS, ordered = NULL) {
  nm <- colnames(df)
  quickdf(stats::setNames(
    lapply(nm, FUN = \(v) if (v %in% ordered) ordinalize(df[[v]], probs = PROBS[[v]])
                          else df[[v]]),
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

  cond1 <- !(lhs == rhs & op == "~~" & !grepl("~", rhs))
  cond2 <- !((isIntTermVariable(lhs) | isIntTermVariable(rhs)) & op == "~~")
  cond3 <- op != "~1"

  out <- parTable[cond1 & cond2 & cond3, , drop = FALSE]
  attr(out, "cond") <- cond1 & cond2 & cond3

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
                                          PROBS,
                                          ordered,
                                          seed = NULL,
                                          clusterSizes = NULL,
                                          clusterName = NULL) {
  sim <- simulateDataParTable(
    parTable     = parTable,
    N            = mc.reps,
    seed         = seed,
    check.hi.ord = model@info$is.high.ord,
    clusterSizes = clusterSizes,
    clusterName  = clusterName
  )

  sim.ord <- ordinalizeDataFrame(
    df      = sim$ov,
    PROBS   = PROBS,
    ordered = ordered
  )

  lvs    <- getLVs(parTable)
  indsLVs <- getIndsLVs(parTable, lVs = lvs)

  SC   <- cov2cor(Rfast::cova(as.matrix(sim$all)))
  ovs  <- colnames(model@matrices$S)
  lvsc <- colnames(model@matrices$C)

  model@matrices$S.ord.expected <- cov2cor(Rfast::cova(as.matrix(sim.ord)[, ovs]))
  model@matrices$S.ord.observed <- cov2cor(Rfast::cova(model@data[, ovs]))
  model@matrices$sim.ov.cont    <- sim$ov
  model@matrices$sim.ov.ord     <- sim.ord
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

  for (lv in colnames(fitMeasurement)) for (ov in rownames(fitMeasurement)) {
    if (!selectLambda[ov, lv]) next
    par <- getpar(lhs = lv, op = "=~", rhs = ov)
    if (is.na(par)) par <- tryCatchNA(SC[ov, lv])
    fitMeasurement[ov, lv] <- par

    if (selectTheta[ov, ov])
      fitTheta[ov, ov] <- max(0, 1 - par)
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
    parTable     = parTable,
    model        = model,
    mc.reps      = mc.reps,
    PROBS        = PROBS,
    ordered      = ordered,
    seed         = seed,
    clusterSizes = clusterSizes,
    clusterName  = clusterName
  )

  refreshModelParams(model, update.names = TRUE)
}


resampleMCPLS_Fit <- function(model, ...) {
  args         <- model@status$mcpls.update.args
  new.args     <- list(...)
  args[names(new.args)] <- new.args
  do.call(updateModelFromFreeParTableMC, args)
}


getPROBS <- function(data, ordered) {
  PROBS <- stats::setNames(vector("list", length(ordered)), nm = ordered)
  for (ord in ordered) {
    freq       <- table(data[, ord])
    pct        <- cumsum(freq) / sum(freq)
    PROBS[[ord]] <- pct[-length(pct)]
  }
  PROBS
}
