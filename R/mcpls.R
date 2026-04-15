mcpls <- function(
  fit0,
  p.start         = fit0$info$mc.args$p.start,
  min.iter        = fit0$info$mc.args$min.iter,
  max.iter        = fit0$info$mc.args$max.iter,
  mc.reps         = fit0$info$mc.args$mc.reps,
  rng.seed        = fit0$info$mc.args$rng.seed,
  tol             = fit0$info$mc.args$tol,
  fixed.seed      = fit0$info$mc.args$fixed.seed,
  verbose         = fit0$info$verbose,
  polyak.juditsky = fit0$info$mc.args$polyak.juditsky,
  fn.args         = fit0$info$mc.args$fn.args,
  ...
) {

  fit0.base <- fit0
  data <- fit0.base$data
  vars <- colnames(data)
  ordered <- fit0$info$ordered

  PROBS <- stats::setNames(vector("list", length(ordered)), nm = ordered)
  for (ord in ordered) {
    freq <- table(data[,ord])
    pct  <- cumsum(freq) / sum(freq)
    PROBS[[ord]] <- pct[-length(pct)]
  }

  par0 <- getFreeParamsTable(fit0)
  par1 <- par0[c("lhs", "op", "rhs", "est", "is.free")]

  if (fixed.seed && is.null(rng.seed)) {
    rng.seed <- floor(stats::runif(1L, min = 0, max = 9999999))
    if (verbose) printf("Using fixed seed %i...\n", rng.seed)
  }

  .f <- function(p) {
    par1[par1$is.free, "est"] <- p

    sim <- simulateDataParTable(par1, N = mc.reps, seed = rng.seed)
    sim.ov <- ordinalizeDataFrame(sim$ov, PROBS = PROBS, ordered = ordered)

    fit.sim <- fit0.base
    fit.sim$data <- Rfast::standardise(as.matrix(sim.ov[vars]))
    fit.sim$matrices$S <- Rfast::cova(fit.sim$data)

    fit2 <- estimatePLS_Inner(fit.sim)
    par2 <- getFreeParamsTable(fit2)

    eps <- par2$est - par0$est + sim$penalty
    eps[par0$is.free]
  }

  p <- par1[par1$is.free, "est"]


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
    warning2("Maximum number of iterations reached!\n",
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

  par1[par1$is.free, "est"] <- as.vector(mcfit$root)
  fit1 <- updateModelFromFreeParTableMC(
    parTable = par1,
    model    = fit0.base,
    mc.reps  = mc.reps,
    PROBS    = PROBS,
    ordered  = ordered, 
    seed     = rng.seed
  )

  fit1$status$iterations <- mcfit$iter
  fit1$info$mc.args$p.start <- as.vector(mcfit$root)

  fit1
}


ordinalize <- function(x, probs) {
  probs <- sort(probs[probs < 1]) # skip 100%, we analytically know it to be Inf
  breaks <- collapse::fquantile(x, probs = probs)

  findInterval(x, vec = breaks)
}


ordinalizeDataFrame <- function(df, PROBS, ordered = NULL) {
  nm <- colnames(df)
  quickdf(stats::setNames(
    lapply(nm, FUN = \(v) if (v %in% ordered) ordinalize(df[[v]], probs = PROBS[[v]])
                          else df[[v]]),
    nm = nm))
}


getFreeParamsTable <- function(model) {
  parTable <- getParTableEstimates(model, rm.tmp = FALSE)

  lhs <- parTable$lhs
  op  <- parTable$op
  rhs <- parTable$rhs

  cond1 <- !(lhs == rhs & op == "~~")
  cond2 <- !((grepl(":", lhs) | grepl(":", rhs)) & op == "~~")

  out <- parTable[cond1 & cond2, , drop = FALSE]
  attr(out, "cond") <- cond1 & cond2

  out$is.free <- out$op != "<~" # Can be made more comprehensive later

  out
}


updateModelFromFreeParTableMC <- function(parTable, model, mc.reps,
                                          PROBS, ordered, seed = NULL) {
  sim <- simulateDataParTable(parTable, N = mc.reps, seed = seed)
  sim.ord <- ordinalizeDataFrame(sim$ov, PROBS = PROBS, ordered = ordered)

  lvs <- getLVs(parTable)
  indsLVs <- getIndsLVs(parTable, lVs = lvs)

  SC <- Rfast::cova(as.matrix(sim$all))

  ovs <- colnames(model$matrices$S)
  lvs <- colnames(model$matrices$C)

  model$matrices$S.ord.expected <- Rfast::cova(as.matrix(sim.ord))
  model$matrices$S.ord.observed <- Rfast::cova(model$data)

  model$matrices$S  <- SC[ovs, ovs, drop = FALSE]
  model$matrices$C  <- SC[lvs, lvs, drop = FALSE]
  model$matrices$SC <- SC[c(ovs, lvs), c(ovs, lvs), drop = FALSE]

  fit <- model$fit
  fitMeasurement <- fit$fitMeasurement
  fitStructural  <- fit$fitStructural
  fitCov         <- fit$fitCov
  fitTheta       <- fit$fitTheta

  select <- model$matrices$select
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

  k <- NCOL(fitCov)
  C <- SC[colnames(fitStructural), colnames(fitStructural)]
  projCov.mc <- t(fitStructural) %*% C %*% fitStructural
  diag(C) <- diag(C) - diag(projCov.mc)

  for (i in seq_len(k)) for (j in seq_len(i)) {
    lhs <- colnames(fitCov)[[i]]
    rhs <- rownames(fitCov)[[j]]

    if (!selectCov[lhs, rhs]) next

    par <- getpar(lhs = lhs, op = "~~", rhs = rhs)
    if (is.na(par)) par <- tryCatchNA(C[lhs, rhs])

    fitCov[i, j] <- fitCov[j, i] <- par
  }

  model$fit$fitMeasurement <- fitMeasurement
  model$fit$fitStructural  <- fitStructural
  model$fit$fitCov         <- fitCov
  model$fit$fitTheta       <- fitTheta
  model$status$is.admissible <- sim$is.admissible

  model$status$mcpls.update.args <- list(
    parTable = parTable,
    model    = model,
    mc.reps  = mc.reps,
    PROBS    = PROBS,
    ordered  = ordered,
    seed     = seed
  )

  estimatePLS_Step8(model)
}


resampleMCPLS_Fit <- function(model, ...) {
  args <- model$status$mcpls.update.args
  new.args <- list(...)
  args[names(new.args)] <- new.args

  do.call(updateModelFromFreeParTableMC, args)
}


robbinsMonro1951 <- function(p, f, tol, min.iter, max.iter, verbose,
                             polyak.juditsky, fn.args, ...) {

  args.required <- list(
      p               = p,
      f               = f,
      tol             = tol,
      miniter         = min.iter,
      maxiter         = max.iter,
      verbose         = verbose,
      Polyak_Juditsky = polyak.juditsky
  )

  args <- c(args.required, fn.args, list(...))

  mcfit <- do.call(SimDesign::RobbinsMonro, args)
  if (verbose) cat("\n")

  mcfit
}
