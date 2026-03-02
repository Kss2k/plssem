mcpls <- function(
  fit0,
  max.iter.mc = 100,
  mc.reps = 20000,
  rng.seed = NULL,
  tol = 1e-3,
  miniter = 5,
  maxiter = 250,
  fixed.seed = FALSE,
  verbose = TRUE,
  ...
) {

  data <- fit0$data
  vars <- colnames(data)
  ordered <- fit0$info$ordered

  PROBS <- stats::setNames(vector("list", length(ordered)), nm = ordered)
  for (ord in ordered) {
    freq <- table(data[,ord])
    pct  <- cumsum(freq) / sum(freq)
    PROBS[[ord]] <- pct[-length(pct)]
  }

  par0 <- getFreeParamsTable(fit0)
  par1 <- par0[c("lhs", "op", "rhs", "est")]

  if (fixed.seed && is.null(rng.seed)) {
    rng.seed <- floor(runif(1L, min = 0, max = 9999999))
    if (verbose) printf("Using fixed seed %i...\n", rng.seed)
  }

  .f <- function(p) {
    par1$est <- p

    sim <- simulateDataParTable(par1, N = mc.reps, seed = rng.seed)
    sim.ov <- ordinalizeDataFrame(sim$ov, PROBS = PROBS, ordered = ordered)

    fit0$data <- Rfast::standardise(as.matrix(sim.ov[vars]))
    fit0$matrices$S <- Rfast::cova(fit0$data)

    fit2 <- estimatePLS_Inner(fit0)
    par2 <- getFreeParamsTable(fit2) 

    par2$est - par0$est
  }

  mcfit <- SimDesign::RobbinsMonro(
    p = par1$est,
    f = .f,
    tol = tol,
    miniter = miniter,
    maxiter = maxiter,
    verbose = verbose
  )

  if (verbose) cat("\n")

  iter <- mcfit$iter
  if (iter >= maxiter && !fixed.seed) {
    rng.seed <- floor(runif(1L, min = 0, max = 9999999))
    warning2("Maximum number of iterations reached!\n",
             sprintf("Attempting to use fixed seed %i...", rng.seed))

    mcfit <- SimDesign::RobbinsMonro(
      p = mcfit$root,
      f = .f,
      tol = tol,
      miniter = miniter,
      maxiter = maxiter,
      verbose = verbose
    )

    if (verbose) cat("\n")

    iter <- iter + mcfit$iter
  }

  par1$est <- mcfit$root
  fit1 <- updateModelFromFreeParTableMC(
    parTable = par1,
    model    = fit0,
    mc.reps  = mc.reps,
    seed     = rng.seed
  )

  fit1$status$iterations <- mcfit$iter

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
  parTable <- getParTableEstimates(model)

  lhs <- parTable$lhs
  op  <- parTable$op
  rhs <- parTable$rhs

  cond1 <- !(lhs == rhs & op == "~~")
  cond2 <- !((grepl(":", lhs) | grepl(":", rhs)) & op == "~~")

  out <- parTable[cond1 & cond2, , drop = FALSE]
  attr(out, "cond") <- cond1 & cond2

  out
}


updateModelFromFreeParTableMC <- function(parTable, model, mc.reps, seed = NULL) {
  sim <- simulateDataParTable(parTable, N = mc.reps, seed = seed)

  lvs <- getLVs(parTable)
  indsLVs <- getIndsLVs(parTable, lVs = lvs)

  SC <- Rfast::cova(as.matrix(sim$all))

  ovs <- colnames(model$matrices$S)
  lvs <- colnames(model$matrices$C)

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
  for (i in seq_len(k)) for (j in seq_len(i - 1)) {
    lhs <- colnames(fitCov)[[i]]
    rhs <- rownames(fitCov)[[j]]

    if (!selectCov[lhs, rhs]) next

    par <- getpar(lhs = lhs, op = "~~", rhs = rhs)
    if (is.na(par)) par <- tryCatchNA(SC[lhs, rhs])

    fitCov[i, j] <- fitCov[j, i] <- par
  }

  model$fit$fitMeasurement <- fitMeasurement
  model$fit$fitStructural  <- fitStructural
  model$fit$fitCov         <- fitCov
  model$fit$fitTheta       <- fitTheta

  estimatePLS_Step8(model)
}
