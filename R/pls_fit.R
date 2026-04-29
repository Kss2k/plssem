getFitPLSModel <- function(model, consistent = TRUE) {
  lambda  <- model@matrices$lambda
  gamma   <- model@matrices$gamma
  preds   <- model@matrices$preds
  etas    <- model@info$etas
  xis     <- model@info$xis
  lvs     <- model@info$lvs
  lvs.lin <- model@info$lvs.linear
  inds    <- model@info$allInds
  inds.a  <- model@info$inds.a
  inds.b  <- model@info$inds.b
  indsLvs <- model@info$indsLvs
  modes   <- model@info$modes
  mode.a  <- model@info$mode.a
  ptl     <- model@parTableInput
  SC      <- model@matrices$SC

  fitMeasurement <- fitLambda <- fitWeights <- lambda
  fitMeasurement[TRUE] <- fitLambda[TRUE] <- fitWeights[TRUE] <- 0

  C <- model@matrices$C

  for (lv in lvs.lin) {
    inds.lv <- indsLvs[[lv]]
    mode.lv <- modes[[lv]]

    wq <- lambda[inds.lv, lv]
    lq <- SC[inds.lv, lv]
    pq <- switch(mode.lv, A = lq, B = wq, NA_real_)

    fitMeasurement[inds.lv, lv] <- pq
    fitWeights[inds.lv, lv]     <- wq
    fitLambda[inds.lv, lv]      <- lq
  }

  if (consistent) {
    Q                   <- getConstructQualities(model)
    fitMeasurement      <- getConsistentLoadings(model, Q = Q)
    fitLambda[, mode.a] <- fitMeasurement[, mode.a]
    C                   <- getConsistentCorrMat(model, Q = Q)

  } else {
    Q <- numeric(0)
    attr(Q, "admissible") <- TRUE

  }

  fitStructural       <- gamma
  fitStructural[TRUE] <- 0
  for (lv in lvs) {
    predsLv <- lvs[preds[, lv, drop = TRUE]]
    if (length(predsLv))
      fitStructural[predsLv, lv] <- getPathCoefs(lv, predsLv, C)
  }

  fitCov     <- C
  fitCovProj <- t(fitStructural) %*% C %*% fitStructural
  fitCovRes  <- diag2(fitCov) - diag2(fitCovProj)
  fitCov[etas, etas]  <- fitCovRes[etas, etas]
  fitCov[etas, xis]   <- fitCov[xis, etas] <- 0

  k        <- length(inds)
  fitTheta <- matrix(0, nrow = k, ncol = k, dimnames = list(inds, inds))
  crossLoaded <- apply(
    X      = fitMeasurement,
    MARGIN = 1L,
    FUN    = \(x) sum(abs(x) > .Machine$double.xmin) > 1L
  )

  fitThetaFull <- model@matrices$SC[inds, inds]
  is.formative <- inds %in% inds.b
  mask         <- outer(is.formative, is.formative, FUN = "&")
  fitTheta[mask] <- fitThetaFull[mask]

  warnif(any(crossLoaded),
         "Did not expect any cross loaded indicators,\n",
         "when calculating indicator residuals!")

  for (ind in inds.a) {                                  # Guard for NaN in fitMeasurement
    j   <- max(which.max(abs(fitMeasurement[ind, ])), 1) # max(numeric(0), 1) = 1
    r   <- fitMeasurement[ind, j]
    v   <- SC[ind, ind]

    fitTheta[ind, ind] <- v - r^2
  }

  list(
    fitMeasurement = plssemMatrix(fitMeasurement),
    fitStructural  = plssemMatrix(fitStructural),
    fitCov         = plssemMatrix(fitCov, symmetric = TRUE),
    fitTheta       = plssemMatrix(fitTheta, symmetric = TRUE),
    fitWeights     = plssemMatrix(fitWeights),
    fitLambda      = plssemMatrix(fitLambda),
    fitC           = plssemMatrix(C),
    Q              = plssemVector(Q)
  )
}


modelFitIsAdmissible <- function(fit) {
  # Simple check to see if model fit is (in)admissible

  Q.admissible <- (
    is.null(attr(fit$Q, "admissible")) ||
    isTRUE(attr(fit$Q, "admissible"))
  )

  (
    !anyNA(fit$fitWeights)       &&
    !anyNA(fit$fitLambda)        &&
    !anyNA(fit$fitStructural)    &&
    !anyNA(fit$fitTheta)         &&
    !anyNA(fit$fitCov)           &&
    !anyNA(fit$fitC)             &&
    isPositiveDefinite(fit$fitC) &&
    all(diag(fit$fitTheta) >= 0) &&
    all(diag(fit$fitCov) >= 0)   &&
    Q.admissible
  )
}


getParamVecNames <- function(model) {
  selectLambda <- model@matrices$select$lambda
  modes        <- model@info$modes
  lvs.linear   <- model@info$lvs.linear
  lambda       <- selectLambda

  for (j in lvs.linear) {
    op <- switch(modes[[j]], A = "=~", B = "<~", "=~")
    for (i in rownames(lambda))
      lambda[i, j] <- paste0(j, op, i)
  }

  selectGamma <- model@matrices$select$gamma
  gamma       <- selectGamma
  for (j in colnames(gamma)) for (i in rownames(gamma))
    gamma[i, j] <- paste0(j, "~", i)

  selectCov <- model@matrices$select$cov
  psi       <- selectCov
  for (j in colnames(psi)) for (i in rownames(psi))
    psi[i, j] <- paste0(j, "~~", i)

  selectTheta <- model@matrices$select$theta
  theta       <- selectTheta
  for (j in colnames(theta)) for (i in rownames(theta))
    theta[i, j] <- paste0(j, "~~", i)

  c(lambda[selectLambda], gamma[selectGamma], psi[selectCov], theta[selectTheta])
}


extractCoefs <- function(model) {
  fit <- model@fit

  lambda       <- fit$fitMeasurement
  selectLambda <- model@matrices$select$lambda

  gamma       <- fit$fitStructural
  selectGamma <- model@matrices$select$gamma

  fitCov    <- fit$fitCov
  selectCov <- model@matrices$select$cov

  fitTheta    <- fit$fitTheta
  selectTheta <- model@matrices$select$theta

  out <- c(
    lambda[selectLambda],
    gamma[selectGamma],
    fitCov[selectCov],
    fitTheta[selectTheta]
  )

  names(out) <- model@params$names
  plssemVector(out)
}


computeFactorScores <- function(model) {
  W <- model@matrices$lambda
  X <- model@data
  F <- X %*% W

  if (!model@info$standardize)
    F <- Rfast::standardise(F)

  F
}


getEstimatorFromInfo <- function(info) {
  consistent <- info$consistent
  is.mcpls   <- info$is.mcpls
  is.mlm     <- info$is.mlm
  is.ord     <- info$is.probit || (info$is.mcpls && length(info$ordered))

  estimator <- "PLS"
  if (consistent || is.mcpls) estimator <- paste0(estimator, "c")
  if (is.mlm)                 estimator <- paste0(estimator, "-MLM")
  if (is.ord)                 estimator <- paste0("Ord", estimator)
  if (is.mcpls)               estimator <- paste0("MC", estimator)

  estimator
}


refreshModelParams <- function(model, update.names = TRUE) {
  # Should we update names?
  if (update.names)
    model@params$names <- getParamVecNames(model)

  # Single level params
  model@params$values <- extractCoefs(model)
  model@params$se     <- rep(NA_real_, length(model@params$values))

  # Multilevel/Mixed-Effect params
  if (isMLM(model))
    model <- refreshLmerParams(model)

  model
}


refreshLmerParams <- function(model) {
  lmerFit <- modelFitLmer(model)

  if (!isMLM(model) || is.null(lmerFit))
    return(model)

  coefs.x <- model@params$values
  coefs.y <- lmerFit$values

  common  <- intersect(names(coefs.x), names(coefs.y))
  new     <- setdiff(names(coefs.y),   names(coefs.x))

  coefs.x[common] <- coefs.y[common]
  coefs.all       <- c(coefs.x, coefs.y[new])

  model@params$values <- plssemVector(coefs.all)
  model@params$se     <- rep(NA_real_, length(coefs.all))

  model
}

