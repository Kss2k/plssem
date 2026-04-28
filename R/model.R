OPERATORS <- c("<~", "~~", "=~", "~1", "~", "|")
MOPS      <- c("<~", "=~")
MAT_STRUCT <- matrix(0, nrow = 0, ncol = 0)
DF_STRUCT  <- data.frame(NULL)
VEC_STRUCT <- numeric(0)


specifyModel <- function(syntax, data, ...) {
  parTable <- modsem::modsemify(syntax, parentheses.as.string = TRUE)
  hiOrd    <- getHigherOrderLVs(parTable)
  specifyModelParTable(parTable = parTable, data = data, higherOrderLVs = hiOrd, ...)
}


specifyModelParTable <- function(parTable, data, higherOrderLVs = NULL, ...) {
  split      <- splitHigherOrderParTable(parTable)
  parTableO2 <- split$parTableO2
  parTableO1 <- split$parTableO1

  firstOrder <- specifySubModel(
    parTable       = parTableO1,
    data           = data,
    is.lower.order = NROW(parTableO2) > 0,
    higherOrderLVs = higherOrderLVs,
    ...
  )

  if (NROW(parTableO2) > 0) {

    higherOrder <- specifyModelParTable(
      parTable       = parTableO2,
      data           = getSecondOrderInputData(firstOrder),
      higherOrderLVs = higherOrderLVs,
      ...
    )

  } else {
    higherOrder <- NULL

  }

  higherOrderModel(firstOrder) <- higherOrder
  firstOrder@info$is.high.ord <- !is.null(higherOrder)
  firstOrder
}


specifySubModel <- function(parTable,
                            data,
                            is.lower.order     = FALSE,
                            consistent         = TRUE,
                            missing            = "listwise",
                            standardize        = TRUE,
                            ordered            = NULL,
                            probit             = NULL,
                            mcpls              = NULL,
                            tolerance          = 1e-5,
                            max.iter.0_5       = 100,
                            mc.max.iter        = 250,
                            mc.min.iter        = 5,
                            mc.reps            = 20000,
                            mc.tol             = 1e-3,
                            mc.fixed.seed      = FALSE,
                            mc.polyak.juditsky = FALSE,
                            mc.fn.args         = list(),
                            verbose            = interactive(),
                            bootstrap          = FALSE,
                            boot.ncpus         = 1L,
                            boot.parallel      = "no",
                            boot.R             = 50L,
                            boot.iseed         = NULL,
                            boot.optimize      = FALSE,
                            mc.boot.control    = list(),
                            knn.k              = 5,
                            reliabilities      = NULL,
                            higherOrderLVs     = NULL) {
  if (is.null(parTable))
    return(NULL)

  parsed <- parseModelArguments(
    parTable       = parTable,
    data           = data,
    ordered        = ordered,
    probit         = probit,
    mcpls          = mcpls,
    consistent     = consistent,
    is.lower.order = is.lower.order
  )

  syntax       <- parsed$syntax
  data         <- parsed$data
  pt           <- parsed$parTable.pls
  cluster      <- parsed$cluster
  lme4.syntax  <- parsed$lme4.syntax
  is.mlm       <- parsed$is.mlm
  is.mcpls     <- parsed$is.mcpls
  intTermElems <- parsed$intTermElems
  intTermNames <- parsed$intTermNames
  is.nlin      <- parsed$is.nlin
  ordered      <- parsed$ordered
  is.probit    <- parsed$is.probit
  consistent   <- parsed$consistent

  matricesAndInfo <- initMatrices(pt, higherOrderLVs = higherOrderLVs)
  matrices        <- matricesAndInfo$matrices
  info            <- matricesAndInfo$info

  preppedData <- getPLS_Data(
    data        = data,
    indicators  = matricesAndInfo$info$allInds,
    consistent  = consistent,
    cluster     = cluster,
    standardize = standardize,
    ordered     = ordered,
    is.probit   = is.probit,
    missing     = missing,
    knn.k       = knn.k
  )

  inds.x    <- info$inds.x
  inds.y    <- info$inds.y
  ordered.x <- intersect(inds.x, ordered)
  ordered.y <- intersect(inds.y, ordered)

  info$lme4.syntax   <- lme4.syntax
  info$is.mlm        <- is.mlm
  info$is.mcpls      <- is.mcpls
  info$is.probit     <- is.probit
  info$cluster       <- cluster
  info$consistent    <- consistent
  info$ordered       <- ordered
  info$ordered.x     <- ordered.x
  info$ordered.y     <- ordered.y
  info$intTermElems  <- intTermElems
  info$intTermNames  <- intTermNames
  info$is.nlin       <- is.nlin
  info$rng.seed      <- floor(stats::runif(1L, min = 0, max = 9999999))
  info$n             <- NROW(preppedData$X)
  info$estimator     <- getEstimatorFromInfo(info)
  info$verbose       <- verbose
  info$standardized  <- standardize
  info$reliabilities <- reliabilities
  info$is.high.ord   <- FALSE
  info$is.lower.order <- isTRUE(is.lower.order)

  info$mc.args <- list(
    min.iter        = mc.min.iter,
    max.iter        = mc.max.iter,
    mc.reps         = mc.reps,
    tol             = mc.tol,
    fixed.seed      = mc.fixed.seed,
    polyak.juditsky = mc.polyak.juditsky,
    fn.args         = mc.fn.args,
    rng.seed        = NULL,
    p.start         = NULL
  )

  info$boot <- list(
    bootstrap       = bootstrap,
    ncpus           = boot.ncpus,
    parallel        = boot.parallel,
    R               = boot.R,
    iseed           = boot.iseed,
    optimize        = boot.optimize,
    mc.boot.control = mc.boot.control
  )

  matrices$S  <- preppedData$S
  matrices$SC <- diagPartitioned(matrices$S, matrices$C)

  model <- PlsModel(
    parTableInput = pt,
    matrices      = matrices,
    data          = preppedData$X,
    info          = info,
    status        = list(
      convergence    = FALSE,
      iterations     = 0L,
      iterations.0_5 = 0L,
      tolerance      = tolerance,
      max.iter.0_5   = max.iter.0_5,
      is.admissible  = TRUE
    )
  )

  parnames <- getParamVecNames(model)
  k        <- length(parnames)

  model@params <- list(
    names      = parnames,
    values     = rep(NA_real_, k),
    values.old = NULL,
    se         = rep(NA_real_, k),
    vcov       = NULL
  )

  model
}


# badly named function -- since it also returns some useful info
initMatrices <- function(pt, higherOrderLVs = NULL) {
  lvs.linear <- getLVs(pt)

  mode.a <- getReflectiveLVs(pt)
  mode.b <- getFormativeLVs(pt)

  getmode <- function(x)
    ifelse(x %in% mode.a, yes = "A", no = ifelse(x %in% mode.b, yes = "B", no = NA))

  modes <- stats::setNames(
    vapply(lvs.linear, FUN.VALUE = character(1L), FUN = getmode),
    nm = lvs.linear
  )

  etas <- unique(pt[pt$op == "~", "lhs"])
  lvs  <- c(lvs.linear, getIntTerms(pt))
  xis  <- lvs[!lvs %in% etas]

  indsLvs        <- vector("list", length(lvs))
  names(indsLvs) <- lvs
  for (lv in lvs) {
    indsLvs[[lv]] <- pt[pt$lhs == lv & pt$op %in% c("=~", "<~"), "rhs"]
  }

  allInds <- unique(unlist(indsLvs))
  inds.x  <- unique(unlist(indsLvs[xis]))
  inds.y  <- unique(unlist(indsLvs[etas]))
  inds.a  <- unique(unlist(indsLvs[mode.a]))
  inds.b  <- unique(unlist(indsLvs[mode.b]))

  # Lambda -----------------------------------------------------------------
  lambda       <- matrix(0,     nrow = length(allInds), ncol = length(lvs),
                         dimnames = list(allInds, lvs))
  selectLambda <- matrix(FALSE, nrow = length(allInds), ncol = length(lvs),
                         dimnames = list(allInds, lvs))
  for (lv in lvs) {
    lambda[indsLvs[[lv]], lv]       <- 1
    selectLambda[indsLvs[[lv]], lv] <- TRUE
  }

  # Gamma ------------------------------------------------------------------
  gamma       <- matrix(0,     nrow = length(lvs), ncol = length(lvs),
                        dimnames = list(lvs, lvs))
  selectGamma <- matrix(FALSE, nrow = length(lvs), ncol = length(lvs),
                        dimnames = list(lvs, lvs))
  preds <- succs <- matrix(FALSE, nrow = length(lvs), ncol = length(lvs),
                           dimnames = list(lvs, lvs))
  for (lv in lvs) {
    predsLv <- pt[pt$lhs == lv & pt$op == "~", "rhs"]
    succsLv <- pt[pt$rhs == lv & pt$op == "~", "lhs"]
    preds[predsLv, lv]       <- TRUE
    succs[succsLv, lv]       <- TRUE
    selectGamma[predsLv, lv] <- TRUE
  }

  is.nlin      <- grepl(":", lvs)
  preds.linear <- preds
  succs.linear <- succs
  preds.linear[is.nlin, ] <- FALSE
  preds.linear[, is.nlin] <- FALSE
  succs.linear[is.nlin, ] <- FALSE
  succs.linear[, is.nlin] <- FALSE

  is.cfa <- !any(preds.linear) && !any(succs.linear)

  preds.cfa              <- preds
  preds.cfa[TRUE]        <- FALSE
  preds.cfa[upper.tri(preds.cfa)] <- TRUE
  preds.cfa[is.nlin, ]   <- FALSE
  preds.cfa[, is.nlin]   <- FALSE
  succs.cfa              <- t(preds.cfa)

  # Covariance xis ---------------------------------------------------------
  xis       <- lvs[!lvs %in% pt[pt$op == "~", "lhs"]]
  selectCov <- matrix(FALSE, nrow = length(lvs), ncol = length(lvs),
                      dimnames = list(lvs, lvs))
  diag(selectCov) <- TRUE
  selectCov[xis, xis][lower.tri(selectCov[xis, xis])] <- TRUE

  # Residual covariances ---------------------------------------------------
  k           <- length(allInds)
  selectTheta <- matrix(FALSE, nrow = k, ncol = k,
                        dimnames = list(allInds, allInds))
  diag(selectTheta)     <- TRUE
  is.formative          <- allInds %in% inds.b
  selectTheta[outer(is.formative, is.formative, "&")] <- TRUE
  selectTheta[upper.tri(selectTheta, diag = FALSE)]   <- FALSE

  Ip <- diag(nrow = nrow(lambda))
  colnames(Ip) <- rownames(Ip) <- rownames(lambda)

  nlinSelectFrom          <- outer(lvs, lvs, Vectorize(f1))
  dimnames(nlinSelectFrom) <- list(lvs, lvs)

  C <- diag(nrow(gamma))
  dimnames(C) <- dimnames(gamma)

  isHigherOrderComposite <- stats::setNames(
    lvs.linear %in% mode.b & lvs.linear %in% higherOrderLVs,
    nm = lvs.linear
  )
  higherOrderComposites <- lvs.linear[isHigherOrderComposite]

  matrices <- list(
    lambda       = lambda,
    gamma        = gamma,
    preds        = preds,
    succs        = succs,
    succs.linear = succs.linear,
    preds.linear = preds.linear,
    preds.cfa    = preds.cfa,
    succs.cfa    = succs.cfa,
    outerWeights = getNonZeroElems(lambda),
    Ip           = Ip,
    C            = C,
    S            = NULL,
    SC           = NULL,
    select       = list(
      lambda   = selectLambda,
      gamma    = selectGamma,
      cov      = selectCov,
      theta    = selectTheta,
      nlinFrom = nlinSelectFrom
    )
  )

  info <- list(
    indsLvs                = indsLvs,
    allInds                = allInds,
    lvs                    = lvs,
    lvs.linear             = lvs.linear,
    etas                   = etas,
    xis                    = xis,
    inds.x                 = inds.x,
    inds.y                 = inds.y,
    inds.a                 = inds.a,
    inds.b                 = inds.b,
    mode.a                 = mode.a,
    mode.b                 = mode.b,
    modes                  = modes,
    is.cfa                 = is.cfa,
    higherOrderLVs         = higherOrderLVs,
    isHigherOrderComposite = isHigherOrderComposite,
    higherOrderComposites  = higherOrderComposites
  )

  list(matrices = matrices, info = info)
}


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
    model@matrices$C    <- getConsistentCorrMat(model, Q = Q)

  } else {
    Q <- NULL
    attr(Q, "admissible") <- TRUE

  }

  fitStructural       <- gamma
  fitStructural[TRUE] <- 0
  for (lv in lvs) {
    predsLv <- lvs[preds[, lv, drop = TRUE]]
    if (length(predsLv))
      fitStructural[predsLv, lv] <- getPathCoefs(lv, predsLv, model@matrices$C)
  }

  fitCov     <- model@matrices$C
  fitCovProj <- t(fitStructural) %*% model@matrices$C %*% fitStructural
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
    fitC           = plssemMatrix(model@matrices$C),
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
  Rfast::standardise(X %*% W)
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
