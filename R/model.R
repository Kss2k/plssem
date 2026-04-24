OPERATORS <- c("<~", "~~", "=~", "~1", "~", "|")
MOPS <- c("<~", "=~")
MAT_STRUCT <- matrix(0, nrow = 0, ncol = 0)
DF_STRUCT <- data.frame(NULL)
VEC_STRUCT <- numeric(0)


specifyModel <- function(syntax, data, ...) {
  split      <- splitHigherOrderModel(syntax)
  parTableO2 <- split$parTableO2
  parTableO1 <- split$parTableO1
  hiOrdLVs   <- split$higherOrderLVs

  firstOrder <- specifySubModel(
    parTable       = parTableO1,
    data           = data,
    is.lower.order = NROW(parTableO2) > 0,
    ...
  )

  secondOrder <- specifySubModel(
    parTable       = parTableO2,
    data           = firstOrder$data %*% firstOrder$matrices$lambda, # placeholder
    higherOrderLVs = hiOrdLVs,
    is.lower.order = FALSE,
    ...
  )

  info1 <- firstOrder$info
  info2 <- secondOrder$info

  if (is.null(secondOrder)) {
    inds.x <- info1$inds.x
    inds.y <- info1$inds.y

  } else {
    lvs     <- info1$lvs
    etas    <- info2$etas
    ovInds  <- info1$allInds
    lvInds  <- intersect(info2$allInds, lvs)
    indsMap <- c(info1$indsLvs, info2$indsLvs)

    lowInds <- unique(unlist(indsMap[lvInds])) # inds which are downstream of HOLVs
    etaInds <- unique(unlist(indsMap[etas]))
    yInds   <- c(lowInds, etaInds)

    inds.y <- intersect(ovInds, yInds)
    inds.x <- setdiff(ovInds, yInds)
  }

  is.high.ord <- !is.null(secondOrder)
  ordered <- union(info1$ordered, info2$ordered)

  is.mcpls <- (
    isTRUE(info1$is.mcpls) || isTRUE(info2$is.mcpls) || # check if eiter sub model is
    (length(ordered) > 0 && is.high.ord)                # We also switch to MC-PLS for
  )                                                     # higher order probit models
    
  if (is.mcpls && isTRUE(firstOrder$info$is.probit)) {
    # If the full model uses MC-PLS we should disable probit estimation
    # for the first order model
    firstOrder$info$is.probit <- FALSE
    firstOrder$matrices$S <- getCorrMat(firstOrder$data, probit = FALSE)
  }

  # Class Fields
  out <- list(
    submodels = list(
      firstOrder  = firstOrder,
      secondOrder = secondOrder
    ),

    matrices = list(
      S           = MAT_STRUCT,
      C           = MAT_STRUCT,
      firstOrder  = firstOrder$matrices,
      secondOrder = secondOrder$matrices
    ),

    data = firstOrder$data,

    info = list(
      lvs.linear   = union(info1$lvs.linear, info2$lvs.linear),
      lvs          = union(info1$lvs, info2$lvs),
      lvs.hi.ord   = hiOrdLVs,

      mode.a       = union(info1$mode.a, info2$mode.a),
      mode.b       = union(info1$mode.b, info2$mode.b),
      modes        = namedListUnion(info1$modes, info2$modes),

      inds.x       = inds.x, # Observed indicators
      inds.y       = inds.y, # Observed indicators

      cluster      = info1$cluster,
      ordered      = info1$ordered,

      is.mlm       = isTRUE(info1$is.mlm) || isTRUE(info2$is.mlm),
      is.mcpls     = is.mcpls,
      is.probit    = isTRUE(info1$is.probit) || isTRUE(info2$is.probit),
      is.cfa       = isTRUE(info1$is.cfa) && (is.null(info1$is.cfa) || info1$is.cfa),
      is.high.ord  = is.high.ord,

      n            = info1$n,
      estimator    = info1$estimator,
      standardized = info1$standardized,
      verbose      = isTRUE(info1$verbose) || isTRUE(info2$verbose),

      mc.args      = info1$mc.args,
      boot         = info1$boot,
      args         = list(...)
    ),

    status = firstOrder$status,

    parTable = DF_STRUCT,

    params = list(
      names  = VEC_STRUCT,
      values = VEC_STRUCT,
      se     = VEC_STRUCT,
      vcov   = NULL
    ),

    fit = list(
      fitMeasurement = MAT_STRUCT,
      fitLambda      = MAT_STRUCT,
      fitWeights     = MAT_STRUCT,
      fitTheta       = MAT_STRUCT,
      fitC           = MAT_STRUCT,
      fitCov         = MAT_STRUCT,
      fitStructural  = MAT_STRUCT,
      Q              = VEC_STRUCT
    )
  )

  class(out) <- "plssem"
  out
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

  syntax        <- parsed$syntax
  data          <- parsed$data
  pt            <- parsed$parTable.pls
  cluster       <- parsed$cluster
  lme4.syntax   <- parsed$lme4.syntax
  is.mlm        <- parsed$is.mlm
  is.mcpls      <- parsed$is.mcpls
  intTermElems  <- parsed$intTermElems
  intTermNames  <- parsed$intTermNames
  is.nlin       <- parsed$is.nlin
  ordered       <- parsed$ordered
  is.probit     <- parsed$is.probit
  consistent    <- parsed$consistent

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
  info$rng.seed      <- floor(stats::runif(1L, min=0, max=9999999))
  info$n             <- NROW(preppedData$X)
  info$estimator     <- getEstimatorFromInfo(info)
  info$verbose       <- verbose
  info$standardized  <- standardize
  info$reliabilities <- reliabilities

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
    bootstrap = bootstrap,
    ncpus     = boot.ncpus,
    parallel  = boot.parallel,
    R         = boot.R,
    iseed     = boot.iseed,
    optimize  = boot.optimize,
    mc.boot.control = mc.boot.control
  )

  matrices$S  <- preppedData$S
  matrices$SC <- diagPartitioned(matrices$S, matrices$C)

  model <- list(
    parTable.input = pt,
    parTable       = NULL,
    matrices       = matrices,
    data           = preppedData$X,
    factorScores   = NULL,
    info           = info,
    params         = NULL,
    fit            = NULL,
    fit.c          = NULL,
    fit.u          = NULL,
    params         = NULL,
    status         = list(
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

  model$params <- list(
    names      = parnames,
    values     = rep(NA_real_, k),
    values.old = NULL,
    se         = rep(NA_real_, k),
    vcov       = NULL
  )

  class(model) <- "plssemSubModel"
  model
}


# badly named function -- since if also returns some useful info
initMatrices <- function(pt, higherOrderLVs = NULL) {
  lvs.linear <- getLVs(pt)

  k      <- length(lvs.linear)
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

  indsLvs <- vector("list", length(lvs))

  names(indsLvs) <- lvs
  for (lv in lvs) {
    indsLv <- pt[pt$lhs == lv & pt$op %in% c("=~", "<~"), "rhs"]
    indsLvs[[lv]] <- indsLv
  }

  allInds <- unique(unlist(indsLvs))
  inds.x <- unique(unlist(indsLvs[xis]))
  inds.y <- unique(unlist(indsLvs[etas]))
  inds.a <- unique(unlist(indsLvs[mode.a]))
  inds.b <- unique(unlist(indsLvs[mode.b]))

  # Lamdba ---------------------------------------------------------------------
  lambda <- matrix(0, nrow = length(allInds), ncol = length(lvs),
                   dimnames = list(allInds, lvs))
  selectLambda <- matrix(FALSE, nrow = length(allInds), ncol = length(lvs),
                   dimnames = list(allInds, lvs))
  for (lv in lvs) {
    lambda[indsLvs[[lv]], lv] <- 1
    selectLambda[indsLvs[[lv]], lv] <- TRUE
  }

  # Gamma ----------------------------------------------------------------------
  gamma <- matrix(0, nrow = length(lvs), ncol = length(lvs),
                  dimnames = list(lvs, lvs))
  # Selection Matrix
  selectGamma <- matrix(FALSE, nrow = length(lvs), ncol = length(lvs),
                  dimnames = list(lvs, lvs))
  # Predecessors and successors
  preds <- succs <- matrix(FALSE, nrow = length(lvs), ncol = length(lvs),
                           dimnames = list(lvs, lvs))
  for (lv in lvs) {
    predsLv <- pt[pt$lhs == lv & pt$op == "~", "rhs"]
    succsLv <- pt[pt$rhs == lv & pt$op == "~", "lhs"]
    preds[predsLv, lv] <- TRUE
    succs[succsLv, lv] <- TRUE
    # selectionmatrix
    selectGamma[predsLv, lv] <- TRUE
  }

  is.nlin      <- grepl(":", lvs)
  preds.linear <- preds
  succs.linear <- succs
  preds.linear[is.nlin,] <- FALSE
  preds.linear[,is.nlin] <- FALSE
  succs.linear[is.nlin,] <- FALSE
  succs.linear[,is.nlin] <- FALSE

  is.cfa <- !any(preds.linear) && !any(succs.linear)

  preds.cfa                       <- preds
  preds.cfa[TRUE]                 <- FALSE
  preds.cfa[upper.tri(preds.cfa)] <- TRUE

  preds.cfa[is.nlin,] <- FALSE
  preds.cfa[,is.nlin] <- FALSE
  succs.cfa           <- t(preds.cfa)

  # Covariance Matrix xis ------------------------------------------------------
  xis <- lvs[!lvs %in% pt[pt$op == "~", "lhs"]]
  selectCov <- matrix(
    FALSE, nrow = length(lvs), ncol = length(lvs),
    dimnames = list(lvs, lvs)
  )

  diag(selectCov) <- TRUE
  selectCov[xis, xis][lower.tri(selectCov[xis, xis])] <- TRUE

  # Residual Covariance Matrix indicators --------------------------------------
  k <- length(allInds)
  selectTheta <- matrix(FALSE, nrow = k, ncol = k,
                        dimnames = list(allInds, allInds))
  diag(selectTheta) <- TRUE
  is.formative <- allInds %in% inds.b
  selectTheta[outer(is.formative, is.formative, "&")] <- TRUE
  selectTheta[upper.tri(selectTheta, diag = FALSE)] <- FALSE

  # Selection Matric Indicators ------------------------------------------------
  Ip <- diag(nrow = nrow(lambda))
  colnames(Ip) <- rownames(Ip) <- rownames(lambda)

  nlinSelectFrom <- outer(lvs, lvs, Vectorize(f1))
  dimnames(nlinSelectFrom) <- list(lvs, lvs)

  C <- diag(nrow(gamma))
  dimnames(C) <- dimnames(gamma)

  # Higher Order Composites ----------------------------------------------------
  isHigherOrderComposite <- stats::setNames(
    lvs.linear %in% mode.b & lvs.linear %in% higherOrderLVs, nm = lvs.linear
  )
  higherOrderComposites <- lvs.linear[isHigherOrderComposite]


  # Create output lists --------------------------------------------------------
  matrices <- list(
    lambda = lambda,
    gamma = gamma,
    preds = preds,
    succs = succs,
    succs.linear = succs.linear,
    preds.linear = preds.linear,
    preds.cfa    = preds.cfa,
    succs.cfa    = succs.cfa,
    outerWeights = getNonZeroElems(lambda),
    Ip = Ip,
    C  = C,
    S  = NULL,
    SC = NULL,
    select = list(
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

  list(
    matrices = matrices,
    info = info
  )
}


getFitPLSModel <- function(model, consistent = TRUE) {
  lambda  <- model$matrices$lambda
  gamma   <- model$matrices$gamma
  preds   <- model$matrices$preds
  etas    <- model$info$etas
  xis     <- model$info$xis
  lvs     <- model$info$lvs
  lvs.lin <- model$info$lvs.linear
  inds    <- model$info$allInds
  inds.a  <- model$info$inds.a
  inds.b  <- model$info$inds.b
  indsLvs <- model$info$indsLvs
  modes   <- model$info$modes
  mode.a  <- model$info$mode.a
  ptl     <- model$parTable.input
  SC      <- model$matrices$SC
  k       <- length(lvs.lin)

  # measurement model
  fitMeasurement <- fitLambda <- fitWeights <- lambda
  fitMeasurement[TRUE] <- fitLambda[TRUE] <- fitWeights[TRUE] <- 0

  for (lv in lvs.lin) {
    inds.lv <- indsLvs[[lv]]
    mode.lv <- modes[[lv]]

    wq <- lambda[inds.lv, lv]
    lq <- SC[inds.lv, lv]
    pq <- switch(mode.lv,
      A = lq,
      B = wq,
      NA_real_
    )

    fitMeasurement[inds.lv, lv] <- pq
    fitWeights[inds.lv, lv]     <- wq
    fitLambda[inds.lv, lv]      <- lq
  }

  # Caluculate consistent weights and correlations
  if (consistent) {
    # We want to correct both for the errors causes by using the CEXP
    # estimator, compared to the probit estimator. As well as the bias
    # caused by ignoring measurement error.

    Q <- getConstructQualities(model)

    fitMeasurement <- getConsistentLoadings(model, Q = Q)
    fitLambda[,mode.a] <- fitMeasurement[,mode.a]
    model$matrices$C <- getConsistentCorrMat(model, Q = Q)

  } else {
    Q <- NULL

  }

  # structural model
  fitStructural       <- gamma
  fitStructural[TRUE] <- 0
  for (lv in lvs) {
    predsLv <- lvs[preds[ , lv, drop = TRUE]]

    if (length(predsLv))
      fitStructural[predsLv, lv] <- getPathCoefs(lv, predsLv,  model$matrices$C)
  }

  # Covariance matrix
  fitCovXi  <- model$matrices$C[xis, xis, drop = FALSE]
  fitCov    <- model$matrices$C
  fitCovProj <- t(fitStructural) %*% model$matrices$C %*% fitStructural
  fitCovRes  <- diag2(fitCov) - diag2(fitCovProj) # Removing diag2() might be better here
                                                  # But I'm not sure the residual covariances
                                                  # will be consistent estimates
  fitCov[etas, etas] <- fitCovRes[etas, etas]
  fitCov[etas, xis] <- fitCov[xis, etas] <- 0

  # Indicator Residuals
  k <- length(inds)
  fitTheta <- matrix(0, nrow = k, ncol = k,
                     dimnames = list(inds, inds))
  crossLoaded <- apply(
    X      = fitMeasurement,
    MARGIN = 1L,
    FUN    = \(x) sum(abs(x) > .Machine$double.xmin) > 1L
  )

  fitThetaFull <- model$matrices$SC[inds, inds]
  is.formative <- inds %in% inds.b
  mask <- outer(is.formative, is.formative, FUN = "&")
  fitTheta[mask] <- fitThetaFull[mask]

  warnif(any(crossLoaded),
         "Did not expect any cross loaded indicators,\n",
         "when calculating indicator residuals!")

  for (ind in inds.a) {
    j  <- which.max(abs(fitMeasurement[ind, ]))
    r  <- fitMeasurement[ind, j]
    v  <- SC[ind, ind]
    res <- v - r^2

    fitTheta[ind, ind] <- res
  }

  list(
    fitMeasurement = plssemMatrix(fitMeasurement),
    fitStructural  = plssemMatrix(fitStructural),
    fitCov         = plssemMatrix(fitCov, symmetric = TRUE),
    fitTheta       = plssemMatrix(fitTheta, symmetric = TRUE),
    fitWeights     = plssemMatrix(fitWeights),
    fitLambda      = plssemMatrix(fitLambda),
    fitC           = plssemMatrix(model$matrices$C),
    Q              = plssemVector(Q)
  )
}


getParamVecNames <- function(model) {
  selectLambda <- model$matrices$select$lambda
  modes        <- model$info$modes
  lvs.linear   <- model$info$lvs.linear
  lambda       <- selectLambda

  for (j in lvs.linear) {
    mode <- modes[[j]]
    op <- switch(mode, A = "=~", B = "<~", "=~")

    for (i in rownames(lambda))
      lambda[i, j] <- paste0(j, op, i)

  }

  selectGamma <- model$matrices$select$gamma
  gamma       <- selectGamma

  for (j in colnames(gamma)) for (i in rownames(gamma))
    gamma[i, j] <- paste0(j, "~", i)

  selectCov <- model$matrices$select$cov
  psi       <- selectCov

  for (j in colnames(psi)) for (i in rownames(psi))
    psi[i, j] <- paste0(j, "~~", i)

  selectTheta <- model$matrices$select$theta
  theta       <- selectTheta

  for (j in colnames(theta)) for (i in rownames(theta))
    theta[i, j] <- paste0(j, "~~", i)

  c(lambda[selectLambda], gamma[selectGamma], psi[selectCov], theta[selectTheta])
}


extractCoefs <- function(model) {
  fit <- model$fit

  lambda       <- fit$fitMeasurement
  selectLambda <- model$matrices$select$lambda

  gamma       <- fit$fitStructural
  selectGamma <- model$matrices$select$gamma

  fitCov    <- fit$fitCov
  selectCov <- model$matrices$select$cov

  fitTheta    <- fit$fitTheta
  selectTheta <- model$matrices$select$theta

  out <- c(
    lambda[selectLambda],
    gamma[selectGamma],
    fitCov[selectCov],
    fitTheta[selectTheta]
  )

  names(out) <- model$params$names

  plssemVector(out)
}


getFactorScores <- function(model) {
  matrices <- model$matrices
  W <- model$matrices$lambda
  X <- model$data

  Rfast::standardise(X %*% W)
}


getEstimatorFromInfo <- function(info) {
  consistent <- info$consistent
  is.mcpls   <- info$is.mcpls
  is.mlm     <- info$is.mlm
  is.ord     <- info$is.probit || (info$is.mcpls && length(info$ordered))

  estimator  <- "PLS"

  if (consistent || is.mcpls) estimator <- paste0(estimator, "c")
  if (is.mlm)                 estimator <- paste0(estimator, "-MLM")
  if (is.ord)                 estimator <- paste0("Ord", estimator)
  if (is.mcpls)               estimator <- paste0("MC", estimator)

  estimator
}
