OPERATORS <- c("~~", "=~", "~1", "~", "|", "<~")

specifyModel <- function(syntax,
                         data,
                         consistent         = TRUE,
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
                         verbose            = interactive()) {

  parsed <- parseModelArguments(
    syntax     = syntax,
    data       = data,
    ordered    = ordered,
    probit     = probit,
    mcpls      = mcpls,
    consistent = consistent
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

  matricesAndInfo <- initMatrices(pt)
  matrices        <- matricesAndInfo$matrices
  info            <- matricesAndInfo$info

  preppedData <- getPLS_Data(
    data        = data,
    indicators  = matricesAndInfo$info$allInds,
    consistent  = consistent,
    cluster     = cluster,
    standardize = standardize,
    ordered     = ordered,
    is.probit   = is.probit
  )

  inds.x    <- info$inds.x
  inds.y    <- info$inds.y
  ordered.x <- intersect(inds.x, ordered)
  ordered.y <- intersect(inds.y, ordered)

  info$lme4.syntax  <- lme4.syntax
  info$is.mlm       <- is.mlm
  info$is.mcpls     <- is.mcpls
  info$is.probit    <- is.probit
  info$cluster      <- cluster
  info$consistent   <- consistent
  info$ordered      <- ordered
  info$ordered.x    <- ordered.x
  info$ordered.y    <- ordered.y
  info$intTermElems <- intTermElems
  info$intTermNames <- intTermNames
  info$is.nlin      <- is.nlin
  info$rng.seed     <- floor(stats::runif(1L, min=0, max=9999999))
  info$n            <- NROW(preppedData$X)
  info$estimator    <- getEstimatorFromInfo(info)
  info$verbose      <- verbose

  info$mc.args <- list(
    min.iter        = mc.min.iter,
    max.iter        = mc.max.iter,
    mc.reps         = mc.reps,
    tol             = mc.tol,
    fixed.seed      = mc.fixed.seed,
    polyak.juditsky = mc.polyak.juditsky,
    fn.args         = mc.fn.args,
    rng.seed        = NULL
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
      finished       = FALSE,
      convergence    = FALSE,
      iterations     = 0L,
      iterations.0_5 = 0L,
      iterations.0_9 = 0L,
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
    se         = rep(NA_real_, k)
  )

  model
}


# badly named function -- since if also returns some useful info
initMatrices <- function(pt) {
  lvs.linear <- unique(pt[pt$op == "=~", "lhs"])
  etas <- unique(pt[pt$op == "~", "lhs"])

  lvs  <- c(lvs.linear, pt[grepl(":", pt$rhs), "rhs"])
  xis  <- lvs[!lvs %in% etas]

  allInds <- vector("character", 0)
  indsLvs <- vector("list", length(lvs))

  names(indsLvs) <- lvs
  for (lv in lvs) {
    indsLv <- pt[pt$lhs == lv & pt$op == "=~", "rhs"]
    allInds <- c(allInds, indsLv)
    indsLvs[[lv]] <- indsLv
  }

  inds.x <- unique(unlist(indsLvs[xis]))
  inds.y <- unique(unlist(indsLvs[etas]))

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

  # Selection Matric Indicators ------------------------------------------------
  Ip <- diag(nrow = nrow(lambda))
  colnames(Ip) <- rownames(Ip) <- rownames(lambda)

  nlinSelectFrom <- outer(lvs, lvs, Vectorize(f1))
  dimnames(nlinSelectFrom) <- list(lvs, lvs)


  C <- diag(nrow(gamma))
  dimnames(C) <- dimnames(gamma)

  # Create output lists --------------------------------------------------------
  matrices <- list(
    lambda = lambda,
    gamma = gamma,
    preds = preds,
    succs = succs,
    succs.linear = succs.linear,
    preds.linear = preds.linear,
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
    indsLvs    = indsLvs,
    allInds    = allInds,
    lvs        = lvs,
    lvs.linear = lvs.linear,
    etas       = etas,
    xis        = xis,
    inds.x     = inds.x,
    inds.y     = inds.y
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
  indsLvs <- model$info$indsLvs
  ptl     <- model$parTable.input
  SC      <- model$matrices$SC
  k       <- length(lvs.lin)

  # measurement model
  fitMeasurement <- lambda
  fitMeasurement[TRUE] <- 0
  for (lv in lvs) for (indsLv in indsLvs[[lv]])
    fitMeasurement[indsLv, lv] <- SC[indsLv, lv]

  # Caluculate consistent weights and correlations
  if (consistent) {
    # We want to correct both for the errors causes by using the CEXP
    # estimator, compared to the probit estimator. As well as the bias
    # caused by ignoring measurement error.

    if (consistent) Q <- getConstructQualities(model)
    else            Q <- stats::setNames(rep(1L, k), nm = lvs.lin) # ignore measurement error

    fitMeasurement <- getConsistentLoadings(model, Q = Q)
    model$matrices$C <- getConsistentCorrMat(model, Q = Q)
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

  # Indicator Residuals
  indicators <- rownames(lambda)
  k <- length(indicators)

  fitTheta <- matrix(0, nrow = k, ncol = k,
                     dimnames = list(indicators, indicators))
  crossLoaded <- apply(
    X      = fitMeasurement,
    MARGIN = 1L,
    FUN    = \(x) sum(abs(x) > .Machine$double.xmin) > 1L
  )

  warnif(any(crossLoaded),
         "Did not expect any cross loaded indicators,\n",
         "when calculating indicator residuals!")

  for (ind in indicators) {
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
    fitTheta       = plssemMatrix(fitTheta, symmetric = TRUE)
  )
}


getParamVecNames <- function(model) {
  selectLambda <- model$matrices$select$lambda
  lambda       <- selectLambda

  for (j in colnames(lambda)) for (i in rownames(lambda))
    lambda[i, j] <- paste0(j, "=~", i)

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
