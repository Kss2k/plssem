OPERATORS <- c("~~", "=~", "~1", "~", "|", "<~")

specifyModel <- function(syntax,
                         data,
                         consistent = TRUE,
                         standardize = TRUE,
                         ordered = NULL,
                         probit = NULL,
                         consistent.probit = NULL,
                         mc.reps = 1e5,
                         tolerance = 1e-5,
                         max.iter.0_5 = 100,
                         max.iter.0_9 = 50) {
  parsed <- parseModelArguments(
    syntax = syntax,
    data = data,
    ordered = ordered,
    probit = probit
  )

  syntax               <- parsed$syntax
  data                 <- parsed$data
  pt                   <- parsed$parTable.pls
  cluster              <- parsed$cluster
  lme4.syntax          <- parsed$lme4.syntax
  intTermElems         <- parsed$intTermElems
  intTermNames         <- parsed$intTermNames
  is.nlin              <- parsed$is.nlin
  ordered              <- parsed$ordered
  is.probit            <- parsed$is.probit
  is.cexp              <- parsed$is.cexp

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
    is.probit   = is.probit,
    is.cexp     = is.cexp
  )

  inds.x    <- info$inds.x
  inds.y    <- info$inds.y
  ordered.x <- intersect(inds.x, ordered)
  ordered.y <- intersect(inds.y, ordered)
 
  info$lme4.syntax       <- lme4.syntax
  info$is.multilevel     <- !is.null(lme4.syntax)
  info$cluster           <- cluster
  info$consistent        <- consistent
  info$is.probit         <- is.probit 
  info$is.cexp           <- is.cexp
  info$ordered           <- ordered
  info$ordered.x         <- ordered.x
  info$ordered.y         <- ordered.y
  info$intTermElems      <- intTermElems
  info$intTermNames      <- intTermNames
  info$is.nlin           <- is.nlin
  info$mc.reps           <- mc.reps
  info$rng.seed          <- floor(stats::runif(1L) * 1e6)
  info$consistent.probit <- consistent.probit

  matrices$S <- preppedData$S
  matrices$C <- diag(nrow(matrices$gamma))

  dimnames(matrices$C) <- dimnames(matrices$gamma)
  matrices$SC <- rbind(
    cbind(matrices$S, matrix(0, nrow = nrow(matrices$S), ncol = nrow(matrices$C))),
    cbind(matrix(0, nrow = nrow(matrices$C), ncol = nrow(matrices$S)), matrices$C)
  )

  colnames(matrices$SC) <- rownames(matrices$SC) <- c(
    colnames(matrices$S), 
    colnames(matrices$C)
  )

  if (consistent.probit && info$is.cexp) {
    matrices$probit2cont <- getCorrMatsProbit2cont(
      data         = preppedData$X,
      selectLambda = matrices$selectLambda,
      ordered      = ordered,
      lvs          = info$lVs
    )
  }

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
      max.iter.0_9   = max.iter.0_9,
      max.iter.0_5   = max.iter.0_5
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
  lVs.linear <- unique(pt[pt$op == "=~", "lhs"])
  etas <- unique(pt[pt$op == "~", "lhs"])

  lVs  <- c(lVs.linear, pt[grepl(":", pt$rhs), "rhs"])
  xis  <- lVs[!lVs %in% etas]

  allInds <- vector("character", 0)
  indsLvs <- vector("list", length(lVs))

  names(indsLvs) <- lVs
  for (lV in lVs) {
    indsLv <- pt[pt$lhs == lV & pt$op == "=~", "rhs"]
    allInds <- c(allInds, indsLv)
    indsLvs[[lV]] <- indsLv
  }

  inds.x <- unique(unlist(indsLvs[xis]))
  inds.y <- unique(unlist(indsLvs[etas]))

  # Lamdba ---------------------------------------------------------------------
  lambda <- matrix(0, nrow = length(allInds), ncol = length(lVs),
                   dimnames = list(allInds, lVs))
  selectLambda <- matrix(FALSE, nrow = length(allInds), ncol = length(lVs),
                   dimnames = list(allInds, lVs))
  for (lV in lVs) {
    lambda[indsLvs[[lV]], lV] <- 1
    selectLambda[indsLvs[[lV]], lV] <- TRUE
  }

  # Gamma ----------------------------------------------------------------------
  gamma <- matrix(0, nrow = length(lVs), ncol = length(lVs),
                  dimnames = list(lVs, lVs))
  # Selection Matrix
  selectGamma <- matrix(FALSE, nrow = length(lVs), ncol = length(lVs),
                  dimnames = list(lVs, lVs))
  # Predecessors and successors
  preds <- succs <- matrix(FALSE, nrow = length(lVs), ncol = length(lVs),
                           dimnames = list(lVs, lVs))
  for (lV in lVs) {
    predsLv <- pt[pt$lhs == lV & pt$op == "~", "rhs"]
    succsLv <- pt[pt$rhs == lV & pt$op == "~", "lhs"]
    preds[predsLv, lV] <- TRUE
    succs[succsLv, lV] <- TRUE
    # selectionmatrix 
    selectGamma[predsLv, lV] <- TRUE
  }

  is.nlin      <- grepl(":", lVs)
  preds.linear <- preds
  succs.linear <- succs
  preds.linear[is.nlin,] <- FALSE
  preds.linear[,is.nlin] <- FALSE
  succs.linear[is.nlin,] <- FALSE
  succs.linear[,is.nlin] <- FALSE

  # Covariance Matrix xis ------------------------------------------------------
  xis <- lVs[!lVs %in% pt[pt$op == "~", "lhs"]]
  selectCov <- matrix(
    FALSE, nrow = length(lVs), ncol = length(lVs),
    dimnames = list(lVs, lVs)
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

  nlinSelectFrom <- outer(lVs, lVs, Vectorize(f1))
  dimnames(nlinSelectFrom) <- list(lVs, lVs)

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
    selectLambda = selectLambda,
    selectGamma = selectGamma,
    selectCov = selectCov,
    selectTheta = selectTheta,
    nlinSelectFrom = nlinSelectFrom,
    probit2cont = NULL
  )

  info <- list(
    indsLvs = indsLvs,
    allInds = allInds,
    lVs = lVs,
    lVs.linear = lVs.linear,
    etas = etas,
    xis = xis,
    inds.x = inds.x,
    inds.y = inds.y
  )

  list(
    matrices = matrices,
    info = info
  )
}


getFitPLSModel <- function(model, consistent = TRUE) {
  model$matrices$C <- model$matrices$C
   
  lambda  <- model$matrices$lambda
  gamma   <- model$matrices$gamma
  preds   <- model$matrices$preds
  etas    <- model$info$etas
  xis     <- model$info$xis
  lVs     <- model$info$lVs
  lVs.lin <- model$info$lVs.linear
  indsLvs <- model$info$indsLvs
  is.cexp <- model$info$is.cexp
  ptl     <- model$parTable.input
  SC      <- model$matrices$SC
  k       <- length(lVs.lin)
  
  # measurement model 
  fitMeasurement <- lambda
  fitMeasurement[TRUE] <- 0
  for (lV in lVs) for (indsLv in indsLvs[[lV]])
    fitMeasurement[indsLv, lV] <- SC[indsLv, lV]

  # Caluculate consistent weights and correlations
  if (consistent || is.cexp) {
    # We want to correct both for the errors causes by using the CEXP
    # estimator, compared to the probit estimator. As well as the bias
    # caused by ignoring measurement error.

    if (consistent) Q <- getConstructQualities(model)
    else            Q <- stats::setNames(rep(1L, k), nm = lVs.lin) # ignore measurement error

    fitMeasurement <- getConsistentLoadings(model, Q = Q)
    model$matrices$C <- getConsistentCorrMat(model, Q = Q)
  }
  
  # structural model
  fitStructural       <- gamma
  fitStructural[TRUE] <- 0
  for (lV in lVs) {
    predsLv <- lVs[preds[ , lV, drop = TRUE]]

    if (length(predsLv))
      fitStructural[predsLv, lV] <- getPathCoefs(lV, predsLv,  model$matrices$C)
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
    r  <- fitMeasurement[ind, which.max(fitMeasurement[ind, ])]
    v  <- SC[ind, ind]
    fitTheta[ind, ind] <- v - r^2
  }

  list(
    fitMeasurement = plssemMatrix(fitMeasurement),
    fitStructural  = plssemMatrix(fitStructural),
    fitCov         = plssemMatrix(fitCov, symmetric = TRUE),
    fitTheta       = plssemMatrix(fitTheta, symmetric = TRUE)
  )
}


getParamVecNames <- function(model) {
  selectLambda <- model$matrices$selectLambda
  lambda       <- selectLambda

  for (j in colnames(lambda)) for (i in rownames(lambda))
    lambda[i, j] <- paste0(j, "=~", i)

  selectGamma <- model$matrices$selectGamma
  gamma       <- selectGamma

  for (j in colnames(gamma)) for (i in rownames(gamma))
    gamma[i, j] <- paste0(j, "~", i)
  
  selectCov <- model$matrices$selectCov
  psi       <- selectCov

  for (j in colnames(psi)) for (i in rownames(psi))
    psi[i, j] <- paste0(j, "~~", i)

  selectTheta <- model$matrices$selectTheta
  theta       <- selectTheta
  
  for (j in colnames(theta)) for (i in rownames(theta))
    theta[i, j] <- paste0(j, "~~", i)

  c(lambda[selectLambda], gamma[selectGamma], psi[selectCov], theta[selectTheta]) 
}


extractCoefs <- function(model) {
  fit <- model$fit 

  lambda       <- fit$fitMeasurement
  selectLambda <- model$matrices$selectLambda

  gamma       <- fit$fitStructural 
  selectGamma <- model$matrices$selectGamma

  fitCov    <- fit$fitCov
  selectCov <- model$matrices$selectCov

  fitTheta    <- fit$fitTheta
  selectTheta <- model$matrices$selectTheta

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

  scale(X %*% W)
}
