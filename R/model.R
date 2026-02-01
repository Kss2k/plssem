OPERATORS <- c("=~", "~~", "~", "|", "<~")

specifyModel <- function(syntax,
                         data,
                         consistent = TRUE,
                         standardize = TRUE,
                         ordered = NULL,
                         probit = NULL) {
  parsed      <- parseModelArguments(syntax = syntax, data = data)
  syntax      <- parsed$syntax
  data        <- parsed$data
  pt          <- parsed$parTable.pls
  cluster     <- parsed$cluster
  lme4.syntax <- parsed$lme4.syntax

  matricesAndInfo <- initMatrices(pt)
  matrices        <- matricesAndInfo$matrices
  info            <- matricesAndInfo$info

  preppedData <- prepData(
    data        = data,
    indicators  = matricesAndInfo$info$allInds,
    consistent  = consistent,
    cluster     = cluster,
    standardize = standardize,
    ordered     = ordered,
    probit      = probit
  ) 
  
  info$lme4.syntax   <- lme4.syntax
  info$is.multilevel <- !is.null(lme4.syntax)
  info$cluster       <- cluster
  info$consistent    <- consistent
  info$is.probit     <- preppedData$probit 
  info$ordered       <- preppedData$ordered

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

  model <- list(
    parTable.input = pt,
    parTable       = NULL,
    matrices       = matrices,
    data           = preppedData$X,
    factorScores   = NULL,
    info           = info,
    params         = NULL,
    fit            = NULL
  )

  parnames <- getParamVecNames(model)
  k        <- length(parnames)

  model$params <- list(
    names  = parnames, 
    values = rep(NA_real_, k),
    se     = rep(NA_real_, k)
  )

  model
}


# badly named function -- since if also returns some useful info
initMatrices <- function(pt) {
  lVs <- unique(pt[pt$op == "=~", "lhs"])
  # add interaction terms at the end
  lVs <- c(lVs, pt[grepl(":", pt$rhs), "rhs", drop = TRUE] |> unique())
  etas <- unique(pt[pt$op == "~", "lhs"])
  xis  <- lVs[!lVs %in% etas]

  allInds <- vector("character", 0)
  indsLvs <- vector("list", length(lVs))
  names(indsLvs) <- lVs
  for (lV in lVs) {
    indsLv <- pt[pt$lhs == lV & pt$op == "=~", "rhs"]
    allInds <- c(allInds, indsLv)
    indsLvs[[lV]] <- indsLv
  }

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

  # Interaction effects 
  interactionPairs <- getInteractionPairs(pt)

  # Covariance Matrix xis ------------------------------------------------------
  xis <- lVs[!lVs %in% pt[pt$op == "~", "lhs"]]
  selectCov <- matrix(
    FALSE, nrow = length(lVs), ncol = length(lVs),
    dimnames = list(lVs, lVs)
  )

  diag(selectCov) <- TRUE
  selectCov[xis, xis][lower.tri(selectCov[xis, xis])] <- TRUE
  
  # Selection Matric Indicators ------------------------------------------------
  Ip <- diag(nrow = nrow(lambda))
  colnames(Ip) <- rownames(Ip) <- rownames(lambda)

  # Create output lists --------------------------------------------------------
  matrices <- list(lambda = lambda, gamma = gamma, 
                   preds = preds, succs = succs, 
                   outerWeights = getNonZeroElems(lambda), 
                   Ip = Ip, 
                   selectLambda = selectLambda, selectGamma = selectGamma,
                   selectCov = selectCov)
  info <- list(indsLvs = indsLvs, allInds = allInds, lVs = lVs,
               interactionPairs = interactionPairs, etas = etas, xis = xis)
  list(matrices = matrices, info = info)
}


getFitPLSModel <- function(model, consistent = TRUE) {
  model$matrices$C <- model$matrices$C
   
  lambda  <- model$matrices$lambda
  gamma   <- model$matrices$gamma
  preds   <- model$matrices$preds
  etas    <- model$info$etas
  xis     <- model$info$xis
  lVs     <- model$info$lVs
  indsLvs <- model$info$indsLvs
  ptl     <- model$parTable.input
  SC      <- model$matrices$SC
  
  # measurement model 
  fitMeasurement <- lambda
  fitMeasurement[TRUE] <- 0
  for (lV in lVs) for (indsLv in indsLvs[[lV]])
    fitMeasurement[indsLv, lV] <- SC[indsLv, lV]

  # Caluculate consistent weights, based on measurement model
  if (consistent) {
    P <- getReliabilityCoefs(model)
    fitMeasurement <- getConsistentLoadings(model, P = P)
    model$matrices$C <- getConsistenCorrMat(model, P = P)
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
  fitCovRes <- diag2(t(fitStructural) %*% model$matrices$C %*% fitStructural)
  fitCov[etas, etas] <- fitCovRes[etas, etas]
  
  list(
    fitMeasurement = plssemMatrix(fitMeasurement),
    fitStructural  = plssemMatrix(fitStructural),
    fitCov         = plssemMatrix(fitCov, symmetric = TRUE)
  )
}


getInteractionPairs <- function(pt) {
  lVs <- pt[pt$op == "=~", "lhs"] |> unique()
  allInds <- vector("character", 0)
  indsLvs <- vector("list", length(lVs))
  names(indsLvs) <- lVs
  for (lV in lVs) {
    indsLv <- pt[pt$lhs == lV & pt$op == "=~", "rhs"]
    allInds <- c(allInds, indsLv)
    indsLvs[[lV]] <- indsLv
  }
  
  intTerms <- pt[grepl(":", pt$rhs), "rhs", drop = TRUE] |> unique()
  interactionPairs <- vector("list", length(intTerms))
  names(interactionPairs) <- intTerms 
  for (int in intTerms) {
    interactionPairs[[int]] <- 
      as.vector(stringr::str_split_fixed(int, ":", n = 2))
  }
  interactionPairs
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

  c(lambda[selectLambda], gamma[selectGamma], psi[selectCov]) 
}


extractCoefs <- function(model) {
  fit <- model$fit 

  lambda       <- fit$fitMeasurement
  selectLambda <- model$matrices$selectLambda

  gamma       <- fit$fitStructural 
  selectGamma <- model$matrices$selectGamma

  fitCov    <- fit$fitCov
  selectCov <- model$matrices$selectCov

  out <- c(lambda[selectLambda], gamma[selectGamma], fitCov[selectCov])
  names(out) <- model$params$names

  plssemVector(out)
}


getFactorScores <- function(model) {
  matrices <- model$matrices
  W <- model$matrices$lambda
  X <- model$data

  X %*% W
}
