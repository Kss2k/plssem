specifyModel <- function(syntax, data) {
  pt <- modsem::modsemify(syntax)
  matricesAndInfo <- initMatrices(pt)
  matrices <- matricesAndInfo$matrices
  sortedData <- sortData(data, matricesAndInfo$info$allInds) 
  matrices$S <- cov(as.data.frame(sortedData))
  matrices$C <- diag(nrow(matrices$gamma))
  dimnames(matrices$C) <- dimnames(matrices$gamma)
  matrices$SC <- rbind(cbind(matrices$S, matrix(0, nrow = nrow(matrices$S), 
                                          ncol = nrow(matrices$C))),
                       cbind(matrix(0, nrow = nrow(matrices$C), 
                                   ncol = nrow(matrices$S)), matrices$C))
  colnames(matrices$SC) <- rownames(matrices$SC) <- c(colnames(matrices$S), 
                                                      colnames(matrices$C))
                

  model <- list(pt = pt, matrices = matrices, data = sortedData, 
                factorScores = NULL,
                info = matricesAndInfo$info,
                params = NULL, fit = NULL)
  model$params <- list(names = getParamVecNames(model), 
                       values = rep(NA, length(getParamVecNames(model))))
  model
}


# badly named function -- since if also returns some useful info
initMatrices <- function(pt) {
  lVs <- pt[pt$op == "=~", "lhs"] |> unique()
  # add interaction terms at the end
  lVs <- c(lVs, pt[grepl(":", pt$rhs), "rhs", drop = TRUE] |> unique())

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
  selectCovXis <- matrix(FALSE, nrow = length(xis), ncol = length(xis),
                   dimnames = list(xis, xis))
  selectCovXis[lower.tri(selectCovXis)] <- TRUE
  
  # Selection Matric Indicators ------------------------------------------------
  Ip <- diag(nrow = nrow(lambda))
  colnames(Ip) <- rownames(Ip) <- rownames(lambda)

  # Create output lists --------------------------------------------------------
  matrices <- list(lambda = lambda, gamma = gamma, 
                   preds = preds, succs = succs, 
                   outerWeights = getNonZeroElems(lambda), 
                   Ip = Ip, 
                   selectLambda = selectLambda, selectGamma = selectGamma,
                   selectCovXis = selectCovXis)
  info <- list(indsLvs = indsLvs, allInds = allInds, lVs = lVs,
               interactionPairs = interactionPairs)
  list(matrices = matrices, info = info)
}


getFit <- function(model, consistent = FALSE) {
  model$matrices$C <- model$matrices$C
   
  lambda <- model$matrices$lambda
  gamma <- model$matrices$gamma
  preds <- model$matrices$preds
  lVs <- model$info$lVs
  indsLvs <- model$info$indsLvs
  pt <- model$pt
  SC <- model$matrices$SC
  etas <- pt[pt$op == "~", "lhs"] |> unique()
  xis <- lVs[!lVs %in% etas]
  
  # measurement model 
  fitMeasurement <- lambda
  fitMeasurement[TRUE] <- 0
  for (lV in model$info$lVs) {
    for (indsLv in indsLvs[[lV]]) {
      fitMeasurement[indsLv, lV] <- SC[indsLv, lV]
    }
  }

  # Caluculate consistent weights, based on measurement model
  if (consistent) {
    P <- getReliabilityCoefs(model)
    fitMeasurement <- getConsistentLoadings(model, P = P)
    model$matrices$C <- getConsistenCorrMat(model, P = P)
  }
  
  # structural model
  fitStructural <- gamma
  fitStructural[TRUE] <- 0
  for (lV in model$info$lVs) {
    predsLv <- model$info$lVs[preds[ , lV, drop = TRUE]]
    if (length(predsLv) > 0) {
      fitStructural[predsLv, lV] <- getPathCoefs(lV, predsLv,  model$matrices$C)
    }
  }

  # Covariance matrix 
  fitCov <- model$matrices$C[xis, xis]
  
  list(fitMeasurement = fitMeasurement, fitStructural = fitStructural, fitCov = fitCov)
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
  lambda <- model$matrices$lambda 
  for (j in colnames(lambda)) {
    for (i in rownames(lambda)) {
      lambda[i, j] <- paste0(j, " =~ ", i)
    }
  }

  selectGamma <- model$matrices$selectGamma
  gamma <- model$matrices$gamma
  for (j in colnames(gamma)) {
    for (i in rownames(gamma)) {
      gamma[i, j] <- paste0(j, " ~ ", i)
    }
  }
  
  selectCov <- model$matrices$selectCovXis
  covXis <- selectCov
  for (j in colnames(covXis)) {
    for (i in rownames(covXis)) {
      covXis[i, j] <- paste0(j, " ~~ ", i)
    }
  }

  c(lambda[selectLambda], gamma[selectGamma], covXis[selectCov]) 
}


extractCoefs <- function(model) {
  fit <- model$fit 
  lambda <- fit$fitMeasurement
  selectLambda <- model$matrices$selectLambda
  gamma <- fit$fitStructural 
  selectGamma <- model$matrices$selectGamma
  fitCov <- fit$fitCov
  selectCov <- model$matrices$selectCovXis

  c(lambda[selectLambda], gamma[selectGamma], fitCov[selectCov])
}
