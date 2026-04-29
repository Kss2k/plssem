splitHigherOrderParTable <- function(parTable) {
  hiOrd     <- getHigherOrderLVs(parTable)
  isHiOrd   <- length(hiOrd) > 0
  structOVs <- getStructOVs(parTable)

  if (isHiOrd) {
    parMsrO1 <- parTable[
      parTable$op %in% MOPS & !parTable$lhs %in% hiOrd, , drop = FALSE
    ]

    if (length(structOVs)) {
      parStrO1 <- data.frame(lhs = structOVs, op = "~~", rhs = structOVs, mod = "")
    } else {
      parStrO1 <- NULL
    }

    parMsrO2 <- parTable[
      parTable$op %in% MOPS & parTable$lhs %in% hiOrd, , drop = FALSE
    ]
    parStrO2  <- parTable[parTable$op == "~", , drop = FALSE]

    parTableO2 <- rbind(parMsrO2, parStrO2)
    parTableO1 <- rbind(parMsrO1, parStrO1)

  } else {
    parTableO2 <- NULL
    parTableO1 <- parTable
  }

  list(
    parTableO2     = parTableO2,
    parTableO1     = parTableO1,
    higherOrderLVs = hiOrd
  )
}


splitHigherOrderModel <- function(syntax) {
  stopif(length(syntax) > 1L || !is.character(syntax),
         "`syntax` must be a string of length 1!")

  parTable <- modsem::modsemify(syntax, parentheses.as.string = TRUE)
  splitHigherOrderParTable(parTable)
}


combineModelResultsFirstSecondOrder <- function(model) {
  fo <- model
  so <- higherOrderModel(model)

  if (is.null(so))
    return(model)

  info1 <- fo@info
  info2 <- so@info

  # Higher-order model - combine first-order and (possibly recursive)
  # higher-order results.

  S   <- fo@matrices$S
  C   <- so@matrices$C
  SC1 <- fo@matrices$SC
  SC2 <- so@matrices$SC

  cn1 <- colnames(SC1)
  cn2 <- colnames(SC2)[!grepl(TEMP_OV_PREFIX, colnames(SC2))]
  add <- setdiff(cn2, cn1)

  if (length(add)) SC <- diagPartitioned(SC1, SC2[add, add, drop = FALSE])
  else             SC <- SC1
  SC[cn2, cn2] <- SC2[cn2, cn2]

  # Fit field
  f1 <- fo@fit
  f2 <- so@fit

  L1 <- f1$fitLambda
  L2 <- f2$fitLambda
  L2 <- L2[!grepl(TEMP_OV_PREFIX, rownames(L2)),
            !colnames(L2) %in% colnames(L1), drop = FALSE]

  W1 <- f1$fitWeights
  W2 <- f2$fitWeights
  W2 <- W2[!grepl(TEMP_OV_PREFIX, rownames(W2)),
            !colnames(W2) %in% colnames(W1), drop = FALSE]

  M1 <- f1$fitMeasurement
  M2 <- f2$fitMeasurement
  M2 <- M2[!grepl(TEMP_OV_PREFIX, rownames(M2)),
            !colnames(M2) %in% colnames(M1), drop = FALSE]

  T1   <- f1$fitTheta
  T2   <- f2$fitTheta
  keep <- !grepl(TEMP_OV_PREFIX, colnames(T2))
  T2   <- T2[keep, keep, drop = FALSE]

  fit <- list(
    fitMeasurement = diagPartitioned(M1, M2),
    fitLambda      = diagPartitioned(L1, L2),
    fitWeights     = diagPartitioned(W1, W2),
    fitTheta       = diagPartitioned(T1, T2),
    fitC           = f2$fitC,
    fitCov         = f2$fitCov,
    fitStructural  = f2$fitStructural,
    Q              = f2$Q
  )

  # Select matrices
  slc1 <- fo@matrices$select
  slc2 <- so@matrices$select

  selectLambda1 <- slc1$lambda
  selectLambda2 <- slc2$lambda
  selectLambda2 <- selectLambda2[
    !grepl(TEMP_OV_PREFIX, rownames(selectLambda2)),
    !colnames(selectLambda2) %in% colnames(selectLambda1), drop = FALSE
  ]

  selectTheta1 <- slc1$theta
  selectTheta2 <- slc2$theta
  keep         <- !grepl(TEMP_OV_PREFIX, colnames(selectTheta2))
  selectTheta2 <- selectTheta2[keep, keep, drop = FALSE]

  select <- list(
    lambda = diagPartitioned(selectLambda1, selectLambda2) != 0,
    gamma  = slc2$gamma,
    cov    = slc2$cov,
    theta  = diagPartitioned(selectTheta1, selectTheta2) != 0
  )

  # status
  s1 <- fo@status
  s2 <- so@status
  status <- list(
    convergence    = s1$convergence && s2$convergence,
    iterations     = s1$iterations + s2$iterations,
    iterations.0_5 = s1$iterations + s2$iterations,
    tolerance      = c(s1$tolerance, s2$tolerance),
    max.iter.0_5   = c(s1$max.iter.0_5, s2$max.iter.0_5),
    is.admissible  = s1$is.admissible && s2$is.admissible,
    mcpls.update.args = NULL
  )

  model@matrices <- list(
    S           = S,
    C           = C,
    SC          = SC,
    firstOrder  = fo@matrices,
    secondOrder = so@matrices,
    select      = select
  )

  model@data     <- fo@data
  model@status   <- status
  model@fit      <- fit

  # Update parameters
  refreshModelParams(model, update.names = TRUE)
}


computeCombinedModel <- function(model, lowerOrderAsEta = FALSE) {
  stopif(!is(model, "PlsModel"), "Expected a PlsModel")

  if (is.null(higherOrderModel(model)))
    return(model)

  so <- computeCombinedModel(higherOrderModel(model))
  fo <- model

  combined <- fo
  combined@higherOrderModel <- NULL
  combined@combinedModel <- NULL
  combined@boot <- list()
  combined@parTable <- NULL

  higherOrderModel(combined) <- so
  combined <- combineModelResultsFirstSecondOrder(combined)

  # Combined model is a terminal representation.
  combined@higherOrderModel <- NULL

  # Merge model metadata so downstream methods see full structure.
  info1 <- fo@info
  info2 <- so@info

  lvs     <- info1$lvs
  etas    <- info2$etas
  ovInds  <- info1$allInds

  allInds2 <- info2$allInds
  if (is.null(allInds2))
    allInds2 <- unique(unlist(info2$indsLvs))

  lvInds  <- intersect(allInds2, lvs)
  indsMap <- c(info1$indsLvs, info2$indsLvs)

  lowInds <- unique(unlist(indsMap[lvInds]))
  etaInds <- unique(unlist(indsMap[etas]))

  etaIndsOv <- intersect(etaInds, ovInds)
  etaIndsLv <- intersect(etaInds, lvs)

  # Technically lower order (latent) variables/indicators are endogenous variables
  # in the model, but this is seldom how people think of these variables. From
  # a PLS sentered standpoint they are predictor variables, not dependent variables.
  # For now we allow lower order variables to be treated as exogenous variables,
  # if their parent is exogenous.
  if (lowerOrderAsEta) yInds <- c(lowInds, etaInds)
  else yInds <- unique(c(etaIndsOv, unlist(indsMap[etaIndsLv])))

  inds.y <- intersect(ovInds, yInds)
  inds.x <- setdiff(ovInds, yInds)

  etas.all <- union(info1$etas, info2$etas)
  xis.all  <- setdiff(union(info1$xis, info2$xis), etas.all)
  ordered  <- union(info1$ordered, info2$ordered)
  ordered.base <- info1$ordered

  # MC-PLSc is an estimator choice, not merely the presence of ordered
  # indicators. A combined model should only be flagged MC-PLSc if any level
  # was actually fitted with MC-PLSc.
  is.mcpls <- isTRUE(info1$is.mcpls) || isTRUE(info2$is.mcpls)

  is.probit <- (
    (isTRUE(info1$is.probit) || isTRUE(info2$is.probit)) && !is.mcpls
  )

  combined@info <- list(
    lvs.linear   = union(info1$lvs.linear, info2$lvs.linear),
    lvs          = union(info1$lvs, info2$lvs),
    lvs.hi.ord   = union(info1$higherOrderLVs, info2$higherOrderLVs),
    allInds      = union(info1$allInds, allInds2),
    xis          = xis.all,
    etas         = etas.all,
    mode.a       = union(info1$mode.a, info2$mode.a),
    mode.b       = union(info1$mode.b, info2$mode.b),
    modes        = namedListUnion(info1$modes, info2$modes),
    inds.x       = inds.x,
    inds.y       = inds.y,
    indsLvs      = namedListUnion(info1$indsLvs, info2$indsLvs),
    cluster      = info1$cluster,
    ordered      = ordered.base,
    is.mlm       = isTRUE(info1$is.mlm) || isTRUE(info2$is.mlm),
    is.mcpls     = is.mcpls,
    is.probit    = is.probit,
    is.cfa       = isTRUE(info1$is.cfa) && (is.null(info2$is.cfa) || info2$is.cfa),
    is.high.ord  = TRUE,
    n            = info1$n,
    estimator    = info1$estimator,
    standardized = info1$standardized,
    verbose      = isTRUE(info1$verbose) || isTRUE(info2$verbose),
    mc.args      = info1$mc.args,
    boot         = info1$boot,
    # Preserve additional fields used elsewhere.
    ordered.x     = intersect(inds.x, ordered.base),
    ordered.y     = intersect(inds.y, ordered.base),
    intTermElems  = namedListUnion(info1$intTermElems, info2$intTermElems),
    intTermNames  = union(info1$intTermNames, info2$intTermNames),
    is.nlin       = isTRUE(info1$is.nlin) || isTRUE(info2$is.nlin),
    lme4.syntax   = info1$lme4.syntax,
    consistent    = info1$consistent,
    reliabilities = info1$reliabilities,
    rng.seed      = info1$rng.seed,
    is.lower.order = FALSE
  )

  # Ensure downstream methods see combined state.
  combined@info$is.high.ord <- TRUE

  # Recompute parameter names using combined metadata.
  combined <- refreshModelParams(combined, update.names = TRUE)
  combined
}


getKeepFromLowerOrderParNames <- function(parnames, lvs) {
  split <- getParTableFromParNames(parnames)
  lhs   <- split$lhs
  op    <- split$op
  rhs   <- split$rhs
  !(lhs %in% lvs & rhs %in% lvs & op == "~~")
}


getKeepFromLowerOrderParTable <- function(parTable, lvs) {
  lhs <- parTable$lhs
  op  <- parTable$op
  rhs <- parTable$rhs
  !(lhs %in% lvs & rhs %in% lvs & op == "~~")
}


highOrdMeasrAsStructParTable <- function(parTable) {
  higherOrder <- getHigherOrderLVs(parTable)

  if (!length(higherOrder))
    return(parTable)

  lvs <- getLVs(parTable)
  lhs <- parTable$lhs
  op  <- parTable$op
  rhs <- parTable$rhs

  isHiOrdMsr <- lhs %in% higherOrder & op %in% MOPS & rhs %in% lvs

  parTable[isHiOrdMsr, "lhs"] <- rhs[isHiOrdMsr]
  parTable[isHiOrdMsr, "op"]  <- "~"
  parTable[isHiOrdMsr, "rhs"] <- lhs[isHiOrdMsr]

  parTable
}


correctLoadingsAndWeightsSecondOrder <- function(firstOrder, secondOrder) {
  if (is.null(secondOrder))
    return(NULL)

  Lambda  <- secondOrder@fit$fitLambda
  Weights <- secondOrder@fit$fitWeights

  is.mode.b <- colnames(Lambda) %in% secondOrder@info$mode.b

  Q <- firstOrder@fit$Q
  Q[setdiff(rownames(Lambda), names(Q))] <- 1

  RelCorrection <- 1 / matrix(Q[removeTempAffixes(rownames(Lambda))],
                              nrow = NROW(Lambda), ncol = NCOL(Lambda),
                              byrow = FALSE)
  Lambda[, !is.mode.b] <- (Lambda * RelCorrection)[, !is.mode.b]

  for (lv in secondOrder@info$higherOrderComposites) {
    indsLv2nd <- secondOrder@info$indsLvs[[lv]] # Second order names
    indsLv1st <- removeTempAffixes(indsLv2nd)   # First order (lv) names

    lambda  <- Lambda[indsLv2nd, lv]
    q       <- lambda / Q[indsLv1st]
    S.ii    <- firstOrder@fit$fitC[indsLv1st, indsLv1st]
    v       <- solve(S.ii) %*% q
    w       <- c(v / sqrt(c(t(v) %*% S.ii %*% v)))
    lambda  <- c(S.ii %*% t(t(w)))

    Weights[indsLv2nd, lv] <- w
    Lambda[indsLv2nd, lv]  <- lambda
  }

  Measurement                  <- Lambda
  Measurement[, is.mode.b]     <- Weights[, is.mode.b]
  Measurement[, !is.mode.b]    <- Lambda[, !is.mode.b]

  secondOrder@fit$fitLambda      <- Lambda
  secondOrder@fit$fitWeights     <- Weights
  secondOrder@fit$fitMeasurement <- Measurement

  refreshModelParams(secondOrder)
}


getSecondOrderInputData <- function(firstOrder) {
  if (is.null(firstOrder))
    return(data.frame())

  secOrdData  <- as.data.frame(firstOrder@data %*% firstOrder@matrices$lambda)
  clusterVals <- attr(firstOrder@data, "cluster")
  clusterName <- firstOrder@info$cluster

  if (!is.null(clusterName) && !is.null(clusterVals))
    secOrdData[, clusterName] <- clusterVals

  secOrdData
}


getSecondOrderDataMatrix <- function(firstOrder, secondOrder) {
  olddata <- modelData(secondOrder)
  Scores <- computeFactorScores(firstOrder)

  want <- colnames(olddata)
  have <- colnames(Scores)

  # Variables with TEMP_OV_PREFIX
  isTmpOvCol <- hasTempOvPrefix(want)

  if (any(isTmpOvCol)) {
    tmpCols <- want[isTmpOvCol]
    clnCols <- removeTempOvPrefix(tmpCols)

    checkMissingConstructScores(have = have, want = clnCols)

    TmpOv <- Scores[,clnCols, drop = FALSE]
    colnames(TmpOv) <- tmpCols

    Scores <- cbind(Scores, TmpOv)
  }

  # Variables with TEMP_IND_SUFFIX
  isTmpIndCol <- hasTempIndSuffix(want)

  if (any(isTmpIndCol)) {
    tmpCols <- want[isTmpIndCol]
    clnCols <- removeTempIndSuffix(tmpCols)

    checkMissingConstructScores(have = have, want = clnCols)

    TmpInd <- Scores[,clnCols, drop = FALSE]
    colnames(TmpInd) <- tmpCols

    Scores <- cbind(Scores, TmpInd)
  }

  # Finalize
  checkMissingConstructScores(have = colnames(Scores), want = want)

  newdata <- Scores[, want, drop = FALSE]
  attr(newdata, "cluster") <- attr(modelData(firstOrder), "cluster")

  newdata
}


checkMissingConstructScores <- function(have, want) {
  if (!all(want %in% have)) {
    stop2("Missing construct scores for: ",
          paste0(setdiff(want, have), collapse = ", "))
  }
}


getRepeatedIndicatorWeights <- function(model) {
  if (!hasHigherOrderModel(model)) {
    fit <- modelFit(model)
    return(fit$fitWeights)
  }

  W1 <- modelFit(model)$fitWeights
  W2 <- getRepeatedIndicatorWeights(higherOrderModel(model))
  W2 <- plssemMatrix(W2, is.public = TRUE)

  want <- colnames(W1)
  have <- rownames(W2)
  both <- intersect(want, have)

  W3 <- W1[,both, drop=FALSE] %*% W2[both, , drop=FALSE]
  cbind(W1, W3[,setdiff(colnames(W3), colnames(W1)), drop=FALSE])
}
