splitHigherOrderModel <- function(syntax) {
  stopif(length(syntax) > 1L || !is.character(syntax),
         "`syntax` must be a string of length 1!")

  parTable  <- modsem::modsemify(syntax, parentheses.as.string = TRUE)
  inds      <- getIndicators(parTable, observed = FALSE)
  hiOrd     <- getHigherOrderLVs(parTable)
  isHiOrd   <- length(hiOrd) > 0
  structOVs <- getStructOVs(parTable)

  if (isHiOrd) {
    stopif(any(hiOrd %in% inds),
      "Higher order models with third order (or higher) constructs are not\n",
      "supported (yet)!"
    )

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


# Combine first- and second-order sub-model results into the parent Plssem
# object.  Returns the updated Plssem object (not a plain list).
combineModelResultsFirstSecondOrder <- function(model) {
  fo <- firstOrder(model)
  so <- secondOrder(model)

  if (is.null(so)) {
    # No higher-order model — promote first-order results directly.
    model@matrices <- c(
      list(
        S           = fo@matrices$S,
        C           = fo@matrices$C,
        SC          = fo@matrices$SC,
        firstOrder  = fo@matrices,
        secondOrder = NULL,
        select      = fo@matrices$select
      )
    )
    model@data     <- fo@data
    model@status   <- fo@status
    model@params   <- fo@params
    model@fit      <- fo@fit

  } else {
    info1 <- fo@info
    info2 <- so@info
    # Higher-order model — combine first-order and second-order results.

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
      is.admissible  = s1$is.admissible && s2$is.admissible
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
    model <- refreshModelParams(model, update.names = TRUE)
  }

  model
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

  RelCorrection <- 1 / matrix(Q[rownames(Lambda)], nrow = NROW(Lambda),
                              ncol = NCOL(Lambda), byrow = FALSE)
  Lambda[, !is.mode.b] <- (Lambda * RelCorrection)[, !is.mode.b]

  for (lv in secondOrder@info$higherOrderComposites) {
    inds.lv <- secondOrder@info$indsLvs[[lv]]
    lambda  <- Lambda[inds.lv, lv]
    q       <- lambda / Q[inds.lv]
    S.ii    <- firstOrder@fit$fitC[inds.lv, inds.lv]
    v       <- solve(S.ii) %*% q
    w       <- c(v / sqrt(c(t(v) %*% S.ii %*% v)))
    lambda  <- c(S.ii %*% t(t(w)))

    Weights[inds.lv, lv] <- w
    Lambda[inds.lv, lv]  <- lambda
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
  scores <- computeFactorScores(firstOrder)

  want <- colnames(olddata)
  have <- colnames(scores)

  have[!have %in% want] <- paste0(TEMP_OV_PREFIX, have[!have %in% want])
  colnames(scores) <- have

  stopif(!all(have %in% want), "Missing construct scores for: ",
         paste0(setdiff(want, have), collapse = ", "))

  newdata <- scores[, want]
  attr(newdata, "cluster") <- attr(olddata, "cluster")

  newdata
}
