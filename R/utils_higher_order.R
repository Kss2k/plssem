splitHigherOrderModel <- function(syntax) {
  stopif(length(syntax) > 1L || !is.character(syntax),
         "`syntax` must be a string of length 1!")

  parTable <- modsem::modsemify(syntax, parentheses.as.string = TRUE)

  inds <- getIndicators(parTable, observed = FALSE)
  hiOrd <- getHigherOrderLVs(parTable)
  isHigherOrder <- length(hiOrd) > 0
  structOVs <- getStructOVs(parTable)

  if (isHigherOrder) {
    stopif(any(hiOrd %in% inds),
      "Higher order models with third order (or higher) constructs are not\n",
      "supported (yet)!"
    )

    parMsrO1 <- parTable[
      parTable$op %in% MOPS & !parTable$lhs %in% hiOrd, , drop = FALSE
    ]

    if (length(structOVs)) {
      # Make sure structOVs are passed on to the model parsing for the lower
      # order model, not just the higher order model. This also appends multilevel
      # cluster exprs. Since they are counted as structural variables, not
      # observed in the measurement model.
      parStrO1 <- data.frame(lhs = structOVs, op = "~~", rhs = structOVs, mod = "")

    } else {
      parStrO1 <- NULL

    }

    parMsrO2 <- parTable[
      parTable$op %in% MOPS & parTable$lhs %in% hiOrd, , drop = FALSE
    ]

    parStrO2 <- parTable[parTable$op == "~", , drop = FALSE]

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


combineModelResultsFirstSecondOrder <- function(model) {
  submodels   <- model$submodels
  firstOrder  <- submodels$firstOrder
  secondOrder <- submodels$secondOrder

  info1 <- firstOrder$info
  info2 <- secondOrder$info

  if (is.null(secondOrder)) {
    S        <- firstOrder$matrices$S
    C        <- firstOrder$matrices$C
    SC       <- firstOrder$matrices$SC
    fit      <- firstOrder$fit
    params   <- firstOrder$params
    parTable <- firstOrder$parTable
    status   <- firstOrder$status
    select   <- firstOrder$matrices$select

  } else {

    # Matrices -----------------------------------------------------------------
    S   <- firstOrder$matrices$S
    C   <- secondOrder$matrices$C
    SC1 <- firstOrder$matrices$SC
    SC2 <- secondOrder$matrices$SC
  
    cn1 <- colnames(SC1)
    cn2 <- colnames(SC2)[!grepl(TEMP_OV_PREFIX, colnames(SC2))]
    add <- setdiff(cn2, cn1)

    # The name order might be wrong here...
    if (length(add)) SC <- diagPartitioned(SC1, SC2[add, add, drop = FALSE])
    else             SC <- SC1

    SC[cn2, cn2] <- SC2[cn2, cn2]


    # Fit field ----------------------------------------------------------------
    f1 <- firstOrder$fit
    f2 <- secondOrder$fit

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
   
    T1 <- f1$fitTheta
    T2 <- f2$fitTheta
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

    # Matrices -> Select field -------------------------------------------------
    slc1 <- firstOrder$matrices$select
    slc2 <- secondOrder$matrices$select

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
      lambda = diagPartitioned(selectLambda1, selectLambda2) != 0, # diagPartitioned converts to numeric
      gamma  = slc2$gamma,
      cov    = slc2$cov,
      theta  = diagPartitioned(selectTheta1, selectTheta2) != 0
    )

    # params field -------------------------------------------------------------
    p1 <- firstOrder$params
    p2 <- secondOrder$params
    keep1 <- getKeepFromLowerOrderParNames(p1$names, lvs = info1$lvs)

    params <- list(
      names  = c(p1$names[keep1], p2$names),
      values = c(p1$values[keep1],p2$values),
      se     = c(p1$se[keep1], p2$se),
      vcov   = NULL
    )

    # parTable field -----------------------------------------------------------
    par1 <- firstOrder$parTable
    par2 <- secondOrder$parTable

    if (!is.null(par1)) {
      keep1 <- getKeepFromLowerOrderParTable(par1, lvs = info1$lvs)
      par1 <- par1[keep1, , drop=FALSE]
    }

    parTable <- rbind(par1, par2)

    # status field -------------------------------------------------------------
    s1 <- firstOrder$status
    s2 <- secondOrder$status

    status <- list(
      convergence    = s1$convergence && s2$convergence,
      iterations     = s1$iterations + s2$iterations,
      iterations.0_5 = s1$iterations + s2$iterations,
      tolerance      = c(s1$tolerance, s2$tolerance),
      max.iter.0_5   = c(s1$max.iter.0_5, s2$max.iter.0_5),
      is.admissible  = s1$is.admissible && s2$is.admissible
    )
  }

  list(
    submodels = list(
      firstOrder  = firstOrder,
      secondOrder = secondOrder
    ),
    matrices = list(
      S           = S,
      C           = C,
      SC          = SC,
      firstOrder  = firstOrder$matrices,
      secondOrder = secondOrder$matrices,
      select      = select
    ),
    data     = firstOrder$data, # Any updates in firstOrder data should be propagated
    info     = model$info,
    status   = status,
    parTable = parTable,
    params   = params,
    fit      = fit
  )
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
  # Redifine higher order measurement model as structural paths
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

  Lambda  <- secondOrder$fit$fitLambda
  Weights <- secondOrder$fit$fitWeights
  
  is.mode.b <- colnames(Lambda) %in% secondOrder$info$mode.b

  # Correct Loadings (Mode A)
  Q <- firstOrder$fit$Q
  Q[setdiff(rownames(Lambda), names(Q))] <- 1

  RelCorrection <- 1 / matrix(Q[rownames(Lambda)], nrow = NROW(Lambda),
                              ncol = NCOL(Lambda), byrow = FALSE)
  Lambda[,!is.mode.b] <- (Lambda * RelCorrection)[,!is.mode.b]


  # Correct Weights + Loadings (Mode B)
  for (lv in secondOrder$info$higherOrderComposites) {

    inds.lv <- secondOrder$info$indsLvs[[lv]]
   
    # Consistent weights from loadings
    lambda <- Lambda[inds.lv, lv]
    q      <- lambda/Q[inds.lv]
    S.ii   <- firstOrder$fit$fitC[inds.lv, inds.lv]
    v      <- solve(S.ii) %*% q

    w      <- c(v / sqrt(c(t(v) %*% S.ii %*% v))) # Standardize weights
    lambda <- c(S.ii %*% t(t(w)))                 # Consistent loadings

    Weights[inds.lv, lv] <- w 
    Lambda[inds.lv, lv]  <- lambda
  }

  # Collect in corrections in fitMeasurement
  Measurement <- Lambda

  Measurement[,is.mode.b]  <- Weights[,is.mode.b]
  Measurement[,!is.mode.b] <- Lambda[,!is.mode.b]

  # Update fit matrices
  secondOrder$fit$fitLambda      <- Lambda
  secondOrder$fit$fitWeights     <- Weights 
  secondOrder$fit$fitMeasurement <- Measurement

  # Update coefficients
  estimatePLS_Step8(secondOrder)  
}


getSecondOrderData <- function(firstOrder) {
  secOrdData <- as.data.frame(
    firstOrder$data %*% firstOrder$matrices$lambda
  )

  clusterVals <- attr(firstOrder$data, "cluster")
  clusterName <- firstOrder$info$cluster

  if (!is.null(clusterName) && !is.null(clusterVals))
    secOrdData[,clusterName] <- clusterVals

  secOrdData
}
