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
                            mc.fast.lmer       = NULL,
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
    mc.fast.lmer   = mc.fast.lmer,
    consistent     = consistent,
    is.lower.order = is.lower.order
  )

  pt           <- parsed$parTable.pls
  cluster      <- parsed$cluster
  consistent   <- parsed$consistent
  ordered      <- parsed$ordered

  matricesAndInfo <- initMatrices(pt, higherOrderLVs = higherOrderLVs)
  matrices        <- matricesAndInfo$matrices
  info            <- matricesAndInfo$info

  preppedData <- getPLS_Data(
    data        = parsed$data,
    indicators  = matricesAndInfo$info$allInds,
    consistent  = consistent,
    cluster     = cluster,
    standardize = standardize,
    ordered     = ordered,
    is.probit   = parsed$is.probit,
    missing     = missing,
    knn.k       = knn.k
  )

  n <- NROW(preppedData$X)

  mc.args <- initModelMcArgs(
    min.iter        = mc.min.iter,
    max.iter        = mc.max.iter,
    mc.reps         = mc.reps,
    tol             = mc.tol,
    fixed.seed      = mc.fixed.seed,
    polyak.juditsky = mc.polyak.juditsky,
    fn.args         = mc.fn.args
  )

  boot.info <- initModelBootInfo(
    bootstrap       = bootstrap,
    ncpus           = boot.ncpus,
    parallel        = boot.parallel,
    R               = boot.R,
    iseed           = boot.iseed,
    optimize        = boot.optimize,
    mc.boot.control = mc.boot.control
  )

  info <- initModelInfo(
    baseInfo       = info,
    parsed         = parsed,
    n              = n,
    ordered        = ordered,
    consistent     = consistent,
    verbose        = verbose,
    standardize    = standardize,
    reliabilities  = reliabilities,
    is.lower.order = is.lower.order,
    mc.args        = mc.args,
    boot           = boot.info
  )

  if (info$is.mlm && !info$is.mcpls) {
    message(
      "Multilevel/Mixed-Effects PLSc models are currently under development!\n",
      "Consider passing `mcpls=TRUE` to yield more consistent results."
    )
  }

  matrices$S  <- preppedData$S
  matrices$SC <- diagPartitioned(matrices$S, matrices$C)

  model <- PlsModel(
    parTableInput = pt,
    matrices      = matrices,
    data          = preppedData$X,
    info          = info,
    status        = initModelStatus(
      tolerance = tolerance,
      max.iter.0_5 = max.iter.0_5
    )
  )

  initModelParams(model)
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

