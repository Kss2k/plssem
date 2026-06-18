OPERATORS <- c("<~", "~~", "=~", "~1", "~", "|")
MOPS      <- c("<~", "=~")
# Operator used to declare a structural variable that has no regression paths,
# e.g. `x ~ 1`. We never model the mean structure ourselves, so the `~1`
# operator is free to repurpose for this. It decouples "is a structural node"
# from the `~~` operator: `~~` now *only* frees a (co)variance among
# already-declared structural nodes, which (among other things) lets us specify
# covariances within the measurement model without accidentally promoting
# indicators to structural variables.
ADDITIONAL_STRUCT_VAR_OP <- "~1"
MAT_STRUCT <- matrix(0, nrow = 0, ncol = 0)
DF_STRUCT  <- data.frame(NULL)
VEC_STRUCT <- numeric(0)


specifyModel <- function(syntax, data, ..., strict = TRUE) {
  parTable <- modsem::modsemify(syntax, parentheses.as.string = TRUE)

  if (strict) {
    # check if we have any user-supplied names which use reserved patterns
    nm <- union(parTable$lhs, parTable$rhs)
    hasTmp <- hasTempAffixes(nm)

    pls_stopif(any(hasTmp),
      "Some variables have reserved keywords/patterns!",
      "Variables:", paste0(nm[hasTmp], collapse = ", ")
    )

    # check if we have user specified intercepts (which aren't allowed)
    hasIntr <- parTable$op == "~1"
    intrPar <- paste0(parTable$lhs[hasIntr], parTable$op[hasIntr])
    pls_stopif(any(hasIntr),
      "Estimation of intercepts is not available!",
      "Intercepts:", paste0(intrPar, collapse = ", ")
    )
  }

  hiOrd <- getHigherOrderLVs(parTable)
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
                            is.lower.order         = FALSE,
                            consistent             = TRUE,
                            missing                = "listwise",
                            standardize            = TRUE,
                            ordered                = NULL,
                            probit                 = NULL,
                            mcpls                  = NULL,
                            mc.fast.lmer           = NULL,
                            tolerance              = 1e-5,
                            max.iter.0_5           = 100,
                            mc.max.iter            = 250,
                            mc.min.iter            = 5,
                            mc.reps                = 20000,
                            mc.tol                 = 0.0005,
                            mc.fixed.seed          = FALSE,
                            mc.polyak.juditsky     = FALSE,
                            mc.pj.extrapolate      = TRUE,
                            mc.delta.se            = TRUE,
                            mc.delta.jacobian.k    = NULL,
                            mc.fn.args             = list(),
                            mc.rescov              = "auto",
                            verbose                = interactive(),
                            bootstrap              = FALSE,
                            boot.ncores            = 1L,
                            boot.parallel          = "no",
                            boot.R                 = 50L,
                            boot.iseed             = NULL,
                            boot.optimize          = FALSE,
                            mc.boot.control        = list(),
                            knn.k                  = 5,
                            reliabilities          = NULL,
                            default.path.estimator = "ols",
                            higherOrderLVs         = NULL) {
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

  pt         <- parsed$parTable.pls
  cluster    <- parsed$cluster
  consistent <- parsed$consistent
  ordered    <- parsed$ordered

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
    min.iter         = mc.min.iter,
    max.iter         = mc.max.iter,
    mc.reps          = mc.reps,
    tol              = mc.tol,
    fixed.seed       = mc.fixed.seed,
    polyak.juditsky  = mc.polyak.juditsky,
    pj.extrapolate   = mc.pj.extrapolate,
    delta.se         = mc.delta.se,
    delta.jacobian.k = mc.delta.jacobian.k,
    fn.args          = mc.fn.args,
    rescov           = match.arg(mc.rescov, c("auto", "reduced", "full"))
  )

  boot.info <- initModelBootInfo(
    bootstrap       = bootstrap,
    ncores          = boot.ncores,
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
    boot           = boot.info,
    scale          = preppedData$scale
  )

  thresholdStruct <- ThresholdStruct( # holds information about thresholds and
    data = preppedData$X,             # category proportions
    ordered = ordered
  )

  # GLS is used automatically whenever the model contains residual covariances
  # (which the OLS path estimator cannot handle). The user may also force GLS via
  # `default.path.estimator = "gls"`, even when OLS would otherwise suffice.
  gls.default <- tolower(default.path.estimator) == "gls"

  if (gls.default || hasResidualCovariances(parTable)) {
    info$path.estimator <- "gls"
    glsPathModel <- GlsPathModel(
      parTable = parTable,
      data.cov = NULL
    )
  } else {
    info$path.estimator <- "ols"
    glsPathModel <- GlsPathModel()
  }

  if (info$is.mlm && !info$is.mcpls) {
    pls_msg_note(
      "Multilevel/Mixed-Effects PLSc models are currently under development!\n",
      "Consider passing `mcpls=TRUE` to yield more consistent results."
    )
  }

  matrices$S  <- preppedData$S
  matrices$SC <- diagPartitioned(matrices$S, matrices$C)

  model <- PlsModel(
    parTableInput    = pt,
    matrices         = matrices,
    data             = preppedData$X,
    thresholdStruct  = thresholdStruct,
    glsPathModel     = glsPathModel,
    info             = info,
    status           = initModelStatus(
      tolerance      = tolerance,
      max.iter.0_5   = max.iter.0_5
    )
  )

  initModelParams(model)
}


# badly named function -- since it also returns some useful info
initMatrices <- function(pt, higherOrderLVs = NULL) {
  lvs.linear <- getLVs(pt)

  mode.a <- getReflectiveLVs(pt)
  mode.b <- getFormativeLVs(pt)
  mode.c <- intersect(mode.a, mode.b)

  # Should be handled by parseModelArguments() but we check
  # anyways...
  pls_stopif(length(mode.c),
    "Constructs must be either of mode A or mode B!",
    "Constructs with both mode A and B:",
    paste0(mode.c, collapse = ", ")
  )

  getmode <- function(x) {
    ifelse(x %in% mode.a,
      yes = "A",
      no = ifelse(x %in% mode.b,
        yes = "B",
        no = pls_msg_stop(
          sprintf("Unrecognized mode for variable `%s`!", x)
        )
      )
    )
  }

  modes <- stats::setNames(
    vapply(lvs.linear, FUN.VALUE = character(1L), FUN = getmode),
    nm = lvs.linear
  )

  checkLhsIntTerms(pt)
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

  # Covariances (lvs) ----------------------------------------------------------
  xis <- lvs[!lvs %in% pt[pt$op == "~", "lhs"]]
  selectCov <- matrix(
    FALSE, nrow = length(lvs), ncol = length(lvs),
    dimnames = list(lvs, lvs)
  )

  # defaults
  diag(selectCov) <- TRUE
  selectCov[xis, xis][lower.tri(selectCov[xis, xis])] <- TRUE

  # additional (residual) covariances
  rescov <- pt[
    pt$op == "~~" & pt$lhs != pt$rhs & pt$lhs %in% lvs & pt$rhs %in% lvs
    , , drop = FALSE
  ]

  for (i in seq_len(NROW(rescov))) {
    lhs <- rescov[i, "lhs"]
    rhs <- rescov[i, "rhs"]
    selectCov[lhs, rhs] <- selectCov[rhs, lhs] <- TRUE
  }

  # make sure upper diagonal is not selected
  selectCov[upper.tri(selectCov)] <- FALSE

  # Residual covariances (inds) ------------------------------------------------
  k           <- length(allInds)
  selectTheta <- matrix(
    FALSE, nrow = k, ncol = k,
    dimnames = list(allInds, allInds)
  )

  # always keep diagonal
  diag(selectTheta) <- TRUE

  # keep formative blocks
  for (b in mode.b) {
    idx <- indsLvs[[b]]
    selectTheta[idx, idx] <- TRUE
  }

  # only lower.tri
  selectTheta[upper.tri(selectTheta, diag = FALSE)] <- FALSE

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
