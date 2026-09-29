mpls <- function(syntax,
                 data,
                 cluster,
                 diag.secant = FALSE,
                 polyak.juditsky = TRUE,
                 pj.extrapolate = TRUE,
                 tol = if (polyak.juditsky) 0.0001 else 0.001,
                 fn.args = list(),
                 max.iter = 1000L,
                 min.iter = 50L,
                 verbose = interactive(),
                 ...,
                 consistent = NULL # capture
                 ) {
  pls_stopif(length(cluster) != 1 || !is.character(cluster),
    "cluster must be a character string of length 1!"
  )

  parsed <- parseMultilevelModelArguments(syntax, data)
  
  data <- as.data.frame(data)
  all.vars <- c(cluster, parsed$ovs.all)
  missing <- setdiff(all.vars, colnames(data))
  pls_stopif(length(missing),
    "Missing variables in `data`:", paste0(missing, ", ")
  )

  parTableL1 <- parsed$level.1 
  parTableL2 <- parsed$level.2

  data[[cluster]] <- as.integer(as.factor(data[[cluster]]))
  clusterIdx <- data[[cluster]]

  if (anyNA(clusterIdx)) {
    pls_stopif(all(is.na(clusterIdx), "cluster is all NA!"))
    pls_msg_warn("removing missing values in `cluster`!")

    data <- data[!is.na(clusterIdx),, drop = FALSE]
    clusterIdx <- data[[cluster]]
  }

  # data must be sorted by the clusters
  data <- as.matrix(data[order(clusterIdx), all.vars,drop=FALSE])
  data[,parsed$ovs.all] <- Rfast::standardise(data[,parsed$ovs.all])
  clusterIdx <- data[,cluster, drop=TRUE]
  n <- NROW(data)

  # fit auxiliary models
  baseFits <- fitAuxiliaryMLM_PLS(
    parsed     = parsed,
    data       = data,
    clusterIdx = clusterIdx,
    ...
  )

  # for now
  is.hi.ord.l1 <- isTRUE(combinedModel(baseFits$level.1)@info$is.high.ord)
  is.hi.ord.l2 <- isTRUE(combinedModel(baseFits$level.2)@info$is.high.ord)
  use.full.rescov.l1 <- combinedModel(baseFits$level.1)@info$path.estimator == "gls"
  use.full.rescov.l2 <- combinedModel(baseFits$level.2)@info$path.estimator == "gls"
  mc.reps <- 20000
  rng.seed <- NULL

  # calibrate mc.reps
  clusterSizes <- table(clusterIdx)
  mc.reps.l1 <- max(mc.reps - mc.reps %% n, n) # must be a multiple of n
  times <- max(floor(mc.reps.l1 / n), 1)
  mc.reps.l2 <- length(clusterSizes) * times

  clusterSizes.sim <- rep(clusterSizes, times)
  clusterIdx.sim <- rep(seq_along(clusterSizes.sim), clusterSizes.sim)

  # get target parameters
  par0L1 <- getFreeParamsTable(combinedModel(baseFits$level.1))
  par0L2 <- getFreeParamsTable(combinedModel(baseFits$level.2))
  par1L1 <- par0L1[c("lhs", "op", "rhs", "est", "is.free")]
  par1L2 <- par0L2[c("lhs", "op", "rhs", "est", "is.free")]
  freeL1 <- par0L1$is.free
  freeL2 <- par0L2$is.free

  # icc
  # `icc0` is cor(x, cluster mean), which is not the ICC itself. It is only used
  # as the calibration target. With a (mean) cluster size of nbar we have
  # cor^2 ~= icc + (1 - icc) / nbar, which we invert to get a starting value.
  icc0 <- baseFits$icc
  nbar <- n / length(clusterSizes)
  lower.icc <- rep(0, length(icc0))
  upper.icc <- rep(1, length(icc0))
  start.icc <- pmin(pmax((icc0^2 - 1 / nbar) / (1 - 1 / nbar), 0.01), 0.99)

  # starting parameters
  start <- c(par1L1[par1L1$is.free, "est"], par1L2[par1L2$is.free, "est"], start.icc)
  lower <- c(getMcLowerBounds(par1L1), getMcLowerBounds(par1L2), lower.icc)
  upper <- c(getMcUpperBounds(par1L1), getMcUpperBounds(par1L2), upper.icc)

  .parStruct <- function(p) {
    parxL1 <- par1L1
    parxL2 <- par1L2
    iccx   <- icc0

    n0 <- sum(parxL1$is.free)
    n1 <- sum(parxL2$is.free)

    parxL1[parxL1$is.free, "est"] <- p[seq_len(n0)]
    parxL2[parxL2$is.free, "est"] <- p[n0 + seq_len(n1)]
    iccx[] <- p[n0 + n1 + seq_along(icc0)]

    list(level.1 = parxL1, level.2 = parxL2, icc = iccx)
  }

  .f <- function(p, sim.ov.cont = NULL) {
    
    if (is.null(sim.ov.cont)) {
      parStruct <- .parStruct(p)
      icc <- parStruct$icc

      simL1 <- simulateDataParTable(
        parTable     = parStruct$level.1,
        N            = mc.reps.l1,
        seed         = rng.seed,
        check.hi.ord = is.hi.ord.l1,
        full         = use.full.rescov.l1
      )

      simL2 <- simulateDataParTable(
        parTable     = parStruct$level.2,
        N            = mc.reps.l2,
        seed         = rng.seed,
        check.hi.ord = is.hi.ord.l2,
        full         = use.full.rescov.l2
      )

      sim.ov.l1 <- simL1$ov
      sim.ov.l2 <- simL2$ov[clusterIdx.sim,,drop=FALSE]
      mix <- parsed$ovs.both

      sim.ov.cont <- cbind(
        sim.ov.l1[,parsed$ovs.only.1,drop=FALSE],
        sim.ov.l2[,parsed$ovs.only.2,drop=FALSE],
        sweep(sim.ov.l1[,mix,drop=FALSE], MARGIN = 2, STATS = sqrt(1 - icc), FUN = "*") +
        sweep(sim.ov.l2[,mix,drop=FALSE], MARGIN = 2, STATS = sqrt(icc), FUN = "*")
      )

    }

    # sim.ov  <- ordinalizeDataFrame(
    #   df = sim$ov[idx,,drop=FALSE], thresholdStruct = thresholdStruct
    # )
    sim.ov <- sim.ov.cont

    refit <- refitAuxiliaryMLM_PLS(
      fits = baseFits,
      parsed = parsed,
      data.sim = sim.ov,
      clusterIdx.sim = clusterIdx.sim
    )

    par2L1 <- getFreeParamsTable(combinedModel(refit$level.1))
    par2L2 <- getFreeParamsTable(combinedModel(refit$level.2))

    out <- c(
      par2L1[freeL1, "est"] - par0L1[freeL1, "est"],
      par2L2[freeL2, "est"] - par0L2[freeL2, "est"],
      refit$icc - icc0
    )

    attr(out, "lower") <- c(simL1$lower[freeL1], simL2$lower[freeL2], lower.icc)
    attr(out, "upper") <- c(simL1$upper[freeL1], simL2$upper[freeL2], upper.icc)

    out
  }
  
  mcfit <- robbinsMonro1951(
    p                = start,
    f                = .f,
    tol              = tol,
    min.iter         = min.iter,
    max.iter         = max.iter,
    verbose          = verbose,
    polyak.juditsky  = polyak.juditsky,
    fn.args          = fn.args,
    pj.extrapolate   = pj.extrapolate,
    lower            = lower,
    upper            = upper,
    diag.secant      = diag.secant,
    ...
  )

  .parStruct(c(mcfit$root))
}


parseMultilevelModelArguments <- function(syntax, data) {
  lines <- stringr::str_split_1(syntax, pattern = "\n|;")
  lines <- stringr::str_trim(lines)
  lines <- lines[lines != ""]

  # Check input formatting...
  idxl2 <- which(grepl("^level\\s*:\\s*2$", lines))
  pls_stopif(!length(lines), "model syntax is empty!")
  pls_stopif(!grepl("^level\\s*:\\s*1", lines[[1]]), "first line must be `level: 1`")
  pls_stopif(!length(idxl2), "`level: 2` must be specified!")
  pls_stopif(idxl2 <= 2, "level: 1 must have more than one line")
  pls_stopif(idxl2 == length(lines), "level: 2 must have more than one line")

  s1 <- paste0(lines[2:(idxl2-1)], collapse = "\n")
  s2 <- paste0(lines[(idxl2+1):length(lines)], collapse = "\n")

  parTableL1 <- modsem::modsemify(s1)
  parTableL2 <- modsem::modsemify(s2)

  ovsL1 <- getOVs(parTableL1)
  ovsL2 <- getOVs(parTableL2)

  list(
    level.1    = parTableL1,
    level.2    = parTableL2,
    syntax.1   = s1,
    syntax.2   = s2,
    ovs.1      = ovsL1,
    ovs.2      = ovsL2,
    ovs.only.1 = setdiff(ovsL1, ovsL2),
    ovs.only.2 = setdiff(ovsL2, ovsL1),
    ovs.both   = intersect(ovsL1, ovsL2),
    ovs.all    = union(ovsL1, ovsL2)
  )
}


groupMean <- function(x, g) {
  drop(rowsum(x, g, reorder = FALSE)) / tabulate(g)
}


decompData <- function(data, ovs.1, ovs.2, clusterIdx) {
  k <- length(unique(clusterIdx))
  n <- NROW(data)

  dataL1 <- matrix(
    NA_real_, nrow = n, ncol = length(ovs.1),
    dimnames = list(NULL, ovs.1)
  )

  dataL2 <- matrix(
    NA_real_, nrow = k, ncol = length(ovs.2),
    dimnames = list(NULL, ovs.2)
  )

  vars <- union(ovs.1, ovs.2)
  ovs.both <- intersect(ovs.1, ovs.2)
  icc <- stats::setNames(numeric(length(ovs.both)), ovs.both)

  for (nm in vars) {
    # x is standardized
    x <- data[,nm,drop=TRUE]

    is.l1 <- nm %in% ovs.1
    is.l2 <- nm %in% ovs.2

    if (is.l2) {
      x.l2 <- groupMean(x, g = clusterIdx)
      dataL2[,nm] <- x.l2
    }

    if (is.l1 && is.l2) {
      x.l2.full <- x.l2[clusterIdx]
      x.l1 <- x - x.l2.full # currently not alligned
      dataL1[,nm] <- x.l1
      icc[[nm]] <- cor(x, x.l2.full)

    } else if (is.l1) {
      # exists only at level 1
      dataL1[,nm] <- x
    }
  }

  list(level.1 = dataL1, level.2 = dataL2, icc = icc)
}


fitAuxiliaryMLM_PLS <- function(parsed, data, clusterIdx, ...) {
  decomp <- decompData(
    data = data,
    ovs.1 = parsed$ovs.1,
    ovs.2 = parsed$ovs.2,
    clusterIdx = clusterIdx
  )

  dataL1 <- decomp$level.1
  dataL2 <- decomp$level.2

  fit0L1 <- pls(parsed$syntax.1, data = dataL1, consistent = FALSE, ...)
  fit0L2 <- pls(parsed$syntax.2, data = dataL2, consistent = FALSE, ...)

  list(level.1 = fit0L1, level.2 = fit0L2, icc = decomp$icc)
}


refitAuxiliaryMLM_PLS <- function(fits, parsed, data.sim, clusterIdx.sim) {
  decomp <- decompData(
    data = data.sim,
    ovs.1 = parsed$ovs.1,
    ovs.2 = parsed$ovs.2,
    clusterIdx = clusterIdx.sim
  )

  varsL1 <- colnames(fits$level.1@data)
  varsL2 <- colnames(fits$level.2@data)

  dataL1 <- decomp$level.1
  dataL2 <- decomp$level.2

  fit1 <- fits$level.1
  fit2 <- fits$level.2

  X1 <- Rfast::standardise(dataL1[,varsL1])
  X2 <- Rfast::standardise(dataL2[,varsL2])

  S1 <- Rfast::cova(X1)
  S2 <- Rfast::cova(X2)
    
  # Update observed-data (lowest-order) model input
  modelData(fit1) <- X1
  modelData(fit2) <- X2

  indCorrMatrix(fit1) <- S1
  indCorrMatrix(fit2) <- S2

  # Update fits
  fit1 <- estimatePLS_Inner(fit1)
  fit2 <- estimatePLS_Inner(fit2)

  list(level.1 = fit1, level.2 = fit2, icc = decomp$icc)
}
