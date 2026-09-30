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
                 ordered = NULL,
                 consistent = FALSE,
                 small.sample = TRUE,
                 small.sample.max.k = 100L,
                 small.sample.point.estimate = c("mean", "median"),
                 mc.reps = 50000,
                 rng.seed = NULL,
                 ...) {
  pls_stopif(length(cluster) != 1 || !is.character(cluster),
    "cluster must be a character string of length 1!"
  )

  small.sample.point.estimate <- match.arg(small.sample.point.estimate)

  parsed <- parseMultilevelModelArguments(
    syntax  = syntax,
    data    = data,
    cluster = cluster
  )

  data       <- parsed$data
  vars.all   <- parsed$vars.all
  parTableL1 <- parsed$level.1 
  parTableL2 <- parsed$level.2

  data[[cluster]] <- as.integer(as.factor(data[[cluster]]))
  clusterIdx <- data[[cluster]]

  if (anyNA(clusterIdx)) {
    pls_stopif(all(is.na(clusterIdx)), "cluster is all NA!")
    pls_msg_warn("removing missing values in `cluster`!")

    data <- data[!is.na(clusterIdx),, drop = FALSE]
    clusterIdx <- data[[cluster]]
  }

  # ordered variables
  is.ord  <- vapply(data[parsed$ovs.all], FUN.VALUE = logical(1L), FUN = is.ordered)
  ordered <- intersect(union(ordered, parsed$ovs.all[is.ord]), parsed$ovs.all)

  for (ord in ordered)
    data[[ord]] <- as.integer(as.ordered(data[[ord]]))

  # data must be sorted by the clusters
  data <- as.matrix(data[order(clusterIdx), vars.all,drop=FALSE])
  thresholdStruct0 <- ThresholdStruct(data, ordered = ordered)
  data[,parsed$ovs.all] <- Rfast::standardise(data[,parsed$ovs.all])
  clusterIdx <- data[,cluster, drop=TRUE]
  n <- NROW(data)

  # fit auxiliary models
  baseFits <- fitAuxiliaryMLM_PLS(
    parsed     = parsed,
    data       = data,
    rpar       = parsed$rpar,
    cluster    = cluster,
    clusterIdx = clusterIdx,
    consistent = consistent,
    ...
  )

  is.hi.ord.l1 <- isTRUE(combinedModel(baseFits$level.1)@info$is.high.ord)
  is.hi.ord.l2 <- isTRUE(combinedModel(baseFits$level.2)@info$is.high.ord)
  use.full.rescov.l1 <- combinedModel(baseFits$level.1)@info$path.estimator == "gls"
  use.full.rescov.l2 <- combinedModel(baseFits$level.2)@info$path.estimator == "gls"

  # The simulated data consists of `times` replicates of the
  # observed cluster structure, stacked on top of each other. With
  # `small.sample = TRUE`, each replicate is refitted separately, with the
  # same sample size and number of clusters as the observed data
  clusterSizes <- table(clusterIdx)
  nclusters    <- length(clusterSizes)
  mc.reps.l1   <- max(mc.reps - mc.reps %% n, n) # must be a multiple of n
  times        <- max(floor(mc.reps.l1 / n), 1)

  if (small.sample) {
    times      <- min(times, max(round(small.sample.max.k), 1))
    mc.reps.l1 <- times * n
  }

  mc.reps.l2 <- nclusters * times

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

  # random slopes (standard deviations)
  rslopes   <- parsed$rslopes
  rsd0      <- baseFits$rsd
  lower.rsd <- rep(0, length(rsd0))
  upper.rsd <- rep(1, length(rsd0))

  target <- c(
    par0L1[freeL1, "est"],
    par0L2[freeL2, "est"],
    icc0,
    rsd0
  )

  # starting parameters
  start <- c(par1L1[par1L1$is.free, "est"], par1L2[par1L2$is.free, "est"], start.icc, rsd0)
  lower <- c(getMcLowerBounds(par1L1), getMcLowerBounds(par1L2), lower.icc, lower.rsd)
  upper <- c(getMcUpperBounds(par1L1), getMcUpperBounds(par1L2), upper.icc, upper.rsd)

  .parStruct <- function(p) {
    parxL1 <- par1L1
    parxL2 <- par1L2
    iccx   <- icc0
    rsdx   <- rsd0

    n0 <- sum(parxL1$is.free)
    n1 <- sum(parxL2$is.free)
    n2 <- length(icc0)

    parxL1[parxL1$is.free, "est"] <- p[seq_len(n0)]
    parxL2[parxL2$is.free, "est"] <- p[n0 + seq_len(n1)]
    iccx[] <- p[n0 + n1 + seq_along(icc0)]
    rsdx[] <- p[n0 + n1 + n2 + seq_along(rsd0)]

    list(level.1 = parxL1, level.2 = parxL2, icc = iccx, rsd = rsdx)
  }

  # The levels must use different seeds, otherwise a fixed seed yields
  # identical draws (and thus correlated components) at both levels.
  rng.seed.l1 <- if (is.null(rng.seed)) NULL else rng.seed + 1L

  .simulate <- function(p) {
    parStruct <- .parStruct(p)
    icc <- parStruct$icc

    simL2 <- simulateDataParTable(
      parTable     = parStruct$level.2,
      N            = mc.reps.l2,
      seed         = rng.seed,
      check.hi.ord = is.hi.ord.l2,
      full         = use.full.rescov.l2
    )

    parTableSimL1 <- parStruct$level.1
    exogenous     <- NULL

    if (NROW(rslopes)) {

      # treat the random effect as an interaction term, where the coefficient
      # is the standard deviation of the random effect
      parTableSimL1 <- rbind(parTableSimL1, data.frame(
        lhs     = rslopes$lhs,
        op      = "~",
        rhs     = paste0(rslopes$name, ":", rslopes$rhs),
        est     = unname(parStruct$rsd[rslopes$name]),
        is.free = FALSE
      ))

      exogenous <- as.data.frame(
        simL2$all[clusterIdx.sim, rslopes$name, drop = FALSE]
      )
    }

    simL1 <- simulateDataParTable(
      parTable     = parTableSimL1,
      N            = mc.reps.l1,
      seed         = rng.seed.l1,
      check.hi.ord = is.hi.ord.l1,
      full         = use.full.rescov.l1,
      exogenous    = exogenous
    )

    sim.ov.l1 <- toOriginalNames(simL1$ov)
    sim.ov.l2 <- toOriginalNames(simL2$ov)[clusterIdx.sim,,drop=FALSE]
    mix <- parsed$ovs.both

    ov <- cbind(
      sim.ov.l1[,parsed$ovs.only.1,drop=FALSE],
      sim.ov.l2[,parsed$ovs.only.2,drop=FALSE],
      sweep(sim.ov.l1[,mix,drop=FALSE], MARGIN = 2, STATS = sqrt(1 - icc), FUN = "*") +
      sweep(sim.ov.l2[,mix,drop=FALSE], MARGIN = 2, STATS = sqrt(icc), FUN = "*")
    )

    # the random slope rows are appended after the level 1 parameters
    nL1   <- NROW(parStruct$level.1)
    idxL1 <- seq_len(nL1)
    idxR  <- nL1 + seq_len(NROW(rslopes))

    lower <- c(
      simL1$lower[idxL1][freeL1], simL2$lower[freeL2], lower.icc,
      pmax(simL1$lower[idxR], lower.rsd)
    )

    upper <- c(
      simL1$upper[idxL1][freeL1], simL2$upper[freeL2], upper.icc,
      pmin(simL1$upper[idxR], upper.rsd)
    )

    list(ov = ov, lower = lower, upper = upper)
  }

  .f <- function(p, sim = NULL) {
    if (is.null(sim))
      sim <- .simulate(p)

    if (length(ordered)) {
      sim.ov <- ordinalizeDataFrame(
        df = as.data.frame(sim$ov), thresholdStruct = thresholdStruct0
      )
    } else {
      sim.ov <- sim$ov
    }

    .estimates <- function(rows = seq_len(NROW(sim.ov)), offset = 0L) {
      refit <- refitAuxiliaryMLM_PLS(
        fits           = baseFits,
        parsed         = parsed,
        data.sim       = sim.ov[rows, , drop = FALSE],
        rpar           = parsed$rpar,
        cluster        = cluster,
        clusterIdx.sim = clusterIdx.sim[rows] - offset
      )

      par2L1 <- getFreeParamsTable(combinedModel(refit$level.1))
      par2L2 <- getFreeParamsTable(combinedModel(refit$level.2))

      c(par2L1[freeL1, "est"], par2L2[freeL2, "est"], refit$icc, refit$rsd)
    }

    if (small.sample) {
      THETA <- matrix(NA_real_, nrow = times, ncol = length(target))

      for (i in seq_len(times)) {
        THETA[i, ] <- tryCatch(
          .estimates(rows = (i - 1L) * n + seq_len(n), offset = (i - 1L) * nclusters),
          error = \(e) NA_real_
        )
      }

      est <- switch(small.sample.point.estimate,
        mean   = colMeans(THETA, na.rm = TRUE),
        median = colMedians(THETA, na.rm = TRUE)
      )

    } else {
      est <- .estimates()
    }

    out <- est - target

    attr(out, "lower") <- sim$lower
    attr(out, "upper") <- sim$upper

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
    diag.secant      = diag.secant
  )

  root <- c(mcfit$root)
  out  <- .parStruct(root)

  # thresholds of the latent response (total) scores
  out$thresholds <- updateThresholds(
    thr = thresholdStruct0, sim.cont = .simulate(root)$ov
  )@thresholds

  out
}


parseMultilevelModelArguments <- function(syntax, data, cluster) {
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

  # The level 2 model is fitted with `strict = FALSE` (see below), so we do
  # the checks of `specifyModel()` here instead.
  nm <- unique(c(parTableL1$lhs, parTableL1$rhs, parTableL2$lhs, parTableL2$rhs))
  pls_stopif(any(hasTempAffixes(nm)),
    "Some variables have reserved keywords/patterns!",
    "Variables:", paste0(nm[hasTempAffixes(nm)], collapse = ", ")
  )

  hasIntr <- c(parTableL1$op, parTableL2$op) == "~1"
  pls_stopif(any(hasIntr),
    "Estimation of intercepts is not available!",
    "Intercepts:", paste0(c(parTableL1$lhs, parTableL2$lhs)[hasIntr], "~1", collapse = ", ")
  )

  # Random slopes are defined at level 1 (e.g., fw ~ rv(s1)*x1), and enter
  # the level 2 model as (observed) variables (approximated in level 1)
  isRand  <- isRandomEffectMod(parTableL1$mod)
  rslopes <- data.frame(
    name = extractRandomEffectName(parTableL1$mod[isRand]),
    lhs  = parTableL1$lhs[isRand],
    rhs  = parTableL1$rhs[isRand]
  )
  rpar <- rslopes$name

  pls_stopif(any(rpar %in% c(parTableL1$lhs, parTableL1$rhs)),
    "Random slopes cannot share names with variables at level 1!",
    "Random slopes:", paste0(rpar, collapse = ", ")
  )

  # Random slopes which don't appear in the level 2 model are declared as
  # (observed) exogenous variables (`s ~ 1`), such that they covary with the
  # other exogenous variables at level 2.
  missingL2 <- setdiff(rpar, c(parTableL2$lhs, parTableL2$rhs))

  if (length(missingL2)) {
    s2 <- paste0(c(s2, paste0(missingL2, " ~ 1")), collapse = "\n")
    parTableL2 <- modsem::modsemify(s2)
  }

  ovsL1 <- getOVs(parTableL1)
  ovsL2 <- setdiff(getOVs(parTableL2), rpar)

  data <- as.data.frame(data)

  vars.all <- c(cluster, union(ovsL1, ovsL2))
  missing <- setdiff(vars.all, colnames(data))
  pls_stopif(length(missing),
    "Missing variables in `data`:", paste0(missing, collapse = ", ")
  )

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
    ovs.all    = union(ovsL1, ovsL2),
    vars.all   = vars.all,
    rpar       = rpar,
    rslopes    = rslopes,
    data       = data
  )
}


groupMean <- function(x, g) {
  drop(rowsum(x, g, reorder = FALSE)) / tabulate(g)
}


decompData <- function(data, ovs.1, ovs.2, cluster, clusterIdx) {
  k <- length(unique(clusterIdx))
  n <- NROW(data)

  dataL1 <- matrix(
    NA_real_, nrow = n, ncol = length(ovs.1) + 1L,
    dimnames = list(NULL, c(ovs.1, cluster))
  )
  dataL1[,cluster] <- clusterIdx

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

  list(
    level.1 = dataL1,
    level.2 = dataL2,
    icc     = icc
  )
}


fitAuxiliaryMLM_PLS <- function(parsed,
                                data,
                                cluster,
                                clusterIdx,
                                rpar,
                                consistent = FALSE,
                                ...) {
  decomp <- decompData(
    data       = data,
    ovs.1      = parsed$ovs.1,
    ovs.2      = parsed$ovs.2,
    cluster    = cluster,
    clusterIdx = clusterIdx
  )

  dataL1 <- decomp$level.1
  dataL2 <- decomp$level.2

  fit0L1 <- pls(
    syntax     = parsed$syntax.1,
    data       = dataL1,
    consistent = consistent,
    cluster    = cluster,
    ...
  )

  slopes <- getRandomSlopes(fit0L1, rpar = rpar, k = NROW(dataL2))
  dataL2 <- cbind(dataL2, slopes)

  # `strict = FALSE` allows the `s ~ 1` declarations of random slopes
  fit0L2 <- pls(
    syntax     = parsed$syntax.2,
    data       = dataL2,
    consistent = consistent,
    strict     = FALSE,
    ...
  )

  list(
    level.1 = fit0L1,
    level.2 = fit0L2,
    icc     = decomp$icc,
    rsd     = slopeSDs(slopes, rpar = rpar)
  )
}


refitAuxiliaryMLM_PLS <- function(fits, parsed, data.sim, rpar, cluster, clusterIdx.sim) {
  decomp <- decompData(
    data       = data.sim,
    ovs.1      = parsed$ovs.1,
    ovs.2      = parsed$ovs.2,
    cluster    = cluster,
    clusterIdx = clusterIdx.sim
  )

  varsL1 <- colnames(fits$level.1@data)
  varsL2 <- colnames(fits$level.2@data)

  dataL1 <- decomp$level.1
  dataL2 <- decomp$level.2

  fit1 <- fits$level.1
  fit2 <- fits$level.2

  # Level 1. The cluster is needed for the random slopes (lmer)
  X1 <- Rfast::standardise(toInternalNames(dataL1, vars = varsL1))
  colnames(X1) <- varsL1
  attr(X1, "cluster") <- as.data.frame(dataL1[, cluster, drop = FALSE])

  modelData(fit1)     <- X1
  indCorrMatrix(fit1) <- Rfast::cova(X1)
  fit1 <- estimatePLS_Inner(fit1)

  # Level 2. The (estimated) random slopes from level 1 enter as variables
  slopes <- getRandomSlopes(fit1, rpar = rpar, k = NROW(dataL2))
  dataL2 <- cbind(dataL2, slopes)

  X2 <- Rfast::standardise(toInternalNames(dataL2, vars = varsL2))
  colnames(X2) <- varsL2

  modelData(fit2)     <- X2
  indCorrMatrix(fit2) <- Rfast::cova(X2)
  fit2 <- estimatePLS_Inner(fit2)

  list(
    level.1 = fit1,
    level.2 = fit2,
    icc     = decomp$icc,
    rsd     = slopeSDs(slopes, rpar = rpar)
  )
}


# cluster specific slopes (fixed + random effects), ordered by cluster index
getRandomSlopes <- function(fit, rpar, k) {
  if (!length(rpar))
    return(NULL)

  randef <- modelFit(fit)$randef
  pls_stopif(is.null(randef), "Random slopes were not estimated at level 1!")

  randef[as.character(seq_len(k)), rpar, drop = FALSE]
}


slopeSDs <- function(slopes, rpar) {
  if (is.null(slopes))
    return(stats::setNames(numeric(0), character(0)))

  apply(slopes[, rpar, drop = FALSE], MARGIN = 2L, FUN = stats::sd)
}


toOriginalNames <- function(X) {
  original <- removeTempAffixes(colnames(X))
  keep     <- !duplicated(original)

  X <- X[, keep, drop = FALSE]
  colnames(X) <- original[keep]
  X
}


toInternalNames <- function(X, vars) {
  X <- X[, removeTempAffixes(vars), drop = FALSE]
  colnames(X) <- vars
  X
}
