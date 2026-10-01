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
                 bootstrap = FALSE,
                 boot.R = 500L,
                 boot.parallel = c("no", "multicore", "multisession"),
                 boot.ncores = 1L,
                 boot.iseed = NULL,
                 delta.jacobian.k = 1L,
                 delta.eps = 5e-3,
                 level2.cov = c("means", "muml"),
                 ...) {
  pls_stopif(length(cluster) != 1 || !is.character(cluster),
    "cluster must be a character string of length 1!"
  )

  small.sample.point.estimate <- match.arg(small.sample.point.estimate)
  boot.parallel <- match.arg(boot.parallel)
  level2.cov    <- match.arg(level2.cov)

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
    level2.cov = level2.cov,
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

  # (naive) statistics of the auxiliary fits, which are matched to `target`
  .stats <- function(refit) {
    par2L1 <- getFreeParamsTable(combinedModel(refit$level.1))
    par2L2 <- getFreeParamsTable(combinedModel(refit$level.2))
    c(par2L1[freeL1, "est"], par2L2[freeL2, "est"], refit$icc, refit$rsd)
  }

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
  .simulate <- function(p, seed = rng.seed) {
    seed.l1   <- if (is.null(seed)) NULL else seed + 1L
    parStruct <- .parStruct(p)
    icc <- parStruct$icc

    simL2 <- simulateDataParTable(
      parTable     = parStruct$level.2,
      N            = mc.reps.l2,
      seed         = seed,
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
      seed         = seed.l1,
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

    list(ov = ov, lower = lower, upper = upper, sim.l1 = simL1, sim.l2 = simL2)
  }

  .f <- function(p, sim = NULL, seed = rng.seed) {
    if (is.null(sim))
      sim <- .simulate(p, seed = seed)

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
        clusterIdx.sim = clusterIdx.sim[rows] - offset,
        level2.cov     = level2.cov
      )

      .stats(refit)
    }

    if (small.sample) {
      est <- averageMcReplicates(
        k              = times,
        point.estimate = small.sample.point.estimate,
        catch          = TRUE,
        fun            = \(i) .estimates(
          rows   = (i - 1L) * n + seq_len(n),
          offset = (i - 1L) * nclusters
        )
      )

    } else {
      est <- .estimates()
    }

    out <- est - target

    attr(out, "lower") <- sim$lower
    attr(out, "upper") <- sim$upper

    out
  }
  
  mcfit <- solveMcRoot(
    p               = start,
    f               = .f,
    lower           = lower,
    upper           = upper,
    tol             = tol,
    min.iter        = min.iter,
    max.iter        = max.iter,
    verbose         = verbose,
    polyak.juditsky = polyak.juditsky,
    pj.extrapolate  = pj.extrapolate,
    fn.args         = fn.args,
    diag.secant     = diag.secant
  )

  root      <- c(mcfit$root)
  parStruct <- .parStruct(root)
  simRoot   <- .simulate(root)

  level.1 <- finalizeMultilevelLevel(
    model = baseFits$level.1, parTable = parStruct$level.1,
    sim = simRoot$sim.l1, iterations = mcfit$iter
  )

  level.2 <- finalizeMultilevelLevel(
    model = baseFits$level.2, parTable = parStruct$level.2,
    sim = simRoot$sim.l2, iterations = mcfit$iter
  )

  # thresholds of the latent response (total) scores
  thresholdStruct <- updateThresholds(thr = thresholdStruct0, sim.cont = simRoot$ov)

  parTableInput <- rbind(
    cbind(parsed$level.1, level = 1L),
    cbind(parsed$level.2, level = 2L)
  )

  is.ord    <- length(ordered) > 0L
  estimator <- paste0("MC", if (is.ord) "Ord" else "", "PLSc-MLM")

  converged  <- mcfit$ok
  admissible <- converged && isAdmissible(level.1) && isAdmissible(level.2)

  model <- PlsMultilevelModel(
    level.1         = level.1,
    level.2         = level.2,
    info            = list(
      estimator    = estimator,
      cluster      = cluster,
      n            = n,
      nclusters    = nclusters,
      ordered      = ordered,
      rslopes      = rslopes,
      small.sample = small.sample,
      level2.cov   = level2.cov,
      mc.reps      = mc.reps,
      rng.seed     = rng.seed
    ),
    data            = data,
    thresholdStruct = thresholdStruct,
    status          = list(
      iterations    = mcfit$iter,
      converged     = converged,
      diverged      = mcfit$diverged,
      is.admissible = admissible
    ),
    params          = list(
      icc        = parStruct$icc,
      rsd        = parStruct$rsd,
      thresholds = thresholdStruct@thresholds
    ),
    fit             = list(mcfit = mcfit),
    factorScores    = list(
      level.1 = modelFactorScores(combinedModel(level.1)),
      level.2 = modelFactorScores(combinedModel(level.2))
    ),
    parTableInput   = parTableInput
  )

  model@parTable <- getParTableMultilevel(model)
  model@params$values <- getCoefsMultilevel(model)

  if (!bootstrap)
    return(model)

  # Standard errors (delta method)
  .values <- function(p, sim) {
    parStruct <- .parStruct(p)

    l1 <- finalizeMultilevelLevel(baseFits$level.1, parStruct$level.1, sim$sim.l1, iterations = 0L)
    l2 <- finalizeMultilevelLevel(baseFits$level.2, parStruct$level.2, sim$sim.l2, iterations = 0L)

    c(
      getLevelParamValues(l1, l2),
      prefixNames(parStruct$icc, prefix = "icc."),
      prefixNames(parStruct$rsd, prefix = "rsd.")
    )
  }

  .bootstrap <- function(b) {
    tryCatch({
      idx  <- sample(nclusters, size = nclusters, replace = TRUE)
      rows <- unlist(split(seq_len(n), clusterIdx)[idx], use.names = FALSE)
      cl   <- rep(seq_along(idx), times = clusterSizes[idx])

      Xb <- data[rows, , drop = FALSE]
      Xb[, cluster] <- cl
      Xb[, parsed$ovs.all] <- Rfast::standardise(Xb[, parsed$ovs.all, drop = FALSE])

      .stats(refitAuxiliaryMLM_PLS(
        fits           = baseFits,
        parsed         = parsed,
        data.sim       = Xb,
        rpar           = parsed$rpar,
        cluster        = cluster,
        clusterIdx.sim = cl,
        level2.cov     = level2.cov
      ))

    }, error = \(e) structure(rep(NA_real_, length(target)), error = conditionMessage(e)))
  }

  .fg <- function(p, seed) {
    sim <- .simulate(p, seed = seed)
    list(f = c(.f(p, sim = sim)), g = .values(p, sim = sim))
  }

  if (is.null(boot.iseed)) boot.iseed <- floor(stats::runif(1L, min = 0, max = 999999999))

  if (verbose) pls_msg_note("Bootstrapping auxiliary models...")
  results <- runMcReplicates(
    R        = boot.R,
    fun      = .bootstrap,
    parallel = boot.parallel,
    ncores   = boot.ncores,
    verbose  = verbose,
    iseed    = boot.iseed,
    label    = "Bootstrap"
  )

  errors <- unlist(lapply(results, FUN = attr, which = "error"))
  T0     <- do.call(rbind, results)

  n.failed <- sum(!stats::complete.cases(T0))
  pls_stopif(n.failed >= boot.R - 1L, # at least two replicates are needed
    sprintf("%d (out of %d) bootstrap replicate(s) failed!", n.failed, boot.R),
    if (length(errors)) paste("First error:", errors[[1L]])
  )

  pls_warnif(n.failed > 0L,
    sprintf("%d (out of %d) bootstrap replicate(s) failed!", n.failed, boot.R),
    if (length(errors)) paste("First error:", errors[[1L]])
  )

  V0 <- stats::cov(T0, use = "complete.obs")

  if (verbose) pls_msg_note("Calculating Jacobian...")
  seeds  <- floor(stats::runif(delta.jacobian.k, min = 0, max = 9999999))
  values <- .values(root, sim = simRoot)

  JAC <- calcMcJacobians(
    .fg      = .fg,
    seeds    = seeds,
    p0       = unname(root), # the names of the root are empty
    p1       = values,
    lower    = lower,
    upper    = upper,
    eps      = delta.eps,
    parallel = boot.parallel,
    ncores   = boot.ncores,
    verbose  = verbose,
    iseed    = boot.iseed
  )

  J0 <- JAC$J0
  J1 <- JAC$J1

  vcov <- deltaMcVcov(invertMcJacobian(J0), V = V0, J1 = J1)
  se   <- sqrt(pmax(diag(vcov), 0))
  se[se <= 1e-10] <- NA_real_ # fixed/constant parameters

  isIcc <- startsWith(names(se), "icc.")
  isRsd <- startsWith(names(se), "rsd.")
  isPar <- !isIcc & !isRsd

  model@params$vcov   <- plssemMatrix(vcov[isPar, isPar, drop = FALSE], is.public = TRUE)
  model@params$se     <- se[isPar]
  model@params$icc.se <- stats::setNames(se[isIcc], names(parStruct$icc))
  model@params$rsd.se <- stats::setNames(se[isRsd], names(parStruct$rsd))
  model@parTable      <- setParTableMultilevelSE(model@parTable, se = se[isPar])

  model@boot <- list(
    R         = boot.R,
    iseed     = boot.iseed,
    naive     = T0,
    vcov.naive = V0,
    J0        = J0,
    J1        = J1
  )

  model
}


prefixNames <- function(x, prefix) {
  if (!length(x)) return(numeric(0L))
  stats::setNames(x, paste0(prefix, names(x)))
}


# Named values of the parameters at both levels (level 2 names are suffixed by
# `.l2`), in the same order as `getCoefsMultilevel()`.
getLevelParamValues <- function(level.1, level.2) {
  pt1 <- parameter_estimates(level.1)
  pt2 <- parameter_estimates(level.2)

  c(
    stats::setNames(pt1$est, paste0(pt1$lhs, pt1$op, pt1$rhs)),
    stats::setNames(pt2$est, paste0(pt2$lhs, pt2$op, pt2$rhs, ".l2"))
  )
}


setParTableMultilevelSE <- function(parTable, se) {
  nm <- paste0(parTable$lhs, parTable$op, parTable$rhs)
  l2 <- which(parTable$level == 2L)
  nm[l2] <- paste0(nm[l2], ".l2")

  parTable$se       <- unname(se[nm]) # NA for thresholds (not included)
  parTable$z        <- parTable$est / parTable$se
  parTable$pvalue   <- 2 * stats::pnorm(-abs(parTable$z))
  parTable$ci.lower <- parTable$est - CI_QUANTILE * parTable$se
  parTable$ci.upper <- parTable$est + CI_QUANTILE * parTable$se

  plssemParTable(parTable)
}


# Update a level model with the calibrated parameters (similar to `mcpls()`),
# and store its parameter table (similar to `pls()`).
finalizeMultilevelLevel <- function(model, parTable, sim, iterations) {
  combined <- combinedModel(model)

  combined <- updateModelFromFreeParTableMC(
    parTable        = parTable,
    model           = combined,
    mc.reps         = NROW(sim$all),
    thresholdStruct = combined@thresholdStruct,
    ordered         = NULL,
    sim             = sim,
    params.only     = TRUE
  )

  combined@status$iterations <- iterations
  combined@parTable <- getParTableEstimates(combined)

  if (hasCombinedModel(model) || hasHigherOrderModel(model)) {
    model@combinedModel <- combined
    return(model)
  }

  combined
}


getParTableMultilevel <- function(model) {
  pt1 <- parameter_estimates(model@level.1)
  pt2 <- parameter_estimates(model@level.2)

  pt1$level <- rep(1L, NROW(pt1))
  pt2$level <- rep(2L, NROW(pt2))

  parTable <- rbind(as.data.frame(pt1), as.data.frame(pt2))

  # The thresholds are defined for the total (latent response) scores, and
  # are thus level agnostic (`level = NA`), similar to custom parameters.
  thr <- model@params$thresholds

  if (length(thr)) {
    split <- stringr::str_split_fixed(names(thr), pattern = "\\|", n = 2L)

    rows <- parTable[rep(NA_integer_, length(thr)), , drop = FALSE]
    rows$lhs   <- split[, 1L]
    rows$op    <- "|"
    rows$rhs   <- split[, 2L]
    rows$label <- ""
    rows$est   <- unname(thr)
    rows$level <- NA_integer_

    parTable <- rbind(parTable, rows)
  }

  plssemParTable(parTable)
}


# Parameters at level 2 get a `.l2` suffix, as the same parameter
# (e.g., `y1~~y1`) can appear at both levels.
getCoefsMultilevel <- function(model) {
  pt <- model@parTable
  nm <- paste0(pt$lhs, pt$op, pt$rhs)
  l2 <- which(pt$level == 2L)
  nm[l2] <- paste0(nm[l2], ".l2")

  stats::setNames(pt$est, nm)
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
                                level2.cov = "means",
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

  if (level2.cov != "means") {
    indCorrMatrix(fit0L2) <- getLevel2CorMat(
      dataL1     = dataL1,
      dataL2     = dataL2,
      clusterIdx = clusterIdx,
      ovs.both   = intersect(parsed$ovs.1, parsed$ovs.2),
      vars       = colnames(fit0L2@data),
      level2.cov = level2.cov
    )

    fit0L2 <- estimatePLS_Inner(fit0L2)
  }

  list(
    level.1 = fit0L1,
    level.2 = fit0L2,
    icc     = decomp$icc,
    rsd     = slopeSDs(slopes, rpar = rpar)
  )
}


refitAuxiliaryMLM_PLS <- function(fits, parsed, data.sim, rpar, cluster, clusterIdx.sim,
                                  level2.cov = "means") {
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

  modelData(fit2) <- X2

  if (level2.cov == "means") {
    indCorrMatrix(fit2) <- Rfast::cova(X2)
  } else {
    indCorrMatrix(fit2) <- getLevel2CorMat(
      dataL1     = dataL1,
      dataL2     = dataL2,
      clusterIdx = clusterIdx.sim,
      ovs.both   = intersect(parsed$ovs.1, parsed$ovs.2),
      vars       = varsL2,
      level2.cov = level2.cov
    )
  }

  fit2 <- estimatePLS_Inner(fit2)

  list(
    level.1 = fit1,
    level.2 = fit2,
    icc     = decomp$icc,
    rsd     = slopeSDs(slopes, rpar = rpar)
  )
}


getLevel2CorMat <- function(dataL1, dataL2, clusterIdx, ovs.both, vars,
                            level2.cov = "muml") {
  S <- stats::cov(dataL2)

  if (level2.cov == "muml" && length(ovs.both)) {
    # Muthen, 1994. Also used for the starting values of two-level models in lavaan
    #  S_PW    = sum_j sum_i (y_ij - ybar_j)(y_ij - ybar_j)' / (N - G)
    #  S_B     = sum_j n_j (ybar_j - ybar)(ybar_j - ybar)' / (G - 1)
    #  Sigma_B = (S_B - S_PW) / c,   c = (N - sum_j n_j^2 / N) / (G - 1)
    n  <- NROW(dataL1)
    G  <- NROW(dataL2)
    nj <- tabulate(clusterIdx, nbins = G)

    W <- dataL1[, ovs.both, drop = FALSE] # within deviations
    M <- dataL2[, ovs.both, drop = FALSE] # cluster means
    D <- sweep(M, MARGIN = 2L, STATS = colSums(M * nj) / n) # deviations from the grand mean

    S.PW <- crossprod(W) / (n - G)
    S.B  <- crossprod(D * sqrt(nj)) / (G - 1)
    c    <- (n - sum(nj^2) / n) / (G - 1)

    S[ovs.both, ovs.both] <- (S.B - S.PW) / c
    S <- clipEigenvalues(S)
  }

  original <- removeTempAffixes(vars)
  S <- stats::cov2cor(S[original, original, drop = FALSE])
  dimnames(S) <- list(vars, vars)

  S
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
