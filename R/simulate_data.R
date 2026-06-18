simulateDataParTable <- function(parTable,
                                 N            = 1e5,
                                 seed         = NULL,
                                 .cortol      = .95,
                                 check.hi.ord = FALSE,
                                 clusterSizes = NULL,
                                 clusterName  = NULL,
                                 standardize  = FALSE,
                                 full         = FALSE) {
  if (!is.null(seed) && exists(".Random.seed")) .Random.seed.orig <- .Random.seed
  else                                          .Random.seed.orig <- NULL

  on.exit({
    if (!is.null(.Random.seed.orig)) .Random.seed <<- .Random.seed.orig
  })

  if (!is.null(seed))
    set.seed(seed)

  is.admissible <- TRUE
  checkFixVar <- function(v) {
    if (is.na(v) || !is.finite(v)) {
      is.admissible <<- FALSE
      v <- 0
      attr(v, "ok") <- FALSE
      return(v)
    }

    if (v < 0) {
      is.admissible <<- FALSE
      v <- 0
      attr(v, "ok") <- FALSE
      return(v)
    }

    attr(v, "ok") <- TRUE
    v
  }

  parTable$penalty <- 0 # parameter penalties for inadmissible solutions

  penalty.cfg <- list(
    corr = list(
      limit = abs(.cortol),
      guard = 0,
      beta  = 10,
      scale = 1,
      max_penalty = 2
    ),
    loading = list( # shouldn't be necessary, as we used bounded optimization
      limit = 1,
      guard = 0.000005,
      clamp = 0.000001,
      beta  = 3,
      scale = 0.02,
      tol   = sqrt(.Machine$double.eps),
      max_penalty = 1
    )
  )

  if (check.hi.ord)
    parTable <- highOrdMeasrAsStructParTable(parTable)

  xis     <- getXis(parTable, isLV = !check.hi.ord)
  etas    <- getSortedEtas(parTable)
  mode.a  <- getReflectiveLVs(parTable)
  mode.b  <- getFormativeLVs(parTable)
  lvs     <- unique(c(mode.a, mode.b))
  indsLVs <- getIndsLVs(parTable, lVs = lvs)
  ovs     <- getOVs(parTable)
  mixed   <- !is.null(clusterSizes) && !is.null(clusterName)

  if (mixed) {
    randeff <- getRandomEffectLabels(parTable)

    ovs  <- setdiff(ovs, randeff)
    lvs  <- setdiff(lvs, randeff)
    xis  <- setdiff(xis, randeff)
    etas <- setdiff(etas, randeff)

    if (N < sum(clusterSizes)) {
      N <- sum(clusterSizes)

    } else {
      K <- floor(N / sum(clusterSizes))
      clusterSizes <- rep(clusterSizes, K)
      N <- sum(clusterSizes)

    }

    ncluster <- length(clusterSizes)
    cluster <- rep(seq_along(clusterSizes), clusterSizes)
    clusterMat <- matrix(cluster, nrow = N, dimnames = list(NULL, clusterName))

  } else {
    randeff    <- NULL
    ncluster   <- 0
    cluster    <- NULL
    clusterMat <- NULL

  }

  res <- buildCovMat(
    vars          = xis,
    parTable      = parTable,
    .cortol       = .cortol,
    penalty.cfg   = penalty.cfg,
    unitVariances = TRUE
  )

  parTable <- res$parTable
  Xi <- as.data.frame(Rfast::standardise(rmvnSafe(N, res$mat)))
  colnames(Xi) <- xis

  # Full mode: track the realised disturbances (including exogenous lvs) and,
  # as they are drawn, so each disturbance can be drawn conditional on the
  # prior noise with the specified residual covariances (eta~~eta and xi~~eta).
  # `rescov(v, w)` returns the residual covariance between two nodes
  if (full) {
    disturbances <- as.matrix(Xi[, xis, drop = FALSE])
    dnames <- xis

    rescovRows <- parTable[
      parTable$op == "~~" &
      parTable$lhs != parTable$rhs, , drop = FALSE
    ]

    rescov <- function(v, w) {
      idx <- (
        (rescovRows$lhs == v & rescovRows$rhs == w) |
        (rescovRows$lhs == w & rescovRows$rhs == v)
      )

      if (any(idx)) rescovRows$est[which(idx)[1L]] else 0
    }
  }

  undefIntTerms <- getIntTerms(parTable)
  elemsIntTerms <- stringr::str_split(undefIntTerms, pattern = ":")
  names(elemsIntTerms) <- undefIntTerms

  for (eta in etas) {

    # forward declare
    U          <- NULL
    U.expanded <- NULL

    if (mixed) {
      randeff.eta <- randeff[startsWith(randeff, paste0(eta, "~"))]

      if (length(randeff.eta)) {

        res <- buildCovMat(
          vars          = randeff.eta,
          parTable      = parTable,
          .cortol       = .cortol,
          penalty.cfg   = penalty.cfg,
          unitVariances = FALSE
        )

        parTable   <- res$parTable
        U          <- rmvnSafe(ncluster, res$mat)
        colnames(U) <- randeff.eta
        U.expanded  <- U[cluster, , drop = FALSE]
      }

    }

    for (intTerm in undefIntTerms) {
      elems <- elemsIntTerms[[intTerm]]

      if (all(elems %in% colnames(Xi))) {
        Xi[[intTerm]] <- multiplyIndicatorsCpp(Xi[elems])
        Xi[[intTerm]] <- Xi[[intTerm]] - mean(Xi[[intTerm]])

        undefIntTerms <- setdiff(undefIntTerms, intTerm)
      }
    }

    cond <- parTable$lhs == eta & parTable$op == "~"
    predRows <- parTable[cond, , drop = FALSE]

    vals <- numeric(N)

    # Random Intercept
    par <- paste0(eta, "~1")
    if (par %in% colnames(U))
      vals <- vals + U.expanded[,par]

    for (i in seq_len(NROW(predRows))) {
      row  <- predRows[i, ]
      beta <- row$est
      pred <- row$rhs

      # Fixed effect
      vals <- vals + beta * Xi[[pred]]

      # Random Effect
      par <- paste0(eta, "~", pred)
      if (par %in% colnames(U))
        vals <- vals + U.expanded[,par] * Xi[[pred]]
    }

    resvar <- checkFixVar(1 - stats::var(vals))

    if (!attr(resvar, "ok")) {
      # Recalc coefficients and penalize
      formulaString <- paste(
        eta, "~",
        paste0(predRows$rhs, collapse = " + ")
      )

      Xi[[eta]] <- standardizeAtomic(vals)
      beta.hat <- betacoef(formulaString, data = Xi)

      beta.y <- beta.hat[parTable[cond, "rhs"]]
      beta.x <- parTable[cond, "est"]

      parTable[cond, "penalty"] <- parTable[cond, "penalty"] + (beta.x - beta.y)
    }

    # Disturbance. In `reduced` mode (or when this eta has no residual
    # covariance) the disturbance is independent with variance `resvar`. In
    # `full` mode it is drawn conditional on the prior disturbances so that its
    # covariance with each prior node equals the specified residual covariance:
    # with target cross-covariances `a` and realised noise covariance `M`, the
    # regression `beta = M^-1 a` gives realised Cov(zeta, noise) = a exactly;
    # the conditional variance `resvar - a' beta` carries the fresh part.
    if (full) a <- vapply(dnames, FUN.VALUE = numeric(1L), FUN = \(v) rescov(v, eta))
    else      a <- 0

    if (full && any(a != 0)) {
      M    <- Rfast::cova(disturbances)
      beta <- tryCatch(as.vector(solve(M, a)), error = \(...) numeric(length(a)))

      cmean   <- as.vector(disturbances %*% beta)
      condvar <- checkFixVar(resvar - sum(a * beta))

      if (!attr(condvar, "ok")) {
        # The requested residual covariances exceed what the variances allow
        # (the conditional variance went negative). Penalise the offending rows
        # back toward zero.
        idx <- which(
          parTable$op == "~~" &
          parTable$lhs != parTable$rhs &
          (parTable$lhs == eta | parTable$rhs == eta)
        )

        parTable[idx, "penalty"] <- parTable[idx, "penalty"] +
          sign(parTable[idx, "est"]) * (sum(a * beta) - resvar)
      }

      zeta <- cmean + Rfast::Rnorm(N, m = 0, s = sqrt(condvar), seed = rfast.seed())

    } else {
      zeta <- Rfast::Rnorm(N, m = 0, s = sqrt(resvar), seed = rfast.seed())
    }

    vals <- vals + zeta

    if (full) {
      disturbances <- cbind(disturbances, zeta)
      dnames       <- c(dnames, eta)
    }

    if (standardize)
      vals <- (vals - mean(vals)) / stats::sd(vals)

    Xi[[eta]] <- vals
  }

  Inds <- list()

  for (lv in mode.a) {
    for (ind in indsLVs[[lv]]) {
      cond <- (
        parTable$lhs == lv &
        parTable$op == "=~" &
        parTable$rhs == ind
      )

      lambda.raw <- parTable[cond, "est"]

      lambda.penalty <- smoothBoundaryPenalty(
        lambda.raw,
        limit = penalty.cfg$loading$limit,
        guard = penalty.cfg$loading$guard,
        beta  = penalty.cfg$loading$beta,
        scale = penalty.cfg$loading$scale,
        penalty.max = penalty.cfg$loading$max_penalty
      )

      if (lambda.penalty != 0) {
        parTable[cond, "penalty"] <- parTable[cond, "penalty"] + lambda.penalty
        is.admissible <- FALSE
      }

      clamp.margin <- max(penalty.cfg$loading$clamp, penalty.cfg$loading$tol)
      clamp.limit <- max(0, penalty.cfg$loading$limit - clamp.margin)
      lambda <- clampAbs(lambda.raw, clamp.limit)
      epsilon <- checkFixVar(1 - lambda^2)

      # vals <- lambda * Xi[[lv]] + rnorm(N, mean = 0, sd = sqrt(epsilon))
      vals <- lambda * Xi[[lv]] + Rfast::Rnorm(N, m = 0, s = sqrt(epsilon), seed = rfast.seed())

      if (standardize)
        vals <- (vals - mean(vals)) / stats::sd(vals)

      Inds[[ind]] <- vals
    }
  }

  for (lv in mode.b) {
    inds.lv <- indsLVs[[lv]]
    nind <- length(inds.lv)

    pls_stopif(nind != 1, "Mode B is not available in MC-OrdPLSc (yet)!")
    Inds[[inds.lv]] <- Xi[[lv]]
  }

  Inds <- as.data.frame(Inds)
  All  <- cbind(Xi, Inds)
  Lv   <- All[lvs]
  Ov   <- All[ovs]

  list(
    all           = All,
    ov            = Ov,
    lv            = Lv,
    is.admissible = is.admissible,
    penalty       = parTable$penalty,
    parTable      = parTable,
    cluster       = clusterMat
  )
}


buildCovMat <- function(vars, parTable, .cortol, penalty.cfg, unitVariances = FALSE) {
  mat <- diag(if (unitVariances) 1 else 0, length(vars))
  dimnames(mat) <- list(vars, vars)

  if (!unitVariances) for (v in vars) {
    cond.v <- parTable$rhs == v & parTable$lhs == v & parTable$op == "~~"
    v.i    <- parTable[cond.v, "est"]

    if (v.i <= 0) {
      penalty <- smoothBoundaryPenalty(
        -v.i,
        limit       = 0,
        guard       = 0,
        beta        = penalty.cfg$corr$beta,
        scale       = penalty.cfg$corr$scale,
        penalty.max = penalty.cfg$corr$max_penalty
      )
      if (penalty != 0)
        parTable[cond.v, "penalty"] <- parTable[cond.v, "penalty"] + penalty
      v.i <- .Machine$double.eps
    }

    mat[v, v] <- v.i
  }

  for (i in seq_along(vars)) {
    var.i <- vars[[i]]

    for (j in seq_len(i - 1)) {
      var.j <- vars[[j]]

      cond <- (
        (parTable$rhs == var.i & parTable$op == "~~" & parTable$lhs == var.j) |
        (parTable$rhs == var.j & parTable$op == "~~" & parTable$lhs == var.i)
      )

      denom <- sqrt(mat[i,i] * mat[j,j])
      p.ij  <- parTable[cond, "est"][1] # cov
      if (is.na(p.ij)) p.ij <- 0        # no `~~` row specified -> uncorrelated
      r.ij  <- p.ij / denom             # cor

      if (r.ij <= -.cortol || r.ij >= .cortol) {

        # Get standardized penalty
        penalty.r <- smoothBoundaryPenalty(
          r.ij,
          limit       = penalty.cfg$corr$limit,
          guard       = penalty.cfg$corr$guard,
          beta        = penalty.cfg$corr$beta,
          scale       = penalty.cfg$corr$scale,
          penalty.max = penalty.cfg$corr$max_penalty
        )

        penalty <- denom * penalty.r
        p.ij    <- denom * sign(p.ij) * abs(.cortol)

        if (penalty != 0)
          parTable[cond, "penalty"] <- parTable[cond, "penalty"] + penalty
      }

      mat[i, j] <- mat[j, i] <- p.ij
    }
  }

  list(mat = mat, parTable = parTable)
}


rmvnSafe <- function(n, mat) {
  decomp <- tryCatch(
    chol(mat),
    error = \(e) tryCatch({
      diag(mat) <- diag(mat) + 0.01
      chol(mat)
    }, error = \(e) diag2(mat))
  )
  mvnfast::rmvn(n = n, mu = rep(0, NCOL(mat)), sigma = decomp, isChol = TRUE)
}


betacoef <- function(formulaString, data) {
  # Here we assume that all product terms have been formed
  # explicitly in the data. The reason we do it this way
  # is because the lm() function doesn't treat X:X as a quadratic term

  INTR_OP <- "__INTR__"
  formula <- stats::formula(
    stringr::str_replace_all(formulaString, pattern = ":", replacement = INTR_OP)
  )

  names(data) <- stringr::str_replace_all(
    names(data), pattern = ":", replacement = INTR_OP
  )

  fit <- stats::lm(formula, data = data)
  beta <- stats::coef(fit)

  names(beta) <- stringr::str_replace_all(names(beta), pattern = INTR_OP, replacement = ":")
  beta
}


smoothBoundaryPenalty <- function(value, limit, guard = 0, beta = 4, scale = 1,
                                  penalty.max = Inf) {
  if (!length(value) || scale == 0)
    return(numeric(length(value)))

  guard <- max(guard, 0)
  limit <- abs(limit)
  start <- max(limit - guard, 0)
  delta <- abs(value) - start
  penalty <- numeric(length(value))
  idx <- delta > 0

  if (!any(idx))
    return(penalty)

  if (guard > 0)
    norm <- delta[idx] / guard
  else
    norm <- delta[idx]

  z <- beta * norm
  z <- pmin(z, 700) # avoid overflow inside expm1

  penalty[idx] <- sign(value[idx]) * scale * (expm1(z) - z)

  if (!is.finite(penalty.max) || penalty.max <= 0)
    return(penalty)

  cap <- abs(penalty.max)

  penalty <- cap * tanh(penalty / cap)
  penalty
}


clampAbs <- function(value, limit) {
  if (limit <= 0)
    return(rep(0, length(value)))

  pmin(pmax(value, -limit), limit)
}


rfast.seed <- function() {
  # Rfast doensn't work correctly with set.seed()
  # Instead we have to pass a seed to Rfast.
  # Here we generate a random seed, yielding deterministic
  # results is set.seed() has been used.
  floor(stats::runif(1L, min = 0, max = 9999999))
}
