simulateDataParTable <- function(parTable, N = 1e5, seed = NULL, .covtol = .95,
                                 check.hi.ord = FALSE) {
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
      limit = abs(.covtol),
      guard = 0,
      beta  = 10,
      scale = 1,
      max_penalty = 2
    ),
    loading = list(
      limit = 1,
      guard = 0.005,
      clamp = 0.001,
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

  Psi.x <- diag(1, length(xis))
  dimnames(Psi.x) <- list(xis, xis)

  for (i in seq_along(xis)) {
    xi.i <- xis[[i]]

    for (j in seq_len(i - 1)) {
      xi.j <- xis[[j]]

      cond <- (
        (parTable$rhs == xi.i & parTable$op == "~~" & parTable$lhs == xi.j) |
        (parTable$rhs == xi.j & parTable$op == "~~" & parTable$lhs == xi.i)
      )

      p.ij <- parTable[cond, "est"][1]

      if (p.ij <= -.covtol || p.ij >= .covtol) {
        penalty <- smoothBoundaryPenalty(
          p.ij,
          limit = penalty.cfg$corr$limit,
          guard = penalty.cfg$corr$guard,
          beta  = penalty.cfg$corr$beta,
          scale = penalty.cfg$corr$scale,
          penalty.max = penalty.cfg$corr$max_penalty
        )

        if (penalty != 0)
          parTable[cond, "penalty"] <- parTable[cond, "penalty"] + penalty

        p.ij <- sign(p.ij) * abs(.covtol)
      }

      Psi.x[i, j] <- Psi.x[j, i] <- p.ij
    }
  }

  # Xi <- mvtnorm::rmvnorm(n = N, mean = rep(0, length(xis)), sigma = Psi.x)
  Psi.x.decomp <- tryCatch(
    chol(Psi.x),
    error = \(e) {
      tryCatch({
        diag(Psi.x) <- 1.01
        chol(Psi.x)
      }, error = \(e) diag2(Psi.x))
    }
  )

  Xi <- mvnfast::rmvn(n = N, mu = rep(0, length(xis)), sigma = Psi.x.decomp, isChol = TRUE)
  Xi <- as.data.frame(Rfast::standardise(Xi))
  colnames(Xi) <- xis

  undefIntTerms <- getIntTerms(parTable)
  elemsIntTerms <- stringr::str_split(undefIntTerms, pattern = ":")
  names(elemsIntTerms) <- undefIntTerms

  for (eta in etas) {

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

    for (i in seq_len(NROW(predRows))) {
      row  <- predRows[i, ]
      beta <- row$est
      pred <- row$rhs

      vals <- vals + beta * Xi[[pred]]
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

    vals <- vals + Rfast::Rnorm(N, m = 0, s = sqrt(resvar), seed = rfast.seed())
    # vals <- vals + rnorm(N, mean = 0, sd = sqrt(resvar))

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

      if (lambda.penalty != 0)
        parTable[cond, "penalty"] <- parTable[cond, "penalty"] + lambda.penalty

      clamp.margin <- max(penalty.cfg$loading$clamp, penalty.cfg$loading$tol)
      clamp.limit <- max(0, penalty.cfg$loading$limit - clamp.margin)
      lambda <- clampAbs(lambda.raw, clamp.limit)
      epsilon <- checkFixVar(1 - lambda^2)

      if (!attr(epsilon, "ok")) {
        penalty <- sign(lambda.raw) * (abs(lambda.raw) - penalty.cfg$loading$limit)
        parTable[cond, "penalty"] <- parTable[cond, "penalty"] + penalty
      }

      # vals <- lambda * Xi[[lv]] + rnorm(N, mean = 0, sd = sqrt(epsilon))
      vals <- lambda * Xi[[lv]] + Rfast::Rnorm(N, m = 0, s = sqrt(epsilon), seed = rfast.seed())

      Inds[[ind]] <- vals
    }
  }

  for (lv in mode.b) {
    inds.lv <- indsLVs[[lv]]
    nind <- length(inds.lv)

    stopif(nind != 1, "Mode B is not available in MC-OrdPLSc (yet)!")
    Inds[[inds.lv]] <- Xi[[lv]]
  }

  Inds <- as.data.frame(Inds)
  All  <- cbind(Xi, Inds)
  Lv   <- All[lvs]
  Ov   <- All[ovs]

  list(
    all = All,
    ov  = Ov,
    lv  = Lv,
    is.admissible = is.admissible,
    penalty = parTable$penalty,
    parTable = parTable
  )
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


simulateDataMixedEffectsParTable <- function(parTable, N = 1e5,
                                             clusterSizes, clusterName,
                                             seed = NULL, .covtol = .95,
                                             check.hi.ord = FALSE) {
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
      limit = abs(.covtol),
      guard = 0,
      beta  = 10,
      scale = 1,
      max_penalty = 2
    ),
    loading = list(
      limit = 1,
      guard = 0.005,
      clamp = 0.001,
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
  randeff <- getRandomEffectLabels(parTable)

  ovs  <- setdiff(ovs, randeff)
  lvs  <- setdiff(lvs, randeff)
  xis  <- setdiff(xis, randeff)
  etas <- setdiff(etas, randeff)

  Psi.x <- diag(1, length(xis))
  dimnames(Psi.x) <- list(xis, xis)

  if (N < sum(clusterSizes)) {
    N <- sum(clusterSizes)
  } else {
    K <- floor(N / sum(clusterSizes))
    clusterSizes <- rep(clusterSizes, K)
    N <- sum(clusterSizes)
  }

  ncluster <- length(clusterSizes)
  cluster <- rep(seq_along(clusterSizes), clusterSizes)

  for (i in seq_along(xis)) {
    xi.i <- xis[[i]]

    for (j in seq_len(i - 1)) {
      xi.j <- xis[[j]]

      cond <- (
        (parTable$rhs == xi.i & parTable$op == "~~" & parTable$lhs == xi.j) |
        (parTable$rhs == xi.j & parTable$op == "~~" & parTable$lhs == xi.i)
      )

      p.ij <- parTable[cond, "est"][1]

      if (p.ij <= -.covtol || p.ij >= .covtol) {
        penalty <- smoothBoundaryPenalty(
          p.ij,
          limit = penalty.cfg$corr$limit,
          guard = penalty.cfg$corr$guard,
          beta  = penalty.cfg$corr$beta,
          scale = penalty.cfg$corr$scale,
          penalty.max = penalty.cfg$corr$max_penalty
        )

        if (penalty != 0)
          parTable[cond, "penalty"] <- parTable[cond, "penalty"] + penalty

        p.ij <- sign(p.ij) * abs(.covtol)
      }

      Psi.x[i, j] <- Psi.x[j, i] <- p.ij
    }
  }

  # Xi <- mvtnorm::rmvnorm(n = N, mean = rep(0, length(xis)), sigma = Psi.x)
  Psi.x.decomp <- tryCatch(
    chol(Psi.x),
    error = \(e) {
      tryCatch({
        diag(Psi.x) <- 1.01
        chol(Psi.x)
      }, error = \(e) diag2(Psi.x))
    }
  )

  Xi <- mvnfast::rmvn(n = N, mu = rep(0, length(xis)), sigma = Psi.x.decomp, isChol = TRUE)
  Xi <- as.data.frame(Rfast::standardise(Xi))
  colnames(Xi) <- xis


  undefIntTerms <- getIntTerms(parTable)
  elemsIntTerms <- stringr::str_split(undefIntTerms, pattern = ":")
  names(elemsIntTerms) <- undefIntTerms

  for (eta in etas) {

    randeff.eta <- randeff[startsWith(randeff, paste0(eta, "~"))]
    G <- diag(0, length(randeff.eta))
    dimnames(G) <- list(randeff.eta, randeff.eta)

    for (i in seq_along(randeff.eta)) {
      rand.i <- randeff.eta[[i]]

      for (j in seq_len(i)) {
        rand.j <- randeff.eta[[j]]

        cond <- (
          (parTable$rhs == rand.i & parTable$op == "~~" & parTable$lhs == rand.j) |
          (parTable$rhs == rand.j & parTable$op == "~~" & parTable$lhs == rand.i)
        )

        p.ij <- parTable[cond, "est"][1]

        if (p.ij <= -.covtol || p.ij >= .covtol) {
          penalty <- smoothBoundaryPenalty(
             p.ij,
             limit = penalty.cfg$corr$limit,
             guard = penalty.cfg$corr$guard,
             beta  = penalty.cfg$corr$beta,
             scale = penalty.cfg$corr$scale,
             penalty.max = penalty.cfg$corr$max_penalty
          )

          if (penalty != 0)
            parTable[cond, "penalty"] <- parTable[cond, "penalty"] + penalty

          p.ij <- sign(p.ij) * abs(.covtol)
        }

        G[i, j] <- G[j, i] <- p.ij
      }
    }

    G.decomp <- tryCatch(
      chol(G),
      error = \(e) {
        tryCatch({
          diag(G) <- diag(G) + 0.01
          chol(G)
        }, error = \(e) diag2(G))
      }
    )

    U <- mvnfast::rmvn(n = ncluster, mu = rep(0, NCOL(G)), sigma = G.decomp, isChol = TRUE)
    colnames(U) <- randeff.eta
    U.expanded <- U[cluster, , drop = FALSE]

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

    for (i in seq_len(NROW(predRows))) {
      row  <- predRows[i, ]
      beta <- row$est
      pred <- row$rhs

      # Fixed effect
      vals <- vals + beta * Xi[[pred]]

      # Random Effect
      par <- paste0(eta, "~", pred)
      if (par %in% colnames(U)) {
        u <- U.expanded[,par]
        vals <- vals + u * Xi[[pred]]
      }
    }

    # Random Intercept
    par <- paste0(eta, "~1")
    if (par %in% colnames(U))
      vals <- vals + U.expanded[,par]

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

    vals <- vals + Rfast::Rnorm(N, m = 0, s = sqrt(resvar), seed = rfast.seed())
    # vals <- vals + rnorm(N, mean = 0, sd = sqrt(resvar))

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

      if (lambda.penalty != 0)
        parTable[cond, "penalty"] <- parTable[cond, "penalty"] + lambda.penalty

      clamp.margin <- max(penalty.cfg$loading$clamp, penalty.cfg$loading$tol)
      clamp.limit <- max(0, penalty.cfg$loading$limit - clamp.margin)
      lambda <- clampAbs(lambda.raw, clamp.limit)
      epsilon <- checkFixVar(1 - lambda^2)

      if (!attr(epsilon, "ok")) {
        penalty <- sign(lambda.raw) * (abs(lambda.raw) - penalty.cfg$loading$limit)
        parTable[cond, "penalty"] <- parTable[cond, "penalty"] + penalty
      }

      # vals <- lambda * Xi[[lv]] + rnorm(N, mean = 0, sd = sqrt(epsilon))
      vals <- lambda * Xi[[lv]] + Rfast::Rnorm(N, m = 0, s = sqrt(epsilon), seed = rfast.seed())

      Inds[[ind]] <- vals
    }
  }

  for (lv in mode.b) {
    inds.lv <- indsLVs[[lv]]
    nind <- length(inds.lv)

    stopif(nind != 1, "Mode B is not available in MC-OrdPLSc (yet)!")
    Inds[[inds.lv]] <- Xi[[lv]]
  }

  Inds <- as.data.frame(Inds)
  All  <- cbind(Xi, Inds)
  Lv   <- All[lvs]
  Ov   <- All[ovs]

  list(
    all = All,
    ov  = Ov,
    lv  = Lv,
    is.admissible = is.admissible,
    penalty = parTable$penalty,
    parTable = parTable,
    cluster = matrix(cluster, nrow = N, dimnames = list(NULL, clusterName))
  )
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
