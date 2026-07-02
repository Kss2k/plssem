simulateDataParTable <- function(parTable,
                                 N            = 1e5,
                                 seed         = NULL,
                                 .cortol      = .95,
                                 tol          = 1e-4,
                                 check.hi.ord = FALSE,
                                 clusterSizes = NULL,
                                 clusterName  = NULL,
                                 standardize  = FALSE,
                                 full         = FALSE,
                                 cut          = FALSE) {
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

  if (check.hi.ord)
    parTable <- highOrdMeasrAsStructParTable(parTable)

  # bounds
  parTable$lower <- -Inf
  parTable$upper <- +Inf

  # info
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
    U           <- NULL
    U.expanded  <- NULL
    randeff.eta <- character(0L)

    if (mixed) {
      randeff.eta <- randeff[startsWith(randeff, paste0(eta, "~"))]

      if (length(randeff.eta)) {

        res <- buildCovMat(
          vars          = randeff.eta,
          parTable      = parTable,
          .cortol       = .cortol,
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

    vals.fixed <- numeric(N)
    vals.random <- numeric(N)

    # Random Intercept
    par <- paste0(eta, "~1")
    if (par %in% colnames(U))
      vals.random <- vals.random + U.expanded[,par]

    for (i in seq_len(NROW(predRows))) {
      row  <- predRows[i, ]
      beta <- row$est
      pred <- row$rhs

      # Fixed effect
      vals.fixed <- vals.fixed + beta * Xi[[pred]]

      # Random Effect
      par <- paste0(eta, "~", pred)
      if (par %in% colnames(U))
        vals.random <- vals.random + U.expanded[,par] * Xi[[pred]]
    }

    vals <- vals.fixed + vals.random
    projvar <- stats::var(vals)
    resvar  <- checkFixVar(1 - projvar)

    if (is.finite(projvar) && projvar > 1 - tol) {
      # Get bounds for (fixed) beta (i.e., resvar=0)
      beta.x <- predRows[,"est"]
      beta.y <- beta.x / sqrt(max(projvar, tol))

      parTable[cond, "lower"] <- pmin(-abs(beta.y) + tol, -tol)
      parTable[cond, "upper"] <- pmax(+abs(beta.y) - tol, +tol)
    }

    if (is.finite(projvar) && projvar > 1 - tol && length(randeff.eta)) {
      v0 <- stats::var(vals.fixed)
      v1 <- stats::var(vals.random)
      vc <- stats::cov(vals.fixed, vals.random)

      if (is.finite(v1) && v1 > 0) {
        limit <- 1 - tol
        roots <- polyroot(c(v0 - limit, 2 * vc, v1))
        roots <- Re(roots[abs(Im(roots)) < 1e-7])
        roots <- roots[is.finite(roots) & roots >= 0 & roots <= 1]
        scale <- if (length(roots)) max(roots) else 0
        scale2 <- scale^2

        re.rows <- which(
          parTable$op == "~~" &
          parTable$lhs %in% randeff.eta &
          parTable$rhs %in% randeff.eta
        )

        for (row in re.rows) {
          est <- parTable[row, "est"]
          if (!is.finite(est)) next

          if (parTable$lhs[[row]] == parTable$rhs[[row]]) {
            lim <- max(.Machine$double.eps, abs(est) * scale2 - tol)
            parTable[row, "lower"] <- pmax(parTable[row, "lower"], .Machine$double.eps)
            parTable[row, "upper"] <- pmin(parTable[row, "upper"], lim)

          } else {
            lim <- max(0, abs(est) * scale2 - tol)
            parTable[row, "lower"] <- pmax(parTable[row, "lower"], -lim)
            parTable[row, "upper"] <- pmin(parTable[row, "upper"],  lim)
          }
        }
      }
    }


    # Disturbance. In `reduced` mode (or when this eta has no residual
    # covariance) the disturbance is independent with variance `resvar`. In
    # `full` mode it is drawn conditional on the prior disturbances so that its
    # covariance with each prior node equals the specified residual covariance:
    # with target cross-covariances `a` and realised noise covariance `M`, the
    # regression `beta = M^-1 a` gives realised Cov(zeta, noise) = a exactly.
    # The fresh part is then sized so that `vals + zeta` has unit variance --
    # `Var(fresh) = 1 - Var(vals + cmean)` -- which keeps the latent variable
    # standardized even when the residual covaries with one of its own
    # predictors (then `cmean` is correlated with `vals`). When it does not,
    # this reduces to `resvar - a' beta`.
    if (full) a <- vapply(dnames, FUN.VALUE = numeric(1L), FUN = \(v) rescov(v, eta))
    else      a <- 0

    if (full && any(a != 0)) {
      M    <- Rfast::cova(disturbances)
      beta <- tryCatch(as.vector(solve(M, a)), error = \(...) numeric(length(a)))

      cmean   <- as.vector(disturbances %*% beta)
      vcmean  <- stats::var(vals + cmean)
      condvar <- checkFixVar(1 - vcmean)

      if (is.finite(vcmean) && vcmean > 1 - tol) {
        scale <- sqrt(max(0, (1 - tol) / vcmean))
        a.bound <- abs(a) * scale

        for (v in names(a.bound)) {
          if (!is.finite(a.bound[[v]]) || a.bound[[v]] == 0) next

          lim <- max(0, a.bound[[v]] - tol)
          idx <- which(
            parTable$op == "~~" &
            parTable$lhs != parTable$rhs &
            (
              (parTable$lhs == eta & parTable$rhs == v) |
              (parTable$lhs == v   & parTable$rhs == eta)
            )
          )

          if (!length(idx)) next

          parTable[idx, "lower"] <- pmax(parTable[idx, "lower"], -lim)
          parTable[idx, "upper"] <- pmin(parTable[idx, "upper"],  lim)
        }
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

      lambda <- parTable[cond, "est"]
      epsilon <- checkFixVar(1 - lambda^2)

      parTable[cond, "lower"] <- -1 + tol
      parTable[cond, "upper"] <-  1 - tol

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

  if (cut) {
    # Create ordinal variables from thresholds. This is set to FALSE in the
    # MC-PLS estimation, as we cut by the observed proportions instead.
    # This is however useful in other circumstances (e.g., mcpls_loglik())

    thrvars <- unique(parTable[parTable$op == "|", "lhs"])
    for (thrvar in thrvars) {
      # get thresholds
      tau <- parTable[parTable$op == "|" & parTable$lhs == thrvar, "est"]

      # cut
      cont <- All[[thrvar]]
      ord  <- cut(cont, breaks = c(-Inf, sort(tau), Inf), labels = FALSE)

      # replace continous values with ordinal categories
      All[[thrvar]] <- ord
    }
  }

  Lv <- All[lvs]
  Ov <- All[ovs]

  list(
    all           = All,
    ov            = Ov,
    lv            = Lv,
    is.admissible = is.admissible,
    lower         = parTable$lower,
    upper         = parTable$upper,
    parTable      = parTable,
    cluster       = clusterMat
  )
}


buildCovMat <- function(vars, parTable, .cortol, unitVariances = FALSE) {
  mat <- diag(if (unitVariances) 1 else 0, length(vars))
  dimnames(mat) <- list(vars, vars)

  if (!unitVariances) for (v in vars) {
    cond.v <- parTable$rhs == v & parTable$lhs == v & parTable$op == "~~"
    v.i    <- parTable[cond.v, "est"]

    if (v.i <= 0)
      v.i <- .Machine$double.eps
    
    parTable[cond.v, "lower"] <- .Machine$double.eps
    parTable[cond.v, "upper"] <- Inf

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

      if (r.ij <= -.cortol || r.ij >= .cortol)
        p.ij <- denom * sign(p.ij) * abs(.cortol)

      parTable[cond, "upper"] <-  abs(.cortol) * denom
      parTable[cond, "lower"] <- -abs(.cortol) * denom

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
