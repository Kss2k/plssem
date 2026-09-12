simulateDataParTable <- function(parTable,
                                 N                       = 1e5,
                                 seed                    = NULL,
                                 tol                     = 1e-3,
                                 .cortol                 = .95,
                                 .varguard               = 5 * tol,
                                 check.hi.ord            = FALSE,
                                 clusterSizes            = NULL,
                                 clusterName             = NULL,
                                 standardize             = FALSE,
                                 full                    = FALSE,
                                 cut                     = FALSE,
                                 collect.empirical.vpars = FALSE,
                                 innovations             = NULL,
                                 return.innovations      = FALSE,
                                 compiled.info           = NULL) {

  if (!is.null(seed) && exists(".Random.seed")) .Random.seed.orig <- .Random.seed
  else                                          .Random.seed.orig <- NULL

  on.exit({
    if (!is.null(.Random.seed.orig)) .Random.seed <<- .Random.seed.orig
  })

  if (!is.null(seed))
    set.seed(seed)

  use.innovations <- return.innovations || !is.null(innovations)
  supplied.innovations <- !is.null(innovations)

  pls_stopif(supplied.innovations &&
    (!inherits(innovations, "Innovations") || !innovations$initialized),
    "`innovations` must be an initialized `Innovations` object!"
  )

  if (is.null(innovations)) innovations.out <- innovationStruct()
  else innovations.out <- innovations

  drawInnovations <- function(key, nrow, ncol = 1) {

    if (supplied.innovations) {

      z <- innovations$blocks[[key]]

      pls_stopif(is.null(z),
        "Innovation block `", key, "` is missing."
      )

      pls_stopif(NROW(z) != nrow || NCOL(z) != ncol,
        "Innovation block `", key, "` has incompatible dimensions."
      )

    } else {
      if (ncol > 1) z <- matrix(frnorm(nrow * ncol), nrow = nrow, ncol = ncol)
      else          z <- frnorm(nrow)

      innovations.out$blocks[[key]] <<- z
    }

    z
  }

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
  if (is.null(compiled.info)) {
    xis           <- getXis(parTable, isLV = !check.hi.ord)
    etas          <- getSortedEtas(parTable)
    mode.a        <- getReflectiveLVs(parTable)
    mode.b        <- getFormativeLVs(parTable)
    lvs           <- unique(c(mode.a, mode.b))
    indsLVs       <- getIndsLVs(parTable, lVs = lvs)
    ovs           <- getOVs(parTable)
    mixed         <- !is.null(clusterSizes) && !is.null(clusterName)
    randeff       <- NULL
    intTerms      <- getIntTerms(parTable)
    undefIntTerms <- intTerms
    elemsIntTerms <- stats::setNames(
      stringr::str_split(intTerms, pattern = ":"),
      nm = intTerms 
    )

  } else {
    xis           <- compiled.info$xis
    etas          <- compiled.info$etas
    mode.a        <- compiled.info$mode.a
    mode.b        <- compiled.info$mode.b
    lvs           <- compiled.info$lvs
    indsLVs       <- compiled.info$indsLVs
    ovs           <- compiled.info$ovs
    mixed         <- compiled.info$mixed
    randeff       <- compiled.info$randeff
    intTerms      <- compiled.info$intTerms
    undefIntTerms <- compiled.info$intTerms
    elemsIntTerms <- compiled.info$elemsIntTerms
  }

  empirical.vpars <- numeric(0L)

  if (mixed) {

    if (is.null(compiled.info)) {
      randeff <- getRandomEffectLabels(parTable)

      ovs  <- setdiff(ovs, randeff)
      lvs  <- setdiff(lvs, randeff)
      xis  <- setdiff(xis, randeff)
      etas <- setdiff(etas, randeff)
    }

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

  if (use.innovations) z.xi <- drawInnovations("exogenous", N, length(xis))
  else z.xi <- NULL

  xiDraw <- rmvnSafe(N, res$mat, innovations = z.xi)
  is.admissible <- is.admissible && xiDraw$is.admissible

  Xi <- as.data.frame(Rfast::standardise(xiDraw$x))
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

        parTable <- res$parTable

        if (use.innovations) {
          z.random <- drawInnovations(
            key = paste0("random-effects:", eta),
            nrow = ncluster,
            ncol = length(randeff.eta)
          )
        } else {
          z.random <- NULL
        }

        uDraw <- rmvnSafe(ncluster, res$mat, innovations = z.random)
        is.admissible <- is.admissible && uDraw$is.admissible

        U <- uDraw$x
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

    if (is.finite(projvar) && projvar >= 1 - .varguard && NROW(predRows) > 0) {
      # Get bounds for the fixed-effect `beta` such that its own implied
      # variance stays within the guard.
      preds  <- predRows$rhs
      beta.x <- predRows[, "est"]

      Sigma <- Rfast::cova(as.matrix(Xi[preds]))

      # Only fall back to `.varguard` when the actual budget isn't positive
      rawMaxvar <- (1 - .varguard) - stats::var(vals.random)
      if (!full && is.finite(rawMaxvar) && rawMaxvar > 0) maxvar <- rawMaxvar
      else maxvar <- .varguard

      beta.y <- projectBetaOntoConstrainedEllipsoid(
        beta   = beta.x,
        Sigma  = Sigma,
        maxvar = maxvar
      )

      for (i in seq_along(beta.x)) {
        if (!is.finite(beta.y[[i]]) || beta.y[[i]] == beta.x[[i]]) next

        # Tighten whichever side `beta.x` needs to move toward
        # to reach the projected boundary point.
        if (beta.x[[i]] > beta.y[[i]]) {
          lim <- beta.y[[i]] - tol
          parTable[cond, "upper"][i] <- min(parTable[cond, "upper"][i], lim)

        } else {
          lim <- beta.y[[i]] + tol
          parTable[cond, "lower"][i] <- max(parTable[cond, "lower"][i], lim)
        }
      }
    }

    if (is.finite(projvar) && projvar > 1 - .varguard && length(randeff.eta)) {
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

    # In `reduced` mode (or when eta has no residual covariance) the residual is
    # independent with var(eta) = resvar. In full mode it's drawn conditional on
    # the prior disturbances such that its covariance with each prior node
    # equals the specified residual covariance: with target cross-covariances
    # `a` and realised noise covariance `M`, the regression `beta = M^-1 a`
    # gives realised Cov(zeta, noise) = a exactly. The fresh part is then sized
    # so that `vals + zeta` has unit variance (`Var(fresh) = 1 - Var(vals + cmean)`),
    # which keeps the latent variable standardized even when the residual covaries
    # with one of its own predictors (then `cmean` is correlated with `vals`).
    # When it does not, this reduces to `resvar - a' beta`.
    if (full) a <- vapply(dnames, FUN.VALUE = numeric(1L), FUN = \(v) rescov(v, eta))
    else      a <- 0

    if (full && any(a != 0)) {
      M    <- Rfast::cova(disturbances)
      beta <- tryCatch(as.vector(solve(M, a)), error = function(...) {
        # fails, so set it to inadmissible
        is.admissible <<- FALSE
        numeric(length(a))
      })

      cmean   <- as.vector(disturbances %*% beta)
      vcmean  <- stats::var(vals + cmean)
      condvar <- checkFixVar(1 - vcmean)

      if (is.finite(vcmean) && vcmean > 1 - tol) {
        # `Var(vals + cmean(a))` is quadratic in `a`, and - unlike the
        # `beta` case above - the resulting ellipsoid is centred at `-c`
        # Projecting `a` under Euclidean distance is the same
        # `projectBetaOntoConstrainedEllipsoid()` problem applied to the
        # shifted point `a + c`, then shifted back.
        Minv <- tryCatch(solve(M), error = \(...) NULL)

        # `a.bound == 0` is overloaded below to also mean "no update needed",
        # so the two cases that deliberately require `a` to be exactly zero
        # need their own flag rather than relying on the value alone.
        force.zero <- FALSE

        if (is.null(Minv)) {
          is.admissible <- FALSE
          force.zero <- TRUE
          a.bound <- rep(0, length(a))

        } else {
          cVec  <- as.vector(stats::cov(disturbances, vals))
          limit <- 1 - tol
          shiftMaxvar <- limit - stats::var(vals) + c(t(cVec) %*% Minv %*% cVec)

          if (is.finite(shiftMaxvar) && shiftMaxvar > 0) {
            aShifted.proj <- projectBetaOntoConstrainedEllipsoid(
              beta   = a + cVec,
              Sigma  = Minv,
              maxvar = shiftMaxvar
            )

            # Keep the sign: the feasible ellipsoid for `a` is centred at
            # `-cVec`, not at 0, so this boundary point is generally not
            # symmetric around zero. Taking `abs()` here would discard that
            # asymmetry and let the (potentially infeasible) mirror-image
            # side back in.
            a.bound <- aShifted.proj - cVec

          } else {
            # If even `a = 0` (resvar=0) would violate the variance constraint,
            # we constrain a to 0 going forward.
            force.zero <- TRUE
            a.bound <- rep(0, length(a))
          }
        }

        names(a.bound) <- names(a)

        for (v in names(a.bound)) {
          idx <- which(
            parTable$op == "~~" &
            parTable$lhs != parTable$rhs & (
              (parTable$lhs == eta & parTable$rhs == v) |
              (parTable$lhs == v   & parTable$rhs == eta)
            )
          )

          if (!length(idx)) next

          if (force.zero) {
            parTable[idx, "lower"] <- 0
            parTable[idx, "upper"] <- 0
            next
          }

          if (!is.finite(a.bound[[v]]) || a.bound[[v]] == a[[v]]) next

          # which side to tighten is determined by which direction the current
          # value needs to move to reach the projected boundary.
          if (a[[v]] > a.bound[[v]]) {
            lim <- a.bound[[v]] - tol
            parTable[idx, "upper"] <- pmin(parTable[idx, "upper"], lim)

          } else {
            lim <- a.bound[[v]] + tol
            parTable[idx, "lower"] <- pmax(parTable[idx, "lower"], lim)
          }
        }
      }

      if (use.innovations)
        z <- drawInnovations(key = paste0("structural:", eta), nrow = N, ncol = 1)
      else
        z <- Rfast::Rnorm(N, m = 0, s = 1, seed = rfast.seed())

      zeta <- cmean + sqrt(condvar) * z

    } else {
      if (use.innovations)
        z <- drawInnovations(key = paste0("structural:", eta), nrow = N, ncol = 1)
      else
        z <- Rfast::Rnorm(N, m = 0, s = 1, seed = rfast.seed())

      zeta <- sqrt(resvar) * z
    }

    vals <- vals + zeta

    if (collect.empirical.vpars && !full)
      empirical.vpars[[paste0(eta, "~~", eta)]] <- resvar

    if (full) {
      disturbances <- cbind(disturbances, zeta)
      dnames <- c(dnames, eta)
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
      if (use.innovations)
        z <- drawInnovations(key = paste0("indicator:", ind), nrow = N, ncol = 1)
      else
        z <- Rfast::Rnorm(N, m = 0, s = 1, seed = rfast.seed())

      eps <- sqrt(epsilon) * z
      vals <- lambda * Xi[[lv]] + eps

      if (standardize)
        vals <- (vals - mean(vals)) / stats::sd(vals)

      Inds[[ind]] <- vals

      if (collect.empirical.vpars) {
        empirical.vpars <- c(
          empirical.vpars,
          stats::setNames(epsilon, nm = paste0(ind, "~~", ind))
        )
      }
    }
  }
  
  if (collect.empirical.vpars) {
    if (full && length(etas)) {
      eta.disturbances <- disturbances[, match(etas, dnames), drop = FALSE]
      eta.variances <- Rfast::colVars(eta.disturbances)
      names(eta.variances) <- paste0(etas, "~~", etas)
      empirical.vpars[names(eta.variances)] <- eta.variances

      rows <- which(
        rescovRows$lhs %in% dnames &
        rescovRows$rhs %in% dnames &
        (rescovRows$lhs %in% etas | rescovRows$rhs %in% etas)
      )

      if (length(rows)) {
        values <- vapply(rows, FUN.VALUE = numeric(1L), FUN = \(row) {
          lhs <- disturbances[, match(rescovRows$lhs[[row]], dnames)]
          rhs <- disturbances[, match(rescovRows$rhs[[row]], dnames)]
          stats::cov(lhs, rhs)
        })
        names(values) <- paste0(rescovRows$lhs[rows], "~~", rescovRows$rhs[rows])
        empirical.vpars[names(values)] <- values
      }
    }

    if (length(intTerms)) {
      X <- as.matrix(Xi[, xis, drop = FALSE])
      I <- as.matrix(Xi[, intTerms, drop = FALSE])
      Sigma.ix <- crossprod(I, X) / (N - 1L)

      int.names <- rep(intTerms, times = length(xis))
      xi.names <- rep(xis, each = length(intTerms))
      values <- c(Sigma.ix)

      empirical.vpars[paste0(int.names, "~~", xi.names)] <- values
      empirical.vpars[paste0(xi.names, "~~", int.names)] <- values
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

  if (use.innovations)
    innovations.out$initialized <- TRUE

  if (is.null(compiled.info)) {
    compiled.info <- list(
      xis           = xis,  
      etas          = etas,
      mode.a        = mode.a,
      mode.b        = mode.b,
      lvs           = lvs,
      indsLVs       = indsLVs,
      ovs           = ovs,  
      mixed         = mixed,
      randeff       = randeff,
      intTerms      = intTerms,
      elemsIntTerms = elemsIntTerms
    )
  }

  list(
    all             = All,
    ov              = Ov,
    lv              = Lv,
    is.admissible   = is.admissible,
    lower           = parTable$lower,
    upper           = parTable$upper,
    parTable        = parTable,
    cluster         = clusterMat,
    empirical.vpars = empirical.vpars,
    innovations     = if (use.innovations) innovations.out else NULL,
    compiled.info   = compiled.info
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


rmvnSafe <- function(n, mat, innovations = NULL) {
  is.admissible <- TRUE

  decomp <- tryCatch(
    chol(mat),
    error = function(e) {
      is.admissible <<- FALSE

      tryCatch({
        # ridge solve?
        diag(mat) <- diag(mat) + 0.01
        chol(mat)
      }, error = \(e) diag2(mat))
    }
  )

  if (is.null(innovations)) {
    x <- mvnfast::rmvn(
      n  = n,
      mu = rep(0, NCOL(mat)),
      sigma = decomp,
      isChol = TRUE
    )

  } else {
    pls_stopif(
      NROW(innovations) != n || NCOL(innovations) != NCOL(mat),
      "Multivariate-normal innovations have incompatible dimensions."
    )

    x <- innovations %*% decomp
  }

  list(
    x             = x,
    is.admissible = is.admissible
  )
}


innovationStruct <- function() {
  out <- list(
    blocks = list(),
    initialized = FALSE
  )

  class(out) <- "Innovations"
  out 
}


perturbInnovations <- function(innovations, scale) {
  pls_stopif(
    !inherits(innovations, "Innovations") || !innovations$initialized,
    "`innovations` must be an initialized `Innovations` object!"
  )

  pls_stopif(
    length(scale) != 1L || !is.finite(scale) || scale <= 0 || scale > 1,
    "`scale` must be in (0, 1]."
  )

  persistence <- sqrt(1 - scale^2)
  out <- innovations

  out$blocks <- lapply(
    X = out$blocks,
    FUN = \(X) X * persistence + scale * frnorm(length(X))
  )

  out
}


frnorm <- function(n, mean = 0, sd = 1) {
  Rfast::Rnorm(n, m = mean, s = sd, seed = rfast.seed())
}


projectBetaOntoConstrainedEllipsoid <- function(beta, Sigma, maxvar, tol = 1e-8, max.delta = 1e12) {
  currentvar <- c(t(beta) %*% Sigma %*% beta)
  if (currentvar <= maxvar)
    return(beta)

  eig <- eigen(Sigma, symmetric = TRUE)
  d   <- pmax(eig$values, 0) # guard against numerical noise below zero
  y   <- as.vector(t(eig$vectors) %*% beta)

  g <- \(delta) sum(d * y^2 / (1 + delta * d)^2) - maxvar

  delta.high <- 1
  while (g(delta.high) > 0 && delta.high < max.delta)
    delta.high <- 10 * delta.high

  # solve for delta
  solved <- stats::uniroot(
    f     = g,
    lower = 0,
    upper = delta.high,
    tol = tol
  )

  delta <- solved$root

  as.vector(eig$vectors %*% (y / (1 + delta * d)))
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
