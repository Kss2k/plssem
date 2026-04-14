simulateDataParTable <- function(parTable, N = 1e5, seed = NULL, .covtol = .95) {
  if (!is.null(seed) && exists(".Random.seed")) .Random.seed.orig <- .Random.seed
  else                                          .Random.seed.orig <- NULL

  on.exit({
    if (!is.null(.Random.seed.orig)) .Random.seed <<- .Random.seed.orig
  })

  if (!is.null(seed))
    set.seed(seed)

  is.admissible <- TRUE
  checkFixVar <- function(v) {
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

  # Generate seed passed to Rfast::Rnorm. Passing seed=NULL does not work
  # If the user has used set.seed()
  rfast.seed <- floor(stats::runif(1L, min = 0, max = 9999999))
  rfast.seed <- NULL

  xis     <- getXis(parTable)
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
        parTable[cond, "penalty"] <- sign(p.ij) * (abs(p.ij) - abs(.covtol))
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
      formula <- formula(paste(
        eta, "~",
        paste0(predRows$rhs, collapse = " + ")
      ))

      Xi[[eta]] <- standardizeAtomic(vals)
      beta.hat <- coef(lm(formula, data = Xi))

      beta.y <- beta.hat[parTable[cond, "rhs"]]
      beta.x <- parTable[cond, "est"]

      parTable[cond, "penalty"] <- beta.x - beta.y
    }

    vals <- vals + Rfast::Rnorm(N, m = 0, s = sqrt(resvar), seed = rfast.seed)
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

      lambda <- parTable[cond, "est"]
      epsilon <- checkFixVar(1 - lambda^2)

      if (!attr(epsilon, "ok")) {
        penalty <- sign(lambda) * (abs(lambda) - 1)
        parTable[cond, "penalty"] <- penalty
      }

      # vals <- lambda * Xi[[lv]] + rnorm(N, mean = 0, sd = sqrt(epsilon))
      vals <- lambda * Xi[[lv]] + Rfast::Rnorm(N, m = 0, s = sqrt(epsilon), seed = rfast.seed)

      Inds[[ind]] <- vals
    }
  }

  for (lv in mode.b) {
    stop("Mode B is not available in MC-OrdPLSc (yet)!")
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
