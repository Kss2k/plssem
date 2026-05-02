# Module-level helpers shared across plslmer, fastMixedModel, meanDiagZGZt.
# The INTR_OP substitution avoids lme4 forming its own interaction terms from ':'.

lmerSafeIntr <- function(x) {
  stringr::str_replace_all(x, stringr::fixed(":"), "__INTR__")
}

lmerRestoreIntr <- function(x) {
  x <- stringr::str_replace_all(x, stringr::fixed("__INTR__"), ":")
  stringr::str_replace_all(x, stringr::fixed("`"), "")
}

lmerGetNames <- function(lhs, nm) {
  rhs <- stringr::str_replace_all(nm, "\\(Intercept\\)", "1")
  rhs <- lmerRestoreIntr(rhs)
  paste0(lhs, "~", rhs)
}

lmerFixVecNames <- function(vec, dep) {
  if (!is.null(names(vec)))
    names(vec) <- lmerGetNames(dep, names(vec))
  vec
}

lmerFixMatNames <- function(mat, dep, cols = TRUE, rows = TRUE) {
  if (!is.null(rownames(mat)) && rows)
    rownames(mat) <- lmerGetNames(dep, rownames(mat))
  if (!is.null(colnames(mat)) && cols)
    colnames(mat) <- lmerGetNames(dep, colnames(mat))
  mat
}


plslmer <- function(plsModel, fast = FALSE) {
  lme4.syntax <- plsModel@info$lme4.syntax
  cluster     <- plsModel@info$cluster
  consistent  <- plsModel@info$consistent

  stopif(!is.character(lme4.syntax), "`lme4.syntax` must be a character vector!")

  stopif(length(cluster) != 1L || !is.character(cluster),
         "`cluster` must be a character string of length 1. If lme4.syntax is provided!")

  fit.c <- plsModel@fitConsistent
  fit.u <- plsModel@fitUncorrected

  selectGamma   <- plsModel@matrices$select$gamma
  # Full (corrected) covariance/correlation matrix for the structural model.
  # Unlike `fitCov`, this retains the total covariances between endogenous and
  # exogenous variables.
  Cov.lv        <- fit.c$fitC
  Q             <- fit.c$Q
  Correction    <- fit.c$fitStructural / fit.u$fitStructural
  CorrectionCov <- fit.c$fitCov / fit.u$fitCov

  Xf <- as.data.frame(plsModel@factorScores)
  Xx <- as.data.frame(plsModel@data)
  Xc <- as.data.frame(attr(plsModel@data, "cluster"))
  X  <- cbind(Xf, Xx, Xc)

  # Mean-center interaction/quadratic terms (e.g., X:Z and X:X).
  intTerms <- plsModel@info$intTermNames
  if (!is.null(intTerms) && length(intTerms)) {
    intTerms <- intersect(intTerms, colnames(X))

    for (term in intTerms) {
      if (is.numeric(X[[term]]))
        X[[term]] <- X[[term]] - mean(X[[term]], na.rm = TRUE)
    }
  }

  # We don't want lmer to form the interaction terms on its own, since we want
  # them to be mean centered, matching the assumptions of our correction procedure.
  # We therefore replace the ':' operator and use explicitly formed interaction
  # terms, renaming any columns containing ':' before fitting.
  colsWithColon     <- colnames(X)[grepl(":", colnames(X), fixed = TRUE)]
  colsWithColonSafe <- lmerSafeIntr(colsWithColon)

  X.safe <- X
  if (length(colsWithColon)) {
    colnames(X.safe)[match(colsWithColon, colnames(X.safe))] <- colsWithColonSafe
  }

  makeSafeFormula <- function(line) {
    if (!length(colsWithColon))
      return(line)

    ord <- order(nchar(colsWithColon), decreasing = TRUE)
    out <- line

    for (i in ord) {
      key  <- colsWithColon[[i]]
      safe <- colsWithColonSafe[[i]]

      out <- stringr::str_replace_all(out, stringr::fixed(paste0("`", key, "`")), safe)
      out <- stringr::str_replace_all(out, stringr::fixed(key), safe)
    }

    out
  }

  FITS    <- list()
  FIXEF   <- list()
  COEF    <- list()
  VCOV    <- list()
  VARCORR <- list()
  SIGMA   <- list()

  for (line in lme4.syntax) {
    line.safe     <- makeSafeFormula(line)
    formula.full  <- stats::formula(line.safe)
    formula.fixed <- reformulas::nobars(formula.full)
    bars          <- reformulas::findbars(formula.full)

    fterms <- stats::terms(formula.fixed)
    vars   <- attr(fterms, "variables")
    dep    <- lmerRestoreIntr(as.character(vars[[2L]]))
    indep  <- colnames(selectGamma)[selectGamma[, dep]]
    c.safe <- lmerSafeIntr(cluster)

    # --- Estimation (branches on fast) -----------------------------------

    if (!fast) {
      lmerFit <- lme4::lmer(
        formula = line.safe,
        data    = X.safe,
        control = lme4::lmerControl(calc.derivs = FALSE)
      )

      fixefRaw   <- lme4::fixef(lmerFit)
      vcovRaw    <- as.matrix(stats::vcov(lmerFit))
      coefFit    <- stats::coef(lmerFit)
      varCorrFit <- as.list(lme4::VarCorr(lmerFit))

    } else {
      lmerFit <- NULL

      y          <- as.numeric(X.safe[[dep]])
      Xmat       <- stats::model.matrix(formula.fixed, data = X.safe)
      grp_factor <- factor(X.safe[[c.safe]])
      grp        <- as.integer(grp_factor)
      J          <- nlevels(grp_factor)
      grp_levels <- levels(grp_factor)

      bar_idx <- which(vapply(bars, function(b) {
        lmerRestoreIntr(deparse(b[[3L]])) == cluster
      }, logical(1L)))
      if (!length(bar_idx)) bar_idx <- 1L
      bar    <- bars[[bar_idx[[1L]]]]
      f_re   <- stats::as.formula(paste0("~", deparse(bar[[2L]])))
      Zsmall <- stats::model.matrix(f_re, data = X.safe)
      k      <- ncol(Zsmall)
      p      <- ncol(Xmat)

      result <- fastMixedModel(y = y, X = Xmat, Zsmall = Zsmall,
                               grp = grp, J = J, grp_levels = grp_levels)

      names(result$beta)         <- colnames(Xmat)
      rownames(result$Sigma_u)   <- colnames(Zsmall)
      colnames(result$Sigma_u)   <- colnames(Zsmall)
      rownames(result$vcov_beta) <- colnames(Xmat)
      colnames(result$vcov_beta) <- colnames(Xmat)

      # Build cluster-level coef matrix (rows = clusters, cols = fixed effects)
      coef_mat <- matrix(rep(result$beta, J), nrow = J, ncol = p, byrow = TRUE,
                         dimnames = list(grp_levels, colnames(Xmat)))
      for (ri in seq_len(k)) {
        cn <- colnames(Zsmall)[[ri]]
        if (cn %in% colnames(coef_mat))
          coef_mat[, cn] <- coef_mat[, cn] + result$u[, ri]
      }

      fixefRaw              <- result$beta
      vcovRaw               <- result$vcov_beta
      coefFit               <- list()
      coefFit[[c.safe]]     <- coef_mat
      varCorrFit            <- list()
      varCorrFit[[c.safe]]  <- result$Sigma_u
    }

    # --- Shared: apply naming, compute corrections, rename coef/vcorr ----

    fixefFit <- lmerFixVecNames(fixefRaw, dep)
    vcovFit  <- lmerFixMatNames(vcovRaw, dep)

    params <- names(fixefFit)
    split  <- stringr::str_split_fixed(params, pattern = "~", n = 2L)
    lhs    <- split[, 1L]
    rhs    <- split[, 2L]

    correctionTerms <- stats::setNames(rep(1, length(fixefFit)), nm = params)
    for (i in seq_along(lhs)) {
      lhs.i <- lhs[[i]]
      rhs.i <- rhs[[i]]

      if (rhs.i == "1") {
        if (!is.null(Q) && lhs.i %in% names(Q)) {
          term <- 1 / Q[[lhs.i]]

          if (is.finite(term) && !is.na(term) && !is.nan(term))
            correctionTerms[[i]] <- term
        }

        next

      } else if (!lhs.i %in% colnames(Correction) || !rhs.i %in% rownames(Correction)) {
        warning("Unable to identify correction term for ", params[[i]], "!")
        next
      }

      term <- Correction[rhs.i, lhs.i]

      if (is.na(term) || is.nan(term)) {
        warning("Correction term for ", params[[i]], " is NaN!")
        next
      }

      correctionTerms[[i]] <- term
    }

    DCorrectionTerms <- diag(correctionTerms)
    dimnames(DCorrectionTerms) <- list(params, params)

    if (consistent) {
      fixefFit <- correctionTerms * fixefFit
      vcovFit  <- DCorrectionTerms %*% vcovFit %*% DCorrectionTerms
    }

    for (c in cluster) {
      c.safe.c <- lmerSafeIntr(c)

      coefFit[[c]]    <- lmerFixMatNames(as.matrix(coefFit[[c.safe.c]]), dep = dep, rows = FALSE)
      varCorrFit[[c]] <- lmerFixMatNames(varCorrFit[[c.safe.c]], dep = dep)

      if (c.safe.c != c) {
        coefFit[[c.safe.c]]    <- NULL
        varCorrFit[[c.safe.c]] <- NULL
      }

      if (consistent) {
        vcPars       <- rownames(varCorrFit[[c]])
        vcCorrection <- DCorrectionTerms[vcPars, vcPars, drop = FALSE]

        coefFit[[c]]    <- coefFit[[c]] %*% diag(correctionTerms[colnames(coefFit[[c]])])
        varCorrFit[[c]] <- vcCorrection %*% varCorrFit[[c]] %*% vcCorrection
      }

      attr(varCorrFit[[c]], "stddev")      <- sqrt(pmax(0, diag(varCorrFit[[c]])))
      attr(varCorrFit[[c]], "correlation") <- tryCatchNA(cov2cor(varCorrFit[[c]]))
    }

    FITS[[dep]]    <- lmerFit
    COEF[[dep]]    <- coefFit
    VCOV[[dep]]    <- vcovFit
    FIXEF[[dep]]   <- fixefFit
    VARCORR[[dep]] <- varCorrFit
    SIGMA[[dep]]   <- getSigmaFromVarCorr(
      fit           = lmerFit,
      beta          = fixefFit,
      varCorr       = varCorrFit,
      dep           = dep,
      indep         = indep,
      CorrectionCov = CorrectionCov,
      Cov.lv        = Cov.lv
    )
  }

  valuesFixef <- unlist(unname(FIXEF))
  valuesSigma <- unlist(unname(SIGMA))
  values      <- c(valuesFixef, valuesSigma)

  list(
    pls    = plsModel,
    fits   = FITS,
    coef   = COEF,
    vcov   = VCOV,
    fixef  = FIXEF,
    vcorr  = VARCORR,
    values = values
  )
}


# Iterative GLS estimator for the linear mixed model y = Xb + Zu + e.
# Uses per-group Woodbury updates (O(J*k^3) per iteration) rather than
# a full n x n matrix factorization. Variance components are updated via
# a moment estimator (biased but fast); n_iter = 2 is sufficient for
# MC-PLS since the outer SA loop corrects the resulting bias.
fastMixedModel <- function(y, X, Zsmall, grp, J, grp_levels = NULL, n_iter = 2L) {
  n <- length(y)
  p <- ncol(X)
  k <- ncol(Zsmall)

  # OLS initialisation
  XtX      <- crossprod(X)
  Xty      <- drop(crossprod(X, y))
  beta_hat <- tryCatch(
    solve(XtX, Xty),
    error = function(e) solve(XtX + diag(1e-8, p), Xty)
  )

  resid     <- y - drop(X %*% beta_hat)
  var_resid <- max(mean(resid^2), 1e-8)

  # Split initial variance between within- and between-group components
  sigma2_e <- var_resid * 0.5
  Sigma_u  <- diag(var_resid * 0.5 / k, k)

  u_hat   <- matrix(0, J, k)
  XtVinvX <- XtX
  XtVinvy <- Xty

  for (iter in seq_len(n_iter)) {
    Sigma_u_inv <- tryCatch(
      solve(Sigma_u),
      error = function(e) diag(1 / pmax(diag(Sigma_u), 1e-10), k)
    )

    XtVinvX <- matrix(0, p, p)
    XtVinvy <- numeric(p)

    # Accumulate X'V^{-1}X and X'V^{-1}y over groups using Woodbury:
    # V_j^{-1} v = (v - Zj M^{-1} Zj'v / sigma2_e) / sigma2_e
    # where M = Sigma_u^{-1} + Zj'Zj / sigma2_e
    for (j in seq_len(J)) {
      idx <- which(grp == j)
      Xj  <- X[idx, , drop = FALSE]
      Zj  <- Zsmall[idx, , drop = FALSE]
      yj  <- y[idx]

      M    <- Sigma_u_inv + crossprod(Zj) / sigma2_e
      Minv <- tryCatch(solve(M), error = function(e) diag(1 / pmax(diag(M), 1e-10), k))

      VinvFun <- function(v) {
        (v - Zj %*% (Minv %*% drop(crossprod(Zj, v))) / sigma2_e) / sigma2_e
      }

      XtVinvX <- XtVinvX + crossprod(Xj, VinvFun(Xj))
      XtVinvy <- XtVinvy + drop(crossprod(Xj, VinvFun(yj)))
    }

    beta_hat <- tryCatch(
      solve(XtVinvX, XtVinvy),
      error = function(e) tryCatch(
        solve(XtVinvX + diag(1e-8, p), XtVinvy),
        error = function(e2) beta_hat
      )
    )

    # BLUPs: u_j = M^{-1} Zj' resj / sigma2_e  (derivation: Sigma_u Zj' Vj^{-1} resj)
    ss_resid <- 0
    for (j in seq_len(J)) {
      idx  <- which(grp == j)
      Xj   <- X[idx, , drop = FALSE]
      Zj   <- Zsmall[idx, , drop = FALSE]
      resj <- y[idx] - drop(Xj %*% beta_hat)
      M    <- Sigma_u_inv + crossprod(Zj) / sigma2_e
      Minv <- tryCatch(solve(M), error = function(e) diag(1 / pmax(diag(M), 1e-10), k))
      u_hat[j, ] <- drop(Minv %*% drop(crossprod(Zj, resj))) / sigma2_e
      ss_resid   <- ss_resid + sum((resj - drop(Zj %*% u_hat[j, ]))^2)
    }

    # Moment update for Sigma_u: (1/J) sum_j u_j u_j'
    Sigma_u_new <- crossprod(u_hat) / J
    ev          <- eigen(Sigma_u_new, symmetric = TRUE)
    ev$values   <- pmax(ev$values, 1e-8)
    Sigma_u     <- ev$vectors %*% diag(ev$values, k) %*% t(ev$vectors)

    sigma2_e <- max(ss_resid / n, 1e-8)
  }

  vcov_beta <- tryCatch(solve(XtVinvX), error = function(e) diag(1e-6, p))

  if (!is.null(grp_levels))
    rownames(u_hat) <- grp_levels

  list(
    beta      = as.vector(beta_hat),
    u         = u_hat,
    Sigma_u   = Sigma_u,
    sigma2_e  = sigma2_e,
    vcov_beta = vcov_beta
  )
}


meanDiagZGZt <- function(fit, varCorr.c, dep = NULL) {
  mf   <- stats::model.frame(fit)
  bars <- reformulas::findbars(stats::formula(fit))

  n <- nrow(mf)
  v <- numeric(n)

  for (b in bars) {
    expr  <- b[[2]]
    gname <- lmerRestoreIntr(deparse(b[[3]]))

    Sigma <- varCorr.c[[gname]]
    if (is.null(Sigma)) {
      stop("varCorr.c is missing grouping factor '", gname, "'.")
    }
    Sigma <- as.matrix(Sigma)

    f_rhs  <- stats::as.formula(paste0("~", deparse(expr)))
    Zsmall <- stats::model.matrix(f_rhs, mf)

    if (!is.null(dep)) {
      colnames(Zsmall) <- lmerGetNames(dep, colnames(Zsmall))
    }

    rn <- rownames(Sigma)
    if (is.null(rn)) rn <- colnames(Sigma)
    if (is.null(rn)) {
      stop("Sigma for grouping factor '", gname, "' must have row/col names.")
    }

    missing_cols <- setdiff(rn, colnames(Zsmall))
    if (length(missing_cols) > 0) {
      stop("For grouping factor '", gname, "', Z columns missing: ",
           paste(missing_cols, collapse = ", "))
    }

    Zsmall <- Zsmall[, rn, drop = FALSE]

    # Rowwise quadratic form: diag(Zsmall %*% Sigma %*% t(Zsmall))
    # computed as rowSums((Zsmall %*% Sigma) * Zsmall)
    A <- Zsmall %*% Sigma
    v <- v + rowSums(A * Zsmall)
  }

  mean(v)
}


sigma2FromUnitVarY <- function(fit, beta.c, varCorr.c, targetVarY = 1,
                               Cov.lv, dep, indep) {
  if (length(indep)) {
    beta.par <- paste0(dep, "~", indep)
    beta.sub <- beta.c[beta.par]
    Cov.x    <- Cov.lv[indep, indep, drop = FALSE]
    v_fi     <- drop(t(beta.sub) %*% Cov.x %*% beta.sub)
  } else {
    v_fi <- 0
  }

  # Expected random-effect variance contribution on the consistent (latent)
  # scale: E[z' Sigma_u z] = tr(Sigma_u E[zz']).
  #
  # Here E[zz'] is built from the corrected covariance matrix of the structural
  # variables (Cov.lv). Since all variables are standardized and centered,
  # intercept cross-moments are assumed to be 0.
  v_re <- 0

  for (VC in varCorr.c) {
    Sigma_u <- as.matrix(VC)
    re_pars <- rownames(Sigma_u)
    stopif(is.null(re_pars), "Random effect covariance matrix must have rownames.")

    dep_regex <- stringr::str_replace_all(
      dep,
      pattern = "([\\.\\^\\$\\|\\(\\)\\[\\]\\*\\+\\?\\{\\}\\\\])",
      replacement = "\\\\\\1"
    )
    rhs <- stringr::str_remove(re_pars, pattern = paste0("^", dep_regex, "~"))
    Ezz <- matrix(0, nrow = length(re_pars), ncol = length(re_pars),
                  dimnames = list(re_pars, re_pars))

    is_int <- rhs == "1"
    Ezz[is_int, is_int] <- 1

    slope_rhs <- rhs[!is_int]
    if (length(slope_rhs)) {
      missing <- setdiff(slope_rhs, colnames(Cov.lv))
      stopif(length(missing),
             "Missing variables in covariance matrix: ",
             paste0(missing, collapse = ", "))

      slope_pars <- re_pars[!is_int]
      Ezz[slope_pars, slope_pars] <- Cov.lv[slope_rhs, slope_rhs, drop = FALSE]
    }

    v_re <- v_re + sum(Sigma_u * Ezz)
  }

  max(targetVarY - v_re - v_fi, 0)
}


getSigmaFromVarCorr <- function(fit, beta, varCorr, dep, indep, CorrectionCov, Cov.lv) {
  rvdep <- sprintf("%s~~%s", dep, dep)

  # Match the simulation parameterization (zeta) on the consistent scale.
  targetVarY <- if (!is.null(Cov.lv) && dep %in% colnames(Cov.lv)) Cov.lv[dep, dep] else 1

  sigma2 <- sigma2FromUnitVarY(fit, beta.c = beta, varCorr.c = varCorr,
                               targetVarY = targetVarY, Cov.lv = Cov.lv,
                               dep = dep, indep = indep)
  sigma  <- stats::setNames(sigma2, nm = rvdep)

  for (VC in varCorr) {
    namesVC <- matrix("", nrow = NROW(VC), ncol = NCOL(VC))
    for (i in seq_len(NROW(VC))) for (j in seq_len(i))
      namesVC[i, j] <- sprintf("%s~~%s", rownames(VC)[i], colnames(VC)[j])

    namesSigma  <- namesVC[lower.tri(namesVC, diag = TRUE)]
    valuesSigma <- VC[lower.tri(VC, diag = TRUE)]
    names(valuesSigma) <- namesSigma

    sigma <- c(sigma, valuesSigma)
  }

  sigma
}
