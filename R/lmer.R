plslmer <- function(plsModel) {
  INTR_OP <- "__INTR__"

  safeIntr <- function(x) {
    stringr::str_replace_all(x, stringr::fixed(":"), INTR_OP)
  }

  restoreIntr <- function(x) {
    x <- stringr::str_replace_all(x, stringr::fixed(INTR_OP), ":")
    stringr::str_replace_all(x, stringr::fixed("`"), "")
  }

  lme4.syntax <- plsModel@info$lme4.syntax
  cluster     <- plsModel@info$cluster
  consistent  <- plsModel@info$consistent

  stopif(!is.character(lme4.syntax), "`lme4.syntax` must be a character vector!")

  stopif(length(cluster) != 1L || !is.character(cluster),
         "`cluster` must be a character string of length 1. If lme4.syntax is provided!")

  fit.c <- plsModel@fitConsistent
  fit.u <- plsModel@fitUncorrected

  selectGamma <- plsModel@matrices$select$gamma
  # Full (corrected) covariance/correlation matrix for the structural model.
  # Unlike `fitCov`, this retains the total covariances between endogenous and
  # exogenous variables.
  Cov.lv <- fit.c$fitC
  Q <- fit.c$Q
  Correction <- fit.c$fitStructural / fit.u$fitStructural
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

  # We don't want lmer to form the interaction terms on its' own, since we want
  # them to be mean centered, matching the assumptions of our correction procedure.
  # We threefore replace the ':' operator, and use explicitly formed interaction
  # terms, also enaming any columns containing ':' before fitting, and map the
  # lme4 syntax accordingly.
  colsWithColon <- colnames(X)[grepl(":", colnames(X), fixed = TRUE)]
  colsWithColonSafe <- safeIntr(colsWithColon)

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

  getNames <- function(lhs, nm) {
    rhs <- stringr::str_replace_all(nm, pattern = "\\(Intercept\\)", replacement = "1")
    rhs <- restoreIntr(rhs)
    paste0(lhs, "~", rhs)
  }

  fixVecNames <- function(vec, dep) {
    if (!is.null(names(vec)))
      names(vec) <- getNames(dep, names(vec))
    vec
  }

  fixMatNames <- function(mat, dep, cols = TRUE, rows = TRUE) {
    if (!is.null(rownames(mat)) && rows)
      rownames(mat) <- getNames(dep, rownames(mat))
    if (!is.null(colnames(mat)) && cols)
      colnames(mat) <- getNames(dep, colnames(mat))
    mat
  }

  FITS    <- list()
  FIXEF   <- list()
  COEF    <- list()
  VCOV    <- list()
  VARCORR <- list()
  SIGMA   <- list()

  for (line in lme4.syntax) {
    line.safe <- makeSafeFormula(line)

    lmerFit <- lme4::lmer(
      formula = line.safe,
      data    = X.safe,
      control = lme4::lmerControl(calc.derivs = FALSE)
    )

    fterms  <- stats::terms(stats::formula(line.safe))
    vars    <- attr(fterms, "variables")
    dep     <- restoreIntr(as.character(vars[[2L]]))
    indep   <- colnames(selectGamma)[selectGamma[,dep]]

    fixefFit   <- fixVecNames(lme4::fixef(lmerFit), dep = dep)
    vcovFit    <- fixMatNames(stats::vcov(lmerFit), dep = dep)
    coefFit    <- stats::coef(lmerFit)
    varCorrFit <- lme4::VarCorr(lmerFit)

    params <- names(fixefFit)
    split  <- stringr::str_split_fixed(params, pattern = "~", n = 2L)
    lhs    <- split[,1L]
    rhs    <- split[,2L]

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
      c.safe <- safeIntr(c)

      coefFit[[c]] <- fixMatNames(as.matrix(coefFit[[c.safe]]), dep = dep, rows = FALSE)
      varCorrFit[[c]] <- fixMatNames(varCorrFit[[c.safe]], dep = dep)

      if (c.safe != c) {
        coefFit[[c.safe]] <- NULL
        varCorrFit[[c.safe]] <- NULL
      }

      if (consistent) {
        vcPars <- rownames(varCorrFit[[c]])
        vcCorrection <- DCorrectionTerms[vcPars, vcPars, drop = FALSE]

        coefFit[[c]] <- coefFit[[c]] %*% diag(correctionTerms[colnames(coefFit[[c]])])
        varCorrFit[[c]] <- vcCorrection %*% varCorrFit[[c]] %*% vcCorrection
      }

      attr(varCorrFit[[c]], "stddev") <- sqrt(diag(varCorrFit[[c]]))
      attr(varCorrFit[[c]], "correlation") <- cov2cor(varCorrFit[[c]])
    }

    FITS[[dep]]    <- lmerFit
    COEF[[dep]]    <- coefFit
    VCOV[[dep]]    <- vcovFit
    FIXEF[[dep]]   <- fixefFit
    VARCORR[[dep]] <- varCorrFit
    SIGMA[[dep]]   <- getSigmaFromVarCorr(fit = lmerFit, varCorr = varCorrFit,
                                          beta = fixefFit, CorrectionCov = CorrectionCov,
                                          Cov.lv = Cov.lv, dep = dep, indep = indep)
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


meanDiagZGZt <- function(fit, varCorr.c, dep = NULL) {
  INTR_OP <- "__INTR__"

  restoreIntr <- function(x) {
    x <- stringr::str_replace_all(x, stringr::fixed(INTR_OP), ":")
    stringr::str_replace_all(x, stringr::fixed("`"), "")
  }

  mf   <- stats::model.frame(fit)
  bars <- reformulas::findbars(stats::formula(fit))

  # helper to match your naming convention dep~(Intercept)->dep~1 etc.
  getNames <- function(lhs, nm) {
    rhs <- stringr::str_replace_all(nm, pattern="\\(Intercept\\)", replacement="1")
    rhs <- restoreIntr(rhs)
    paste0(lhs, "~", rhs)
  }

  # Build the "small" random-effects model matrix (n x k) for each grouping factor,
  # then compute rowwise z_i' Sigma z_i and sum across factors.
  n <- nrow(mf)
  v <- numeric(n)

  for (b in bars) {
    expr  <- b[[2]]                 # random part, e.g. 1 + x
    gname <- restoreIntr(deparse(b[[3]]))        # grouping factor name

    Sigma <- varCorr.c[[gname]]
    if (is.null(Sigma)) {
      stop("varCorr.c is missing grouping factor '", gname, "'.")
    }
    Sigma <- as.matrix(Sigma)

    # "Small" Z for this term (no expansion by levels), just n x k:
    f_rhs <- stats::as.formula(paste0("~", deparse(expr)))
    Zsmall <- stats::model.matrix(f_rhs, mf)

    # Name alignment: columns of Zsmall must match row/colnames of Sigma
    if (!is.null(dep)) {
      colnames(Zsmall) <- getNames(dep, colnames(Zsmall))
    }

    # Reorder Zsmall to match Sigma ordering
    rn <- rownames(Sigma)
    if (is.null(rn)) rn <- colnames(Sigma)
    if (is.null(rn)) {
      stop("Sigma for grouping factor '", gname, "' must have row/col names.")
    }

    # Check that all Sigma terms exist in Zsmall
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

  # Keep returning the (corrected) random-effect covariances as you already do
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
