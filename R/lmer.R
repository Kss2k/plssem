plslmer <- function(plsModel) {
  lme4.syntax <- plsModel$info$lme4.syntax
  cluster     <- plsModel$info$cluster
  consistent  <- plsModel$info$consistent

  stopif(!is.character(lme4.syntax), "`lme4.syntax` must be a character vector!")

  stopif(length(cluster) != 1L || !is.character(cluster),
         "`cluster` must be a character string of length 1. If lme4.syntax is provided!")
    
  fit.c <- plsModel$fit.c
  fit.u <- plsModel$fit.u

  selectGamma <- plsModel$matrices$select$gamma
  Cov.lv <- fit.c$fitCov
  Correction <- fit.c$fitStructural / fit.u$fitStructural
  CorrectionCov <- fit.c$fitCov / fit.u$fitCov

  Xf <- as.data.frame(plsModel$factorScores)
  Xx <- as.data.frame(plsModel$data)
  Xc <- as.data.frame(attr(plsModel$data, "cluster"))
  X  <- cbind(Xf, Xx, Xc)

  getNames <- function(lhs, nm) {
    rhs <- stringr::str_replace_all(nm, pattern = "\\(Intercept\\)", replacement = "1")
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
    lmerFit <- lme4::lmer(line, data = X)
    fterms  <- stats::terms(stats::formula(line))
    vars    <- attr(fterms, "variables")
    dep     <- as.character(vars[[2L]])
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
      coefFit[[c]] <- fixMatNames(as.matrix(coefFit[[c]]), dep = dep, rows = FALSE)
      varCorrFit[[c]] <- fixMatNames(varCorrFit[[c]], dep = dep)
     
      if (consistent) {
        vcPars <- rownames(varCorrFit[[c]])
        vcCorrection <- DCorrectionTerms[vcPars, vcPars, drop = FALSE]

        coefFit[[c]] <- coefFit[[c]] %*% correctionTerms
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
  mf   <- stats::model.frame(fit)
  bars <- reformulas::findbars(stats::formula(fit))

  # helper to match your naming convention dep~(Intercept)->dep~1 etc.
  getNames <- function(lhs, nm) {
    rhs <- stringr::str_replace_all(nm, pattern="\\(Intercept\\)", replacement="1")
    paste0(lhs, "~", rhs)
  }

  # Build the "small" random-effects model matrix (n x k) for each grouping factor,
  # then compute rowwise z_i' Sigma z_i and sum across factors.
  n <- nrow(mf)
  v <- numeric(n)

  for (b in bars) {
    expr  <- b[[2]]                 # random part, e.g. 1 + x
    gname <- deparse(b[[3]])        # grouping factor name

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
  beta.par <- paste0(dep, "~", indep)
  beta.sub <- beta.c[beta.par]
  Cov.x    <- Cov.lv[indep, indep]

  v_re <- meanDiagZGZt(fit, varCorr.c, dep = dep)
  v_fi <- t(beta.sub) %*% Cov.x %*% beta.sub

  max(targetVarY - v_re - v_fi, 0)
}


getSigmaFromVarCorr <- function(fit, beta, varCorr, dep, indep, CorrectionCov, Cov.lv) {
  rvdep <- sprintf("%s~~%s", dep, dep)

  # If Var(Y)=1 on your reporting scale, use targetVarY=1
  # If instead you want the "post-hoc scaled" total variance, set targetVarY = CorrectionCov[dep, dep]
  targetVarY <- 1

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
