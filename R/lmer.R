plslmer <- function(plsModel) {
  lme4.syntax <- plsModel$info$lme4.syntax
  cluster     <- plsModel$info$cluster
  consistent  <- plsModel$info$consistent

  stopif(!is.character(lme4.syntax), "`lme4.syntax` must be a character vector!")

  stopif(length(cluster) != 1L || !is.character(cluster),
         "`cluster` must be a character string of length 1. If lme4.syntax is provided!")
    
  fit.c <- plsModel$fit.c
  fit.u <- plsModel$fit.u

  Correction <- fit.c$fitStructural / fit.u$fitStructural

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

    fixefFit   <- fixVecNames(lme4::fixef(lmerFit), dep = dep)
    vcovFit    <- fixMatNames(vcov(lmerFit), dep = dep)
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
    SIGMA[[dep]]   <- getSigmaFromVarCorr(fit = lmerFit, varCorr = varCorrFit, dep = dep)
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


getSigmaFromVarCorr <- function(fit, varCorr, dep) {
  rvdep  <- sprintf("%s~~%s", dep, dep)
  sigma <- stats::setNames(lme4::getME(fit, "sigma")^2, nm = rvdep)

  for (VC in varCorr) {
    namesVC <- matrix("", nrow = NROW(VC), ncol = NCOL(VC))
    for (i in seq_len(NROW(VC))) for (j in seq_len(i))
      namesVC[i, j] <- sprintf("%s~~%s", rownames(VC)[i], colnames(VC)[j])

    namesSigma <- namesVC[lower.tri(namesVC, diag = TRUE)]
    valuesSigma <- VC[lower.tri(VC, diag = TRUE)]
    names(valuesSigma) <- namesSigma

    sigma <- c(sigma, valuesSigma)
  }

  sigma 
}


# ORIGINAL SKETCH:
#    plscr <- function(.model, .data,
#                      .disattenuate = NA, # capture
#                      .lme4.model = NULL,
#                      cluster = NULL,
#                      mix = FALSE,
#                      correct.scores = FALSE,
#                      ...) {
#      pt <- modsem::modsemify(.model)
#      .model.step1 <- modsem:::parTableToSyntax(pt[!grepl(":", pt$rhs), , drop = FALSE])
#      fit.c <- csem(.model = .model.step1, .data = .data, .disattenuate = TRUE, .id = NULL)
#    
#      if (!correct.scores) {
#        fit.u <- csem(.model = .model.step1, .data = .data, .disattenuate = FALSE, .id = NULL)
#        # Create correction matrix for correction of path coefficients
#        # and do not correct factor scores directly
#    
#        if (!is.null(cluster)) {
#          cfit.u <- csem(.model = .model.step1, .data = .data, .disattenuate = FALSE, .id = cluster)
#          cfit.c <- csem(.model = .model.step1, .data = .data, .disattenuate = TRUE, .id = cluster)
#    
#          correction <- cfit.c[[1L]]$Estimates$Path_estimates
#          correction[TRUE] <- 0
#     
#          k <- length(cfit.c)
#          for (i in seq_along(cfit.c)) {
#            Gamma.c <- cfit.c[[i]]$Estimates$Path_estimates
#            Gamma.u <- cfit.u[[i]]$Estimates$Path_estimates
#            correction <- correction + Gamma.c / (k * Gamma.u)
#          }
#    
#        } else {
#          Gamma.c <- fit.c$Estimates$Path_estimates
#          Gamma.u <- fit.u$Estimates$Path_estimates
#          correction <- Gamma.c / Gamma.u
#        
#        }
#        
#        .X <- getConstructScores(fit.c)[[1L]]
#        .Y <- cbind(.X, .data)
#      } else {
#        .model.step2 <- modsem:::parTableToSyntax(pt[pt$op == "~", ])
#        .X <- getConstructScores(fit.c)[[1L]]
#        .S <- fit.c$Estimates$Construct_VCV
#        .Y <- fitX2Cov(X = .X, S = .S)
#        .Y <- cbind(.Y, .data)
#        correction <- NULL
#      }
#     
#      if (is.null(.lme4.model)) {
#        list(fit = sem(model = .model.step2, data = .Y), correction = correction)
#      } else {
#        fit.lmer <- lmer(.lme4.model, data = .Y)
#        coef <- coef(fit.lmer)
#        VarCorr <- VarCorr(fit.lmer)
#        fixef <- fixef(fit.lmer)
#        vcov  <- vcov(fit.lmer)
#    
#        # for now we just hard code this stuff
#        cx <- correction["Y", "X"]
#        cz <- correction["Y", "Z"]
#        attr(VarCorr$cluster, "stddev")["X"] <- cx * attr(VarCorr$cluster, "stddev")["X"]
#        attr(VarCorr$cluster, "stddev")["Z"] <- cz * attr(VarCorr$cluster, "stddev")["Z"]
#        fixef["X"] <- cx * fixef["X"]
#        fixef["Z"] <- cz * fixef["Z"]
#        coef$cluster[,"X"] <- cx * coef$cluster[,"X"]
#        coef$cluster[,"Z"] <- cz * coef$cluster[,"Z"]
#        se <- sqrt(diag(vcov))
#        se["X"] <- cx * se["X"]
#        se["Z"] <- cz * se["Z"]
#    
#        parTable <- data.frame(
#          lhs = c( "Y", "Y", "Y"),
#          op  = c("~1", "~", "~"),
#          rhs = c(  "", "X", "Z"),
#          est = fixef,
#          se  = se
#        )
#    
#        out <- list(
#          fit = fit.lmer,
#          correction = correction,
#          coef = coef,
#          VarCorr = VarCorr,
#          fixef = fixef,
#          parTable = parTable
#        )
#    
#        class(out) <- "MlPlsSem"
#        out
#      }
#    }
#    
