calcImpliedConstructCorrMat <- function(model, saturated = FALSE) {
  fit <- model$fit
  M   <- model$matrices

  is.hi.ord <- model$info$is.high.ord

  if (is.null(is.hi.ord) || !is.hi.ord) {
    if (saturated)
      return(M$C)

    Gamma  <- t(fit$fitStructural)
    Lambda <- fit$fitLambda
    Phi    <- fit$fitCov

  } else {
    firstOrder  <- model$submodels$firstOrder
    secondOrder <- model$submodels$secondOrder

    C <- calcImpliedIndicatorCorrMat(secondOrder, saturated = saturated)
    C <- plssemMatrix(C, is.public = TRUE) # handle any tmp variables

    Lambda <- firstOrder$fit$fitLambda
    Gamma  <- t(firstOrder$fit$fitStructural)
    cols   <- colnames(Gamma)
    Phi    <- C[cols, cols, drop = FALSE]

  }

  I <- diag(NCOL(Gamma))
  B <- I - Gamma
  B.inv <- solve(B)

  B.inv %*% Phi %*% t(B.inv)
}


calcImpliedIndicatorCorrMat <- function(model, saturated = FALSE) {
  fit    <- model$fit
  M      <- model$matrices
  mode.b <- model$info$mode.b

  S      <- M$S
  Phi    <- calcImpliedConstructCorrMat(model, saturated = saturated)
  Lambda <- fit$fitLambda # Contains Lvs for higher orded models

  ovs <- intersect(colnames(S), rownames(Lambda))
  lvs <- intersect(colnames(Phi), colnames(Lambda))

  mode.b <- intersect(mode.b, lvs)
  Lambda <- Lambda[ovs, lvs, drop = FALSE]
  Phi    <- Phi[lvs, lvs, drop = FALSE]
  S      <- S[ovs, ovs, drop = FALSE]

  Theta  <- diag2(S) - diag2(Lambda %*% t(Lambda))
  Implied <- Lambda %*% Phi %*% t(Lambda) + Theta

  # Correct blocks for formative constructs
  if (length(mode.b)) {
    MeasurementB <- Lambda[ ,mode.b, drop=FALSE] != 0
    BlocksB <- (MeasurementB %*% t(MeasurementB)) == 1
    Implied[BlocksB] <- S[BlocksB]
  }

  Implied
}


fitMeasures <- function(model, saturated = FALSE, mc.reps = 1e6) {
  tryCatch({

    if (model$info$is.mcpls) {
      warning2("Fit measures for MC-PLSc models are under development!\n",
               "Traditional fit criteria will likely be too strict.")
      message(sprintf("Resampling MC-PLSc Model (R = %d)...", mc.reps))

      resampled <- resampleMCPLS_Fit(model, mc.reps = mc.reps)

      Expected <- resampled$matrices$S.ord.expected
      Observed <- resampled$matrices$S.ord.observed

    } else {
      Expected <- calcImpliedIndicatorCorrMat(model, saturated = FALSE)
      Observed <- model$matrices$S
    }

    # Make sure both are standardized
    Expected <- cov2cor(Expected)
    Observed <- cov2cor(Observed)

    N     <- NROW(model$data)
    chisq <- calcChisq(model, saturated = saturated,
                       E = Expected, O = Observed)
    chisq.df <- calcChisqDf(model, saturated = saturated)

    list(
      chisq    = chisq,
      chisq.df = chisq.df,
      rmsea    = calcRMSEA(chisq, df = chisq.df, N = N)$rmsea,
      srmr     = calcSRMR(model, saturated = saturated,
                          Expected = Expected, Observed = Observed)
    )

  }, error = function(e) {

    warning2("Calculation of fit measures failed, message:\n",
             conditionMessage(e))

    list(chisq = NA_real_, chisq.df = NA_real_,
         rmsea = NA_real_, srmr = NA_real_)
  })
}


calcSRMR <- function(model,
                     saturated = FALSE,
                     diagonal = TRUE,
                     Expected = calcImpliedIndicatorCorrMat(model, saturated = saturated),
                     Observed = model$matrices$S) {

  tryCatch({

    Diff <- cov2cor(Expected) - cov2cor(Observed)
    sqrt(mean(Diff[lower.tri(Diff, diag = diagonal)]^2))

  }, error = function(e) {

    warning2("Calculation of SRMR failed! Message:\n",
             conditionMessage(e))
    NA_real_

  })
}


calcChisq <- function(model,
                      saturated = FALSE,
                      O = model$matrices$S,
                      E = calcImpliedIndicatorCorrMat(model, saturated = saturated),
                      N = NROW(model$data),
                      p = NCOL(O)) {
  tryCatch({

    Einv <- solve(E)
    as.vector((N - 1) * (tr(O %*% Einv) - log(det(O %*% Einv)) - p))

  }, error = function(e) {

    warning2("Calculation of chi-square failed! Message:\n",
             conditionMessage(e))
    NA_real_

  })
}


calcChisqDf <- function(model, saturated = FALSE, count.diag = TRUE) {
  tryCatch({

    info <- model$info
    fit  <- model$fit
    M    <- model$matrices

    mode.a <- info$mode.a   # reflective blocks
    mode.b <- info$mode.b   # formative/composite blocks
    lvs    <- info$lvs.linear
    inds   <- info$indsLvs

    p <- NCOL(M$S)
    Gamma <- fit$fitStructural

    # Standardized model: work with correlations only, not variances
    c <- if (count.diag) 1 else -1
    total.df <- p * (p + c) / 2
    df <- total.df

    #  For now we don't distinguish between mode A and mode B
    #
    #  # Reflective blocks:
    #  # standardized case -> count only free loadings
    #  for (lv in mode.a) {
    #    pa <- length(inds[[lv]])
    #    df.block <- pa
    #    df <- df - df.block
    #  }

    #  # Formative/composite blocks:
    #  # within-block correlation matrix is saturated
    #  for (lv in mode.b) {
    #    pb <- length(inds[[lv]])
    #    df.block <- pb * (pb - 1) / 2
    #    df <- df - df.block
    #  }

    for (lv in lvs) {
      pb <- length(inds[[lv]])
      df.block <- pb * (pb + c) / 2
      df <- df - df.block
    }

    if (info$is.cfa || saturated) {
      # Freely correlated constructs/composites, standardized:
      # count only off-diagonal latent correlations
      ps <- NROW(Gamma)
      df.structural <- ps * (ps + c) / 2
      df <- df - df.structural

    } else {
      # Structural model specified by directed paths.
      # In the standardized case, path coefficients still count as free parameters.
      df.structural <- sum(M$select$gamma)
      df <- df - df.structural
    }

    df

  }, error = function(e) {

    warning2("Calculation of chisq degrees of freedom failed! Message:\n",
             conditionMessage(e))
    NA_real_

  })
}


tryCatchUniroot <- function(f, lower, upper, errorVal = NA) {
  tryCatch(
    stats::uniroot(f, lower = lower, upper = upper)$root,
    error = function(e) errorVal
  )
}


calcRMSEA <- function(chi.sq, df, N, ci.level = 0.90, close.h0=0.05) {
  tryCatch({

    alpha  <- 1 - ci.level

    fLower <- \(lambda) stats::pchisq(chi.sq, df, ncp=lambda) - (1 - alpha/2)
    fUpper <- \(lambda) stats::pchisq(chi.sq, df, ncp=lambda) - (    alpha/2)
    fRMSEA <- \(lambda) sqrt(max(lambda, 0) / (df * (N - 1)))

    point <- chi.sq - df
    lower <- tryCatchUniroot(fLower, lower=0, upper=chi.sq, errorVal=0)
    upper <- tryCatchUniroot(fUpper, lower=0, upper=10*chi.sq, errorVal=df * (N - 1)) # i.e., RMSEA.upper = 1

    rmseaLower  <- fRMSEA(lower)
    rmseaUpper  <- fRMSEA(upper)
    rmseaHat    <- fRMSEA(point)
    rmseaPvalue <- 1 - stats::pchisq(chi.sq, df=df, ncp=df*(N-1)*close.h0^2)

    list(rmsea = rmseaHat, lower = rmseaLower, upper = rmseaUpper,
         ci.level = ci.level, pvalue = rmseaPvalue, close.h0 = close.h0)

  }, error = function(e) {

    warning2("Calculation of RMSEA failed! Message:\n",
             conditionMessage(e))

    list(rmsea = NA_real_, lower = NA_real_, upper = NA_real_,
         ci.level = NA_real_, pvalue = NA_real_, close.h0 = NA_real_)

  })
}
