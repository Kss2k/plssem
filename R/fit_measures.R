setGeneric("impliedConstructCorrMat",
           function(object, saturated = FALSE) standardGeneric("impliedConstructCorrMat"))

setMethod("impliedConstructCorrMat", "PlsBaseModel", function(object, saturated = FALSE) {
  fit <- object@fit
  M   <- object@matrices

  if (saturated) return(M$C)

  Gamma  <- t(fit$fitStructural)
  Lambda <- fit$fitLambda
  Phi    <- fit$fitCov

  I     <- diag(NCOL(Gamma))
  B     <- I - Gamma
  B.inv <- solve(B)
  B.inv %*% Phi %*% t(B.inv)
})

setMethod("impliedConstructCorrMat", "PlsModel", function(object, saturated = FALSE) {
  is.hi.ord <- object@info$is.high.ord

  if (is.null(is.hi.ord) || !is.hi.ord)
    return(callNextMethod())

  fo <- firstOrder(object)
  so <- secondOrder(object)

  C      <- impliedIndicatorCorrMat(so, saturated = saturated)
  C      <- plssemMatrix(C, is.public = TRUE)
  Lambda <- fo@fit$fitLambda
  Gamma  <- t(fo@fit$fitStructural)
  cols   <- colnames(Gamma)
  Phi    <- C[cols, cols, drop = FALSE]

  I     <- diag(NCOL(Gamma))
  B     <- I - Gamma
  B.inv <- solve(B)
  B.inv %*% Phi %*% t(B.inv)
})


setGeneric("impliedIndicatorCorrMat",
           function(object, saturated = FALSE) standardGeneric("impliedIndicatorCorrMat"))

setMethod("impliedIndicatorCorrMat", "PlsBaseModel", function(object, saturated = FALSE) {
  fit    <- object@fit
  M      <- object@matrices
  mode.b <- object@info$mode.b

  S      <- M$S
  Phi    <- impliedConstructCorrMat(object, saturated = saturated)
  Lambda <- fit$fitLambda

  ovs    <- intersect(colnames(S), rownames(Lambda))
  lvs    <- intersect(colnames(Phi), colnames(Lambda))
  mode.b <- intersect(mode.b, lvs)

  Lambda  <- Lambda[ovs, lvs, drop = FALSE]
  Phi     <- Phi[lvs, lvs, drop = FALSE]
  S       <- S[ovs, ovs, drop = FALSE]

  Theta   <- diag2(S) - diag2(Lambda %*% t(Lambda))
  Implied <- Lambda %*% Phi %*% t(Lambda) + Theta

  if (length(mode.b)) {
    MeasurementB <- Lambda[, mode.b, drop = FALSE] != 0
    BlocksB      <- (MeasurementB %*% t(MeasurementB)) == 1
    Implied[BlocksB] <- S[BlocksB]
  }

  Implied
})


fitMeasures <- function(model, saturated = FALSE, mc.reps = 1e6) {
  tryCatch({

    if (model@info$is.mcpls) {
      warning2("Fit measures for MC-PLSc models are under development!\n",
               "Traditional fit criteria will likely be too strict.")
      message(sprintf("Resampling MC-PLSc Model (R = %d)...", mc.reps))

      resampled <- resampleMCPLS_Fit(model, mc.reps = mc.reps)
      Expected  <- resampled@matrices$S.ord.expected
      Observed  <- resampled@matrices$S.ord.observed

    } else {
      Expected <- impliedIndicatorCorrMat(model, saturated = FALSE)
      Observed <- model@matrices$S
    }

    Expected <- cov2cor(Expected)
    Observed <- cov2cor(Observed)

    N        <- NROW(model@data)
    chisq    <- calcChisq(model, saturated = saturated, E = Expected, O = Observed)
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
         rmsea = NA_real_, srmr    = NA_real_)
  })
}


calcSRMR <- function(model,
                     saturated = FALSE,
                     diagonal  = TRUE,
                     Expected  = impliedIndicatorCorrMat(model, saturated = saturated),
                     Observed  = model@matrices$S) {
  tryCatch({
    Diff <- cov2cor(Expected) - cov2cor(Observed)
    sqrt(mean(Diff[lower.tri(Diff, diag = diagonal)]^2))
  }, error = function(e) {
    warning2("Calculation of SRMR failed! Message:\n", conditionMessage(e))
    NA_real_
  })
}


calcChisq <- function(model,
                      saturated = FALSE,
                      O         = model@matrices$S,
                      E         = impliedIndicatorCorrMat(model, saturated = saturated),
                      N         = NROW(model@data),
                      p         = NCOL(O)) {
  tryCatch({
    Einv <- solve(E)
    as.vector((N - 1) * (tr(O %*% Einv) - log(det(O %*% Einv)) - p))
  }, error = function(e) {
    warning2("Calculation of chi-square failed! Message:\n", conditionMessage(e))
    NA_real_
  })
}


calcChisqDf <- function(model, saturated = FALSE) {
  tryCatch({

    info   <- model@info
    fit    <- model@fit
    M      <- model@matrices

    mode.a <- info$mode.a
    mode.b <- info$mode.b
    lvs    <- info$lvs.linear
    xis    <- info$xis
    inds   <- info$indsLvs

    p         <- NCOL(M$S)
    Gamma     <- fit$fitStructural
    total.df  <- p * (p - 1) / 2
    df        <- total.df

    for (lv in mode.a) {
      df <- df - length(inds[[lv]])
    }

    for (lv in mode.b) {
      nw <- length(inds[[lv]])
      ic <- nw * (nw - 1) / 2
      df <- df - (nw - 1 + ic)
    }

    if (length(xis)) {
      ps <- length(xis)
      df <- df - ps * (ps - 1) / 2
    }

    if (!info$is.cfa && !saturated)
      df <- df - sum(M$select$gamma)

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


calcRMSEA <- function(chi.sq, df, N, ci.level = 0.90, close.h0 = 0.05) {
  tryCatch({

    alpha  <- 1 - ci.level
    fLower <- \(lambda) stats::pchisq(chi.sq, df, ncp = lambda) - (1 - alpha / 2)
    fUpper <- \(lambda) stats::pchisq(chi.sq, df, ncp = lambda) - (    alpha / 2)
    fRMSEA <- \(lambda) sqrt(max(lambda, 0) / (df * (N - 1)))

    point <- chi.sq - df
    lower <- tryCatchUniroot(fLower, lower = 0, upper = chi.sq,       errorVal = 0)
    upper <- tryCatchUniroot(fUpper, lower = 0, upper = 10 * chi.sq,  errorVal = df * (N - 1))

    list(
      rmsea    = fRMSEA(point),
      lower    = fRMSEA(lower),
      upper    = fRMSEA(upper),
      ci.level = ci.level,
      pvalue   = 1 - stats::pchisq(chi.sq, df = df, ncp = df * (N - 1) * close.h0^2),
      close.h0 = close.h0
    )

  }, error = function(e) {
    warning2("Calculation of RMSEA failed! Message:\n", conditionMessage(e))
    list(rmsea = NA_real_, lower = NA_real_, upper = NA_real_,
         ci.level = NA_real_, pvalue = NA_real_, close.h0 = NA_real_)
  })
}
