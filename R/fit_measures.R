impliedConstructCorrMat <- function(model, saturated = FALSE, mc.reps = 1e6) {
  
  if (is_mcpls(model)) {
    pls_msg_note(sprintf("Resampling MC-PLSc Model (R = %d)...", mc.reps))
    combined <- combinedModel(model)
    resampled <- resampleMCPLS_Fit(combined, mc.reps = mc.reps)
    return(resampled@matrices$C)
  }

  fit  <- model@fit
  info <- model@info
  M    <- model@matrices

  higherOrder <- higherOrderModel(model)

  if (!is.null(higherOrder)) {

    CTmp <- impliedIndicatorCorrMat(higherOrder, saturated = saturated)
    C <- plssemMatrix(CTmp, is.public = TRUE) # handle any tmp variables

    Gamma  <- t(fit$fitStructural)
    cols   <- colnames(Gamma)
    C      <- C[cols, cols, drop = FALSE]

    if (saturated)
      return(C)

    # Whenever there exists a higher order model, we should have Gamma = 0
    # Such that we always return saturated. This is however not garuanteed
    # for future behaviours. Thus we do it properly anyways.
    xis  <- info$xis
    etas <- info$etas

    Cproj <- Gamma %*% C %*% Gamma
    Phi   <- diag2(C) - diag2(Cproj)

    Phi[xis, xis]   <- C[xis, xis]
    Phi[etas, etas] <- Phi[etas, etas]
    Phi[etas, xis]  <- Phi[xis, etas] <- 0

  } else {

    if (saturated)
      return(M$C)

    Gamma  <- t(fit$fitStructural)
    Phi    <- fit$fitCov

  }

  I <- diag(NCOL(Gamma))
  B <- I - Gamma
  B.inv <- solve(B)

  B.inv %*% Phi %*% t(B.inv)
}


impliedIndicatorCorrMat <- function(object, saturated = FALSE, mc.reps = 1e6) {

  if (is_mcpls(object)) {
    pls_msg_note(sprintf("Resampling MC-PLSc Model (R = %d)...", mc.reps))
    combined <- combinedModel(object)
    resampled <- resampleMCPLS_Fit(combined, mc.reps = mc.reps)
    return(resampled@matrices$S.ord.expected)
  }

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
}


# Joint model-implied correlation matrix of observed and latent variables.
# The observed-observed and latent-latent blocks reuse the canonical implied
# helpers; the cross block is Cov(y, eta) = Lambda %*% Phi. The full loading
# matrix is used so that associations running through single-indicator stand-in
# latents (observed structural variables, composite-MIMIC reflective indicators)
# are captured.
impliedJointCorrMat <- function(object, saturated = FALSE, mc.reps = 1e6) {
  Phi    <- impliedConstructCorrMat(object, saturated = saturated, mc.reps = mc.reps)
  SigmaO <- impliedIndicatorCorrMat(object, saturated = saturated, mc.reps = mc.reps)

  ovs <- rownames(SigmaO)
  lvs <- colnames(Phi)

  Lambda    <- matrix(0, nrow = length(ovs), ncol = length(lvs),
                      dimnames = list(ovs, lvs))
  fitLambda <- object@fit$fitLambda
  ovs.f     <- intersect(ovs, rownames(fitLambda))
  lvs.f     <- intersect(lvs, colnames(fitLambda))
  Lambda[ovs.f, lvs.f] <- fitLambda[ovs.f, lvs.f]

  cross <- Lambda %*% Phi # Cov(observed, latent)

  rbind(
    cbind(SigmaO, cross),
    cbind(t(cross), Phi)
  )
}


fitMeasures <- function(model, saturated = FALSE, mc.reps = 1e6) {
  tryCatch({

    if (is_mcpls(model)) {
      pls_msg_note(
        paste0("Fit measures for MC-PLSc models are under development!\n",
               "Traditional fit criteria will likely be too strict.")
      )

      pls_msg_note(sprintf("Resampling MC-PLSc Model (R = %d)...", mc.reps))

      combined  <- combinedModel(model) # combined chained models
      resampled <- resampleMCPLS_Fit(combined, mc.reps = mc.reps)
      Expected  <- resampled@matrices$S.ord.expected
      Observed  <- resampled@matrices$S.ord.observed

    } else {
      Expected <- impliedIndicatorCorrMat(model, saturated = FALSE)
      Observed <- model@matrices$S
    }

    Expected <- cov2cor(Expected)
    Observed <- cov2cor(Observed)

    N        <- NROW(model@data)
    chisq    <- calcChisq(Expected = Expected, Observed = Observed, N = N)
    chisq.df <- calcChisqDf(model, saturated = saturated)

    list(
      chisq    = chisq,
      chisq.df = chisq.df,
      rmsea    = calcRMSEA(chisq, df = chisq.df, N = N)$rmsea,
      srmr     = calcSRMR(Expected = Expected, Observed = Observed,
                          saturated = saturated)
    )

  }, error = function(e) {
    pls_msg_warn(paste0("Calculation of fit measures failed, message:\n",
                 conditionMessage(e)))
    list(chisq = NA_real_, chisq.df = NA_real_,
         rmsea = NA_real_, srmr    = NA_real_)
  })
}


calcSRMR <- function(Expected, Observed, saturated = FALSE, diagonal  = TRUE) {
  tryCatch({
    Diff <- cov2cor(Expected) - cov2cor(Observed)
    sqrt(mean(Diff[lower.tri(Diff, diag = diagonal)]^2))
  }, error = function(e) {
    pls_msg_warn(paste0("Calculation of SRMR failed! Message:\n", conditionMessage(e)))
    NA_real_
  })
}


calcChisq <- function(Observed, Expected, N, p = NCOL(Observed)) {
  tryCatch({
    Einv <- solve(Expected)
    as.vector((N - 1) * (tr(Observed %*% Einv) - log(det(Observed %*% Einv)) - p))
  }, error = function(e) {
    pls_msg_warn(paste0("Calculation of chi-square failed! Message:\n", conditionMessage(e)))
    NA_real_
  })
}


calcChisqDf <- function(model, saturated = FALSE, count.interactions = FALSE) {
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
    freeGamma <- M$select$gamma
    total.df  <- p * (p - 1) / 2
    df        <- total.df

    if (!count.interactions) {
      xis <- intersect(info$xis, lvs)
      freeGamma <- freeGamma[
        !grepl(":", rownames(freeGamma)),
        !grepl(":", colnames(freeGamma)), drop = FALSE
      ]
    }

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
      df <- df - sum(freeGamma)

    if (hasHigherOrderModel(model))
      df <- df + calcChisqDf(higherOrderModel(model))

    df

  }, error = function(e) {
    pls_msg_warn(paste0("Calculation of chisq degrees of freedom failed! Message:\n",
                 conditionMessage(e)))
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
    pls_msg_warn(paste0("Calculation of RMSEA failed! Message:\n", conditionMessage(e)))
    list(rmsea = NA_real_, lower = NA_real_, upper = NA_real_,
         ci.level = NA_real_, pvalue = NA_real_, close.h0 = NA_real_)
  })
}


plsImpliedR2 <- function(object, output = c("all", "lv", "ov")) {
  output <- match.arg(output)

  parTable <- combinedModel(object)@parTable

  getR2 <- function(x, pt = parTable) {
    rvar <- pt[pt$lhs == x & pt$op == "~~" & pt$rhs == x, "est"]
    if (!length(rvar)) 0 else 1 - rvar
  }

  etas <- getEtas(parTable, checkAny = FALSE)
  inds.a <- getReflectiveIndicators(parTable)

  r2.etas <- vapply(etas,   FUN.VALUE = numeric(1L), FUN = getR2)
  r2.inds <- vapply(inds.a, FUN.VALUE = numeric(1L), FUN = getR2)

  names(r2.etas) <- etas
  names(r2.inds) <- inds.a

  switch(output,
    all = c(r2.inds, r2.etas),
    lv  = r2.etas,
    ov  = r2.inds
  )
}


mcplsLoglik <- function(object, boot.R = 500, verbose = interactive()) {
  combined <- combinedModel(object)
  status   <- modelStatus(combined)

  pls_stopif(!is_mcpls(combined),
    "`mcpls_loglik()` is only implemented for MC-PLS models!"
  )

  par0 <- status$par0
  fit0 <- status$fit0

  pls_stopif(is.null(par0) || is.null(fit0),
    "Missing auxiliary estimates! This is likely a bug!"
  )

  parTable <- getParTableEstimates(
    combined, rm.tmp.ov = FALSE, clean.tmp.ind = FALSE
  )

  # Here we get a whole lot of code duplication compared to R/mcpls.R
  # which isn't optimal.. The code is similar, but quite different as well...
  data <- modelData(object)
  vars <- colnames(data)
  n    <- NROW(data)

  if (isMLM(object)) {
    clusterSizes <- as.numeric(table(attr(data, "cluster")))
    clusterName  <- colnames(attr(data, "cluster"))
  } else {
    clusterSizes <- NULL
    clusterName  <- NULL
  }

  estimator <- combined@info$path.estimator
  use.full.rescov <- switch(object@info$mc.args$rescov,
    full    = TRUE,
    reduced = FALSE,
    auto    = estimator == "gls",
    pls_msg_stop("Unrecognized value for `mc.rescov` argument:", mc.rescov)
  )

  is.hi.ord <- isTRUE(combined@info$is.high.ord)

  .f <- function(i = 0, pb = NULL) {

    if (!is.null(pb)) {
      tryCatch(
        utils::setTxtProgressBar(pb, i),
        error = \(e) pls_msg_warn(
          "Unable to update progress bar!\nMessage:", conditionMessage(e)
        )
      )
    }

    sim <- simulateDataParTable(
      parTable     = parTable,
      N            = n,
      check.hi.ord = is.hi.ord,
      clusterSizes = clusterSizes,
      clusterName  = clusterName,
      full         = use.full.rescov,
      cut          = TRUE # cut from estimated thresholds
    )

    fit.sim <- fit0
    X       <- Rfast::standardise(as.matrix(sim$ov[vars]))
    S       <- Rfast::cova(X)

    if (!is.null(sim$cluster))
      attr(X, "cluster") <- sim$cluster

    # Update observed-data (lowest-order) model input
    modelData(fit.sim)  <- X
    indCorrMatrix(fit.sim) <- S

    fit2 <- estimatePLS_Inner(fit.sim)
    par2 <- getFreeParamsTable(combinedModel(fit2))

    par2[par0$is.free, "est"]
  }

  free0 <- par0[par0$is.free, , drop = FALSE]
  nm0   <- getParNamesFromParTable(free0)

  X <- matrix(NA_real_, nrow = boot.R, ncol = NROW(free0))
  colnames(X) <- nm0

  if (verbose) {
    pb <- utils::txtProgressBar(
      min     = 0,
      max     = boot.R,
      initial = 0,
      style   = 3,
      file    = stderr()
    )

    on.exit(close(pb), add = TRUE)

  } else {
    pb <- NULL

  }

  for (i in seq_len(boot.R)) {

    X[i,] <- tryCatch(.f(i = i, pb = pb), error = function(e) {
      pls_msg_warn("The %dth estimate failed! Message:", conditionMessage(e))
      NA_real_
    })

  }

  complete <- stats::complete.cases(X)
  ncomplete <- sum(complete)
  nparams <- NROW(free0)

  pls_stopif(ncomplete <= nparams,
    "Unable to compute `mcpls_loglik()` with a stable covariance matrix.",
    sprintf(
      "Only %d complete replicate(s) remain for %d free parameter(s).",
      ncomplete, nparams
    ),
    "Increase `boot.R` or inspect warnings from failed estimates."
  )

  observed <- stats::setNames(free0$est, nm = nm0)
  expected <- colMeans(X, na.rm = TRUE)
  sigma    <- cov(X, use = "complete.obs")

  loglik <- tryCatch(
    mvnfast::dmvn(
      X     = matrix(observed, nrow = 1L),
      mu    = expected,
      sigma = sigma,
      log   = TRUE
    ),
    error = function(e) {
      pls_msg_stop(
        "Unable to compute `mcpls_loglik()` with a stable covariance matrix.",
        "Increase `boot.R` or inspect warnings from failed estimates.",
        "Message:", conditionMessage(e)
      )
    }
  )

  list(
    observed = plssemVector(observed),
    expected = plssemVector(expected),
    sigma    = plssemMatrix(sigma),
    loglik   = plssemVector(loglik)
  )
}
