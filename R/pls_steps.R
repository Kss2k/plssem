# PLS estimation steps

estimatePLS_Step0_5 <- function(model) {
  force(model)

  matrices <- model@matrices

  if (model@info$is.cfa) {
    succs <- matrices$succs.cfa
    preds <- matrices$preds.cfa
  } else {
    succs <- matrices$succs.linear
    preds <- matrices$preds.linear
  }

  result <- estimatePLS_Step0_5_Cpp(
    lambda       = matrices$lambda,
    gamma        = matrices$gamma,
    S            = matrices$S,
    C            = matrices$C,
    SC           = matrices$SC,
    R_IndsIdxLVs = matrices$cpp$indsIdxLVs,
    lvColIdx     = matrices$cpp$lvColIdx,
    modeB        = matrices$cpp$modeB,
    preds        = preds,
    succs        = succs,
    tolerance    = model@status$tolerance,
    maxiter      = model@status$max.iter.0_5
  )

  if (!result$convergence) {
    pls_msg_warn("Convergence not reached. Stopping.")
    model@status$is.admissible <- FALSE
  }

  # status
  model@status$iterations.0_5 <- model@status$iterations.0_5 + result$iterations
  model@status$iterations     <- model@status$iterations + result$iterations
  model@status$convergence    <- result$convergence

  # fit
  model@matrices$C[]      <- result$C
  model@matrices$SC[]     <- result$SC
  model@matrices$lambda[] <- result$lambda
  model@matrices$gamma[]  <- result$gamma

  model
}


estimatePLS_Step6 <- function(model, cpp = TRUE) {
  force(model)

  if (cpp && !model@info$is.probit) {
    par <- colnames(model@matrices$C)

    result <- estimatePLS_Step6_Cpp(
      X = model@data,
      W = model@matrices$lambda,
      prodElemsIdx = model@matrices$cpp$prodElemsIdx,
      prodColIdx   = model@matrices$cpp$prodColIdx,
      standardize  = !model@info$standardized
    )

    F <- result$F
    C <- result$C

    colnames(F) <- colnames(model@matrices$C)
    dimnames(C) <- dimnames(model@matrices$C)

    model@factorScores          <- F
    model@matrices$C[par, par]  <- C
    model@matrices$SC[par, par] <- C

    return(model)
  }

  model@factorScores <- computeFactorScores(model)

  if (!model@info$is.nlin)
    return(model)

  # Update variance and covariances of interaction terms.
  elems <- model@info$intTermElems
  X     <- model@factorScores

  for (elems.xz in elems) {
    xz       <- paste0(elems.xz, collapse = ":")
    X[, xz]  <- Rfast::rowprods(X[, elems.xz])
  }

  Cxz <- Rfast::cova(X)
  par <- colnames(X)

  model@factorScores    <- X
  model@matrices$C[par, par]  <- Cxz
  model@matrices$SC[par, par] <- Cxz
  model
}


estimatePLS_Step7 <- function(model) {
  force(model)

  is.mlm     <- model@info$is.mlm
  is.mcpls   <- model@info$is.mcpls
  consistent <- model@info$consistent
  is.probit  <- model@info$is.probit

  if (!is.mlm) {
    if (consistent) {
      modelFitConsistent(model)  <- getFitPLSModel(model, consistent = TRUE)
      modelFitUncorrected(model) <- list(NULL)
      modelFit(model)            <- modelFitConsistent(model)
    } else {
      modelFitConsistent(model)  <- list(NULL)
      modelFitUncorrected(model) <- getFitPLSModel(model, consistent = FALSE)
      modelFit(model)            <- modelFitUncorrected(model)
    }

    model@matrices$C <- model@fit$fitC
    return(model)
  }

  model.c <- model
  model.u <- model

  if (is.probit || is.mcpls) {
    model.u <- updateModelInfo(model.u, is.probit = FALSE, is.mcpls = FALSE)
    model.u@matrices$S <- getCorrMat(model.u@data, probit = FALSE)
    model.u <- updateOuterWeights(model.u) |> updateFactorScores()
  }

  modelFitConsistent(model)  <- getFitPLSModel(model.c, consistent = consistent)
  modelFitUncorrected(model) <- getFitPLSModel(model.u, consistent = FALSE)
  modelFit(model)            <- modelFitConsistent(model)

  model@matrices$C <- model@fit$fitC
  model
}


estimatePLS_Step8 <- function(model) {
  force(model)

  model@params$values <- extractCoefs(model)
  model@params$se     <- rep(NA_real_, length(model@params$values))

  if (!isMLM(model))
    return(model)

  modelFitLmer(model) <- plslmer(
    plsModel = model, fast = isTRUE(model@info$mc.fast.lmer)
  )

  refreshLmerParams(model) # Update params with Mixed-Effects coefficients
}
