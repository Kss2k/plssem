# PLS estimation steps

estimatePLS_Step0_5 <- function(model) {
  force(model)

  matrices <- model@matrices

  approach <- model@info$approach.weights
  scheme   <- model@info$inner.weights
  centroid <- FALSE

  if (identical(approach, "pca")) {
    succs <- matrices$succs.pca
    preds <- matrices$preds.pca
  } else if (scheme %in% c("centroid", "factorial")) {
    succs    <- matrices$succs.factorial
    preds    <- matrices$preds.factorial
    centroid <- scheme == "centroid"
  } else if (model@info$is.cfa) {
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
    maxiter      = model@status$max.iter.0_5,
    centroid     = centroid
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

  if (cpp && !model@info$is.probit && model@info$is.nlin) {
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

    attr(F, "cluster") <- attr(model@data, "cluster")

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

  if (model@info$consistent) {
    modelFitConsistent(model)  <- getFitPLSModel(model, consistent = TRUE)
    modelFitUncorrected(model) <- list(NULL)
    modelFit(model)            <- modelFitConsistent(model)
  } else {
    modelFitConsistent(model)  <- list(NULL)
    modelFitUncorrected(model) <- getFitPLSModel(model, consistent = FALSE)
    modelFit(model)            <- modelFitUncorrected(model)
  }

  model@matrices$C <- model@fit$fitC
  model
}


estimatePLS_Step8 <- function(model) {
  force(model)

  model@params$values <- extractCoefs(model)
  model@params$se     <- rep(NA_real_, length(model@params$values))

  model
}
