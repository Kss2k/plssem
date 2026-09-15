# PLS estimation steps

estimatePLS_Step0_5 <- function(model) {
  force(model)

  max.iter.0_5 <- model@status$max.iter.0_5

  # Steps 0-5 are many tiny matrix operations that depend only on `S` and the
  # model graph. In R the interpreter overhead dominates the arithmetic, and
  # MC-PLS runs this loop `mc.reps` times per iteration, so it is done in C++
  # when the model has the standard (linear, mode A/B) structure.
  fast <- tryCatch(estimatePLS_Step0_5_fast(model), error = function(e) NULL)
  if (!is.null(fast)) return(fast)

  model <- estimatePLS_Step0(model)

  for (i in seq_len(max.iter.0_5)) {
    model <- estimatePLS_Step1(model)
    model <- estimatePLS_Step2(model)
    model <- estimatePLS_Step3(model)
    model <- estimatePLS_Step4(model)
    model <- estimatePLS_Step5(model)

    if (model@status$convergence) {
      break
    } else if (i >= max.iter.0_5) {
      pls_msg_warn("Convergence not reached. Stopping.")
      model@status$is.admissible <- FALSE
      break
    }
  }

  model@status$iterations.0_5 <- model@status$iterations.0_5 + i
  model@status$iterations     <- model@status$iterations + i
  model
}


estimatePLS_Step0 <- function(model) {
  force(model)

  lvs     <- model@info$lvs.linear
  indsLvs <- model@info$indsLvs
  lambda  <- model@matrices$lambda
  SC      <- model@matrices$SC

  for (lv in lvs) {
    inds <- indsLvs[[lv]]
    wj   <- rep(1, length(inds))
    Sjj  <- SC[inds, inds]
    wj   <- wj / c(sqrt(t(wj) %*% Sjj %*% wj))
    lambda[inds, lv] <- wj
  }

  partLambda <- cbind(model@matrices$Ip, lambda)
  S          <- model@matrices$S

  model@matrices$C  <- t(lambda) %*% S %*% lambda
  model@matrices$SC <- t(partLambda) %*% S %*% partLambda
  model@matrices$lambda <- lambda
  model
}


estimatePLS_Step1 <- function(model) {
  force(model)

  lvs   <- model@info$lvs.linear
  succs <- model@matrices$succs.linear
  preds <- model@matrices$preds.linear
  gamma <- model@matrices$gamma
  C     <- model@matrices$C
  SC    <- model@matrices$SC

  if (model@info$is.cfa) {
    succs <- model@matrices$succs.cfa
    preds <- model@matrices$preds.cfa
  }

  for (lv in lvs) {
    predsLv <- lvs[preds[, lv, drop = TRUE]]
    succsLv <- lvs[succs[, lv, drop = TRUE]]

    for (succ in succsLv)
      gamma[succ, lv] <- C[lv, succ]

    if (length(predsLv) > 0)
      gamma[predsLv, lv] <- solve(SC[predsLv, predsLv]) %*% SC[predsLv, lv]

    scalef <- c(sqrt(t(gamma[, lv]) %*% C %*% gamma[, lv]))
    if (scalef)
      gamma[, lv] <- gamma[, lv] / scalef
  }

  model@matrices$gamma <- gamma
  model
}


estimatePLS_Step2 <- function(model) {
  force(model)

  Ip         <- model@matrices$Ip
  lambda     <- model@matrices$lambda
  gamma      <- model@matrices$gamma
  C          <- model@matrices$C
  S          <- model@matrices$S
  SC         <- model@matrices$SC

  if (NROW(gamma) <= 1)
    return(model)

  partLambda <- cbind(Ip, lambda)
  partGamma  <- rbind(
    cbind(Ip, matrix(0, nrow = nrow(Ip), ncol = ncol(gamma))),
    cbind(matrix(0, nrow = nrow(gamma), ncol = ncol(Ip)), gamma)
  )

  newC  <- t(gamma) %*% C %*% gamma
  newSC <- t(partGamma) %*% t(partLambda) %*% S %*% partLambda %*% partGamma

  dimnames(newSC) <- dimnames(SC)

  model@matrices$C  <- newC
  model@matrices$SC <- newSC
  model
}


estimatePLS_Step3 <- function(model) {
  force(model)

  lvs     <- model@info$lvs.linear
  indsLvs <- model@info$indsLvs
  lambda  <- model@matrices$lambda
  SC      <- model@matrices$SC
  modes   <- model@info$modes

  for (lv in lvs) {
    mode.lv <- modes[[lv]]
    inds    <- indsLvs[[lv]]

    wj <- switch(mode.lv,
      A = getWeightsModeA(lv = lv, lambda = lambda, SC = SC, inds = inds),
      B = getWeightsModeB(lv = lv, lambda = lambda, SC = SC, inds = inds),
      NA_real_
    )

    Sjj <- SC[inds, inds]
    wj  <- wj / c(sqrt(t(wj) %*% Sjj %*% wj))
    lambda[inds, lv] <- wj
  }

  model@matrices$lambda <- lambda
  model
}


getWeightsModeA <- function(lv, lambda, SC, inds) {
  as.vector(SC[inds, lv])
}


getWeightsModeB <- function(lv, lambda, SC, inds) {
  getOlsPathCoefs(y = lv, X = inds, C = SC)
}


# Step 4 is structurally identical to step 0: recompute C and SC from the
# updated outer weights after step 3.
estimatePLS_Step4 <- function(model) {
  force(model)

  lambda     <- model@matrices$lambda
  partLambda <- cbind(model@matrices$Ip, lambda)
  S          <- model@matrices$S

  model@matrices$C  <- t(lambda) %*% S %*% lambda
  model@matrices$SC <- t(partLambda) %*% S %*% partLambda
  model
}


estimatePLS_Step5 <- function(model) {
  force(model)

  oldWeights <- model@matrices$outerWeights
  newWeights <- getNonZeroElems(model@matrices$lambda)

  weightDiff <- oldWeights - newWeights
  model@status$convergence <- all(abs(weightDiff) < model@status$tolerance)
  model@matrices$outerWeights <- newWeights
  model
}


estimatePLS_Step6 <- function(model) {
  force(model)

  fast <- estimatePLS_Step6_fast(model)
  if (!is.null(fast)) return(fast)

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
    model.u <- updateOuterWeights(model.u)
    model.u <- updateFactorScores(model.u)
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


# C++ fast path for the outer-weight iteration. Returns NULL when the model
# shape is not supported, so the caller falls back to the reference R code.
# The structural arguments to `plsStep05Cpp()` (index sets, modes, the path
# graph) depend only on the model specification, never on the data, so they are
# identical across every Monte Carlo replication. Rebuilding them per call cost
# as much as the C++ itself, so keep the last one and reuse it when the
# structure is unchanged.
.step05Cache <- new.env(parent = emptyenv())

step05Plan <- function(model) {
  M    <- model@matrices
  info <- model@info

  lvs <- info$lvs.linear
  if (!length(lvs)) return(NULL)

  cfa   <- isTRUE(info$is.cfa)
  succs <- if (cfa) M$succs.cfa else M$succs.linear
  preds <- if (cfa) M$preds.cfa else M$preds.linear
  cn    <- colnames(M$lambda)
  ovs   <- rownames(M$lambda)

  key <- list(cn = cn, ovs = ovs, lvs = lvs, modes = info$modes[lvs],
              preds = preds, succs = succs, inds = info$indsLvs[lvs])

  if (!is.null(.step05Cache$key) && identical(.step05Cache$key, key))
    return(.step05Cache$plan)

  modes <- unlist(info$modes[lvs], use.names = FALSE)
  lvCol <- match(lvs, cn)
  if (!all(modes %in% c("A", "B")) || anyNA(lvCol) || is.null(ovs))
    return(NULL)

  inds <- lapply(info$indsLvs[lvs], function(v) as.integer(match(v, ovs) - 1L))
  if (any(vapply(inds, anyNA, TRUE))) return(NULL)

  plan <- list(
    inds  = inds,
    lvCol = as.integer(lvCol - 1L),
    modeB = as.integer(modes == "B"),
    preds = matrix(as.integer(preds[cn, cn, drop = FALSE]), length(cn)),
    succs = matrix(as.integer(succs[cn, cn, drop = FALSE]), length(cn))
  )

  .step05Cache$key  <- key
  .step05Cache$plan <- plan
  plan
}


estimatePLS_Step0_5_fast <- function(model) {
  M <- model@matrices

  plan <- step05Plan(model)
  if (is.null(plan)) return(NULL)

  res <- plsStep05Cpp(
    S            = M$S,
    lambda       = M$lambda,
    gamma        = M$gamma,
    C            = M$C,
    SC           = M$SC,
    indsLvs      = plan$inds,
    lvCol        = plan$lvCol,
    modeB        = plan$modeB,
    preds        = plan$preds,
    succs        = plan$succs,
    outerWeights = as.numeric(M$outerWeights),
    tol          = model@status$tolerance,
    maxIter      = as.integer(model@status$max.iter.0_5)
  )

  dimnames(res$lambda) <- dimnames(M$lambda)
  dimnames(res$gamma)  <- dimnames(M$gamma)
  dimnames(res$C)      <- dimnames(M$C)
  dimnames(res$SC)     <- dimnames(M$SC)

  if (!res$convergence) {
    pls_msg_warn("Convergence not reached. Stopping.")
    model@status$is.admissible <- FALSE
  }

  model@matrices$lambda       <- res$lambda
  model@matrices$gamma        <- res$gamma
  model@matrices$C            <- res$C
  model@matrices$SC           <- res$SC
  model@matrices$outerWeights <- stats::setNames(
    res$outerWeights, names(M$outerWeights)
  )
  model@status$convergence    <- res$convergence
  model@status$iterations.0_5 <- model@status$iterations.0_5 + res$iterations
  model@status$iterations     <- model@status$iterations + res$iterations
  model
}


# C++ fast path for step 6. Only handles the plain (non-probit, non-higher-
# order) case; anything else falls back to the reference R code.
.step6Cache <- new.env(parent = emptyenv())

estimatePLS_Step6_fast <- function(model) {
  info <- model@info
  if (isTRUE(info$is.probit) || !isTRUE(info$is.nlin)) return(NULL)

  W     <- model@matrices$lambda
  cn    <- colnames(W)
  elems <- info$intTermElems
  if (is.null(cn) || !length(elems)) return(NULL)

  key <- list(cn = cn, elems = elems)
  plan <- .step6Cache$plan
  if (is.null(.step6Cache$key) || !identical(.step6Cache$key, key)) {
    tgt <- match(names(elems), cn)
    src <- lapply(elems, function(e) match(e, cn))
    if (anyNA(tgt) || any(vapply(src, anyNA, TRUE)) ||
        any(lengths(src) != 2L)) return(NULL)
    plan <- list(prodElems = matrix(as.integer(unlist(src) - 1L), nrow = 2L),
                 prodCol   = as.integer(tgt - 1L))
    .step6Cache$key  <- key
    .step6Cache$plan <- plan
  }

  res <- plsStep6Cpp(
    X             = model@data,
    W             = W,
    prodElems     = plan$prodElems,
    prodCol       = plan$prodCol,
    doStandardise = !isTRUE(info$standardized) || isTRUE(info$is.probit)
  )

  F <- res$F
  dimnames(F) <- list(rownames(model@data), cn)
  Cxz <- res$C
  dimnames(Cxz) <- list(cn, cn)

  model@factorScores          <- F
  model@matrices$C[cn, cn]    <- Cxz
  model@matrices$SC[cn, cn]   <- Cxz
  model
}
