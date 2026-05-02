# PLS estimation steps

estimatePLS_Step0_5 <- function(model) {
  force(model)

  max.iter.0_5 <- model@status$max.iter.0_5

  model <- estimatePLS_Step0(model)

  for (i in seq_len(max.iter.0_5)) {
    model <- model |>
      estimatePLS_Step1() |>
      estimatePLS_Step2() |>
      estimatePLS_Step3() |>
      estimatePLS_Step4() |>
      estimatePLS_Step5()

    if (model@status$convergence) {
      break
    } else if (i >= max.iter.0_5) {
      warning("Convergence not reached. Stopping.")
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

  lambda <- model@matrices$lambda

  for (i in seq_len(ncol(lambda))) {
    Li <- sum(lambda[, i])
    if (Li <= 0) next
    lambda[, i] <- lambda[, i] / Li
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
  getPathCoefs(y = lv, X = inds, C = SC)
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

  weightDiff <- (oldWeights - newWeights) / oldWeights
  model@status$convergence    <- all(abs(weightDiff) < model@status$tolerance)
  model@matrices$outerWeights <- newWeights
  model
}


estimatePLS_Step6 <- function(model) {
  force(model)

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
    model.u@matrices$S      <- getCorrMat(model.u@data, probit = FALSE)
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

  modelFitLmer(model) <- plslmer(model)

  refreshLmerParams(model) # Update params with Mixed-Effects coefficients
}
