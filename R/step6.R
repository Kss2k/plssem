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
