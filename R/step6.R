estimatePLS_Step6 <- function(model) {
  model$factorScores <- getFactorScores(model)

  if (model$info$is.nlin) {
    # Update variance and coefficients of interaction terms
    elems <- model$info$intTermElems
    names <- names(elems)
    X     <- model$factorScores

    SC <- model$matrices$SC
    C  <- model$matrices$C
    L  <- model$matrices$lambda
    G  <- model$matrices$gamma

    for (elems.xz in elems) {
      xz <- paste0(elems.xz, collapse = ":")
      X[,xz] <- Rfast::rowprods(X[,elems.xz])
    }

    Cxz <- Rfast::cova(X)
    par <- colnames(X)

    model$factorScores <- X
    model$matrices$C[par, par] <- Cxz
    model$matrices$SC[par, par] <- Cxz
  }

  model
}
