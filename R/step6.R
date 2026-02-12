estimatePLS_Step6 <- function(model) {
  model$factorScores <- getFactorScores(model)

  if (model$info$is.nlin) {
    # Update variance and coefficients of interaction terms
    elems <- model$info$intTermElems
    names <- names(elems)
    X     <- as.data.frame(model$factorScores)

    SC <- model$matrices$SC
    C  <- model$matrices$C
    L  <- model$matrices$lambda
    G  <- model$matrices$gamma

    for (elems.xz in elems) {
      xz <- paste0(elems.xz, collapse = ":")
      values <- multiplyIndicatorsCpp(X[elems.xz])
      model$factorScores[,xz] <- X[[xz]] <- values
    }

    Cxz <- stats::cov(X)
    par <- colnames(X)

    model$matrices$C[par, par] <- Cxz
    model$matrices$SC[par, par] <- Cxz
  }

  model
}
