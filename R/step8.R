estimatePLS_Step8 <- function(model) {
  model@params$values <- extractCoefs(model)
  model@params$se     <- rep(NA_real_, length(model@params$values))

  if (!isMLM(model))
    return(model)

  modelFitLmer(model) <- plslmer(model)

  refreshLmerParams(model) # Update params with Mixed-Effects coefficients
}
