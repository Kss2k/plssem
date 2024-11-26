step5 <- function(model, convergence = 1e-5)  {
  oldOuterWeights <- model$matrices$outerWeights
  newOuterWeights <- getNonZeroElems(model$matrices$lambda)
  if (all(abs((oldOuterWeights - newOuterWeights) / oldOuterWeights) < convergence)) {
    model$convergence <- TRUE 
  } else {
    model$convergence <- FALSE
  }
  model$matrices$outerWeights <- newOuterWeights
  model
}
