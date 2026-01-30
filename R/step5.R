step5 <- function(model, convergence = 1e-5)  {
  oldOuterWeights <- model$matrices$outerWeights
  newOuterWeights <- getNonZeroElems(model$matrices$lambda)

  weightDiff <- (oldOuterWeights - newOuterWeights) / oldOuterWeights
  model$info$convergence <- all(abs(weightDiff) < convergence)

  model$matrices$outerWeights <- newOuterWeights
  model
}
