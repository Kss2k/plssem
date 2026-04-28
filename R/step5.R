estimatePLS_Step5 <- function(model)  {
  oldOuterWeights <- model$matrices$outerWeights
  newOuterWeights <- getNonZeroElems(model$matrices$lambda)

  weightDiff <- (oldOuterWeights - newOuterWeights) / oldOuterWeights
  model$status$convergence <- all(abs(weightDiff) < model$status$tolerance)

  model$matrices$outerWeights <- newOuterWeights
  model
}
