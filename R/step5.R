estimatePLS_Step5 <- function(model) {
  force(model)

  oldWeights <- model@matrices$outerWeights
  newWeights <- getNonZeroElems(model@matrices$lambda)

  weightDiff <- (oldWeights - newWeights) / oldWeights
  model@status$convergence    <- all(abs(weightDiff) < model@status$tolerance)
  model@matrices$outerWeights <- newWeights
  model
}
