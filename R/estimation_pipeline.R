updateOuterWeights <- function(model) {
  resetPLS_ModelLowerOrder(model, hard.reset = TRUE) |>
    estimatePLS_Step0_5()
}


updateFactorScores <- function(model) {
  estimatePLS_Step6(model)
}


updateFitObjects <- function(model) {
  estimatePLS_Step7(model)
}


updateParamVector <- function(model) {
  estimatePLS_Step8(model)
}


updateEstimationStatus <- function(model) {
  estimatePLS_Status(model)
}

