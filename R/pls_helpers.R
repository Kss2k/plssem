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


updateParamVector <- function(model, include.thresholds = TRUE) {
  estimatePLS_Step8(model, include.thresholds = include.thresholds)
}


updateEstimationStatus <- function(model) {
  estimatePLS_Status(model)
}


estimateHigherOrderChain <- function(model, include.thresholds = TRUE) {
  force(model)

  model <- estimatePLS_InnerLocal(
    model, include.thresholds = include.thresholds
  )
  model@combinedModel <- NULL

  if (!hasHigherOrderModel(model))
    return(model)

  higherOrder <- higherOrderModel(model)
  pls_stopif(is.null(higherOrder), "Expected a higher-order model")

  newdata <- getSecondOrderDataMatrix(firstOrder = model, secondOrder = higherOrder)
  modelData(higherOrder)  <- newdata
  indCorrMatrix(higherOrder) <- Rfast::cova(newdata)
  inputReliabilities(higherOrder) <- constructReliabilities(model)

  higherOrder <- estimateHigherOrderChain(
    higherOrder, include.thresholds = include.thresholds
  )
  higherOrder <- correctLoadingsAndWeightsSecondOrder(
    firstOrder = model,
    secondOrder = higherOrder
  )

  higherOrderModel(model) <- higherOrder
  model@combinedModel <- combinedModel(model)
  model
}
