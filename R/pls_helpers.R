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


estimateHigherOrderChain <- function(model) {
  force(model)

  model <- estimatePLS_InnerLocal(model)
  model@combinedModel <- NULL

  if (!hasHigherOrderModel(model))
    return(model)

  higherOrder <- higherOrderModel(model)
  stopif(is.null(higherOrder), "Expected a higher-order model")

  newdata <- getSecondOrderDataMatrix(firstOrder = model, secondOrder = higherOrder)
  modelData(higherOrder)  <- newdata
  indCorrMatrix(higherOrder) <- Rfast::cova(newdata)
  inputReliabilities(higherOrder) <- constructReliabilities(model)

  higherOrder <- estimateHigherOrderChain(higherOrder)
  higherOrder <- correctLoadingsAndWeightsSecondOrder(
    firstOrder = model,
    secondOrder = higherOrder
  )

  higherOrderModel(model) <- higherOrder
  model@combinedModel <- combinedModel(model)
  model
}
