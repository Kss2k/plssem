modelData <- function(object) {
  object@data
}


`modelData<-` <- function(object, value) {
  object@data <- value
  object
}


modelMatrices <- function(object) {
  object@matrices
}


`modelMatrices<-` <- function(object, value) {
  object@matrices <- value
  object
}


modelInfo <- function(object) {
  object@info
}


`modelInfo<-` <- function(object, value) {
  object@info <- value
  object
}


modelStatus <- function(object) {
  object@status
}


`modelStatus<-` <- function(object, value) {
  object@status <- value
  object
}


modelParams <- function(object) {
  object@params
}


`modelParams<-` <- function(object, value) {
  object@params <- value
  object
}


modelFit <- function(object) {
  object@fit
}


`modelFit<-` <- function(object, value) {
  object@fit <- value
  object
}


modelFitConsistent <- function(object) {
  object@fitConsistent
}


`modelFitConsistent<-` <- function(object, value) {
  object@fitConsistent <- value
  object
}


modelFitUncorrected <- function(object) {
  object@fitUncorrected
}


`modelFitUncorrected<-` <- function(object, value) {
  object@fitUncorrected <- value
  object
}


modelFitLmer <- function(object) {
  object@fitLmer
}


`modelFitLmer<-` <- function(object, value) {
  object@fitLmer <- value
  object
}


modelFactorScores <- function(object) {
  object@factorScores
}


`modelFactorScores<-` <- function(object, value) {
  object@factorScores <- value
  object
}


parTableInput <- function(object) {
  object@parTableInput
}


`parTableInput<-` <- function(object, value) {
  object@parTableInput <- value
  object
}


higherOrderModel <- function(object) {
  object@higherOrderModel
}


`higherOrderModel<-` <- function(object, value) {
  object@higherOrderModel <- value
  object
}


hasHigherOrderModel <- function(object) {
  !is.null(object@higherOrderModel)
}


combinedModel <- function(object, refresh = FALSE) {
  if (!refresh && is(object@combinedModel, "PlsModel"))
    return(object@combinedModel)

  if (!hasHigherOrderModel(object))
    return(object)

  computeCombinedModel(object)
}


rootModel <- function(object) {
  object@combinedModel <- NULL
  object@higherOrderModel <- NULL
  object
}


# Backwards-compatible accessors (two-level higher-order API).
firstOrder <- function(object) {
  object
}


secondOrder <- function(object) {
  object@higherOrderModel
}


modelParTable <- function(object) {
  object@parTable
}


`modelParTable<-` <- function(object, value) {
  object@parTable <- value
  object
}


modelBoot <- function(object) {
  object@boot
}


`modelBoot<-` <- function(object, value) {
  object@boot <- value

  # Also update se and vcov in model params
  object@params$se   <- value$se
  object@params$vcov <- value$vcov

  if (is(object@combinedModel, "PlsModel")) {
    object@combinedModel@boot <- value
    object@combinedModel@params$se   <- value$se
    object@combinedModel@params$vcov <- value$vcov
  }

  object
}


isAdmissible <- function(object) {
  isTRUE(modelStatus(object)$is.admissible)
}


`isAdmissible<-` <- function(object, value) {
  object@status$is.admissible <- value
  object
}


isMLM <- function(object) {
  isTRUE(object@info$is.mlm)
}


isMCPLS <- function(object) {
  combined <- combinedModel(object)
  higherOrder <- higherOrderModel(object)

  if      (!is.null(combined))    isTRUE(combined@info$is.mcpls)
  else if (!is.null(higherOrder)) isMCPLS(higherOrder)
  else                            isTRUE(object@info$is.mlm)
}


constructReliabilities <- function(object) {
  modelFit(object)$Q^2
}


inputReliabilities <- function(object) {
  modelInfo(object)$reliabilities
}


`inputReliabilities<-` <- function(object, value) {
  object@info$reliabilities <- value
  object
}


corrMatrix <- function(object) {
  object@matrices$S
}


`corrMatrix<-` <- function(object, value) {
  object@matrices$S <- value
  object
}

