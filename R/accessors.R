setGeneric("modelData", function(object) standardGeneric("modelData"))
setGeneric("modelData<-", function(object, value) standardGeneric("modelData<-"))

setMethod("modelData", "PlsBaseModel", function(object) object@data)
setMethod("modelData<-", "PlsBaseModel", function(object, value) {
  object@data <- value
  object
})


setGeneric("modelMatrices", function(object) standardGeneric("modelMatrices"))
setGeneric("modelMatrices<-", function(object, value) standardGeneric("modelMatrices<-"))

setMethod("modelMatrices", "PlsBaseModel", function(object) object@matrices)
setMethod("modelMatrices<-", "PlsBaseModel", function(object, value) {
  object@matrices <- value
  object
})


setGeneric("modelInfo", function(object) standardGeneric("modelInfo"))
setGeneric("modelInfo<-", function(object, value) standardGeneric("modelInfo<-"))

setMethod("modelInfo", "PlsBaseModel", function(object) object@info)
setMethod("modelInfo<-", "PlsBaseModel", function(object, value) {
  object@info <- value
  object
})


setGeneric("modelStatus", function(object) standardGeneric("modelStatus"))
setGeneric("modelStatus<-", function(object, value) standardGeneric("modelStatus<-"))

setMethod("modelStatus", "PlsBaseModel", function(object) object@status)
setMethod("modelStatus<-", "PlsBaseModel", function(object, value) {
  object@status <- value
  object
})


setGeneric("modelParams", function(object) standardGeneric("modelParams"))
setGeneric("modelParams<-", function(object, value) standardGeneric("modelParams<-"))

setMethod("modelParams", "PlsBaseModel", function(object) object@params)
setMethod("modelParams<-", "PlsBaseModel", function(object, value) {
  object@params <- value
  object
})


setGeneric("modelFit", function(object) standardGeneric("modelFit"))
setGeneric("modelFit<-", function(object, value) standardGeneric("modelFit<-"))

setMethod("modelFit", "PlsBaseModel", function(object) object@fit)
setMethod("modelFit<-", "PlsBaseModel", function(object, value) {
  object@fit <- value
  object
})


setGeneric("modelFitConsistent", function(object) standardGeneric("modelFitConsistent"))
setGeneric("modelFitConsistent<-", function(object, value) standardGeneric("modelFitConsistent<-"))

setMethod("modelFitConsistent", "PlsBaseModel", function(object) object@fitConsistent)
setMethod("modelFitConsistent<-", "PlsBaseModel", function(object, value) {
  object@fitConsistent <- value
  object
})


setGeneric("modelFitUncorrected", function(object) standardGeneric("modelFitUncorrected"))
setGeneric("modelFitUncorrected<-", function(object, value) standardGeneric("modelFitUncorrected<-"))

setMethod("modelFitUncorrected", "PlsBaseModel", function(object) object@fitUncorrected)
setMethod("modelFitUncorrected<-", "PlsBaseModel", function(object, value) {
  object@fitUncorrected <- value
  object
})


setGeneric("modelFitLmer", function(object) standardGeneric("modelFitLmer"))
setGeneric("modelFitLmer<-", function(object, value) standardGeneric("modelFitLmer<-"))

setMethod("modelFitLmer", "PlsBaseModel", function(object) object@fitLmer)
setMethod("modelFitLmer<-", "PlsBaseModel", function(object, value) {
  object@fitLmer <- value
  object
})


setGeneric("modelFactorScores", function(object) standardGeneric("modelFactorScores"))
setGeneric("modelFactorScores<-", function(object, value) standardGeneric("modelFactorScores<-"))

setMethod("modelFactorScores", "PlsBaseModel", function(object) object@factorScores)
setMethod("modelFactorScores<-", "PlsBaseModel", function(object, value) {
  object@factorScores <- value
  object
})


setGeneric("parTableInput", function(object) standardGeneric("parTableInput"))
setGeneric("parTableInput<-", function(object, value) standardGeneric("parTableInput<-"))

setMethod("parTableInput", "PlsSubModel", function(object) object@parTableInput)
setMethod("parTableInput<-", "PlsSubModel", function(object, value) {
  object@parTableInput <- value
  object
})


setGeneric("submodels", function(object) standardGeneric("submodels"))
setGeneric("submodels<-", function(object, value) standardGeneric("submodels<-"))

setMethod("submodels", "PlsModel", function(object) object@submodels)
setMethod("submodels<-", "PlsModel", function(object, value) {
  object@submodels <- value
  object
})


setGeneric("firstOrder", function(object) standardGeneric("firstOrder"))
setGeneric("firstOrder<-", function(object, value) standardGeneric("firstOrder<-"))

setMethod("firstOrder", "PlsModel", function(object) object@submodels$firstOrder)
setMethod("firstOrder<-", "PlsModel", function(object, value) {
  object@submodels$firstOrder <- value
  object
})


setGeneric("secondOrder", function(object) standardGeneric("secondOrder"))
setGeneric("secondOrder<-", function(object, value) standardGeneric("secondOrder<-"))

setMethod("secondOrder", "PlsModel", function(object) object@submodels$secondOrder)
setMethod("secondOrder<-", "PlsModel", function(object, value) {
  object@submodels$secondOrder <- value
  object
})


setGeneric("modelParTable", function(object) standardGeneric("modelParTable"))
setGeneric("modelParTable<-", function(object, value) standardGeneric("modelParTable<-"))

setMethod("modelParTable", "PlsModel", function(object) object@parTable)
setMethod("modelParTable<-", "PlsModel", function(object, value) {
  object@parTable <- value
  object
})


setGeneric("modelBoot", function(object) standardGeneric("modelBoot"))
setGeneric("modelBoot<-", function(object, value) standardGeneric("modelBoot<-"))

setMethod("modelBoot", "PlsModel", function(object) object@boot)
setMethod("modelBoot<-", "PlsModel", function(object, value) {
  object@boot <- value

  # Should also update se and vcov in model params
  object@params$se   <- value$se
  object@params$vcov <- value$vcov

  object
})


setGeneric("isAdmissible", function(object) standardGeneric("isAdmissible"))
setGeneric("isAdmissible<-", function(object, value) standardGeneric("isAdmissible<-"))

setMethod("isAdmissible", "PlsBaseModel", function(object) isTRUE(modelStatus(object)$is.admissible))
setMethod("isAdmissible<-", "PlsBaseModel", function(object, value) {
  object@status$is.admissible <- value
  object
})


setGeneric("isMLM", function(object) standardGeneric("isMLM"))
setMethod("isMLM", "PlsBaseModel", function(object) {
  isTRUE(object@info$is.mlm)
})


setGeneric("constructReliabilities", function(object) standardGeneric("constructReliabilities"))
setMethod("constructReliabilities", "PlsBaseModel", function(object) modelFit(object)$Q^2)

setGeneric("inputReliabilities", function(object) standardGeneric("inputReliabilities"))
setGeneric("inputReliabilities<-", function(object, value) standardGeneric("inputReliabilities<-"))

setMethod("inputReliabilities", "PlsBaseModel", function(object) modelInfo(object)$reliabilities)
setMethod("inputReliabilities<-", "PlsBaseModel", function(object, value) {
  object@info$reliabilities <- value
  object
})


setGeneric("corrMatrix", function(object) standardGeneric("corrMatrix")) 
setGeneric("corrMatrix<-", function(object, value) standardGeneric("corrMatrix<-"))

setMethod("corrMatrix", "PlsBaseModel", function(object) object@matrices$S) 
setMethod("corrMatrix<-", "PlsBaseModel", function(object, value) {
  object@matrices$S <- value
  object
}) 
