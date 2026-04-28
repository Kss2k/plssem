PlsBaseModel <- function(matrices       = matrix(),
                         data           = matrix(),
                         info           = list(),
                         status         = list(),
                         params         = list(),
                         fit            = list(),
                         fitConsistent  = list(),
                         fitUncorrected = list(),
                         fitLmer        = list(),
                         factorScores   = NULL) {

  methods::new("PlsBaseModel",
    matrices       = matrices,
    data           = data,
    info           = info,
    status         = status,
    params         = params,
    fit            = fit,
    fitConsistent  = fitConsistent,
    fitUncorrected = fitUncorrected,
    fitLmer        = fitLmer,
    factorScores   = factorScores
  )
}


PlsModel <- function(matrices        = list(),
                     data            = matrix(),
                     info            = list(),
                     status          = list(),
                     params          = list(),
                     fit             = list(),
                     fitConsistent   = list(),
                     fitUncorrected  = list(),
                     fitLmer         = list(),
                     factorScores    = NULL,
                     parTableInput   = data.frame(),
                     higherOrderModel = NULL,
                     combinedModel   = NULL,
                     parTable        = NULL,
                     boot            = list()) {

  methods::new("PlsModel",
    matrices         = matrices,
    data             = data,
    info             = info,
    status           = status,
    params           = params,
    fit              = fit,
    fitConsistent    = fitConsistent,
    fitUncorrected   = fitUncorrected,
    fitLmer          = fitLmer,
    factorScores     = factorScores,
    parTableInput    = parTableInput,
    higherOrderModel = higherOrderModel,
    combinedModel    = combinedModel,
    parTable         = parTable,
    boot             = boot
  )
}
