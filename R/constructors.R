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


PlsSubModel <- function(matrices       = list(),
                        data           = matrix(),
                        info           = list(),
                        status         = list(),
                        params         = list(),
                        fit            = list(),
                        fitConsistent  = list(),
                        fitUncorrected = list(),
                        fitLmer        = list(),
                        factorScores   = NULL,
                        parTableInput  = data.frame()) {

  methods::new("PlsSubModel",
    matrices       = matrices,
    data           = data,
    info           = info,
    status         = status,
    params         = params,
    fit            = fit,
    fitConsistent  = fitConsistent,
    fitUncorrected = fitUncorrected,
    fitLmer        = fitLmer,
    factorScores   = factorScores,
    parTableInput  = parTableInput
  )
}


PlsModel <- function(submodels      = list(),
                     parTable       = NULL,
                     boot           = list(),
                     matrices       = list(),
                     data           = matrix(),
                     info           = list(),
                     status         = list(),
                     params         = list(),
                     fit            = list(),
                     fitConsistent  = list(),
                     fitUncorrected = list(),
                     fitLmer        = list(),
                     factorScores   = NULL) {

  methods::new("PlsModel",
    submodels      = submodels,
    parTable       = parTable,
    boot           = boot,
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
