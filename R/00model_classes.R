#' @import methods
NULL

PLS_MODEL <- "ANY" # NULL or PlsModel
PAR_TABLE <- "ANY" # NULL or data.frame/PlsSemParTable
LIST      <- "ANY" # NULL or list
MATRIX    <- "ANY"

setClass(
  "PlsModel",
  slots     = c(
    matrices         = "list",
    data             = "matrix",
    info             = "list",
    status           = "list",
    params           = "list",
    fit              = "list",
    fitConsistent    = "list",
    fitUncorrected   = "list",
    fitLmer          = LIST,     # NULL or list
    factorScores     = MATRIX,   # NULL or matrix
    parTableInput    = "data.frame",
    higherOrderModel = PLS_MODEL, # NULL or PlsModel
    combinedModel    = PLS_MODEL, # NULL or PlsModel
    parTable         = PAR_TABLE, # NULL or data.frame/PlsSemParTable
    boot             = "list"
  ),
  prototype = list(
    matrices         = list(),
    data             = matrix(numeric(0), 0, 0),
    info             = list(),
    status           = list(),
    params           = list(),
    fit              = list(),
    fitConsistent    = list(),
    fitUncorrected   = list(),
    fitLmer          = NULL,
    factorScores     = NULL,
    parTableInput    = data.frame(),
    higherOrderModel = NULL,
    combinedModel    = NULL,
    parTable         = NULL,
    boot             = list()
  )
)
