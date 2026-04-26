#' @import methods
NULL


setClass(
  "PlsBaseModel",
  contains  = "VIRTUAL",
  slots     = c(
    matrices       = "list",
    data           = "matrix",
    info           = "list",
    status         = "list",
    params         = "list",
    fit            = "list",
    fitConsistent  = "list",
    fitUncorrected = "list",
    fitLmer        = "ANY",  # NULL or list
    factorScores   = "ANY"   # NULL or matrix
  ),
  prototype = list(
    matrices       = list(),
    data           = matrix(numeric(0), 0, 0),
    info           = list(),
    status         = list(),
    params         = list(),
    fit            = list(),
    fitConsistent  = list(),
    fitUncorrected = list(),
    fitLmer        = NULL,
    factorScores   = NULL
  )
)


setClass(
  "PlsSubModel",
  contains  = "PlsBaseModel",
  slots     = c(
    parTableInput = "data.frame"
  ),
  prototype = list(
    parTableInput = data.frame()
  )
)


#' S4 class for fitted PLS models
#'
#' @slot matrices Model matrices.
#' @slot data Model data matrix.
#' @slot info Model information.
#' @slot status Model status. 
#' @slot params Model parameters.
#' @slot fit Model fit matrices.
#' @slot fitConsistent Consistent model fit matrices (if \code{consistent=TRUE}).
#' @slot fitUncorrected Uncorrected model fit matrices.
#' @slot fitLmer Mixed effects model fit for multilevel models.
#' @slot factorScores factorScores matrix.
#' @slot parTable Parameter table.
#' @slot submodels First and (if applicable) second order sub models.
#' @slot boot Bootstrap results.
#' 
#' @export
setClass(
  "PlsModel",
  contains  = "PlsBaseModel",
  slots     = c(
    submodels = "list",
    parTable  = "ANY",   # NULL or data.frame/PlsSemParTable
    boot      = "list"
  ),
  prototype = list(
    submodels = list(),
    parTable  = NULL,
    boot      = list()
  )
)
