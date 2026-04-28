initModelStatus <- function(tolerance, max.iter.0_5) {
  list(
    convergence    = FALSE,
    iterations     = 0L,
    iterations.0_5 = 0L,
    tolerance      = tolerance,
    max.iter.0_5   = max.iter.0_5,
    is.admissible  = TRUE
  )
}


initModelParams <- function(model) {
  parnames <- getParamVecNames(model)
  k <- length(parnames)

  model@params <- list(
    names      = parnames,
    values     = rep(NA_real_, k),
    values.old = NULL,
    se         = rep(NA_real_, k),
    vcov       = NULL
  )

  model
}
