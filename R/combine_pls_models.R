combinePLS_Models <- function(models) {
  stopif(!length(models), "Expected at least one model!")

  model0 <- models[[1L]]

  browser()

}


averageMatVec <- function(matrices) {
  M <- 0
  n  <- length(matrices)

  for (m in matrices)
    M <- M + m

  M / n
}


averageMatVecField <- function(models, field) {
  averageMatVec(lapply(models, FUN = \(mod) mod[[field]]))
}


averageMatVecFieldField <- function(models, field0, field1) {
  averageMatVec(lapply(models, FUN = \(mod) mod[[field0]][[field1]]))
}
