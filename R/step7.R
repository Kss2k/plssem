estimatePLS_Step7 <- function(model) {
  is.mlm     <- model@info$is.mlm
  is.mcpls   <- model@info$is.mcpls
  consistent <- model@info$consistent
  is.probit  <- model@info$is.probit

  if (!is.mlm) {
    if (consistent) {
      modelFitConsistent(model)  <- getFitPLSModel(model, consistent = TRUE)
      modelFitUncorrected(model) <- list(NULL)
      modelFit(model)            <- modelFitConsistent(model)
    } else {
      modelFitConsistent(model)  <- list(NULL)
      modelFitUncorrected(model) <- getFitPLSModel(model, consistent = FALSE)
      modelFit(model)            <- modelFitUncorrected(model)
    }
    return(model)
  }

  model.c <- model
  model.u <- model

  if (is.probit || is.mcpls) {
    model.u@info$is.probit  <- FALSE
    model.u@info$is.mcpls   <- FALSE
    model.u@matrices$S      <- getCorrMat(model.u@data, probit = FALSE)
    model.u <- estimatePLS_Step0_5(model.u) |> estimatePLS_Step6()
  }

  modelFitConsistent(model)  <- getFitPLSModel(model.c, consistent = consistent)
  modelFitUncorrected(model) <- getFitPLSModel(model.u, consistent = FALSE)
  modelFit(model)            <- modelFitConsistent(model)
  model
}
