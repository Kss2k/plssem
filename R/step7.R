estimatePLS_Step7  <- function(model) {
  # Step 7 get baseline fit for correcting path coefficients in
  # multilevel model.

  is.mlm     <- model$info$is.mlm
  is.mcpls   <- model$info$is.mcpls
  consistent <- model$info$consistent
  is.probit  <- model$info$is.probit

  if (!is.mlm) {

    if (consistent) {
      model$fit.c <- getFitPLSModel(model, consistent = TRUE)
      model$fit.u <- list(NULL)
      model$fit   <- model$fit.c
    } else {
      model$fit.c <- list(NULL)
      model$fit.u <- getFitPLSModel(model, consistent = FALSE)
      model$fit   <- model$fit.u
    }

    return(model)
  }

  model.c <- model
  model.u <- model

  # Even if consistent is FALSE we want to correct for not using probit
  # factor scores only.
  if (is.probit || is.mcpls) {
    model.u$info$is.probit <- FALSE
    model.u$info$is.mcpls  <- FALSE

    model.u$matrices$S <- getCorrMat(model.u$data, probit = FALSE)
    model.u <- estimatePLS_Step0_5(model.u) |> estimatePLS_Step6()
  }

  model$fit.c <- getFitPLSModel(model.c, consistent = consistent)
  model$fit.u <- getFitPLSModel(model.u, consistent = FALSE)
  model$fit   <- model$fit.c

  model
}
