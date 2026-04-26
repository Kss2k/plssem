estimatePLS_Step8 <- function(model, update.names = FALSE) {
  model <- updateModelParams(model, update.names = update.names)

  if (!model@info$is.mlm)
    return(model)

  lmer.result         <- plslmer(model)
  modelFitLmer(model) <- lmer.result

  coefs.x <- model@params$values
  coefs.y <- lmer.result$values

  common  <- intersect(names(coefs.x), names(coefs.y))
  new     <- setdiff(names(coefs.y),   names(coefs.x))

  coefs.x[common] <- coefs.y[common]
  coefs.all       <- c(coefs.x, coefs.y[new])

  model@params$values <- plssemVector(coefs.all)
  model@params$se     <- rep(NA_real_, length(coefs.all))
  model
}


updateModelParams <- function(model, update.names = FALSE) {
  if (update.names)
    model@params$names <- getParamVecNames(model)

  model@params$values <- extractCoefs(model)
  model@params$se     <- rep(NA_real_, length(model@params$values))

  model
}
