estimatePLS_Step8 <- function(model) {
  model$params$values <- extractCoefs(model)
  model$params$se     <- rep(NA_real_, length(model$params$values))

  if (model$info$is.multilevel) {
    model$fit.lmer <- plslmer(model)

    coefs.x <- model$params$values
    coefs.y <- model$fit.lmer$values

    common  <- intersect(names(coefs.x), names(coefs.y))
    new     <- setdiff(names(coefs.y), names(coefs.x))

    coefs.x[common] <- coefs.y[common]
    coefs.all       <- c(coefs.x, coefs.y[new])

    model$params$values <- plssemVector(coefs.all)
    model$params$se     <- rep(NA_real_, length(coefs.all))
  }

  model
}
