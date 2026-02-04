step0 <- function(model) {
  lambda <- model$matrices$lambda

  for (i in seq_len(ncol(lambda))) {
    Li <- sum(lambda[, i])

    if (Li <= 0) next
    lambda[, i] <- lambda[, i] / Li
  }

  partLambda <- cbind(model$matrices$Ip, lambda)
  S  <- model$matrices$S
  C  <- model$matrices$C
  SC <- model$matrices$SC

  # get new expected matrices
  model$matrices$C <- t(lambda) %*% S %*% lambda
  model$matrices$SC <- t(partLambda) %*% S %*% partLambda

  model 
}
