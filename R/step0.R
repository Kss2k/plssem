step0 <- function(model) {
  lambda <- model$matrices$lambda
  for (i in 1:ncol(lambda)) {
    lambda[, i] <- lambda[, i] / sum(lambda[, i])
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
