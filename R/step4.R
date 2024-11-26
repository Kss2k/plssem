# basically just step 0
step4 <- function(model) {
  lambda <- model$matrices$lambda
  partLambda <- cbind(model$matrices$Ip, lambda)
  S <- model$matrices$S 
  C <- model$matrices$C 
  SC <- model$matrices$SC
  
  # get new expected matrices
  model$matrices$C <- t(lambda) %*% S %*% lambda
  model$matrices$SC <- t(partLambda) %*% S %*% partLambda
  model 
}


step4_old <- function(model) {
  lambda <- model$matrices$lambda
  data <- model$data 
  factorScores <- data %*% lambda
  for (intTerm in names(model$info$interactionPairs)) {
    intPair <- model$info$interactionPairs[[intTerm]]
    factorScores[, intTerm] <- factorScores[, intPair[[1]]] * 
      factorScores[, intPair[[2]]]
  }
  model$factorScores <- factorScores 
  model
}
