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
