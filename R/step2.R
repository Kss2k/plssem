estimatePLS_Step2 <- function(model) {
  lVs        <- model$info$lVs.linear
  Ip         <- model$matrices$Ip
  lambda     <- model$matrices$lambda
  partLambda <- cbind(Ip, lambda)
  gamma      <- model$matrices$gamma
  partGamma  <- rbind(cbind(Ip, matrix(0, nrow = nrow(Ip), ncol = ncol(gamma))),
                      cbind(matrix(0, nrow = nrow(gamma), ncol = ncol(Ip)), gamma))

  S  <- model$matrices$S
  C  <- model$matrices$C
  SC <- model$matrices$SC

  model$matrices$C <- t(gamma) %*% C %*% gamma
  model$matrices$SC <- t(partGamma) %*% t(partLambda) %*% S %*% partLambda %*% partGamma

  dimnames(model$matrices$SC) <- dimnames(SC)
  model
}
