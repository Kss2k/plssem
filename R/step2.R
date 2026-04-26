estimatePLS_Step2 <- function(model) {
  Ip         <- model@matrices$Ip
  lambda     <- model@matrices$lambda
  gamma      <- model@matrices$gamma
  C          <- model@matrices$C
  S          <- model@matrices$S
  SC         <- model@matrices$SC

  partLambda <- cbind(Ip, lambda)
  partGamma  <- rbind(
    cbind(Ip,                                      matrix(0, nrow = nrow(Ip),    ncol = ncol(gamma))),
    cbind(matrix(0, nrow = nrow(gamma), ncol = ncol(Ip)), gamma)
  )

  newC  <- t(gamma) %*% C %*% gamma
  newSC <- t(partGamma) %*% t(partLambda) %*% S %*% partLambda %*% partGamma

  dimnames(newSC) <- dimnames(SC)

  model@matrices$C  <- newC
  model@matrices$SC <- newSC
  model
}
