# Step 4 is structurally identical to step 0: recompute C and SC from the
# updated outer weights after step 3.
estimatePLS_Step4 <- function(model) {
  force(model)

  lambda     <- model@matrices$lambda
  partLambda <- cbind(model@matrices$Ip, lambda)
  S          <- model@matrices$S

  model@matrices$C  <- t(lambda) %*% S %*% lambda
  model@matrices$SC <- t(partLambda) %*% S %*% partLambda
  model
}
