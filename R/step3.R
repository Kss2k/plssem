# assuming that all lVs are reflective (i.e., using mode A)
estimatePLS_Step3 <- function(model) {
  lVs     <- model$info$lVs.linear
  indsLvs <- model$info$indsLvs
  lambda  <- model$matrices$lambda
  S       <- model$matrices$S
  SC      <- model$matrices$SC

  for (lV in lVs) {
    indsLv <- indsLvs[[lV]]

    for (indLv in indsLv) {
      lambda[indLv, lV] <- SC[indLv, lV]
    }
    lambda[, lV] <- lambda[, lV, drop = TRUE] / c(sqrt(t(lambda[, lV]) %*% S %*% lambda[, lV]))
  }

  model$matrices$lambda <- lambda
  model
}
