# assuming that all lVs are reflective (i.e., using mode A)
step3 <- function(model) {
  lVs <- model$info$lVs
  indsLvs <- model$info$indsLvs
  lambda <- model$matrices$lambda
  S <- model$matrices$S
  SC <- model$matrices$SC 
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


step3_old <- function(model) {
  lVs <- model$info$lVs
  indsLvs <- model$info$indsLvs
  factorScores <- model$factorScores
  data <- model$data
  lambda <- model$matrices$lambda
  S <- model$matrices$S
  for (lV in lVs) {
    indsLv <- indsLvs[[lV]]
    for (indLv in indsLv) {
      est <- cor(factorScores[, lV], data[, indLv])
      lambda[indLv, lV] <- est
    }
    lambda[, lV] <- lambda[, lV, drop = TRUE] / c(sqrt(t(lambda[, lV]) %*% S %*% lambda[, lV]))
  }
  model$matrices$lambda <- lambda
  model
}
