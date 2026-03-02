# assuming that all lvs are reflective (i.e., using mode A)
estimatePLS_Step3 <- function(model) {
  lvs     <- model$info$lvs.linear
  indsLvs <- model$info$indsLvs
  lambda  <- model$matrices$lambda
  S       <- model$matrices$S
  SC      <- model$matrices$SC

  for (lv in lvs) {
    indsLv <- indsLvs[[lv]]

    for (indLv in indsLv) {
      lambda[indLv, lv] <- SC[indLv, lv]
    }
    lambda[, lv] <- lambda[, lv, drop = TRUE] / c(sqrt(t(lambda[, lv]) %*% S %*% lambda[, lv]))
  }

  model$matrices$lambda <- lambda
  model
}
