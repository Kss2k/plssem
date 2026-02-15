estimatePLS_Step1 <- function(model) {
  lVs   <- model$info$lVs.linear
  succs <- model$matrices$succs.linear
  preds <- model$matrices$preds.linear
  gamma <- model$matrices$gamma
  S     <- model$matrices$S 
  C     <- model$matrices$C
  SC    <- model$matrices$SC 

  for (lV in lVs) {
    predsLv <- lVs[preds[ , lV, drop = TRUE]]
    succsLv <- lVs[succs[ , lV, drop = TRUE]]

    for (succ in succsLv)
      gamma[succ, lV] <- C[lV, succ]

    if (length(predsLv) > 0)
      gamma[predsLv, lV] <- solve(SC[predsLv, predsLv]) %*% SC[predsLv, lV]

    # standardize 
    scalef <- c(sqrt(t(gamma[, lV]) %*% C %*% gamma[, lV]))
    
    if (scalef)
      gamma[, lV] <- gamma[, lV] / scalef
  }

  model$matrices$gamma <- gamma
  model
}
