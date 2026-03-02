estimatePLS_Step1 <- function(model) {
  lvs   <- model$info$lvs.linear
  succs <- model$matrices$succs.linear
  preds <- model$matrices$preds.linear
  gamma <- model$matrices$gamma
  S     <- model$matrices$S 
  C     <- model$matrices$C
  SC    <- model$matrices$SC 

  for (lv in lvs) {
    predsLv <- lvs[preds[ , lv, drop = TRUE]]
    succsLv <- lvs[succs[ , lv, drop = TRUE]]

    for (succ in succsLv)
      gamma[succ, lv] <- C[lv, succ]

    if (length(predsLv) > 0)
      gamma[predsLv, lv] <- solve(SC[predsLv, predsLv]) %*% SC[predsLv, lv]

    # standardize 
    scalef <- c(sqrt(t(gamma[, lv]) %*% C %*% gamma[, lv]))
    
    if (scalef)
      gamma[, lv] <- gamma[, lv] / scalef
  }

  model$matrices$gamma <- gamma
  model
}
