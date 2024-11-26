# step 1: using path scheme
step1 <- function(model) {
  lVs <- model$info$lVs
  succs <- model$matrices$succs
  preds <- model$matrices$preds
  gamma <- model$matrices$gamma
  S <- model$matrices$S 
  C <- model$matrices$C
  SC <- model$matrices$SC 

  for (lV in lVs) {
    predsLv <- lVs[preds[ , lV, drop = TRUE]]
    succsLv <- lVs[succs[ , lV, drop = TRUE]]
    for (succ in succsLv) {
      gamma[succ, lV] <- C[lV, succ]
    }
    if (length(predsLv) > 0) {
      gamma[predsLv, lV] <- solve(SC[predsLv, predsLv]) %*% SC[predsLv, lV]
    }
    # standardize 
    gamma[, lV] <- gamma[, lV, drop = TRUE] / c(sqrt(t(gamma[, lV]) %*% C %*% gamma[, lV]))
  }
  model$matrices$gamma <- gamma
  model
}


step1_old <- function(model) {
  lVs <- model$info$lVs
  factorScores <- model$factorScores
  succs <- model$matrices$succs
  preds <- model$matrices$preds
  gamma <- model$matrices$gamma
  for (lV in lVs) {
    predsLv <- lVs[preds[ , lV, drop = TRUE]]
    succsLv <- lVs[succs[ , lV, drop = TRUE]]
    for (succ in succsLv) {
      gamma[succ, lV] <- cor(factorScores[, succ], factorScores[, lV])
    }
    if (length(predsLv) > 0) {
      gamma[predsLv, lV] <- lm(as.data.frame(factorScores[, c(lV, predsLv)]))$coefficients[-1]
    }
    # standardize 
    gamma[, lV] <- gamma[, lV, drop = TRUE] / c(sqrt(t(gamma[, lV]) %*% cov(as.data.frame(factorScores)) %*% gamma[, lV]))
  }
  model$matrices$gamma <- gamma
  model
}
