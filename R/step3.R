estimatePLS_Step3 <- function(model) {
  lvs     <- model@info$lvs.linear
  indsLvs <- model@info$indsLvs
  lambda  <- model@matrices$lambda
  SC      <- model@matrices$SC
  modes   <- model@info$modes

  for (lv in lvs) {
    mode.lv <- modes[[lv]]
    inds    <- indsLvs[[lv]]

    wj <- switch(mode.lv,
      A = getWeightsModeA(lv = lv, lambda = lambda, SC = SC, inds = inds),
      B = getWeightsModeB(lv = lv, lambda = lambda, SC = SC, inds = inds),
      NA_real_
    )

    Sjj <- SC[inds, inds]
    wj  <- wj / c(sqrt(t(wj) %*% Sjj %*% wj))
    lambda[inds, lv] <- wj
  }

  model@matrices$lambda <- lambda
  model
}


getWeightsModeA <- function(lv, lambda, SC, inds) {
  as.vector(SC[inds, lv])
}


getWeightsModeB <- function(lv, lambda, SC, inds) {
  getPathCoefs(y = lv, X = inds, C = SC)
}
