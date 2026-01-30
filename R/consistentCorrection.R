getConsistenCorrMat <- function(model, P) {
  C <- model$matrices$C

  for (i in seq_len(nrow(C))) {
    for (j in seq_len(i - 1)) {
      x <- rownames(C)[[i]]
      y <- colnames(C)[[j]]
      C[x, y] <- C[y, x] <- C[x, y] / sqrt(P[[x]] * P[[y]])
    }
  }

  C
}


getConsistentLoadings <- function(model, P) {
  lVs <- model$info$lVs
  indsLvs <- model$info$indsLvs
  lambda <- model$matrices$lambda

  for (lV in lVs) {
    wq <- lambda[indsLvs[[lV]], lV]
    lambda[indsLvs[[lV]], lV] <- wq %*% (sqrt(P[[lV]]) / t(wq) %*% wq)
  }

  lambda
}


getReliabilityCoefs <- function(model) {
  gamma <- model$matrices$gamma
  lVs <- model$info$lVs
  indsLvs <- model$info$indsLvs
  lambda <- model$matrices$lambda
  # S = sample covariance matrix 
  S <- model$matrices$S
  P <- vector("numeric", length(lVs))
  names(P) <- lVs
  for (lV in lVs) {
    if (grepl(":", lV)) {
      intPair <- model$info$interactionPairs[[lV]]
      subS <- covProdInds(indsLvs[[intPair[[1]]]], 
                          indsLvs[[intPair[[2]]]], model$data)
      wq <- weightsProdInds(lambda[indsLvs[[intPair[[1]]]], intPair[[1]]],
                            lambda[indsLvs[[intPair[[2]]]], intPair[[2]]])
      # standardize wq
      wq <- wq <- wq / as.vector(sqrt(t(wq) %*% subS %*% wq))
    } else {
      wq <- lambda[indsLvs[[lV]], lV]
      subS <- S[indsLvs[[lV]], indsLvs[[lV]]]
    }
    A <- subS
    diag(A) <- 0
    B <- wq %*% t(wq)
    diag(B) <- 0
    D <- t(wq) %*% A %*% wq 
    N <- t(wq) %*% B %*% wq 
    s <- (t(wq) %*% wq) ^ 2
    P[[lV]] <- s * D / N
  }
  # for (intTerm in names(model$info$interactionPairs)) {
  #   intPair <- model$info$interactionPairs[[intTerm]]
  #   P[[intTerm]] <- P[[intPair[[1]]]] * P[[intPair[[2]]]]
  # }
  P
}
