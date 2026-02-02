getConsistenCorrMat <- function(model, P) {
  lVs          <- model$info$lVs.linear
  intTermElems <- model$info$intTermElems
  intTermNames <- model$info$intTermNames
  C            <- model$matrices$C

  for (i in seq_along(lVs)) {
    for (j in seq_len(i - 1)) {
      x <- lVs[[i]]
      y <- lVs[[j]]
      C[x, y] <- C[y, x] <- C[x, y] / sqrt(P[[x]] * P[[y]])
    }
  }

  for (i in seq_along(intTermElems)) {
    elems.xz <- intTermElems[[i]]
    xz       <- intTermNames[[i]]

    stopif(length(elems.xz) > 2, "three order interaction effects are not allowed!")
    x <- elems.xz[[1L]]
    z <- elems.xz[[2L]]

    # X, Z, and Y are the true latent variables
    # Xc, Zc, and Yc are the composite proxies
    # [1] E[Xc^2 * Zc] = Qx^2 * Qz * E[X^2 * Z]
    # [2] E[Xc * Zc * Yc] = Qx * Qz * Qy * E[X * Z * Z]
    # [3] E[Xc^2 * Zc^2] = Qx^2 * Qy^2 * (E[X^2 * Z^2] - 1) + 1
    # [4] E[Xc^2 * Zc * Yc] = Qx^2 * Qz^2 * Qy^2 * E[X^2 * Z * Y] + p(Z,Y)*Qz*Qy*(1-Qx^2)

    # First go through the normal variables (i.e., not products)
    for (y in lVs) {
      # Eq 1 and 2 simply to this
      C[xz, y] <- C[y, xz] <- C[xz, y] / sqrt(P[[x]] * P[[z]] * P[[y]])
    }

    # Get Variance of xz (Eq 3)
    p.x.z <- C[x, z]
    E.xc2.zc2 <- C[xz, xz] + p.x.z^2
    E.x2.z2 <- (E.xc2.zc2 - 1) / (P[[x]] * P[[z]]) + 1
    C[xz, xz] <- E.x2.z2 - p.x.z^2

    # Then go through the other products
    message("DEBUG: Handling of cross variances between product terms is not implemented (yet)!")
    for (j in seq_len(i-1L)) {
      elems.yw <- intTermElems[[j]]
      # yw <- intTermNames[[j]]
      # stopif(length(elems.yw) > 2, "three order interaction effects are not allowed!")
      # y <- elems.yw[[1L]]
      # w <- elems.yw[[2L]]
      #
      # # Distinct interaction terms (XZ vs. YW); cf. eq. 21 in Dijkstra et al. (2014)
      # if (length(unique(c(x, z, y, w))) == 4L) {
      #   denom <- sqrt(P[[x]] * P[[z]] * P[[y]] * P[[w]])
      #   C[xz, yw] <- C[yw, xz] <- C[xz, yw] / denom
      # }
      #
      # # TODO: Handle shared-latent cases (e.g., XZ vs. XY) using eqs. 22–23.
    }
  }

  C
}


getConsistentLoadings <- function(model, P) {
  lVs     <- model$info$lVs.linear
  indsLvs <- model$info$indsLvs
  lambda  <- model$matrices$lambda

  for (lV in lVs) {
    wq <- lambda[indsLvs[[lV]], lV]
    lambda[indsLvs[[lV]], lV] <- wq %*% (sqrt(P[[lV]]) / t(wq) %*% wq)
  }

  lambda
}


getReliabilityCoefs <- function(model) {
  gamma   <- model$matrices$gamma
  lVs     <- model$info$lVs.linear
  indsLvs <- model$info$indsLvs
  lambda  <- model$matrices$lambda
  # S = sample covariance matrix 
  S <- model$matrices$S
  P <- vector("numeric", length(lVs))
  names(P) <- lVs

  for (lV in lVs) {
    wq <- lambda[indsLvs[[lV]], lV, drop = FALSE]
    subS <- S[indsLvs[[lV]], indsLvs[[lV]], drop = FALSE]

    if (length(wq) <= 1L) {
      P[[lV]] <- 1
      next
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

  P
}
