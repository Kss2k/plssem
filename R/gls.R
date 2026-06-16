createGlsModel <- function(parTable, data.cov = NULL) {
  reg  <- parTable[parTable$op == "~", , drop = FALSE]
  etas <- union(reg$lhs, reg$rhs)
  xis  <- setdiff(reg$rhs, reg$lhs) # purely exogenous variables
  k    <- length(etas)

  # covariances
  cov <- parTable[
    parTable$lhs %in% etas &
    parTable$op == "~~" &
    parTable$rhs %in% etas, , drop = FALSE
  ]

  gamma <- psi <- matrix(
    0, nrow = k, ncol = k,
    dimnames = list(etas, etas)
  )

  # paths 
  for (i in seq_len(NROW(reg))) {
    lhs <- reg[i, "lhs"]
    rhs <- reg[i, "rhs"]
    gamma[lhs, rhs] <- NA
  }

  # covariances
  for (i in seq_len(NROW(cov))) {
    lhs <- reg[i, "lhs"]
    rhs <- reg[i, "rhs"]
    psi[lhs, rhs] <- psi[rhs, lhs] <- NA
  }

  diag(psi) <- NA
  psi[upper.tri(psi)] <- 0

  if (!is.null(cov)) {
    S <- data.cov[etas, etas]
    S.inv <- MASS::ginv(S)
  } else {
    S <- psi
    S[] <- NA
    S.inv <- S
  }

  # starting values
  gamma.start <- gamma
  gamma.start[] <- 0

  psi.start <- psi
  psi.start[] <- 0
  diag(psi.start) <- 1

  list(
    matrices = list(
      psi = psi,
      gamma = gamma,
      psi.free = is.na(psi),
      gamma.free = is.na(gamma),
      I = diag(k),
      S = S,
      S.inv = S.inv
    ),
    info = list(
      start = c(psi.start[is.na(psi)], gamma.start[is.na(gamma)]),
      idx.psi = seq_len(sum(is.na(psi))),
      idx.gamma = seq_len(sum(is.na(gamma))) + sum(is.na(psi)),
      k = k,
      etas = etas,
      xis = xis
    )
  )
}


glsFillModel <- function(model, par) {
  M <- model$matrices

  model$matrices$psi[model$matrices$psi.free] <- par[model$info$idx.psi]
  model$matrices$gamma[model$matrices$gamma.free] <- par[model$info$idx.gamma]

  model
}


glsCalcSigmaHat <- function(model) {
  M     <- model$matrices
  gamma <- M$gamma
  psi   <- M$psi
  I     <- M$I

  # eta = gamma * eta + zeta
  # eta - gamma * eta = zeta
  # (I - gamma) * eta = zeta
  # B * eta = zeta
  # eta = B^1 * zeta
  # E[eta*eta'] = E[(B^1 * zeta)*(B^1 * zeta)']
  # E[eta*eta'] = E[(B^1 * zeta)*(zeta' * B^1']
  # E[eta*eta'] = B^1 * E[zeta*zeta'] * B^1'
  # E[eta*eta'] = B^1 * Psi * B^1'

  binv <- solve(I - gamma)
  binv %*% psi %*% t(binv)
}


glsObjective <- function(model) {
  sigma.hat <- glsCalcSigmaHat(model)
  S <- model$matrices$S
  S.inv <- model$matrices$S.inv
  0.5 * sum(S.inv %*% (S - sigma.hat))
}
