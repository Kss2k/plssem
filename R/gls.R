setClass(
  "GlsPathModel",
  slots     = c(
    matrices  = "list",
    info      = "list",
    parTable  = "data.frame"
  ),
  prototype = list(
    matrices  = list(
      psi         = NULL,
      gamma       = NULL,
      psi.free    = NULL,
      gamma.free  = NULL,
      I           = NULL,
      S           = NULL,
      S.inv       = NULL
    ),
    info = list(
      start     = numeric(0L),
      npar      = 0,
      idx.psi   = integer(0L),
      idx.gamma = integer(0L),
      k         = 0,
      etas      = character(0L),
      xis       = character(0L)
    ),
    parTable = data.frame()
  )
)


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
  psi[xis, xis] <- NA
  psi[upper.tri(psi)] <- 0

  if (!is.null(data.cov)) {
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

  # parTable
  pgamma <- gamma
  pgamma[] <- ""
  for (i in etas) for (j in etas) {
    if (!is.na(gamma[i,j])) next
    pgamma[i, j] <- paste0(i, "~", j)
  }
  
  ppsi <- psi
  ppsi[] <- ""
  for (i in etas) for (j in etas) {
    if (!is.na(psi[i,j])) next
    ppsi[i, j] <- paste0(i, "~~", j)
  }

  pars <- c(ppsi[is.na(psi)], pgamma[is.na(gamma)])
  parTable <- as.data.frame(splitParameterNames(pars))
  parTable$est <- NA_real_

  methods::new("GlsPathModel",
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
      npar = sum(is.na(psi)) + sum(is.na(gamma)),
      idx.psi = seq_len(sum(is.na(psi))),
      idx.gamma = seq_len(sum(is.na(gamma))) + sum(is.na(psi)),
      k = k,
      etas = etas,
      xis = xis
    ),
    parTable = parTable
  )
}


`glsModelCovMatrix<-` <- function(model, value) {
  etas <- model@info$etas
  S    <- value[etas, etas]

  model@matrices$S <- S
  model@matrices$S.inv <- MASS::ginv(S)

  model
}


glsFillModel <- function(model, par) {
  M <- model@matrices
  model@parTable$est[] <- par

  psi <- model@matrices$psi
  gamma <- model@matrices$gamma

  psi[model@matrices$psi.free] <- par[model@info$idx.psi]
  psi[upper.tri(psi)] <- t(psi)[upper.tri(psi)]
  gamma[model@matrices$gamma.free] <- par[model@info$idx.gamma]

  model@matrices$psi <- psi
  model@matrices$gamma <- gamma

  model
}


glsCalcSigmaHat <- function(model) {
  M     <- model@matrices
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
  # Full equation
  # Sigma.hat = inv(I - gamma) * Psi * inv(I - gamma)'
  # tmp = inv(S) * (S - Sigma.hat)
  # tmp = inv(S) * (S - inv(I - gamma) * Psi * inv(I - gamma)')
  # objective = 0.5 * sum(tmp * tmp')
  # objective = 0.5 * sum((inv(S) * (S - inv(I - gamma) * Psi * inv(I - gamma)')) .* (inv(S) * (S - inv(I - gamma) * Psi * inv(I - gamma)'))')
  sigma.hat <- glsCalcSigmaHat(model)
  S <- model@matrices$S
  S.inv <- model@matrices$S.inv
  tmp <- S.inv %*% (S - sigma.hat)
  0.5 * sum(tmp * t(tmp))
}


glsGradientPsi <- function(matrices) {
  I <- matrices$I
  S <- matrices$S
  gamma <- matrices$gamma
  psi   <- matrices$psi

  binv <- solve(I-gamma)
  sinv <- matrices$S.inv
  tbsinv <- t(binv) %*% sinv
  eps <- S - binv %*% psi %*% t(binv)
  dv <- diag(rep(0.5, NCOL(S)))

  tmp <- -0.5 * (
    t(sinv %*% binv) %*% (t(S)-binv %*% t(psi) %*% t(binv)) %*% dv %*% t(sinv) %*% binv +
    tbsinv %*% eps %*% sinv %*% dv %*% binv +
    tbsinv %*% dv %*% eps %*% sinv %*% binv +
    t(binv) %*% dv %*% sinv %*% eps %*% t(sinv) %*% binv
  )

  (2 * tmp - diag(diag(tmp)))
}


glsGradientGamma <- function(matrices) {
  I <- matrices$I
  S <- matrices$S
  gamma <- matrices$gamma
  psi   <- matrices$psi

  binv <- solve(I-gamma)
  sinv <- matrices$S.inv
  eps <- S - binv %*% psi %*% t(binv)
  dv <- diag(rep(0.5, NCOL(sinv)))

  -(
    2 * t(binv) %*% sinv %*% dv %*% eps %*% sinv %*% binv %*% psi %*% t(binv) +
    t(binv) %*% dv %*% t(sinv) %*% (t(S)-binv %*% t(psi) %*% t(binv)) %*% sinv %*% binv %*% psi %*% t(binv)+
    t(binv) %*% dv %*% sinv %*% eps %*% t(sinv) %*% binv %*% t(psi) %*% t(binv)
  )
}


glsGradient <- function(model) {
  M <- model@matrices

  ggamma <- glsGradientGamma(M)
  gpsi   <- glsGradientPsi(M)

  grad <- numeric(model@info$npar)
  c(gpsi[M$psi.free], ggamma[M$gamma.free])
}


glsEstimateModel <- function(model, data.cov = NULL, ...) {
  if (!is.null(data.cov)) # update sample covariance matrix?
    glsModelCovMatrix(model) <- data.cov 

  # objective and gradient
  fn <- \(x) glsObjective(glsFillModel(model, x))
  gr <- \(x) glsGradient(glsFillModel(model, x))

  opt <- nlminb(
    start = model@info$start,
    objective = fn,
    gradient = gr,
    ...
  )

  glsFillModel(model, opt$par)
}

