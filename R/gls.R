GlsPathModel <- function(parTable = NULL, data.cov = NULL) {
  if (!NROW(parTable))
    return(methods::new("GlsPathModel"))

  # Measurement
  inds <- union(
    parTable[parTable$op == "<~", "rhs"],
    parTable[parTable$op == "=~", "rhs"]
  )

  # Structural (possibly with a few indicators)
  vars <- unique(c(
    parTable[parTable$op %in% c("=~", "<~"), "lhs"],
    parTable[parTable$op %in% c("~",  "~~"), "lhs"],
    parTable[parTable$op %in% c("~",  "~~"), "rhs"]
  ))

  # Currently we treat any variable with a "~~" as a structural variable
  # This makes sense as we don't allow the user to specify the covariance
  # structure of the measurement model. If this ever changes we will have to
  # to things differently, particularly for higher order models.
  if (PLS_IGNORE_INDCOV) etas <- setdiff(vars, inds) # all structural variables
  else                   etas <- vars

  reg  <- parTable[parTable$op == "~", , drop = FALSE]
  xis  <- setdiff(etas, reg$lhs) # purely exogenous variables
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
    lhs <- cov[i, "lhs"]
    rhs <- cov[i, "rhs"]
    psi[lhs, rhs] <- psi[rhs, lhs] <- NA
  }

  diag(psi) <- NA
  psi[xis, xis] <- NA
  psi[upper.tri(psi)] <- 0

  if (!is.null(data.cov)) {
    pair <- glsGetSampleInformation(
      data.cov = data.cov, etas = etas
    )

    S     <- pair$S
    S.inv <- pair$S.inv

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


glsGetSampleInformation <- function(data.cov, etas) {
  pls_stopif(!is.matrix(data.cov),
    "`data.cov` must be a matrix!"
  )

  rm <- setdiff(etas, rownames(data.cov))
  pls_stopif(length(rm),
    "Missing rownames in `data.cov`:", paste0(rm, collapse = ", ")
  )

  cm <- setdiff(etas, colnames(data.cov))
  pls_stopif(length(cm),
    "Missing rownames in `data.cov`:", paste0(cm, collapse = ", ")
  )

  S <- data.cov[etas, etas]

  tryCatch(
    S.inv <- MASS::ginv(S),
    error = function(e) {
      pls_msg_stop(
        "Could not invert sample covariance matrix!",
        "Message:", condtionMessage(e)
      )
    }
  )

  list(S = S, S.inv = S.inv)
}


`glsModelCovMatrix<-` <- function(model, value) {
  pair <- glsGetSampleInformation(
    data.cov = value, etas = model@info$etas
  )

  model@matrices$S <- pair$S
  model@matrices$S.inv <- pair$S.inv

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

  binv <- solve(I - gamma)
  binv %*% psi %*% t(binv)
}


glsObjective <- function(model) {
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


glsEstimateParameters <- function(model, data.cov = NULL,
                                  control = list(eval.max = 1500, iter.max = 1000),
                                  ...) {
  if (!is.null(data.cov)) # update sample covariance matrix?
    glsModelCovMatrix(model) <- data.cov

  # objective and gradient
  fn <- \(x) glsObjective(glsFillModel(model, x))
  gr <- \(x) glsGradient(glsFillModel(model, x))

  opt <- nlminb(
    start = model@info$start,
    objective = fn,
    gradient = gr,
    control = control,
    ...
  )

  glsFillModel(model, opt$par)
}
