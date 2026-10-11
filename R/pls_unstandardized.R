plsUnstandardizedEstimates <- function(model, unstandardized = "all",
                                       se = c("delta", "none"),
                                       scale.uncertainty = FALSE,
                                       eps = 1e-4, zero.tol = 1e-10,
                                       rm.tmp.ov = TRUE,
                                       clean.tmp.ind = TRUE,
                                       clean.tmp.mimic = TRUE,
                                       vcov = modelParams(combinedModel(model))$vcov,
                                       sds.override = NULL) {
  se <- tolower(se)
  se <- match.arg(se)

  combined <- combinedModel(model)
  fit      <- modelFit(combined)
  info     <- modelInfo(combined)

  lvs      <- info$lvs.linear
  inds.a   <- info$inds.a
  inds.b   <- info$inds.b
  ovs      <- union(inds.a, inds.b)
  ordered  <- intersect(info$ordered, ovs)
  cont     <- setdiff(ovs, ordered)
  scale    <- info$scale
  resvar   <- diag(fit$fitTheta)
  n        <- NROW(modelData(combined))

  if (length(unstandardized) == 1L) {

    uvars <- switch(
      tolower(unstandardized),
      all = c(lvs, ovs),
      ov  = ovs,
      lv  = lvs,
      unstandardized
    )

  } else {
    uvars <- unstandardized

  }


  # check for potential temp_ov/temp_ind
  # relevant for all but `unstandardized="all"`
  if (length(setdiff(ovs, uvars))) {
    tmp0 <- paste0(TEMP_OV_PREFIX, uvars)
    tmp1 <- paste0(uvars, TEMP_IND_SUFFIX)
    tmp2 <- paste0(tmp0, TEMP_IND_SUFFIX)
    tmp  <- c(tmp0, tmp1, tmp2)

    uvars <- union(uvars, intersect(tmp, ovs))
  }

  # validate
  missing <- setdiff(uvars, c(ovs, lvs))
  pls_stopif(length(missing),
    "Unknown variables in `unstandardized`! Variables:",
    paste0(missing, collapse = ", ")
  )

  parTable <- getParTableEstimates(
    combined,
    rm.tmp.ov = FALSE,
    clean.tmp.ind = FALSE,
    clean.tmp.mimic = FALSE
  )

  # get target sds for observed variables, and placeholders for lvs
  sds <- stats::setNames(
    rep(1, length(ovs) + length(lvs)), nm = c(ovs, lvs)
  )

  # for continuous (unstandardized) variables we can just use the observed sds
  for (x in intersect(cont, uvars)) {
    sds[[x]] <- tryCatch(scale[[x]], error = \(...) NA)
  }

  # for ordered variables we try to fix the residual variance to 1
  # (for indicators which have residual variances)
  for (x in intersect(ordered, uvars)) {
    if (x %in% inds.b) sds[[x]] <- 1 # best guess
    else sds[[x]] <- tryCatch(1 / sqrt(resvar[[x]]), error = \(...) NA)
  }

  # if we don't know, we assume they are standard normal
  sds[is.na(sds)] <- 1

  override <- intersect(names(sds.override), names(sds))
  sds[override] <- sds.override[override]

  parTableUstd <- unstandardizedEstimatesInternal(
    parTable      = parTable,
    unstandardize = uvars,
    lvs           = lvs,
    inds.a        = inds.a,
    inds.b        = inds.b,
    sds           = sds
  )

  sds.out <- attr(parTableUstd, "sds")

  # full vcov
  if (se == "delta" && NROW(vcov) && NCOL(vcov)) {
    # V[pars]
    pars <- paste0(parTable$lhs, parTable$op, parTable$rhs)
    Vpar <- expandVcov(vcov, labels = pars)

    # V[vars]
    # use approximated std.errors
    sdsu <- sds[uvars]
    sdss <- sds[setdiff(names(sds), uvars)]

    if (scale.uncertainty) {
      # add uncertainty estimates to the scale of the variables
      Vsds <- diag(sdsu^2 / (2 * (n - 1)), nrow = length(sdsu))
      dimnames(Vsds) <- list(names(sdsu), names(sdsu))

      # V[vars,pars]
      V <- diagPartitioned(Vpar, Vsds)

    } else {
      # else we just consider the uncertainty in the parameters
      V <- Vpar

    }

    # jacobian for delta method standard errors
    J <- matrix(
      0, nrow = length(pars), ncol = NCOL(V),
      dimnames = list(pars, colnames(V))
    )

    # if scale.uncertainty=FALSE i never reaches any of the values in
    # idx.var, and thus it stays constant.
    idx.par <- seq_along(pars)
    idx.var <- seq_along(sdsu) + length(pars)

    sds0  <- c(sdsu, sdss)
    par0  <- parTable
    x0    <- c(parTable$est, sdsu)
    y0    <- parTableUstd$est

    for (i in seq_len(NROW(V))) {
      # get templates
      pari   <- par0
      sdsi   <- sds0
      xi     <- x0
      xi[i]  <- x0[i] + eps

      # fill in templats
      pari$est                 <- xi[idx.par]
      sdsi[seq_along(idx.var)] <- xi[idx.var]

      # get estimated solution
      parTableUstd.i <- unstandardizedEstimatesInternal(
        parTable      = pari,
        unstandardize = uvars,
        lvs           = lvs,
        inds.a        = inds.a,
        inds.b        = inds.b,
        sds           = sdsi
      )

      yi <- parTableUstd.i$est

      J[,i] <- (yi - y0) / eps
    }

    Vstd <- J %*% V %*% t(J)

    se <- sqrt(diag(Vstd))
    se[se <= zero.tol] <- NA_real_

    parTableUstd$se <- se

  } else {
    parTableUstd$se <- NA_real_
    Vstd <- NULL
  }

  # add z-stats and format
  parTableOut <- addZStatsParTable(parTableUstd)

  if (rm.tmp.ov)
    parTableOut <- removeTempOvRowsParTable(parTableOut)

  if (clean.tmp.ind)
    parTableOut <- cleanTempIndRowsParTable(parTableOut)

  if (clean.tmp.mimic)
    parTableOut <- cleanTempMimicRowsParTable(parTableOut)

  parTableOut <- plssemParTable(parTableOut)

  # add vcov?
  if (!is.null(Vstd))
    attr(parTableOut, "vcov") <- plssemMatrix(Vstd, is.public = TRUE)

  attr(parTableOut, "sds") <- sds.out
  parTableOut
}


plsUnstandardizedEstimatesMultilevel <- function(model, unstandardized = "all",
                                                 se = c("delta", "none"),
                                                 icc.default = 0.5,
                                                 ...) {
  se <- match.arg(tolower(se), c("delta", "none"))
  vcov <- model@params$vcov

  parTable <- as.data.frame(model@parTable)
  scale1   <- modelInfo(combinedModel(model@level.1))$scale
  scale2   <- modelInfo(combinedModel(model@level.2))$scale

  # is unstandardized a set of variables?
  ustd.is.vars <- !(
    length(unstandardized) == 1L &&
    tolower(unstandardized) %in% c("all", "lv", "ov")
  )

  # (clean) names of the variables of a level
  getLevelVars <- function(level) {
    info <- modelInfo(combinedModel(level))
    unique(removeTempAffixes(c(info$lvs.linear, info$inds.a, info$inds.b)))
  }

  vars1 <- getLevelVars(model@level.1)
  vars2 <- getLevelVars(model@level.2)

  if (ustd.is.vars) {
    # A variable might only exist at one of the levels, so we validate the
    # selection against both levels, and pass each level its own variables
    missing <- setdiff(unstandardized, union(vars1, vars2))
    pls_stopif(length(missing),
      "Unknown variables in `unstandardized`! Variables:",
      paste0(missing, collapse = ", ")
    )

    unstandardized1 <- intersect(unstandardized, vars1)
    unstandardized2 <- intersect(unstandardized, vars2)

  } else {
    unstandardized1 <- unstandardized2 <- unstandardized
  }

  # ordered (observed) variables which are unstandardized
  if (!ustd.is.vars) {
    ordered.ustd <- switch(tolower(unstandardized),
      all = model@info$ordered,
      ov  = model@info$ordered,
      lv  = character(0L),
      pls_msg_stop("unrecognized value of unstandardized!")
    )

  } else {
    ordered.ustd <- intersect(unstandardized, model@info$ordered)
  }

  # (possibly temporary) names of a variable in the scales of a level
  getScaleNames <- function(v, scale)
    names(scale)[removeTempAffixes(names(scale)) == v]

  # we fix the total residual variances for ordinals to 1, where
  #   theta = (1 - icc) * theta_W + icc * theta_B
  getStdResVar <- function(v, l) {
    cond <- parTable$level %in% l & parTable$op == "~~" &
            parTable$lhs == v & parTable$rhs == v
    theta <- parTable$est[cond]
    if (length(theta) && is.finite(theta[[1L]])) theta[[1L]] else 1
  }

  # non-finite iccs?
  icc.levels <- model@params$icc
  icc.levels[!is.finite(icc.levels)] <- icc.default

  ordered   <- model@info$ordered
  sds.ord   <- stats::setNames(rep(1, length(ordered)), nm = ordered)
  override1 <- NULL
  override2 <- NULL

  for (v in ordered) {
    nm1 <- getScaleNames(v, scale1)
    nm2 <- getScaleNames(v, scale2)

    is.l1 <- length(nm1) > 0
    is.l2 <- length(nm2) > 0

    icc   <- if (!is.l2) 0 else if (!is.l1) 1 else icc.levels[[v]]
    theta <- (1 - icc) * getStdResVar(v, l = 1L) + icc * getStdResVar(v, l = 2L)

    sds.ord[[v]]   <- 1 / sqrt(max(theta, .Machine$double.eps))
    override1[nm1] <- sqrt(1 - icc) * sds.ord[[v]]
    override2[nm2] <- sqrt(icc) * sds.ord[[v]]
  }

  # standardize L1 first (L2 depends on L1 if there are random slopes)
  pt1 <- plsUnstandardizedEstimates(
    model          = model@level.1,
    unstandardized = unstandardized1,
    se             = se,
    vcov           = subsetVcovByLevel(vcov, l = 1L),
    sds.override   = override1,
    ...
  )

  # scales of the random slopes
  sds1    <- attr(pt1, "sds")
  rslopes <- model@info$rslopes
  rsd     <- model@params$rsd

  getStdDevL1 <- function(v) {
    sd <- tryCatch(sds1[[v]], error = \(e) NA_real_)
    if (is.null(sd) || !is.finite(sd)) 1 else sd
  }

  for (k in seq_len(NROW(rslopes))) {
    rs <- rslopes$name[[k]]
    nm <- getScaleNames(rs, scale2)

    if (!length(nm) || !rs %in% names(rsd))
      next

    sdl <- getStdDevL1(rslopes$lhs[[k]])
    sdr <- getStdDevL1(rslopes$rhs[[k]])

    sd.s <- rsd[[rs]] * sdl / sdr
    override2[nm] <- sd.s
  }

  pt2 <- plsUnstandardizedEstimates(
    model          = model@level.2,
    unstandardized = unstandardized2,
    se             = se,
    vcov           = subsetVcovByLevel(vcov, l = 2L),
    sds.override   = override2,
    ...
  )

  pt1 <- as.data.frame(pt1)
  pt2 <- as.data.frame(pt2)
  pt1$level <- rep(1L, NROW(pt1))
  pt2$level <- rep(2L, NROW(pt2))

  # The thresholds are defined for the total latent responses
  thr  <- parTable[is.na(parTable$level), , drop = FALSE]
  cols <- union(colnames(pt1), colnames(pt2))
  is.thr <- thr$op == "|" & thr$lhs %in% intersect(names(sds.ord), ordered.ustd)

  if (se == "none" && NROW(thr)) {
    # clear the (standardized) statistics, as for the levels
    cols.stat <- intersect(
      c("se", "z", "pvalue", "ci.lower", "ci.upper"),
      colnames(thr)
    )

    thr[, cols.stat] <- NA_real_
  }

  if (any(is.thr)) {
    target <- sds.ord[thr$lhs[is.thr]]
    cols.thr <- intersect(
      c("est", "se", "ci.lower", "ci.upper"),
      colnames(thr)
    )

    for (col in cols.thr)
      thr[is.thr, col] <- thr[is.thr, col] * target
  }

  for (col in setdiff(cols, colnames(thr)))
    thr[[col]] <- NA

  out <- rbind(pt1[cols], pt2[cols], thr[cols])
  rownames(out) <- NULL

  out <- plssemParTable(out)
  attr(out, "sds") <- list(
    level.1 = sds1,
    level.2 = attr(pt2, "sds"),
    total = sds.ord
  )

  out
}


unstandardizedEstimatesInternal <- function(parTable,
                                            unstandardize,
                                            lvs, inds.a, inds.b,
                                            sds) {
  # keep original parTable for restoring the higher order measurement model
  uinds.a  <- intersect(inds.a, unstandardize)
  uinds.b  <- intersect(inds.b, unstandardize)
  ulvs     <- intersect(lvs,    unstandardize)

  # Indicators Mode A and Mode B
  for (ind in c(uinds.a, uinds.b)) {
    idxl  <- which(parTable$op == "=~" & parTable$rhs == ind)
    idxc  <- which(parTable$op == "<~" & parTable$rhs == ind)
    idxvl <- which(parTable$op == "~~" & parTable$lhs == ind)
    idxvr <- which(parTable$op == "~~" & parTable$rhs == ind)
    idxt  <- which(parTable$op == "|"  & parTable$lhs == ind)

    target <- sds[[ind]]

    parTable[idxl,  "est"] <- parTable[idxl,  "est"] * target # loadings
    parTable[idxc,  "est"] <- parTable[idxc,  "est"] / target # weights
    parTable[idxvl, "est"] <- parTable[idxvl, "est"] * target
    parTable[idxvr, "est"] <- parTable[idxvr, "est"] * target
    parTable[idxt,  "est"] <- parTable[idxt,  "est"] * target
  }

  # Latent Variables (Mode A) and Composites (Mode B)
  for (lv in ulvs) {

    # Measurement model
    idxma <- which(parTable$op == "=~" & parTable$lhs == lv)
    idxmb <- which(parTable$op == "<~" & parTable$lhs == lv)

    pls_stopif(length(idxma) && length(idxmb),
      "Did not expect latent variable to be both of mode A and B!",
      "Variable:", paste0("`", lv, "`")
    )

    if (length(idxma)) { # MODE A
      # fix first loading to 1
      loadings <- parTable[idxma, "est", drop = TRUE]
      parTable[idxma, "est"] <- loadings / loadings[[1L]]

      first.ind <- parTable[idxma[[1L]], "rhs", drop = TRUE]
      first.res <- parTable[
        parTable$lhs == first.ind &
        parTable$op  == "~~"      &
        parTable$rhs == first.ind, "est", drop = TRUE
      ][[1L]]
      first.var <- sds[[first.ind]]^2

      target <- sqrt(first.var - first.res)

    } else if (length(idxmb)) { # MODE B

      # fix first loading to 1
      weights <- parTable[idxmb, "est", drop = TRUE]

      uweights <- weights / weights[[1L]]
      names(uweights) <- parTable[idxmb, "rhs", drop = TRUE]

      parTable[idxmb, "est"] <- unname(uweights)

      k <- length(uweights)
      S <- matrix(
        0, nrow = k, ncol = k,
        dimnames = list(names(uweights), names(uweights))
      )

      for (i in seq_len(k)) {
        lhs <- names(uweights)[[i]]

        for (j in seq_len(i)) {
          rhs <- names(uweights)[[j]]

          p <- parTable[
            (parTable$lhs == lhs & parTable$op == "~~" & parTable$rhs == rhs) |
            (parTable$lhs == rhs & parTable$op == "~~" & parTable$rhs == lhs),
            "est", drop = TRUE
          ]

          S[i, j] <- S[j, i] <- p
        }
      }

      target <- c(sqrt(t(uweights) %*% S %*% uweights))
    }

    sds[[lv]] <- target

    # (co-)variances
    idxvl <- which(parTable$op == "~~" & parTable$lhs == lv)
    idxvr <- which(parTable$op == "~~" & parTable$rhs == lv)
    parTable[idxvl, "est"] <- parTable[idxvl, "est"] * target
    parTable[idxvr, "est"] <- parTable[idxvr, "est"] * target

    # outgoing paths
    idxp0 <- which(parTable$rhs == lv & parTable$op == "~")
    parTable[idxp0, "est"] <- parTable[idxp0, "est"] / target

    # incoming paths
    idxp1 <- which(parTable$lhs == lv & parTable$op == "~")
    parTable[idxp1, "est"] <- parTable[idxp1, "est"] * target

  }

  checkLhsIntTerms(parTable)
  intTerms <- getIntTerms(parTable)

  for (intTerm in intTerms) {
    # first we rescale the interaction term
    # note that we've already rescaled the lhs side, we just need to
    # update the rhs side, we have
    # std(b3) =     b3  * sd(x) *  sd(z) / sd(y)
    #     b3  = std(b3) * sd(y) / (sd(x) * sd(z))
    # at this point we have b3 = std(b3) * sd(y), and thus divide
    # by sd(x) * sd(z).

    elems <- stringr::str_split_1(intTerm, pattern = ":")
    idxp  <- which(parTable$rhs == intTerm & parTable$op == "~")
    idxvl <- which(parTable$lhs == intTerm & parTable$op == "~~")
    idxvr <- which(parTable$rhs == intTerm & parTable$op == "~~")
    target <- prod(sds[elems])

    parTable[idxp,  "est"] <- parTable[idxp,  "est"] / target
    parTable[idxvl, "est"] <- parTable[idxvl, "est"] * target
    parTable[idxvr, "est"] <- parTable[idxvr, "est"] * target

  }

  # Custom parameters?
  custom <- parTable$op == ":="

  if (any(custom)) {
    idx <- which(custom)
    has.lab <- parTable$label != ""

    # create environment
    envir <- customExpressionEnv(stats::setNames(
      parTable[has.lab, "est"], nm = parTable[has.lab, "label"]
    ))

    for (i in idx) {
      lab <- parTable[i, "label"]
      val <- eval(parse(text=parTable[i, "rhs"]), envir = envir)

      # Remember to update both the parameter table, and environment
      assign(unname(lab), val, envir = envir)
      parTable[i, "est"] <- val
    }
  }

  attr(parTable, "sds") <- sds
  parTable
}
