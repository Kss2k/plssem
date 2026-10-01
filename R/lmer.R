MOD_OP_ALIAS <- "__MOD__"


lmerEstimateParameters <- function(parTable, data, cluster, control = lmerFastControl(),
                                   vtol = 1e-8, ...) {
  pls_stopif(length(cluster) != 1 || !is.character(cluster),
    "cluster must be a character string of length 1!"
  )

  # The variable names must match the (aliased) names in the formulas
  clusterData <- attr(data, "cluster")
  data <- as.data.frame(data)
  colnames(data) <- stringr::str_replace_all(
    colnames(data), pattern = ":", replacement = MOD_OP_ALIAS
  )

  if (!is.null(clusterData))
    data[cluster] <- clusterData[cluster]

  etas <- getEtas(parTable, checkAny = FALSE)

  if (!length(etas)) { # no structural model (e.g., a CFA model)
    return(list(
      pars   = data.frame(lhs = character(0L), op = character(0L), rhs = character(0L),
                          est = numeric(0L), mod = character(0L)),
      randef = NULL,
      resvar = numeric(0L)
    ))
  }

  parTable$lhs <- stringr::str_replace_all(
    parTable$lhs, pattern = ":", replacement = MOD_OP_ALIAS
  )

  parTable$rhs <- stringr::str_replace_all(
    parTable$rhs, pattern = ":", replacement = MOD_OP_ALIAS
  )

  randef <- NULL
  pars   <- NULL
  resvar <- stats::setNames(numeric(length(etas)), nm = etas)

  for (eta in etas) {
    rows <- parTable[parTable$lhs == eta & parTable$op == "~",,drop=FALSE]
    has.randef <- isRandomEffectMod(rows$mod)

    randef.map <- stats::setNames(
      extractRandomEffectName(rows[has.randef, "mod"]),
      nm = rows[has.randef, "rhs"]
    )

    if (any(has.randef)) {
      formula <- paste0(
        eta, "~", paste0(rows$rhs, collapse = " + "), " + (",
        paste0(rows$rhs[has.randef], collapse = " + "), " - 1 | ",
        paste0(cluster, collapse = "/"), ")"
      )

      fit <- lme4::lmer(formula, data = data, control = control, ...)
      coef <- lme4::fixef(fit)
      resvar[[eta]] <- stats::sigma(fit)^2

      # cluster specific slopes (i.e., fixed + random effects)
      randef.eta <- as.matrix(stats::coef(fit)[[cluster]][, names(randef.map), drop = FALSE])
      colnames(randef.eta) <- randef.map[colnames(randef.eta)]

      for (i in seq_len(NCOL(randef.eta))) {
        # we need at least some variablity for the variable passed to level 2
        var.i <- stats::var(randef.eta[,i,drop=TRUE])

        if (is.na(var.i) || var.i <= vtol) {
          nm <- colnames(randef.eta)
          fixed <- tryCatch(coef[[nm]], error = \(e) 0)
          randef.eta[,i] <- stats::rnorm(NROW(randef.eta), mean = fixed, sd = sqrt(vtol*10))
        }
      }

      if (is.null(randef)) randef <- randef.eta
      else randef <- cbind(randef, randef.eta)

    } else {
      formula <- paste0(
        eta, "~", paste0(rows$rhs, collapse = " + ")
      )

      fit <- stats::lm(formula, data = data, ...)
      coef <- stats::coef(fit)
      resvar[[eta]] <- sum(stats::residuals(fit)^2) / (NROW(data) - 1)
    }

    # We keep a fixed intercept, since the interaction terms aren't centered
    coef <- coef[names(coef) != "(Intercept)"]

    pars <- rbind(pars,
      data.frame(lhs = eta, op = "~", rhs = names(coef), est = unname(coef))
    )
  }

  pars <- merge(
    x = pars,
    y = parTable[c("lhs", "op", "rhs", "mod")],
    by = c("lhs", "op", "rhs"),
    all.x = TRUE, all.y = FALSE
  )

  pars$lhs <- stringr::str_replace_all(
    pars$lhs, pattern = MOD_OP_ALIAS, replacement = ":" 
  )

  pars$rhs <- stringr::str_replace_all(
    pars$rhs, pattern = MOD_OP_ALIAS, replacement = ":" 
  )
 
  list(pars = pars, randef = randef, resvar = resvar)
}


# The post-estimation derivatives (gradient and Hessian) are only used for
# convergence checks, and are the main overhead besides the optimization itself.
# The remaining checks are cheap, but noisy when called repeatedly
# (e.g., in bootstrap and MC-PLS iterations).
lmerFastControl <- function() {
  lme4::lmerControl(
    calc.derivs         = FALSE,
    check.nobs.vs.rankZ = "ignore",
    check.nobs.vs.nlev  = "ignore",
    check.nlev.gtreq.5  = "ignore",
    check.nlev.gtr.1    = "ignore",
    check.nobs.vs.nRE   = "ignore",
    check.rankX         = "ignore",
    check.scaleX        = "ignore",
    check.formula.LHS   = "ignore",
    check.conv.grad     = "ignore",
    check.conv.singular = "ignore",
    check.conv.hess     = "ignore"
  )
}
