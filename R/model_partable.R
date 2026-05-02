CI_QUANTILE <- qnorm(0.05)


getParTableEstimates <- function(model, rm.tmp.ov = TRUE, clean.tmp.ind = TRUE) {
  est    <- model@params$values
  se     <- model@params$se
  names  <- names(est)

  split    <- splitParameterNames(names)
  lhs      <- split$lhs
  op       <- split$op
  rhs      <- split$rhs
  z        <- est / se
  pvalue   <- 2 * stats::pnorm(-abs(z))
  ci.lower <- est - CI_QUANTILE * z
  ci.upper <- est + CI_QUANTILE * z

  parTable <- data.frame(
    lhs      = lhs,
    op       = op,
    rhs      = rhs,
    est      = est,
    se       = se,
    z        = z,
    pvalue   = pvalue,
    ci.lower = ci.lower,
    ci.upper = ci.upper
  )

  if (rm.tmp.ov)
    parTable <- removeTempOV_RowsParTable(parTable)

  if (clean.tmp.ind)
    parTable <- cleanTempInd_RowsParTable(parTable)

  plssemParTable(parTable)
}


parTableToParams <- function(parTable) {
  lhs <- parTable$lhs
  op  <- parTable$op
  rhs <- parTable$rhs
  est <- parTable$est
  se  <- parTable$se

  k          <- NROW(parTable)
  names      <- paste0(lhs, op, rhs)
  values     <- stats::setNames(est, nm = names)

  list(
    names      = names,
    values     = values,
    values.old = NULL,
    se         = se,
    vcov       = NULL
  )
}


splitParameterNames <- function(names) {
  hasBeenSplit <- logical(length(names))

  lhs <- rep(NA_character_, length(names))
  op  <- rep(NA_character_, length(names))
  rhs <- rep(NA_character_, length(names))

  for (OP in OPERATORS) { # go by precedence
    split <- stringr::str_split_fixed(names, pattern = stringr::coll(OP), n = 2L)
    success <- stringr::str_detect(names, pattern = stringr::coll(OP))

    replace <- !hasBeenSplit & success
    hasBeenSplit <- hasBeenSplit | success

    lhs[replace] <- split[replace, 1L]
    rhs[replace] <- split[replace, 2L]
    op[replace]  <- OP
  }

  list(lhs = lhs, op = op, rhs = rhs)
}


removeTempOV_RowsParTable <- function(parTable) {
  tmp <- hasTempOvPrefix(parTable$lhs) | hasTempOvPrefix(parTable$rhs)
  parTable[!tmp, , drop = FALSE]
}


cleanTempInd_RowsParTable <- function(parTable) {
  rhs <- unique(parTable$rhs) # Only injected into the rhs column
  tmp <- rhs[hasTempIndSuffix(rhs)]
  cln <- removeTempIndSuffix(tmp)

  # We should remove any (co-)variances which are non-residuals, as it by
  # definition is an endogenous variable in the model
  parTable <- parTable[
    !((parTable$lhs %in% cln | parTable$rhs %in% cln) & parTable$op == "~~"),
    , drop = FALSE
  ]

  parTable$rhs <- removeTempIndSuffix(parTable$rhs)
  parTable$lhs <- removeTempIndSuffix(parTable$lhs)

  parTable
}


addColonPI_ParTable <- function(parTable, model, label.renamed.prod = FALSE) {
  elems <- model@info$intTermElems

  if (length(elems) && !"label" %in% colnames(parTable))
    parTable$label <- ""

  if (label.renamed.prod)
    origLabels <- getParTableLabels(parTable, labelCol = "label")
  else
    origLabels <- parTable$label

  for (xz in names(elems)) {
    xzColon <- paste0(elems[[xz]], collapse = ":")
    rmatch <- parTable$rhs == xz
    lmatch <- parTable$lhs == xz

    parTable[rmatch | lmatch, "label"] <- origLabels[rmatch | lmatch]

    parTable[rmatch, "rhs"] <- xzColon
    parTable[lmatch, "lhs"] <- xzColon # shouldn't be necessary, but just in case
                                       # the user has done something weird...
  }

  parTable
}
