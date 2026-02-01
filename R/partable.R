CI_QUANTILE <- qnorm(0.05)


getParTableEstimates <- function(model) {
  params <- model$params
  est    <- model$params$values
  se     <- model$params$se
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

  plssemParTable(removeTempOV_RowsParTable(parTable))
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
  tmp <- startsWith(parTable$lhs, TEMP_OV_PREFIX) | startsWith(parTable$rhs, TEMP_OV_PREFIX)
  parTable[!tmp, , drop = FALSE]
}
