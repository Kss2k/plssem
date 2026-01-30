CI_QUANTILE <- qnorm(0.05)


getParTableEstimates <- function(model) {
  params <- model$params
  names  <- model$params$names
  est    <- model$params$values
  se     <- model$params$se

  split    <- stringr::str_split_fixed(names, pattern = OP_REGEX, n = 2L)
  lhs      <- split[, 1L]
  rhs      <- split[, 2L]
  op       <- stringr::str_extract(names, pattern = OP_REGEX)
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

  plssemParTable(parTable)
}
