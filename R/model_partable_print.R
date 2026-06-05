plsPrintParTable <- function(parTable,
                             scientific  = FALSE,
                             ci          = FALSE,
                             digits      = 3,
                             loadings    = TRUE,
                             regressions = TRUE,
                             covariances = TRUE,
                             intercepts  = TRUE,
                             variances   = TRUE,
                             thresholds  = TRUE,
                             ...) {

  if (!"label" %in% colnames(parTable))
    parTable$label <- ""

  parTable <- rename(
    .X      = parTable,
    est.std = "est",
    se      = "std.error",
    pvalue  = "p.value",
    p       = "p.value",
    zvalue  = "z.value",
    z       = "z.value"
  )

  modPrintParTable(
    parTable    = parTable,
    scientific  = scientific,
    ci          = ci,
    digits      = digits,
    loadings    = loadings,
    regressions = regressions,
    covariances = covariances,
    intercepts  = intercepts,
    variances   = variances,
    thresholds  = thresholds,
    ...
  )
}


plsGetWidthPrintedParTable <- function(parTable,
                                       scientific  = FALSE,
                                       ci          = FALSE,
                                       digits      = 3,
                                       loadings    = TRUE,
                                       regressions = TRUE,
                                       covariances = TRUE,
                                       intercepts  = TRUE,
                                       variances   = TRUE,
                                       thresholds  = TRUE,
                                       ...) {

  if (!"label" %in% colnames(parTable))
    parTable$label <- ""

  modGetWidthPrintedParTable(
    parTable    = parTable,
    scientific  = scientific,
    ci          = ci,
    digits      = digits,
    loadings    = loadings,
    regressions = regressions,
    covariances = covariances,
    intercepts  = intercepts,
    variances   = variances,
    thresholds  = thresholds
  )
}
