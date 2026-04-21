#' Summarize a fitted \code{plssem} model
#'
#' @param object An object of class \code{plssem}.
#' @param fit Should fit measures be calculated?
#' @param ... Additional arguments passed to or from methods.
#' @return A \code{SummaryPlsSem} object containing formatted parameter estimates.
#'
#' @export
summary.plssem <- function(object, fit = TRUE, ...) {
  parTable <- parameter_estimates(object)

  lvs <- getLVs(parTable)
  ovs <- getOVs(parTable)
  etas <- getEtas(parTable, checkAny=FALSE)
  inds <- getIndicators(parTable)
  inds.a <- getReflectiveIndicators(parTable)

  strParTableLines <- utils::capture.output(modsem::summarize_partable(parTable))
  strParTable <- paste0(paste0(strParTableLines[-(1:6)], collapse = "\n"), "\n") # [-(1:6)] to skip headers

  is.ord <- object$info$is.probit || (length(object$info$ordered) && object$info$is.mcpls)
  if (is.ord) link <- "PROBIT"
  else        link <- "LINEAR"

  getR2 <- function(x, pt = parTable) {
    rvar <- pt[pt$lhs == x & pt$op == "~~" & pt$rhs == x, "est"]
    if (!length(rvar)) 0 else 1 - rvar
  }

  r2.etas <- vapply(etas,   FUN.VALUE = numeric(1L), FUN = getR2)
  r2.inds <- vapply(inds.a, FUN.VALUE = numeric(1L), FUN = getR2)

  out <- list(
    print = list(
      strParTable = strParTable,
      width = max(nchar(strParTableLines))
    ),
    fit   = object,
    info  = list(
      iterations = object$status$iterations,
      estimator  = object$info$estimator,
      n          = object$info$n,
      nlvs       = length(lvs),
      novs       = length(ovs),
      link       = link,
      etas       = etas,
      inds       = inds
    ),
    r2 = list(
      etas = r2.etas,
      inds = r2.inds
    ),
    fit.measures = if (fit) fitMeasures(object) else NULL
  )

  class(out) <- "SummaryPlsSem"
  out
}


#' Print a \code{SummaryPlsSem} object
#'
#' @param x A \code{SummaryPlsSem} object as returned by [summary.plssem()].
#' @param ... Additional arguments for compatibility with the generic.
#' @return The input object, invisibly.
#'
#' @export
print.SummaryPlsSem <- function(x, ...) {
  formatValue <- function(val, digits = 3) {
    if (length(val) <= 0 || all(is.na(val))) return("NA")
    formatNumeric(val, digits = digits, scientific = FALSE)
  }

  printf("plssem (%s) ended normally after %i iterations\n\n",
         PKG_INFO$version, x$info$iterations)

  width.out <- x$print$width

  headerNames <- c(
    "Estimator",
    "Link",
    "",
    "Number of observations",
    "Number of iterations",
    "Number of latent variables",
    "Number of observed variables"
  )

  headerValues <- c(
    x$info$estimator,
    x$info$link,
    "",
    x$info$n,
    x$info$iterations,
    x$info$nlvs,
    x$info$novs
  )

  cat(allignLhsRhs(lhs = headerNames, rhs = headerValues, pad = "  ",
                   width.out = width.out), "\n", sep = "")

  if (!is.null(x$fit.measures)) {
    cat("Fit Measures:\n")
    headerNames <- c(
      "Chi-Square",
      "Degrees of Freedom",
      "SRMR",
      "RMSEA"
    )

    headerValues <- c(
      sprintf("%.3f", x$fit.measures$chisq),
      sprintf("%d", x$fit.measures$chisq.df),
      sprintf("%.3f", x$fit.measures$srmr),
      sprintf("%.3f", x$fit.measures$rmsea)
    )

    cat(allignLhsRhs(lhs = headerNames, rhs = headerValues, pad = "  ",
                     width.out = width.out), "\n", sep = "")
  }

  if (length(x$r2$inds)) {
    cat("R-squared (indicators):\n")

    headerNames <- names(x$r2$inds)
    headerValues <- formatNumeric(x$r2$inds)

    cat(allignLhsRhs(lhs = headerNames, rhs = headerValues, pad = "  ",
                     width.out = width.out), "\n", sep = "")
  }

  if (length(x$r2$etas)) {
    cat("R-squared (latents):\n")

    headerNames <- names(x$r2$etas)
    headerValues <- formatNumeric(x$r2$etas)

    cat(allignLhsRhs(lhs = headerNames, rhs = headerValues, pad = "  ",
                     width.out = width.out), "\n", sep = "")
  }

  cat(x$print$strParTable)
  invisible(x)
}


#' Print a \code{plssem} object
#'
#' @param x An object of class \code{plssem}.
#' @param ... Additional arguments for compatibility with the generic.
#' @return The input object, invisibly.
#'
#' @export
print.plssem <- function(x, ...) {
  printf("plssem (%s) ended normally after %i iterations\n",
         PKG_INFO$version, x$status$iterations)
  print(parameter_estimates(x))
  invisible(x)
}


#' @export
#' @importFrom stats coefficients
coefficients.plssem <- function(object, ...) {
  object$params$values
}


#' @export
#' @importFrom stats coef
coef.plssem <- function(object, ...) {
  coefficients(object, ...)
}


#' @export
#' @importFrom stats vcov
vcov.plssem <- function(object, ...) {
  object$params$vcov
}


#' Parameter estimates for \code{plssem} objects
#'
#' @param object An object of class \code{plssem}.
#' @param colon.pi Logical; whether to replace labels for interaction terms with colon notation.
#' @param label.renamed.prod Logical; whether renamed product labels should be retained when colon expansion occurs.
#' @param ... Additional arguments (not used).
#' @return A parameter table (data frame) describing the fitted model.
#'
#' @export
parameter_estimates.plssem <- function(object,
                                       colon.pi = TRUE,
                                       label.renamed.prod = FALSE,
                                       ...) {
  parTable <- object$parTable

  if (colon.pi)
    parTable <- addColonPI_ParTable(parTable, model = object,
                                    label.renamed.prod = label.renamed.prod)

  parTable
}


#' Generic accessor for model parameter estimates
#'
#' @param object A fitted model object.
#' @param ... Additional arguments passed to methods.
#' @return A parameter table describing the fitted model.
#'
#' @export
parameter_estimates <- function(object, ...) {
  UseMethod("parameter_estimates")
}


#' @export
isMC_PLS.plssem <- function(object) {
  object$info$is.mcpls 
}


#' Check if object is a MC-PLS model
#'
#' @param object A fitted model object.
#' @return \code{TRUE}/\code{FALSE}.
#'
#' @export
isMC_PLS <- function(object) {
  UseMethod("isMC_PLS")
}
