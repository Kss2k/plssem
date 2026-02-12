#' Summarize a fitted `plssem` model
#'
#' @param object An object of class `plssem`.
#' @param ... Additional arguments passed to or from methods.
#' @return A `SummaryPlsSem` object containing formatted parameter estimates.
#'
#' @export
summary.plssem <- function(object, ...) {
  parTable    <- parameter_estimates(object)

  lvs <- getLVs(parTable)
  ovs <- getOVs(parTable)
  etas <- getEtas(parTable)
  inds <- getIndicators(parTable)

  strParTableLines <- utils::capture.output(modsem::summarize_partable(parTable))
  strParTable <- paste0(paste0(strParTableLines[-(1:6)], collapse = "\n"), "\n") # [-(1:6)] to skip headers

  if      (object$info$is.probit) link <- "PROBIT"
  else if (object$info$is.cexp)   link <- "PROBIT-CEXP"
  else                            link <- "LINEAR"

  getR2 <- function(x, pt = parTable) {
    rvar <- pt[pt$lhs == x & pt$op == "~~" & pt$rhs == x, "est"]
    if (!length(rvar)) 0 else 1 - rvar
  }

  r2.etas <- vapply(etas, FUN.VALUE = numeric(1L), FUN = getR2)
  r2.inds <- vapply(inds, FUN.VALUE = numeric(1L), FUN = getR2)

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
    )
  )



  class(out) <- "SummaryPlsSem"
  out
}


#' Print a `SummaryPlsSem` object
#'
#' @param x A `SummaryPlsSem` object as returned by [summary.plssem()].
#' @param ... Additional arguments for compatibility with the generic.
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

  cat("R-squared (indicators):\n")

  headerNames <- names(x$r2$inds)
  headerValues <- formatNumeric(x$r2$inds)

  cat(allignLhsRhs(lhs = headerNames, rhs = headerValues, pad = "  ",
                   width.out = width.out), "\n", sep = "")

  cat("R-squared (latents):\n")

  headerNames <- names(x$r2$etas)
  headerValues <- formatNumeric(x$r2$etas)

  cat(allignLhsRhs(lhs = headerNames, rhs = headerValues, pad = "  ",
                   width.out = width.out), "\n", sep = "")

  cat(x$print$strParTable)
}


#' Print a `plssem` object
#'
#' @param x An object of class `plssem`.
#' @param ... Additional arguments for compatibility with the generic.
#'
#' @export
print.plssem <- function(x, ...) {
  printf("plssem (%s) ended normally after %i iterations\n",
         PKG_INFO$version, x$status$iterations)
  print(parameter_estimates(x))
  invisible(x)
}


#' Parameter estimates for `plssem` objects
#'
#' @param object An object of class `plssem`.
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
#'
#' @export
parameter_estimates <- function(object, ...) {
  UseMethod("parameter_estimates")
}
