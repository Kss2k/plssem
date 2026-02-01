#' Summarize a fitted `plssem` model
#'
#' @param object An object of class `plssem`.
#' @param ... Additional arguments passed to or from methods.
#' @return A `SummaryPlsSem` object containing formatted parameter estimates.
#'
#' @export
summary.plssem <- function(object, ...) {
  parTable    <- parameter_estimates(object)
  strParTable <- utils::capture.output(modsem::summarize_partable(parTable))
  strParTable <- paste0(paste0(strParTable[-1], collapse = "\n"), "\n") # [-1] to skip header

  out <- list(
    print = list(strParTable = strParTable),
    fit   = object,
    info  = list(iterations = object$info$iterations)
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
  printf("plssem (%s) ended normally after %i iterations\n",
         PKG_INFO$version, x$info$iterations)
  cat(x$print$strParTable)
}


#' Print a `plssem` object
#'
#' @param object An object of class `plssem`.
#' @param ... Additional arguments for compatibility with the generic.
#'
#' @export
print.plssem <- function(object, ...) {
  print(parameter_estimates(object))
  invisible(object)
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
  
  if (colon.pi) {
    elems <- object$info$intTermElems

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
  }

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
