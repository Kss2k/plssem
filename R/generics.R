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


#' @export
print.SummaryPlsSem <- function(x, ...) {
  printf("plssem (%s) ended normally after %i iterations\n",
         PKG_INFO$version, x$info$iterations)
  cat(x$print$strParTable)
}


#' @export
print.plssem <- function(object, ...) {
  print(parameter_estimates(object))
  invisible(object)
}


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


#' @export
parameter_estimates <- function(object, ...) {
  UseMethod("parameter_estimates")
}
