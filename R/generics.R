#' @export
summary.plssem <- function(object, ...) {
  strParTable <- utils::capture.output(modsem::summarize_partable(object$parTable))
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
  print(object$parTable)
  invisible(object)
}
