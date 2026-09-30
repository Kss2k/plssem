# S4 methods for the PlsMultilevelModel class (see `mpls()`).


printMultilevelStatusHeader <- function(model) {
  printStatusHeader(isTRUE(model@status$is.admissible), iterations = model@status$iterations)
}


#' Show a \code{PlsMultilevelModel} object
#'
#' @param object A \code{PlsMultilevelModel} object.
#' @return \code{object}, invisibly.
#' @export
setMethod("show", "PlsMultilevelModel", function(object) {
  printMultilevelStatusHeader(object)
  print(parameter_estimates(object))
  invisible(object)
})


#' Summarize a fitted \code{PlsMultilevelModel} model
#'
#' @param object A \code{PlsMultilevelModel} object.
#' @param ... Currently unused.
#' @return A \code{SummaryPlsMultilevel} list with formatted results.
#' @export
setMethod("summary", "PlsMultilevelModel", function(object, ...) {
  parTable <- parameter_estimates(object)
  info     <- object@info

  levels <- lapply(1:2, FUN = \(l) {
    pt <- parTable[parTable$level %in% l, , drop = FALSE]

    list(
      parTable = pt,
      r2.etas  = getR2ParTable(getEtas(pt, checkAny = FALSE), parTable = pt),
      r2.inds  = getR2ParTable(getReflectiveIndicators(pt), parTable = pt)
    )
  })

  rsd <- object@params$rsd
  rslopes <- info$rslopes

  if (length(rsd)) {
    paths <- paste0(rslopes$lhs, " ~ ", rslopes$rhs)
    names(rsd) <- paste0(names(rsd), " (", paths[match(names(rsd), rslopes$name)], ")")
  }

  out <- list(
    fit       = object,
    levels    = levels,
    agnostic  = parTable[is.na(parTable$level), , drop = FALSE], # e.g., thresholds
    print  = list(width = plsGetWidthPrintedParTable(parTable)),
    info   = list(
      estimator  = info$estimator,
      link       = if (length(info$ordered)) "PROBIT" else "LINEAR",
      n          = info$n,
      nclusters  = info$nclusters,
      iterations = object@status$iterations
    ),
    icc        = object@params$icc,
    rsd        = rsd
  )

  class(out) <- "SummaryPlsMultilevel"
  out
})


#' Print a \code{SummaryPlsMultilevel} object
#'
#' @param x A \code{SummaryPlsMultilevel} object as returned by
#'   \code{\link[=summary,PlsMultilevelModel-method]{summary}()}.
#' @param ... Additional arguments for compatibility with the generic.
#' @return The input object, invisibly.
#' @export
print.SummaryPlsMultilevel <- function(x, ...) {
  width.out <- x$print$width

  printMultilevelStatusHeader(x$fit)

  printSummarySection(width.out = width.out, values = stats::setNames(
    c(x$info$estimator, x$info$link, "", x$info$n, x$info$nclusters, x$info$iterations),
    nm = c("Estimator", "Link", "", "Number of observations", "Number of clusters",
           "Number of iterations")
  ))

  printSummarySection(x$icc, title = "Intraclass correlations:", width.out = width.out)
  printSummarySection(x$rsd, title = "Random slopes (standard deviations):", width.out = width.out)

  titles <- c("Level 1 [within]:", "Level 2 [between]:")

  for (l in 1:2) {
    level <- x$levels[[l]]
    cat("\n", titles[[l]], "\n\n", sep = "")

    printSummarySection(level$r2.inds, title = "R-squared (indicators):", width.out = width.out)
    printSummarySection(level$r2.etas, title = "R-squared (latents):",    width.out = width.out)

    pt <- level$parTable
    pt$level <- NULL

    if (l < 2L || !NROW(x$agnostic)) {
      plsPrintParTable(pt)
      next
    }

    # Level agnostic parameters (e.g., thresholds) are printed at the end of
    # level 2 (similar to defined parameters in lavaan). The printer needs the
    # other rows for formatting, so we print the combined table twice.
    agnostic <- x$agnostic
    agnostic$level <- NULL
    pt <- rbind(pt, agnostic)

    plsPrintParTable(pt, thresholds = FALSE)
    plsPrintParTable(pt, loadings = FALSE, regressions = FALSE, covariances = FALSE,
                     intercepts = FALSE, variances = FALSE, thresholds = TRUE)
  }

  invisible(x)
}


#' Extract coefficients from a \code{PlsMultilevelModel} model
#'
#' Parameters at the between level (level 2) are suffixed with \code{.l2}.
#'
#' @param object A \code{PlsMultilevelModel} object.
#' @param ... Currently unused.
#' @return A named \code{PlsSemVector} of parameter estimates.
#' @export
setMethod("coef", "PlsMultilevelModel", function(object, ...) {
  plssemVector(object@params$values, is.public = TRUE)
})


#' @rdname coef-PlsMultilevelModel-method
#' @export
setMethod("coefficients", "PlsMultilevelModel", function(object, ...) {
  plssemVector(object@params$values, is.public = TRUE)
})


#' Parameter estimates for \code{PlsMultilevelModel} objects
#'
#' @param object A \code{PlsMultilevelModel} object.
#' @param ... Currently unused.
#' @return A \code{PlsSemParTable} data frame, with a \code{level} column.
#' @export
setMethod("parameter_estimates", "PlsMultilevelModel", function(object, ...) {
  object@parTable
})


#' @rdname is_admissible
#' @export
setMethod("is_admissible", "PlsMultilevelModel", function(object) {
  isTRUE(object@status$is.admissible)
})
