# S4 generics and methods for the PlsModel class.
# Internal generics follow camelCase, whilst public ones follow snake_case
# Replaces the former S3 methods (summary.plssem, print.plssem, coef.plssem, …).


#' Show a \code{PlsModel} object
#'
#' Called automatically when an object is printed at the prompt.  Displays
#' the package version, iteration count, and the parameter table.
#'
#' @param object A \code{PlsModel} object.
#' @return \code{object}, invisibly.
#' @export
setMethod("show", "PlsModel", function(object) {
  combined <- combinedModel(object)
  admissible <- isAdmissible(combined)

  statusString <- if (admissible) "ended normally" else "did NOT END NORMALLY"

  printf("plssem (%s) %s after %i iterations\n",
         PKG_INFO$version, statusString, combined@status$iterations)

  parTable <- parameter_estimates(combined)

  if (NROW(parTable)) print(parTable)
  else                print(object@parTableInput)

  invisible(object)
})


#' Summarize a fitted \code{PlsModel} model
#'
#' @param object A \code{PlsModel} object.
#' @param fit Logical; whether to compute fit measures.
#' @param ... Currently unused.
#' @return A \code{SummaryPlsSem} list with formatted results.
#' @export
setMethod("summary", "PlsModel", function(object, fit = TRUE, ...) {
  combined <- combinedModel(object)
  parTable <- parameter_estimates(combined)

  lvs    <- getLVs(parTable)
  ovs    <- getOVs(parTable)
  etas   <- getEtas(parTable, checkAny = FALSE)
  inds   <- getIndicators(parTable)
  inds.a <- getReflectiveIndicators(parTable)

  strParTableLines <- utils::capture.output(modsem::summarize_partable(parTable))
  strParTable <- paste0(paste0(strParTableLines[-(1:6)], collapse = "\n"), "\n")

  is.ord <- combined@info$is.probit || (length(combined@info$ordered) && combined@info$is.mcpls)
  link   <- if (is.ord) "PROBIT" else "LINEAR"

  getR2 <- function(x, pt = parTable) {
    rvar <- pt[pt$lhs == x & pt$op == "~~" & pt$rhs == x, "est"]
    if (!length(rvar)) 0 else 1 - rvar
  }

  r2.etas <- vapply(etas,   FUN.VALUE = numeric(1L), FUN = getR2)
  r2.inds <- vapply(inds.a, FUN.VALUE = numeric(1L), FUN = getR2)

  out <- list(
    print = list(
      strParTable = strParTable,
      width       = max(nchar(strParTableLines))
    ),
    fit  = object,
    info = list(
      iterations = combined@status$iterations,
      estimator  = combined@info$estimator,
      n          = combined@info$n,
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
    fit.measures = if (fit) fit_measures(object) else NULL
  )

  class(out) <- "SummaryPlsSem"
  out
})


#' Print a \code{SummaryPlsSem} object
#'
#' @param x A \code{SummaryPlsSem} object as returned by
#'   \code{\link[=summary,PlsModel-method]{summary}()}.
#' @param ... Additional arguments for compatibility with the generic.
#' @return The input object, invisibly.
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
    headerNames <- c("Chi-Square", "Degrees of Freedom", "SRMR", "RMSEA")
    headerValues <- c(
      sprintf("%.3f", x$fit.measures$chisq),
      sprintf("%d",   x$fit.measures$chisq.df),
      sprintf("%.3f", x$fit.measures$srmr),
      sprintf("%.3f", x$fit.measures$rmsea)
    )
    cat(allignLhsRhs(lhs = headerNames, rhs = headerValues, pad = "  ",
                     width.out = width.out), "\n", sep = "")
  }

  if (length(x$r2$inds)) {
    cat("R-squared (indicators):\n")
    cat(allignLhsRhs(lhs = names(x$r2$inds),
                     rhs = formatNumeric(x$r2$inds), pad = "  ",
                     width.out = width.out), "\n", sep = "")
  }

  if (length(x$r2$etas)) {
    cat("R-squared (latents):\n")
    cat(allignLhsRhs(lhs = names(x$r2$etas),
                     rhs = formatNumeric(x$r2$etas), pad = "  ",
                     width.out = width.out), "\n", sep = "")
  }

  cat(x$print$strParTable)
  invisible(x)
}


#' Extract coefficients from a \code{PlsModel} model
#'
#' @param object A \code{PlsModel} object.
#' @param ... Currently unused.
#' @return A named \code{PlsSemVector} of parameter estimates.
#' @importFrom stats coef
#' @export
setMethod("coef", "PlsModel", function(object, ...) {
  combined <- combinedModel(object)
  plssemVector(combined@params$values, is.public = TRUE)
})


#' @rdname coef-PlsModel-method
#' @importFrom stats coefficients
#' @export
setMethod("coefficients", "PlsModel", function(object, ...) {
  coef(object, ...)
})


#' Extract the variance-covariance matrix from a \code{PlsModel} model
#'
#' @param object A \code{PlsModel} object.
#' @param ... Currently unused.
#' @return A \code{PlsSemMatrix} (bootstrap-based vcov, or \code{NULL}).
#' @importFrom stats vcov
#' @export
setMethod("vcov", "PlsModel", function(object, ...) {
  combined <- combinedModel(object)
  plssemMatrix(combined@params$vcov, is.public = TRUE)
})


#' Generic accessor for model parameter estimates
#'
#' @param object A fitted model object.
#' @param ... Additional arguments passed to methods.
#' @return A parameter table describing the fitted model.
#' @export
setGeneric("parameter_estimates",
           function(object, ...) standardGeneric("parameter_estimates"))


#' Parameter estimates for \code{PlsModel} objects
#'
#' @param object A \code{PlsModel} object.
#' @param colon.pi Logical; replace product-indicator labels with colon
#'   notation (\code{X:Z}).
#' @param label.renamed.prod Logical; retain renamed product labels when colon
#'   expansion occurs.
#' @param ... Currently unused.
#' @return A \code{PlsSemParTable} data frame.
#' @export
setMethod("parameter_estimates", "PlsModel",
          function(object, colon.pi = TRUE, label.renamed.prod = FALSE, ...) {
  object <- combinedModel(object)
  parTable <- object@parTable

  if (colon.pi)
    parTable <- addColonPI_ParTable(parTable, model = object,
                                    label.renamed.prod = label.renamed.prod)
  parTable
})


#' Check whether an object uses the MC-OrdPLSc estimator
#'
#' @param object A fitted model object.
#' @return \code{TRUE} or \code{FALSE}.
#' @export
setGeneric("is_mcpls", function(object) standardGeneric("is_mcpls"))

#' @rdname is_mcpls
#' @export
setMethod("is_mcpls", "PlsModel", function(object) {
  object <- combinedModel(object)
  isTRUE(object@info$is.mcpls)
})


#' Retrieve bootstrap coefficient matrix
#'
#' @param object A fitted model object.
#' @return A \code{PlsSemMatrix} of bootstrap replicate parameter vectors
#'   (rows = replicates, cols = parameters).
#'
#' @examples
#' library(modsem)
#' library(plssem)
#'
#' m <- "
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#'   Y ~ X + Z + X:Z
#' "
#'
#' fit <- pls(m, oneInt, bootstrap = TRUE, boot.R = 50)
#' boot(fit)
#'
#' @export
setGeneric("boot", function(object) standardGeneric("boot"))

#' @rdname boot
#' @export
setMethod("boot", "PlsModel", function(object) {
  object <- combinedModel(object)
  plssemMatrix(object@boot$boot, is.public = TRUE)
})



#' Check whether a fitted model has admissible parameter estimates
#'
#' @param object A fitted \code{PlsModel} object.
#' @return A single logical value.
#' @export
setGeneric("is_admissible", function(object) standardGeneric("is_admissible"))


#' @rdname is_admissible
#' @export
setMethod("is_admissible", "PlsModel", function(object) isAdmissible(object))


#' Implied Construct Correlation Matrix
#'
#' Returns the implied construct correlation matrix for a fitted model.
#'
#' For higher-order models, this is computed for the combined model returned by
#' [combinedModel()].
#'
#' @param object A fitted [PlsModel] object.
#' @param saturated Logical; if `TRUE`, return the saturated implied matrix.
#' @param mc.reps Integer; number of Monte Carlo resamples used for MC-PLSc.
#' @param ... Reserved for future extensions.
#' @return A [PlsSemMatrix].
#' @export
setGeneric(
  "implied_construct_corr",
  function(object, saturated = FALSE, mc.reps = 1e6, ...) standardGeneric("implied_construct_corr")
)

#' @rdname implied_construct_corr
#' @export
setMethod("implied_construct_corr", "PlsModel",
          function(object, saturated = FALSE, mc.reps = 1e6, ...) {
  plssemMatrix(
    impliedConstructCorrMat(object, saturated = saturated, mc.reps = mc.reps),
    is.public = TRUE
  )
})


#' Implied Indicator Correlation Matrix
#'
#' Returns the implied indicator correlation matrix for a fitted model.
#'
#' For higher-order models, this is computed for the combined model returned by
#' [combinedModel()].
#'
#' @param object A fitted [PlsModel] object.
#' @param saturated Logical; if `TRUE`, return the saturated implied matrix.
#' @param mc.reps Integer; number of Monte Carlo resamples used for MC-PLSc.
#' @param ... Reserved for future extensions.
#' @return A numeric matrix.
#' @export
setGeneric(
  "implied_indicator_corr",
  function(object, saturated = FALSE, mc.reps = 1e6, ...) standardGeneric("implied_indicator_corr")
)


#' @rdname implied_indicator_corr 
#' @export
setMethod("implied_indicator_corr", "PlsModel",
          function(object, saturated = FALSE, mc.reps = 1e6, ...) {
  plssemMatrix(
    impliedIndicatorCorrMat(object, saturated = saturated, mc.reps = mc.reps),
    is.public = TRUE
  )
})


#' Fit Measures
#'
#' Computes global fit measures (e.g., chi-square, SRMR, RMSEA) for a fitted
#' model.
#'
#' @param object A fitted [PlsModel] object.
#' @param saturated Logical; if `TRUE`, compute the saturated fit.
#' @param mc.reps Integer; number of Monte Carlo resamples used for MC-PLSc fit.
#' @param ... Reserved for future extensions.
#' @return A named list with fit statistics.
#' @export
setGeneric(
  "fit_measures",
  function(object, saturated = FALSE, mc.reps = 1e6, ...) standardGeneric("fit_measures")
)


#' @rdname fit_measures
#' @export
setMethod("fit_measures", "PlsModel", function(object, saturated = FALSE, mc.reps = 1e6, ...) {
  fitMeasures(object, saturated = saturated, mc.reps = mc.reps)
})
