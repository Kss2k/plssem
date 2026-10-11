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
#' @param fit Logical; Whether to compute fit measures (for each level).
#' @param unstandardized Logical; Should unstandardized estimates be included?
#' @param ci Logical; Should confidence intervals for parameter estimates be included?
#' @param ... Arguments passed to \code{\link{unstandardized_estimates}}.
#' @return A \code{SummaryPlsMultilevel} list with formatted results.
#' @export
setMethod("summary", "PlsMultilevelModel", function(object, fit = TRUE, unstandardized = FALSE,
                                                    ci = FALSE, ...) {
  parTable   <- parameter_estimates(object)
  info       <- object@info
  extra.cols <- NULL

  if (unstandardized) {
    tryCatch({

      # we don't need any standard errors, so save the computation
      parTableu <- unstandardized_estimates(object, se = "none", ...)

      idx0 <- paste0(parTable$lhs, parTable$op, parTable$rhs, parTable$level)
      idx1 <- paste0(parTableu$lhs, parTableu$op, parTableu$rhs, parTableu$level)

      parTable$Unstd <- parTableu$est[match(idx0, idx1)]
      extra.cols <- c(extra.cols, "Unstd")

    }, error = function(e) {

      pls_msg_warn(
        "Calculation of unstandardized estimates failed!",
        "Message:", conditionMessage(e)
      )

    })
  }

  if (ci) {
    extra.cols <- c(extra.cols, "ci.lower", "ci.upper")
  }

  if (fit) fit.measures <- fit_measures(object)
  else fit.measures <- NULL

  levels <- lapply(1:2, FUN = \(l) {
    pt <- parTable[parTable$level %in% l, , drop = FALSE]

    list(
      parTable     = pt,
      r2.etas      = getR2ParTable(getEtas(pt, checkAny = FALSE), parTable = pt),
      r2.inds      = getR2ParTable(getReflectiveIndicators(pt), parTable = pt),
      fit.measures = fit.measures[[l]]
    )
  })

  icc <- object@params$icc
  rsd <- object@params$rsd

  icc <- stats::setNames(formatNumeric(icc), nm = names(icc))
  rsd <- stats::setNames(formatNumeric(rsd), nm = names(rsd))
  rslopes <- info$rslopes

  if (length(rsd)) {
    paths <- paste0(rslopes$lhs, " ~ ", rslopes$rhs)
    names(rsd) <- paste0(names(rsd), " (", paths[match(names(rsd), rslopes$name)], ")")
  }

  out <- list(
    fit      = object,
    levels   = levels,
    agnostic = parTable[is.na(parTable$level), , drop = FALSE], # e.g., thresholds
    print  = list(
      width      = plsGetWidthPrintedParTable(parTable),
      extra.cols = extra.cols
    ),
    info   = list(
      estimator  = info$estimator,
      link       = if (length(info$ordered)) "PROBIT" else "LINEAR",
      n          = info$n,
      nclusters  = info$nclusters,
      iterations = object@status$iterations,
      se         = if (length(object@boot)) "Delta (cluster bootstrap)" else "None"
    ),
    icc = icc,
    rsd = rsd,
    fit.measures = fit.measures
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
    c(x$info$estimator, x$info$link, x$info$se, "", x$info$n, x$info$nclusters,
      x$info$iterations),
    nm = c("Estimator", "Link", "Standard errors", "", "Number of observations",
           "Number of clusters", "Number of iterations")
  ))

  printSummarySection(x$icc, title = "Intraclass correlations:", width.out = width.out)
  printSummarySection(x$rsd, title = "Random slopes (standard deviations):", width.out = width.out)

  fm <- x$fit.measures
  if (!is.null(fm)) {

    .levelValues <- function(i) {
      c(
        "  "                 = c("Within", "Between")[[i]],
        "Chi-Square"         = sprintf("%.3f", fm[[i]]$chisq),
        "Degrees of Freedom" = sprintf("%d",   as.integer(fm[[i]]$chisq.df)),
        "SRMR"               = sprintf("%.3f", fm[[i]]$srmr),
        "RMSEA"              = sprintf("%.3f", fm[[i]]$rmsea)
      )
    }

    printSummarySection(
      title        = "Fit Measures:",
      width.out    = width.out,
      values       = .levelValues(1L),
      values.extra = .levelValues(2L)
    )
  }


  titles <- c("Level 1 [within]:", "Level 2 [between]:")

  for (l in 1:2) {
    level <- x$levels[[l]]
    cat("\n", titles[[l]], "\n\n", sep = "")

    printSummarySection(level$r2.inds, title = "R-squared [indicators]:", width.out = width.out)
    printSummarySection(level$r2.etas, title = "R-squared [latents]:",    width.out = width.out)

    pt <- level$parTable
    pt$level <- NULL

    extra.cols <- x$print$extra.cols

    if (l < 2L || !NROW(x$agnostic)) {
      plsPrintParTable(pt, extra.cols = extra.cols)
      next
    }

    # Level agnostic parameters (e.g., thresholds) are printed at the end of
    # level 2 (similar to defined parameters in lavaan). The printer needs the
    # other rows for formatting, so we print the combined table twice.
    agnostic <- x$agnostic
    agnostic$level <- NULL
    pt <- rbind(pt, agnostic)

    plsPrintParTable(pt, thresholds = FALSE, extra.cols = extra.cols)
    plsPrintParTable(pt, loadings = FALSE, regressions = FALSE, covariances = FALSE,
                     intercepts = FALSE, variances = FALSE, thresholds = TRUE,
                     extra.cols = extra.cols)
  }

  invisible(x)
}


#' Fit measures for \code{PlsMultilevelModel} objects
#'
#' Computes the fit measures (chi-square, SRMR, RMSEA) for each level. The
#' observed (auxiliary) correlation matrices of the levels are compared with
#' the ones implied by the model.
#'
#' @param object A \code{PlsMultilevelModel} object.
#' @param saturated Logical; if \code{TRUE}, compute the saturated fit.
#' @param mc.reps Integer; number of observations in the simulated data.
#' @param ... Currently unused.
#' @return A list with the fit measures of each level (\code{level.1} and \code{level.2}).
#' @export
setMethod("fit_measures", "PlsMultilevelModel", function(object, saturated = FALSE, mc.reps = 1e6, ...) {
  fitMeasuresMultilevel(object, saturated = saturated, mc.reps = mc.reps)
})


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


#' Extract the variance-covariance matrix from a \code{PlsMultilevelModel} model
#'
#' Delta-method (co-)variances, see the \code{bootstrap} argument of \code{mpls()}.
#' Parameters at the between level (level 2) are suffixed with \code{.l2}.
#'
#' @param object A \code{PlsMultilevelModel} object.
#' @param ... Currently unused.
#' @return A \code{PlsSemMatrix}, or \code{NULL} if no standard errors were computed.
#' @export
setMethod("vcov", "PlsMultilevelModel", function(object, ...) {
  object@params$vcov
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


#' Unstandardized estimates for \code{PlsMultilevelModel} objects
#'
#' Computes the unstandardized estimates of each level of a two-level model.
#'
#' @details
#' The levels are unstandardized separately, using the (decomposed) scales of
#' the variables at each level, i.e., \code{sqrt(1 - icc) * sd} at level 1, and
#' \code{sqrt(icc) * sd} at level 2 for variables at both levels. The latent
#' variables are scaled by their first indicators. The scales of the random
#' slopes (\code{y ~ rv(s) * x}) are \code{rsd * sd(y) / sd(x)}, where
#' \code{sd(y)} and \code{sd(x)} are the scales of \code{y} and \code{x} at
#' level 1. For ordered variables, the latent responses are scaled such that
#' their (total) residual variances are 1, i.e.,
#' \code{(1 - icc) * theta_W + icc * theta_B = 1}, where \code{theta_W} and
#' \code{theta_B} are the standardized residual variances at each level. The
#' thresholds are rescaled accordingly.
#'
#' @inheritParams unstandardized_estimates
#' @return A \code{PlsSemParTable} data frame, with a \code{level} column.
#' @export
setMethod("unstandardized_estimates", "PlsMultilevelModel",
          function(model,
                   unstandardized    = "all",
                   se                = c("delta", "none"),
                   scale.uncertainty = FALSE,
                   eps               = 1e-4,
                   zero.tol          = 1e-10,
                   rm.tmp.ov         = TRUE,
                   clean.tmp.ind     = TRUE,
                   clean.tmp.mimic   = TRUE) {
  plsUnstandardizedEstimatesMultilevel(
    model             = model,
    unstandardized    = unstandardized,
    se                = se,
    scale.uncertainty = scale.uncertainty,
    eps               = eps,
    zero.tol          = zero.tol,
    rm.tmp.ov         = rm.tmp.ov,
    clean.tmp.ind     = clean.tmp.ind,
    clean.tmp.mimic   = clean.tmp.mimic
  )
})


#' @rdname is_admissible
#' @export
setMethod("is_admissible", "PlsMultilevelModel", function(object) {
  isTRUE(object@status$is.admissible)
})
