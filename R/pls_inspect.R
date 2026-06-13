#' Inspect a fitted PLS-SEM model
#'
#' Extract important information from a fitted \code{PlsModel}. The interface is
#' modelled after \code{lavaan::lavInspect()}: a single \code{what} argument
#' selects which piece of information to return.
#'
#' @param object A fitted \code{PlsModel} object.
#' @param what A single string selecting what to extract (case-insensitive).
#'   One of:
#'   \describe{
#'     \item{\code{"fit"}}{The model fit object, i.e. \code{\link{modelFit}(object)}.}
#'     \item{\code{"info"}}{A list with the number of observations (\code{nobs}),
#'       the number of latent variables (\code{nlv}), the number of observed
#'       variables (\code{nov}), and the estimation \code{modes} (\code{"A"}/\code{"B"})
#'       for each latent variable.}
#'     \item{\code{"status"}}{A list with the number of \code{iterations}, whether
#'       the algorithm \code{converged}, and whether the solution is
#'       \code{admissible}.}
#'     \item{\code{"cov.lv"}}{The model-implied covariance matrix of the latent
#'       variables.}
#'     \item{\code{"cov.ov"}}{The model-implied covariance matrix of the observed
#'       variables.}
#'     \item{\code{"cov.all"}}{The joint model-implied covariance matrix of the
#'       observed and latent variables.}
#'     \item{\code{"data"}}{The (standardized) data matrix used for estimation.}
#'   }
#' @param ... Currently ignored.
#'
#' @return The requested information; the type depends on \code{what} (see above).
#'
#' @examples
#' \dontrun{
#'   fit <- pls(model, data = data)
#'   pls_inspect(fit, "info")
#'   pls_inspect(fit, "cov.lv")
#' }
#'
#' @export
setGeneric("pls_inspect", function(object, what = "fit", ...) {
  standardGeneric("pls_inspect")
})


#' @rdname pls_inspect
#' @export
setMethod("pls_inspect", "PlsModel", function(object, what = "fit", ...) {
  what <- match.arg(tolower(what), choices = c(
    "fit", "info", "status", "cov.lv", "cov.ov", "cov.all", "data"
  ))

  # The cov.* matrices walk the higher-order model chain explicitly (inside the
  # implied_*() methods), so they must operate on the model as passed. Everything
  # else should report the combined higher-order state (for first-order models
  # combinedModel() simply returns the object unchanged).
  if (!what %in% c("cov.lv", "cov.ov", "cov.all"))
    object <- combinedModel(object)

  switch(what,
    fit     = modelFit(object),
    info    = plsInspectInfo(object),
    status  = plsInspectStatus(object),
    cov.lv  = implied_construct_corr(object),
    cov.ov  = implied_indicator_corr(object),
    cov.all = implied_joint_corr(object),
    data    = plsInspectData(object)
  )
})


plsInspectInfo <- function(object) {
  info <- modelInfo(object)

  modes        <- info$modes
  names(modes) <- removeTempAffixes(names(modes))

  list(
    nobs  = info$n,
    nlv   = length(info$lvs),
    nov   = length(unique(removeTempAffixes(info$allInds))),
    modes = modes
  )
}


plsInspectStatus <- function(object) {
  status <- modelStatus(object)

  list(
    iterations = status$iterations,
    converged  = isTRUE(status$convergence),
    admissible = isAdmissible(object)
  )
}


plsInspectData <- function(object) {
  data <- modelData(object)
  colnames(data) <- removeTempAffixes(colnames(data))
  data
}
