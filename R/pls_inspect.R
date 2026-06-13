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
# #export
setGeneric("pls_inspect", function(object, what = "fit", ...) {
  standardGeneric("pls_inspect")
})


#' @rdname pls_inspect
#' @export
setMethod("pls_inspect", "PlsModel", function(object, what = "fit", ...) {
  what <- match.arg(tolower(what), choices = c(
    "fit", "info", "status", "cov.lv", "cov.ov", "cov.all", "data"
  ))

  switch(what,
    fit     = modelFit(object),
    info    = plsInspectInfo(object),
    status  = plsInspectStatus(object),
    cov.lv  = plsInspectCovLv(object),
    cov.ov  = plssemMatrix(impliedIndicatorCorrMat(object), is.public = TRUE),
    cov.all = plsInspectCovAll(object),
    data    = plsInspectData(object)
  )
})


# Genuine latent variables, i.e. excluding the single-indicator stand-ins that
# the parser creates for observed structural variables (e.g. `y ~ x`) and for
# the reflective indicators of a composite-MIMIC block (`A =~ y; A <~ x`). Such
# a stand-in has exactly one indicator whose (cleaned) name equals the (cleaned)
# latent name; a real latent is named differently from its indicator(s).
genuineLatentVars <- function(object) {
  info    <- modelInfo(object)
  lvs     <- info$lvs.linear
  indsLvs <- info$indsLvs

  isStandin <- vapply(lvs, FUN.VALUE = logical(1L), FUN = function(lv) {
    inds <- indsLvs[[lv]]
    length(inds) == 1L && removeTempAffixes(inds) == removeTempAffixes(lv)
  })

  lvs[!isStandin]
}


plsInspectInfo <- function(object) {
  info <- modelInfo(object)
  lvs  <- genuineLatentVars(object)

  modes        <- info$modes[lvs]
  names(modes) <- removeTempAffixes(names(modes))

  list(
    nobs  = info$n,
    nlv   = length(lvs),
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


plsInspectCovLv <- function(object) {
  Phi <- impliedConstructCorrMat(object)
  lvs <- intersect(genuineLatentVars(object), colnames(Phi))

  plssemMatrix(Phi[lvs, lvs, drop = FALSE], is.public = TRUE)
}


# Joint model-implied covariance matrix of observed and latent variables.
# The observed-observed and latent-latent blocks reuse the package's canonical
# implied-covariance helpers; the cross block is Cov(y, eta) = Lambda %*% Phi.
plsInspectCovAll <- function(object) {
  Phi    <- impliedConstructCorrMat(object)
  SigmaO <- impliedIndicatorCorrMat(object)

  ovs    <- rownames(SigmaO)
  allLvs <- colnames(Phi)

  # Full loading matrix [observed x all latents], so that Cov(observed, latent)
  # picks up associations that run through single-indicator stand-in latents.
  Lambda <- matrix(0, nrow = length(ovs), ncol = length(allLvs),
                   dimnames = list(ovs, allLvs))

  fitLambda <- object@fit$fitLambda
  ovs.f     <- intersect(ovs, rownames(fitLambda))
  lvs.f     <- intersect(allLvs, colnames(fitLambda))
  Lambda[ovs.f, lvs.f] <- fitLambda[ovs.f, lvs.f]

  crossFull <- Lambda %*% Phi # Cov(observed, each latent)

  # Keep only genuine latents in the latent block / cross columns.
  lvs   <- intersect(genuineLatentVars(object), allLvs)
  cross <- crossFull[, lvs, drop = FALSE]
  Phi   <- Phi[lvs, lvs, drop = FALSE]

  covall <- rbind(
    cbind(SigmaO, cross),
    cbind(t(cross), Phi)
  )

  plssemMatrix(covall, is.public = TRUE)
}
