#' Inspect a fitted PLS-SEM model
#'
#' Extract important information from a fitted \code{PlsModel}. The interface is
#' modelled after \code{lavaan::lavInspect()}: a single \code{what} argument
#' selects which piece of information to return.
#'
#' @param object A fitted \code{PlsModel} object.
#' @param what A single string selecting what to extract (case-insensitive);
#'   defaults to \code{"estimates"}. Several values accept aliases, given in
#'   parentheses. One of:
#'   \describe{
#'     \item{\code{"estimates"} (aliases \code{"est"}, \code{"x"},
#'       \code{"matrices"})}{A list of the estimated model matrices in
#'       lavaan-style representation (\code{lambda}, \code{wmat}, \code{theta},
#'       \code{psi}, \code{C}, \code{gamma}).}
#'     \item{\code{"lambda"}, \code{"wmat"}, \code{"theta"}, \code{"psi"},
#'       \code{"C"}, \code{"gamma"}}{The corresponding single matrix from the
#'       \code{"estimates"} list.}
#'     \item{\code{"coef"} (alias \code{"coefficients"})}{The model coefficients.}
#'     \item{\code{"par"} (alias \code{"partable"})}{The parameter table.}
#'     \item{\code{"fit"}}{Fit measures.}
#'     \item{\code{"chisq"}}{The model chi-square statistic.}
#'     \item{\code{"chisq.df"} (alias \code{"df"})}{The chi-square degrees of
#'       freedom.}
#'     \item{\code{"srmr"}}{The standardized root mean square residual.}
#'     \item{\code{"rmsea"}}{The root mean square error of approximation.}
#'     \item{\code{"se"}}{The standard errors of the estimates.}
#'     \item{\code{"vcov"}}{The variance-covariance matrix of the estimates.}
#'     \item{\code{"boot"}}{The bootstrap results.}
#'     \item{\code{"info"}}{A list with the number of observations (\code{nobs}),
#'       the number of latent variables (\code{nlv}), the number of observed
#'       variables (\code{nov}), and the estimation \code{modes} (\code{"A"}/\code{"B"})
#'       for each latent variable.}
#'     \item{\code{"status"}}{A list with the number of \code{iterations}, whether
#'       the algorithm \code{converged}, and whether the solution is
#'       \code{admissible}.}
#'     \item{\code{"qualities"}}{The construct qualities (\eqn{Q^2}).}
#'     \item{\code{"reliabilities"} (alias \code{"rel"})}{The construct
#'       reliabilities.}
#'     \item{\code{"cov.lv"}}{The model-implied covariance matrix of the latent
#'       variables.}
#'     \item{\code{"cov.ov"}}{The model-implied covariance matrix of the observed
#'       variables.}
#'     \item{\code{"cov.all"}}{The joint model-implied covariance matrix of the
#'       observed and latent variables.}
#'     \item{\code{"r2.lv"}, \code{"r2.ov"}, \code{"r2.all"}}{The model-implied
#'       \eqn{R^2} for the latent variables, the observed variables, or both.}
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
setGeneric("pls_inspect", function(object, what = "estimates", ...) {
  standardGeneric("pls_inspect")
})


#' @rdname pls_inspect
#' @export
setMethod("pls_inspect", "PlsModel", function(object, what = "estimates", ...) {
  what <- tolower(what)

  pls_stopif(length(what) != 1L, "`what` must be an argument of length 1!")

  # `object` is passed through untouched: the helpers/methods below decide for
  # themselves whether they need the combined higher-order model (most do so
  # internally, the fit/cov.* family walk the chain explicitly). `combined` is
  # provided for the raw attribute selectors that explicitly rely on it.
  combined <- combinedModel(object)

  switch(what,
    # Fit Measures
    fit      = fitMeasures(object),
    chisq    = pls_chisq(object),
    chisq.df = pls_chisq_df(object),
    df       = pls_chisq_df(object),
    srmr     = pls_srmr(object),
    rmsea    = pls_rmsea(object),

    # Estimates/Parameters 
    x         = plsMatricesLavRep(object),
    est       = plsMatricesLavRep(object),
    estimates = plsMatricesLavRep(object),
    matrices  = plsMatricesLavRep(object),
    lambda    = plsMatricesLavRep(object)$lambda,
    wmat      = plsMatricesLavRep(object)$wmat,
    theta     = plsMatricesLavRep(object)$theta,
    psi       = plsMatricesLavRep(object)$psi,
    c         = plsMatricesLavRep(object)$C,
    gamma     = plsMatricesLavRep(object)$gamma,

    par       = combined@parTable,
    partable  = combined@parTable,
    
    coef         = coef(object),
    coefficients = coef(object),

    # Se/boot
    se        = sqrt(diag(vcov(object))),
    vcov      = vcov(object),
    boot      = pls_boot(object),

    # Information
    info          = plsInspectInfo(object),
    status        = plsInspectStatus(object),
    qualities     = plsConstructQualities(object),
    reliabilities = plsConstructReliabilities(object),
    rel           = plsConstructReliabilities(object),

    # Implied
    cov.lv   = pls_implied_construct_corr(object),
    cov.ov   = pls_implied_indicator_corr(object),
    cov.all  = pls_implied_joint_corr(object),
    r2.lv    = plsImpliedR2(object, output = "lv"),
    r2.ov    = plsImpliedR2(object, output = "ov"),
    r2.all   = plsImpliedR2(object, output = "all"),

    # Data
    data     = plsInspectData(object),

    pls_msg_stop("Unrecognized value for `what`:", what)
  )
})


plsInspectInfo <- function(object) {
  combined <- combinedModel(object)

  cinfo <- modelInfo(combined)
  finfo <- modelInfo(object)

  modes        <- cinfo$modes
  names(modes) <- removeTempAffixes(names(modes))

  list(
    nobs  = finfo$n,
    nlv   = length(cinfo$lvs),
    nov   = length(unique(removeTempAffixes(finfo$allInds))),
    modes = modes
  )
}


plsInspectStatus <- function(object) {
  combined <- combinedModel(object)
  status   <- modelStatus(combined)

  list(
    iterations = status$iterations,
    converged  = isTRUE(status$convergence),
    admissible = isAdmissible(combined)
  )
}


plsInspectData <- function(object) {
  data <- modelData(combinedModel(object))
  colnames(data) <- removeTempAffixes(colnames(data))
  data
}
