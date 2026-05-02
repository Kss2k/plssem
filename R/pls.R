USE_NON_LINEAR_PROBIT_CORR_MAT <- FALSE


#' Fit Partial Least Squares Structural Equation Models
#'
#' \code{pls()} estimates Partial Least Squares Structural Equation Models (PLS-SEM)
#' and their consistent (PLSc) variants. The function accepts \code{lavaan}-style
#' syntax, handles ordered indicators through polychoric correlations and probit
#' factor scores, and supports multilevel specifications expressed with
#' \code{lme4}-style random effects terms inside the structural model.
#'
#' @param syntax Character string with \code{lavaan}-style model syntax describing
#'   both measurement (\code{=~}) and structural (\code{~}) relations. Random effects are
#'   specified with \code{(term | cluster)} statements.
#'
#' @param data A \code{data.frame} or coercible object containing the manifest
#'   indicators referenced in \code{syntax}. Ordered factors are automatically
#'   detected, but can also be supplied explicitly through \code{ordered}.
#'
#' @param standardize Logical; if \code{TRUE}, indicators are standardized before
#'   estimation so that factor scores have comparable scales.
#'
#' @param consistent Logical; \code{TRUE} requests PLSc corrections, whereas \code{FALSE}
#'   fits the traditional PLS model.
#'
#' @param bootstrap Logical; if \code{TRUE}, nonparametric bootstrap standard errors
#'   are computed with \code{boot.R} resamples.
#'
#' @param ordered Optional character vector naming manifest indicators that
#'   should be treated as ordered when computing polychoric correlations.
#'
#' @param missing Character string specifying how to handle missing indicator data.
#'   \code{"listwise"} removes rows with missing values (listwise deletion).
#'   \code{"mean"} imputes missing indicator values using simple univariate
#'   imputation: the mean for continuous variables, the median for ordered variables
#'   with more than two categories, and the mode for binary ordered variables (two
#'   categories) or nominal variables.
#'   \code{"kNN"} (or \code{"knn"}) imputes missing indicator values using
#'   k-nearest neighbors imputation (kNN). When \code{missing = "kNN"}, rows with
#'   all indicators missing are removed prior to imputation, and rows with missing
#'   \code{cluster} values are removed for multilevel models.
#'
#' @param knn.k Integer specifying the number of neighbors (\code{k}) used when
#'   \code{missing = "kNN"}.
#'
#' @param mcpls Should a Monte-Carlo consistency correction be applied?
#'
#' @param probit Logical; overrides the automatic choice of probit factor scores
#'   that is based on whether ordered indicators are present.
#'
#' @param tolerance Numeric; Convergence criteria/tolerance.
#'
#' @param max.iter.0_5 Maximum number of PLS iterations performed when estimating
#'   the measurement and structural models.
#'
#' @param boot.ncpus Integer: number of processes to be used in parallel operation.
#'   By default this is the number of cores (as detected by \code{parallel::detectCores()}) minus one.
#'
#' @param boot.parallel The type of parallel operation to be used (if any). If missing,
#'   the default is \code{"no"}.
#'
#' @param boot.R Integer giving the number of bootstrap resamples drawn when
#'   \code{bootstrap = TRUE}.
#'
#' @param boot.iseed An integer to set the bootstrap seed. Or \code{NULL} if no
#'   reproducible results are needed. This works for both serial (non-parallel) and
#'   parallel settings. Internally, \code{RNGkind()} is set to \code{"L'Ecuyer-CMRG"}
#'   if \code{parallel = "multicore"}. If \code{parallel = "snow"} (under windows),
#'   \code{parallel::clusterSetRNGStream()} is called which automatically switches to
#'   \code{"L'Ecuyer-CMRG"}. When iseed is not \code{NULL}, \code{.Random.seed}
#'   (if it exists) in the global environment is left untouched.
#'
#' @param sample DEPRECATED. Integer giving the number of bootstrap resamples drawn when
#'   \code{bootstrap = TRUE}.
#'
#' @param mc.min.iter Minimum number of iterations in MC-PLS algorithm.
#'
#' @param mc.max.iter Maximum number of iterations in MC-PLS algorithm.
#'
#' @param mc.reps Monte-Carlo sample size in MC-PLS algorithm.
#'
#' @param mc.tol Tolerance in MC-PLS algorithm.
#'
#' @param mc.fixed.seed Should a fixed seed be used in the MC-PLS algorithm?
#'   Setting a fixed seed will likely yield less accurate estimates, but can
#'   substantially improve the stability and computational efficiency of the
#'   algorithm.
#'
#' @param mc.polyak.juditsky Should the polyak.juditsky running average method
#'   be applied in the MC-PLS algorithm?
#'
#' @param mc.fn.args Additional arguments to MC-PLS algorithm, mainly for controlling
#'   the step size.
#'
#' @param verbose Should verbose output be printed?
#'
#' @param boot.optimize Logical; if \code{TRUE} and \code{bootstrap = TRUE}, applies
#'   the settings in \code{mc.boot.control} inside each bootstrap replicate (MC-PLS only).
#'
#' @param mc.boot.control List of control parameters passed to the MC-PLS algorithm
#'   inside each bootstrap replicate when \code{boot.optimize = TRUE}.
#'
#' @param reliabilities Optional named numeric vector of user-supplied reliabilities
#'   used for the PLSc consistency correction.
#'
#' @param ... Internal arguments. For advanced users only.
#'
#' @return A \code{Plssem} object containing the estimated parameters, fit measures,
#'   factor scores, and any bootstrap results. Methods such as \code{summary()},
#'   \code{coef()}, and \code{parameter_estimates()} can be applied to inspect the fit.
#'
#' @seealso \code{\link[=summary,PlsModel-method]{summary}},
#'   \code{\link[=show,PlsModel-method]{show}}
#'
#' @examples
#' \donttest{
#' library(plssem)
#' library(modsem)
#'
#' tpb <- '
#'   ATT =~ att1 + att2 + att3 + att4 + att5
#'   SN =~ sn1 + sn2
#'   PBC =~ pbc1 + pbc2 + pbc3
#'   INT =~ int1 + int2 + int3
#'   BEH =~ b1 + b2
#'   INT ~ ATT + SN + PBC
#'   BEH ~ INT + PBC
#' '
#'
#' fit <- pls(tpb, TPB, bootstrap = TRUE)
#' summary(fit)
#' }
#' @export
pls <- function(syntax,
                data,
                standardize = TRUE,
                consistent = TRUE,
                bootstrap = FALSE,
                ordered = NULL,
                missing = c("listwise", "mean", "kNN"),
                knn.k = 5,
                mcpls = NULL,
                probit = NULL,
                tolerance = 1e-5,
                max.iter.0_5 = 100L,
                boot.ncpus = 1L,
                boot.parallel = c("no", "multicore", "snow"),
                boot.R = 50L,
                boot.iseed = NULL,
                sample = NULL,
                mc.min.iter = 5L,
                mc.max.iter = 250L,
                mc.reps = 20000L,
                mc.tol = 1e-3,
                mc.fixed.seed = FALSE,
                mc.polyak.juditsky = FALSE,
                mc.fn.args = list(),
                verbose = interactive(),
                boot.optimize = TRUE,
                mc.boot.control = list(
                  min.iter        = mc.min.iter,
                  max.iter        = mc.max.iter,
                  mc.reps         = floor(0.5 * mc.reps),
                  tol             = 2 * mc.tol,
                  polyak.juditsky = FALSE,
                  verbose         = FALSE,
                  fixed.seed      = TRUE,
                  reuse.p.start   = TRUE
                ),
                reliabilities = NULL,
                ...) {

  missing       <- match.arg(tolower(missing), c("listwise", "mean", "knn"))
  boot.parallel <- match.arg(boot.parallel)

  if (!is.null(sample)) {
    warning("The sample argument is deprecated, please use the boot.R argument instead!")
    boot.R <- sample
  }

  data <- as.data.frame(data)

  model <- specifyModel(
    syntax             = syntax,
    data               = data,
    consistent         = consistent,
    missing            = missing,
    standardize        = standardize,
    ordered            = ordered,
    probit             = probit,
    mcpls              = mcpls,
    tolerance          = tolerance,
    max.iter.0_5       = max.iter.0_5,
    mc.min.iter        = mc.min.iter,
    mc.max.iter        = mc.max.iter,
    mc.reps            = mc.reps,
    mc.tol             = mc.tol,
    mc.fixed.seed      = mc.fixed.seed,
    mc.polyak.juditsky = mc.polyak.juditsky,
    mc.fn.args         = mc.fn.args,
    verbose            = verbose,
    bootstrap          = bootstrap,
    boot.ncpus         = boot.ncpus,
    boot.parallel      = boot.parallel,
    boot.R             = boot.R,
    boot.iseed         = boot.iseed,
    boot.optimize      = boot.optimize,
    mc.boot.control    = mc.boot.control,
    knn.k              = knn.k,
    reliabilities      = reliabilities,
    ...
  )

  model <- estimatePLS(model = model)

  if (isTRUE(modelInfo(model)$boot$bootstrap)) {
    boot <- tryCatch(
      bootstrap(model),
      error = \(e) {
        warning2("Bootstrapping FAILED!\nMessage: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(boot))
      modelBoot(model) <- boot
  }

  cm <- combinedModel(model)
  cm@parTable <- getParTableEstimates(cm)

  if (hasCombinedModel(model) || hasHigherOrderModel(model)) {
    model@combinedModel <- cm
  } else {
    model@parTable <- cm@parTable
  }

  model
}


resetPLS_ModelLowerOrder <- function(model, hard.reset = FALSE) {
  resetModelStatusLowerOrder(model, hard.reset = hard.reset)
}


estimatePLS_InnerLocal <- function(model) {
  model |>
    updateOuterWeights() |>
    updateFactorScores() |>
    updateFitObjects()   |>
    updateParamVector()  |>
    updateEstimationStatus()
}


estimatePLS_Inner <- function(model) {
  estimateHigherOrderChain(model)
}


estimatePLS_Outer <- function(model, ...) {
  force(model)

  if (is_mcpls(model))
    return(mcpls(model, ...))

  model
}


estimatePLS_Status <- function(model, ...) {
  prev    <- isTRUE(model@status$is.admissible)
  current <- modelFitIsAdmissible(model@fit)

  model@status$is.admissible <- prev && current
  model
}


estimatePLS <- function(model, ...) {
  tryCatch({
    model |>
      estimatePLS_Inner() |>
      estimatePLS_Outer(...) |>
      updateEstimationStatus()

  }, error = function(e) {
    stop2("Model estimation FAILED!\n", "Message: ", conditionMessage(e))
  })
}
