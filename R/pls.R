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
#' @param mcpls Should the model be estimated using the Monte-Carlo Consistent
#'   Partial Least Squares (MC-PLSc) algorithm?
#'
#' @param mc.fast.lmer Should a faster (biased) GLS based estimator of the
#'   Mixed-Effects model be used in conjunction with the MC-PLS algorithm?
#'
#' @param probit Logical; overrides the automatic choice of probit factor scores
#'   that is based on whether ordered indicators are present.
#'
#' @param tolerance Numeric; Convergence criteria/tolerance.
#'
#' @param max.iter.0_5 Maximum number of PLS iterations performed when estimating
#'   the measurement and structural models.
#'
#' @param boot.ncores Integer: number of workers to be used for parallel bootstrapping.
#'   Parallel bootstrapping is enabled when \code{boot.ncores > 1}.
#'
#' @param boot.ncpus Deprecated alias for \code{boot.ncores}.
#'
#' @param boot.parallel The type of parallel operation to be used (if any). The
#'   default is \code{"no"}. \code{"multisession"} runs the bootstrap using multiple
#'   background \code{R} sessions (works on all platforms), while \code{"multicore"}
#'   uses forked processes (not available on Windows). \code{"snow"} is kept for
#'   backwards compatibility and is treated as an alias for \code{"multisession"}.
#'   Internally this is implemented using the \code{future} package 
#'
#' @param boot.R Integer giving the number of bootstrap resamples drawn when
#'   \code{bootstrap = TRUE}.
#'
#' @param boot.iseed An integer to set the bootstrap seed. Or \code{NULL} if no
#'   reproducible results are needed. This works for both serial and parallel
#'   settings. When \code{boot.iseed} is not \code{NULL}, \code{.Random.seed} (if it
#'   exists) in the global environment is left untouched.
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
#' @param mc.fixed.seed Should a fixed seed be used in the MC-PLS algorithm?
#'   Setting a fixed seed will likely yield less accurate estimates, but can
#'   substantially improve the stability and computational efficiency of the
#'   algorithm.
#'
#' @param mc.polyak.juditsky Should the polyak.juditsky running average method
#'   be applied in the MC-PLS algorithm?
#'
#' @param mc.pj.extrapolate Logical; if \code{TRUE} (the default), the Polyak-Juditsky
#'   convergence point is estimated via NLS exponential extrapolation (with
#'   Aitken \eqn{\delta^2} as a fallback). If \code{FALSE} a warm start is performed
#'   instead, and the plain Polyak-Juditsky average is used.
#'   Only relevant when \code{mc.polyak.juditsky = TRUE}.
#'
#' @param mc.tol Tolerance in MC-PLS algorithm.
#'
#' @param mc.delta.se Should delta-method standard errors be computed for
#'   MC-PLS estimates?
#'
#' @param mc.delta.jacobian.k Integer number of Monte-Carlo Jacobians to average
#'   when computing delta-method standard errors. Defaults to
#'   one per 100 bootstrap resamples, with a minimum of 1.
#'
#' @param mc.fn.args Additional arguments to MC-PLS algorithm, mainly for controlling
#'   the step size.
#'
#' @param mc.rescov How residual covariances are treated in MC-PLS. One of
#'   \code{"auto"} (the default), \code{"reduced"}, or \code{"full"}. In
#'   \code{"reduced"} mode residual covariances are not treated as free
#'   parameters; they are identified from the rest of the structural model and
#'   the disturbances are simulated independently. In \code{"full"} mode the
#'   endogenous residual covariances (\code{eta ~~ eta} and \code{xi ~~ eta})
#'   are free parameters, explicitly simulated. \code{"auto"} uses \code{"full"}
#'   when the structural model is estimated by GLS (i.e. the model contains
#'   residual covariances) and \code{"reduced"} otherwise.
#'
#' @param verbose Should verbose output be printed?
#'
#' @param boot.optimize Logical; if \code{TRUE} and \code{bootstrap = TRUE}, applies
#'   the settings in \code{mc.boot.control} inside each bootstrap replicate (MC-PLS only).
#'   In general it will lead to slightly larger and less accurate standard errors.
#'
#' @param boot.drop.inadmissible Logical; if \code{TRUE} and \code{bootstrap = TRUE},
#'   bootstrap replicates that yield an inadmissible solution (e.g. a Heywood case)
#'   are discarded before computing standard errors, just like replicates where the
#'   estimation procedure fails outright. Defaults to \code{FALSE}, which keeps
#'   inadmissible replicates. In either case the number of inadmissible replicates
#'   is reported. Note that dropping inadmissible solutions conditions the bootstrap
#'   distribution on well-behaved resamples and may bias the standard errors downward.
#'
#' @param mc.boot.control List of control parameters passed to the MC-PLS algorithm
#'   inside each bootstrap replicate when \code{boot.optimize = TRUE}.
#'
#' @param reliabilities Optional named numeric vector of user-supplied reliabilities
#'   used for the PLSc consistency correction.
#'
#' @param default.path.estimator Character string selecting the estimator used for
#'   the structural (path) model when the model does not require Generalized Least
#'   Squares (GLS). The default \code{"ols"} uses Ordinary Least Squares whenever
#'   possible, falling back to GLS automatically when the model contains residual
#'   covariances. Setting \code{default.path.estimator = "gls"}
#'   forces GLS estimation of the structural model even when OLS would otherwise be
#'   used.
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
                mc.fast.lmer = mcpls,
                probit = NULL,
                tolerance = 1e-5,
                max.iter.0_5 = 100L,
                boot.ncores = 1L,
                boot.ncpus = NULL,
                boot.parallel = c("no", "multicore", "multisession", "snow"),
                boot.R = 500L,
                boot.iseed = NULL,
                sample = NULL,
                mc.min.iter = 50L,
                mc.max.iter = 500L,
                mc.reps = 20000L,
                mc.fixed.seed = FALSE,
                mc.polyak.juditsky = TRUE,
                mc.pj.extrapolate = TRUE,
                mc.tol = if (mc.polyak.juditsky) 0.0001 else 0.001,
                mc.delta.se = TRUE,
                mc.delta.jacobian.k = max(floor(boot.R / 100L), 1),
                mc.fn.args = list(),
                mc.rescov = c("auto", "reduced", "full"),
                verbose = interactive(),
                boot.optimize = TRUE,
                boot.drop.inadmissible = FALSE,
                mc.boot.control = list(
                  min.iter        = mc.min.iter,
                  max.iter        = mc.max.iter,
                  mc.reps         = floor(0.5 * mc.reps), # increase variance
                  tol             = mc.tol,
                  polyak.juditsky = mc.polyak.juditsky,
                  pj.extrapolate  = FALSE,                # decrease variance
                  verbose         = FALSE,
                  fixed.seed      = TRUE,                 # increase variance
                  reuse.p.start   = TRUE
                ),
                reliabilities = NULL,
                default.path.estimator = c("ols", "gls"),
                ...) {

  missing       <- match.arg(tolower(missing), c("listwise", "mean", "knn"))
  boot.parallel <- match.arg(tolower(boot.parallel), c("no", "multicore", "multisession", "snow"))
  default.path.estimator <- match.arg(tolower(default.path.estimator), c("ols", "gls"))

  if (!is.null(boot.ncpus)) {
    pls_msg_warn("The `boot.ncpus` argument is deprecated; please use `boot.ncores` instead.")
    boot.ncores <- boot.ncpus
  }

  if (!is.null(sample)) {
    pls_msg_warn("The sample argument is deprecated, please use the boot.R argument instead!")
    boot.R <- sample
  }

  data <- asDataFrame(data)

  model <- specifyModel(
    syntax                 = syntax,
    data                   = data,
    consistent             = consistent,
    missing                = missing,
    standardize            = standardize,
    ordered                = ordered,
    probit                 = probit,
    mcpls                  = mcpls,
    mc.fast.lmer           = mc.fast.lmer,
    tolerance              = tolerance,
    max.iter.0_5           = max.iter.0_5,
    mc.min.iter            = mc.min.iter,
    mc.max.iter            = mc.max.iter,
    mc.reps                = mc.reps,
    mc.tol                 = mc.tol,
    mc.fixed.seed          = mc.fixed.seed,
    mc.polyak.juditsky     = mc.polyak.juditsky,
    mc.pj.extrapolate      = mc.pj.extrapolate,
    mc.delta.se            = mc.delta.se,
    mc.delta.jacobian.k    = mc.delta.jacobian.k,
    mc.fn.args             = mc.fn.args,
    mc.rescov              = match.arg(mc.rescov, c("auto", "reduced", "full")),
    verbose                = verbose,
    bootstrap              = bootstrap,
    boot.ncores            = boot.ncores,
    boot.parallel          = boot.parallel,
    boot.R                 = boot.R,
    boot.iseed             = boot.iseed,
    boot.optimize          = boot.optimize,
    boot.drop.inadmissible = boot.drop.inadmissible,
    mc.boot.control        = mc.boot.control,
    knn.k                  = knn.k,
    reliabilities          = reliabilities,
    default.path.estimator = default.path.estimator,
    ...
  )

  model <- estimatePLS(model = model)

  if (isTRUE(modelInfo(model)$boot$bootstrap)) {
    boot <- tryCatch(
      bootstrap(model),
      error = \(e) {
        pls_msg_warn(paste0("Bootstrapping FAILED!\nMessage: ", conditionMessage(e)))
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
    updateParamVector() |>
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
    pls_msg_stop(paste0("Model estimation FAILED!\nMessage: ", conditionMessage(e)))
  })
}
