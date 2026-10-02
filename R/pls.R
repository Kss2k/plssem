SE_NON_LINEAR_PROBIT_CORR_MAT <- FALSE


#' Fit Partial Least Squares Structural Equation Models
#'
#' \code{pls()} estimates Partial Least Squares Structural Equation Models (PLS-SEM)
#' and their consistent (PLSc) variants. The function accepts \code{lavaan}-style
#' syntax, handles ordered indicators through polychoric correlations and probit
#' factor scores, and estimates two-level (multilevel) models, with separate
#' models for the within (\code{level: 1}) and between (\code{level: 2}) levels.
#'
#' @param syntax Character string with \code{lavaan}-style model syntax describing
#'   both measurement (\code{=~}) and structural (\code{~}) relations. Two-level
#'   models are specified with \code{level: 1} and \code{level: 2} blocks (see
#'   \code{mlm}), where random slopes are specified with \code{rv()} modifiers at
#'   level 1 (e.g., \code{fw ~ rv(s1)*x1}). The random slopes (e.g., \code{s1})
#'   are variables at level 2.
#'
#' @param data A \code{data.frame} or coercible object containing the manifest
#'   indicators referenced in \code{syntax}. Ordered factors are automatically
#'   detected, but can also be supplied explicitly through \code{ordered}.
#'
#' @param standardize Logical; if \code{TRUE}, indicators are standardized before
#'   estimation so that factor scores have comparable scales.
#'
#' @param consistent Logical; \code{TRUE} requests PLSc corrections, whereas \code{FALSE}
#'   fits the traditional PLS model. If \code{NULL} (default), \code{FALSE} is
#'   used for MC-PLS (including two-level models) and for models using \code{lmer}
#'   as the path estimator, and \code{TRUE} otherwise.
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
#'   all indicators missing are removed prior to imputation. Rows with missing
#'   \code{cluster} values are always removed.
#'
#' @param knn.k Integer specifying the number of neighbors (\code{k}) used when
#'   \code{missing = "kNN"}.
#'
#' @param mcpls Should the model be estimated using the Monte-Carlo Consistent
#'   Partial Least Squares (MC-PLSc) algorithm?
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
#' @param mc.small.sample Logical; if \code{TRUE}, average MC-PLS
#'   estimating equations over simulated samples of the observed sample size.
#'   If \code{FALSE} (default), fit one sample using all simulated observations.
#'   This assumes that the auxiliary estimator (i.e., traditional PLS) has
#'   negligible finite sample bias. I.e., that the bias of the estimator is
#'   not affected by the sample size.
#'
#' @param mc.small.sample.max.k Maximum number of simulated samples to average
#'   when \code{mc.small.sample = TRUE}. Defaults to 100. The number of samples
#'   is also limited by \code{mc.reps}, rounded down to a multiple of the
#'   observed sample size, with at least one sample.
#'
#' @param mc.small.sample.point.estimate Which point estimate of the simulated
#'   auxiliary parameters the root equation matches to the observed ones, when
#'   \code{mc.small.sample = TRUE}? \code{"mean"} solves
#'   \eqn{E[\theta^{*}|\theta] = \hat{\theta}^{*}}. \code{"median"} (the default)
#'   solves \eqn{median[\theta^{*}|\theta] = \hat{\theta}^{*}}.
#'
#'   The median commutes with the (monotone) binding function where the mean
#'   does not, so median-matching targets a median-unbiased estimator. This
#'   removes the finite-sample bias. which can be introduced by the curvature
#    of the inverse binding function. This is most pronounced for small sample
#'   size models, and when the indicators are uninformative (small loadings,
#'   few categories, strongly assymetric thresholds). It makes the estimating
#'   function somewhat noisier for a given number of simulated samples.
#'   Ignored when \code{mc.small.sample = FALSE}.
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
#' @param mc.diag.secant Logical; if \code{TRUE}, the MC-PLS root-finding algorithm
#'   uses a per-coordinate diagonal-secant step (estimating each coordinate's local
#'   slope from consecutive iterates), instead of the standard Robbins-Monro step.
#'   This might be usefull if the binding funtion is non-monotone.
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
#'   Squares (GLS), or Linear Mixed-Effects Regression (LMER).
#'   The default \code{"ols"} uses Ordinary Least Squares whenever
#'   possible, falling back to GLS automatically when the model contains residual
#'   covariances. Setting \code{default.path.estimator = "gls"}
#'   forces GLS estimation of the structural model even when OLS would otherwise be
#'   used.
#'
#' @param cluster Optional character vector naming the cluster variable(s) in
#'   \code{data}. The cluster variables are stored alongside the (standardized)
#'   indicators, and bootstrapping resamples whole clusters instead of rows.
#'   Required for two-level models (a single cluster variable).
#'
#' @param inner.weights Character string selecting the inner weighting scheme
#'   used with \code{approach.weights = "pls"}. One of: \code{"centroid"},
#'   \code{"factorial"}, or \code{"path"}. Defaults to \code{"path"}.
#'
#' @param approach.weights Character string selecting the approach used to
#'   estimate the outer weights. \code{"pls"} (default) uses the PLS algorithm,
#'   with the inner weighting scheme given by \code{inner.weights}. \code{"pca"}
#'   forms the weights from each construct's own indicators only, ignoring the
#'   structural.
#'
#' @param mlm Should the model be estimated as a two-level (multilevel) model?
#'   If \code{NULL} (default), this is detected from \code{syntax}, i.e., whether
#'   it has \code{level: 1} and \code{level: 2} blocks. Two-level models are
#'   estimated using an extension of the MC-PLSc estimator.
#'   The \code{mc.*} arguments (e.g., \code{mc.max.iter}, \code{mc.reps}) are used when they are
#'   specified, whereas arguments like \code{standardize}
#'   are (currently) not supported. Standard errors (\code{bootstrap = TRUE})
#'   are computed using the delta method.
#'
#' @param level2.cov Two-level models only. How the (co-)variances of the
#'   variables at level 2 are computed. \code{"means"} uses the
#'   covariances of the cluster means. \code{"muml"} (default) uses Muthen's (1994)
#'   estimator of the between-cluster covariance matrix, which corrects for the
#'   within-cluster variation in the cluster means.
#'
#' @param level2.approach.weights Two-level models only. The approach used to
#'   estimate the outer weights of the level 2 model (\code{approach.weights}
#'   applies to the level 1 model). Defaults to \code{"pca"}.
#'
#' @param ... Internal arguments. For advanced users only.
#'
#' @return A \code{PlsModel} object containing the estimated parameters, fit measures,
#'   factor scores, and any bootstrap results (a \code{PlsMultilevelModel} object
#'   for two-level models). Methods such as \code{summary()}, \code{coef()}, and
#'   \code{parameter_estimates()} can be applied to inspect the fit.
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
                consistent = NULL,
                bootstrap = FALSE,
                ordered = NULL,
                missing = c("listwise", "mean", "kNN"),
                knn.k = 5,
                mcpls = NULL,
                probit = NULL,
                tolerance = 1e-5,
                max.iter.0_5 = 500L,
                boot.ncores = 1L,
                boot.ncpus = NULL,
                boot.parallel = c("no", "multicore", "multisession", "snow"),
                boot.R = 500L,
                boot.iseed = NULL,
                sample = NULL,
                mc.min.iter = 50L,
                mc.max.iter = 1000L,
                mc.reps = 20000L,
                mc.fixed.seed = FALSE,
                mc.polyak.juditsky = TRUE,
                mc.pj.extrapolate = TRUE,
                mc.tol = if (mc.polyak.juditsky) 0.0001 else 0.001,
                mc.delta.se = TRUE,
                mc.delta.jacobian.k = max(floor(boot.R / 100L), 1),
                mc.fn.args = list(),
                mc.rescov = c("auto", "reduced", "full"),
                mc.diag.secant = FALSE,
                mc.small.sample = FALSE,
                mc.small.sample.max.k = 100L,
                mc.small.sample.point.estimate = c("median", "mean"),
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
                default.path.estimator = c("ols", "gls", "lmer"),
                cluster = NULL,
                inner.weights = c("path", "centroid", "factorial"),
                approach.weights = c("pls", "pca"),
                mlm = NULL,
                level2.cov = c("muml", "means"),
                level2.approach.weights = c("pca", "pls"),
                ...) {

  missing       <- match.arg(tolower(missing), c("listwise", "mean", "knn"))
  boot.parallel <- match.arg(tolower(boot.parallel), c("no", "multicore", "multisession", "snow"))
  default.path.estimator <- match.arg(tolower(default.path.estimator), c("ols", "gls", "lmer"))
  inner.weights    <- match.arg(tolower(inner.weights), c("path", "centroid", "factorial"))
  approach.weights <- match.arg(tolower(approach.weights), c("pls", "pca"))
  level2.cov       <- match.arg(tolower(level2.cov), c("muml", "means"))
  level2.approach.weights <- match.arg(tolower(level2.approach.weights), c("pca", "pls"))
  mc.small.sample.point.estimate <- match.arg(tolower(mc.small.sample.point.estimate), c("median", "mean"))

  if (!is.null(boot.ncpus)) {
    pls_msg_warn("The `boot.ncpus` argument is deprecated; please use `boot.ncores` instead.")
    boot.ncores <- boot.ncpus
  }

  if (!is.null(sample)) {
    pls_msg_warn("The sample argument is deprecated, please use the boot.R argument instead!")
    boot.R <- sample
  }

  # Two-level (multilevel) models are estimated by `mpls()`
  is.mlm.syntax <- isMultilevelSyntax(syntax)
  if (is.null(mlm)) mlm <- is.mlm.syntax

  pls_stopif(isTRUE(mlm) && !is.mlm.syntax,
    "`mlm = TRUE` requires a two-level model syntax, with `level: 1` and",
    "`level: 2` blocks!"
  )

  pls_stopif(!isTRUE(mlm) && is.mlm.syntax,
    "The model syntax has `level:` blocks, which requires `mlm = TRUE`",
    "(or `mlm = NULL`)!"
  )

  if (isTRUE(mlm)) {
    pls_stopif(is.null(cluster),
      "`cluster` must be specified for two-level (multilevel) models!"
    )

    supplied <- names(as.list(match.call()))[-1L]

    unsupported <- intersect(supplied, MLM_UNSUPPORTED_ARGS)
    pls_warnif(length(unsupported),
      "The following arguments are (currently) ignored for two-level models:",
      paste0("`", unsupported, "`", collapse = ", ")
    )

    boot.parallel <- if (boot.parallel == "snow") "multisession" else boot.parallel

    args <- list(
      syntax        = syntax,
      data          = data,
      cluster       = cluster,
      ordered       = ordered,
      verbose       = verbose,
      bootstrap     = bootstrap,
      boot.R        = boot.R,
      boot.parallel = boot.parallel,
      boot.ncores   = boot.ncores,
      boot.iseed    = boot.iseed,
      level2.cov    = level2.cov,
      level2.approach.weights = level2.approach.weights
    )

    for (arg in intersect(supplied, names(MLM_PLS_ARGS)))
      args[[MLM_PLS_ARGS[[arg]]]] <- get(arg)

    if (isTRUE(mc.fixed.seed))
      args$rng.seed <- floor(stats::runif(1L, min = 0, max = 9999999))

    dots <- list(...) # e.g., `level2.cov`, or `rng.seed` (arguments of `mpls()`)
    args[names(dots)] <- dots

    return(do.call("mpls", args)) # by name, for the header of messages (see `pls_msg()`)
  }

  data <- asDataFrame(data)

  model <- specifyModel(
    syntax                         = syntax,
    data                           = data,
    consistent                     = consistent,
    missing                        = missing,
    standardize                    = standardize,
    ordered                        = ordered,
    probit                         = probit,
    mcpls                          = mcpls,
    tolerance                      = tolerance,
    max.iter.0_5                   = max.iter.0_5,
    mc.min.iter                    = mc.min.iter,
    mc.max.iter                    = mc.max.iter,
    mc.reps                        = mc.reps,
    mc.tol                         = mc.tol,
    mc.fixed.seed                  = mc.fixed.seed,
    mc.polyak.juditsky             = mc.polyak.juditsky,
    mc.pj.extrapolate              = mc.pj.extrapolate,
    mc.delta.se                    = mc.delta.se,
    mc.delta.jacobian.k            = mc.delta.jacobian.k,
    mc.fn.args                     = mc.fn.args,
    mc.rescov                      = match.arg(mc.rescov, c("auto", "reduced", "full")),
    mc.diag.secant                 = mc.diag.secant,
    mc.small.sample                = mc.small.sample,
    mc.small.sample.max.k          = mc.small.sample.max.k,
    mc.small.sample.point.estimate = mc.small.sample.point.estimate,
    verbose                        = verbose,
    bootstrap                      = bootstrap,
    boot.ncores                    = boot.ncores,
    boot.parallel                  = boot.parallel,
    boot.R                         = boot.R,
    boot.iseed                     = boot.iseed,
    boot.optimize                  = boot.optimize,
    boot.drop.inadmissible         = boot.drop.inadmissible,
    mc.boot.control                = mc.boot.control,
    knn.k                          = knn.k,
    reliabilities                  = reliabilities,
    default.path.estimator         = default.path.estimator,
    cluster                        = cluster,
    inner.weights                  = inner.weights,
    approach.weights               = approach.weights,
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
  if (model@status$quick)
    return(model)

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
