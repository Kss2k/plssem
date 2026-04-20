USE_NON_LINEAR_PROBIT_CORR_MAT <- FALSE # for now we stick with the linear assumption


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
#' @param sample DEPRECTATED. Integer giving the number of bootstrap resamples drawn when
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
#'
#' @param mc.polyak.juditsky Should the polyak.juditsky running average method
#'   be applied in the MC-PLS algorithm?
#'
#' @param mc.fn.args Additional arguments to MC-PLS algorithm, mainly for controling
#'   the step size.
#'
#' @param verbose Should verbose output be printed?
#'
#' @param boot.optimize Logical; if \code{TRUE} and \code{bootstrap = TRUE}, applies
#'   the settings in \code{mc.boot.control} inside each bootstrap replicate (MC-PLS only).
#'
#' @param missing How to handle missing values in the observed indicators.
#'   \code{"listwise"} removes any rows with missing values in the indicators
#'   (the behavior prior to introducing this argument). \code{"pairwise"} keeps
#'   partially observed rows and uses pairwise complete observations when
#'   computing covariances/correlations and when standardizing, which can be
#'   useful when missingness is limited but will generally change results.
#'
#' @param mc.boot.control List of control parameters passed to the MC-PLS algorithm
#'   inside each bootstrap replicate when \code{boot.optimize = TRUE}. This can be used
#'   to speed up bootstrapping by, for example, increasing the tolerance or reducing
#'   \code{mc.reps}. The element \code{reuse.p.start} controls whether to reuse the
#'   original \code{p.start} for the bootstrap replicates.
#'
#' @param ... Currently unused, reserved for future extensions.
#'
#' @return An object of class \code{plssem} containing the estimated parameters, fit
#'   measures, factor scores, and any bootstrap results. Methods such as
#'   \code{summary()}, \code{print()}, and \code{coef()} can be applied to inspect the fit.
#'
#' @seealso [summary.plssem()], [print.plssem()]
#'
#' @examples
#' # Linear Model with Continuous Data
#' \donttest{
#'
#' library(plssem)
#' library(modsem)
#'
#' tpb <- '
#' # Outer Model (Based on Hagger et al., 2007)
#'   ATT =~ att1 + att2 + att3 + att4 + att5
#'   SN =~ sn1 + sn2
#'   PBC =~ pbc1 + pbc2 + pbc3
#'   INT =~ int1 + int2 + int3
#'   BEH =~ b1 + b2
#'
#' # Inner Model (Based on Steinmetz et al., 2011)
#'   INT ~ ATT + SN + PBC
#'   BEH ~ INT + PBC
#' '
#'
#' fit <- pls(tpb, TPB, bootstrap = TRUE)
#' summary(fit)
#'
#' # Linear Model with Ordered Data
#' tpb <- '
#' # Outer Model (Based on Hagger et al., 2007)
#'   ATT =~ att1 + att2 + att3 + att4 + att5
#'   SN =~ sn1 + sn2
#'   PBC =~ pbc1 + pbc2 + pbc3
#'   INT =~ int1 + int2 + int3
#'   BEH =~ b1 + b2
#'
#' # Inner Model (Based on Steinmetz et al., 2011)
#'   INT ~ ATT + SN + PBC
#'   BEH ~ INT + PBC
#' '
#'
#' fit <- pls(tpb, TPB_Ordered, bootstrap = TRUE)
#' summary(fit)
#'
#' # Multilevel Random Slopes Model with Continuous Data
#' syntax <- "
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#'   W =~ w1 + w2 + w3
#'   Y ~ X + Z + (1 + X + Z | cluster)
#'   W ~ X + Z + (1 + X + Z | cluster)
#' "
#'
#' fit <- pls(syntax, data = randomSlopes, bootstrap = TRUE)
#' summary(fit)
#'
#' # Multilevel Random Slopes Model with Ordered Data
#' syntax <- "
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#'   W =~ w1 + w2 + w3
#'   Y ~ X + Z + (1 + X + Z | cluster)
#'   W ~ X + Z + (1 + X + Z | cluster)
#' "
#'
#' fit <- pls(syntax, data = randomSlopesOrdered, bootstrap = TRUE)
#' summary(fit)
#'
#' # Multilevel Random Intercepts Model with Continuous Data
#' syntax <- '
#'   f =~ y1 + y2 + y3
#'   f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)
#' '
#'
#' fit <- pls(syntax, data = randomIntercepts, bootstrap = TRUE)
#' summary(fit)
#'
#' # Multilevel Random Intercepts Model with Ordered Data
#' syntax <- '
#'   f =~ y1 + y2 + y3
#'   f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)
#' '
#'
#' fit <- pls(syntax, data = randomInterceptsOrdered, bootstrap = TRUE)
#' summary(fit)
#'
#' # Interaction Model with Continuous Data
#' m <- '
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#'
#'   Y ~ X + Z + X:Z
#' '
#'
#' fit <- pls(m, modsem::oneInt, bootstrap = TRUE)
#' summary(fit)
#'
#' # Interaction Model with Ordered Data
#' m <- '
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#'
#'   Y ~ X + Z + X:Z
#' '
#'
#' fit <- pls(m, oneIntOrdered, bootstrap = TRUE)
#' summary(fit)
#'
#' }
#' @export
pls <- function(syntax,
                data,
                standardize = TRUE,
                consistent = TRUE,
                bootstrap = FALSE,
                ordered = NULL,
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
                missing = c("listwise", "pairwise"),
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
                ...) {

  boot.parallel <-  match.arg(tolower(boot.parallel), c("no", "multicore", "snow"))
  missing <- match.arg(tolower(missing), c("listwise", "pairwise"))

  if (!is.null(sample)) {
    warning("The sample argument is deprecated, please use the boot.R argument instead!")
    boot.R <- sample
  }

  # preprocess data
  data <- as.data.frame(data)

  # Define model
  model <- specifyModel(
    syntax             = syntax,
    data               = data,
    consistent         = consistent,
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
    missing            = missing
  )

  # Fit model
  model <- estimatePLS(model = model)

  # Bootstrap
  if (model$info$boot$bootstrap) {
    model$boot <- bootstrap(model)
    model$params$se <- model$boot$se
  }

  model$parTable <- getParTableEstimates(model)
  class(model) <- "plssem"
  model
}


resetPLS_Model <- function(model, hard.reset = FALSE) {
  model$status$finished       <- FALSE
  model$status$convergence    <- FALSE
  model$status$iterations.0_5 <- 0L

  if (hard.reset) {
    model$status$iterations     <- 0L
    model$status$iterations.0_9 <- 0L
    model$params$values.old     <- NULL
  }

  model
}


estimatePLS_Inner <- function(model) {
  model <- resetPLS_Model(model, hard.reset = TRUE)

  resetPLS_Model(model) |>
    estimatePLS_Step0_5() |>
    estimatePLS_Step6() |>
    estimatePLS_Step7() |>
    estimatePLS_Step8()
}


estimatePLS_Outer <- function(model, ...) {
  if (model$info$is.mcpls)
    return(mcpls(model, ...))

  model
}


estimatePLS <- function(model, ...) {
  model |>
    estimatePLS_Inner() |>
    estimatePLS_Outer(...)
}
