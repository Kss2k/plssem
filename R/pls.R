USE_NON_LINEAR_PROBIT_CORR_MAT <- FALSE # for now we stick with the linear assumption


#' Fit Partial Least Squares Structural Equation Models
#'
#' `pls()` estimates Partial Least Squares Structural Equation Models (PLS-SEM)
#' and their consistent (PLSc) variants. The function accepts `lavaan`-style
#' syntax, handles ordered indicators through polychoric correlations and probit
#' factor scores, and supports multilevel specifications expressed with
#' `lme4`-style random effects terms inside the structural model.
#'
#' @param syntax Character string with `lavaan`-style model syntax describing
#'   both measurement (`=~`) and structural (`~`) relations. Random effects are
#'   specified with `(term | cluster)` statements.
#'
#' @param data A `data.frame` or coercible object containing the manifest
#'   indicators referenced in `syntax`. Ordered factors are automatically
#'   detected, but can also be supplied explicitly through `ordered`.
#'
#' @param standardize Logical; if `TRUE`, indicators are standardized before
#'   estimation so that factor scores have comparable scales.
#'
#' @param max.iter.0_5 Maximum number of PLS iterations performed when estimating
#'   the measurement and structural models.
#'
#' @param consistent Logical; `TRUE` requests PLSc corrections, whereas `FALSE`
#'   fits the traditional PLS model.
#'
#' @param bootstrap Logical; if `TRUE`, nonparametric bootstrap standard errors
#'   are computed with `sample` resamples.
#'
#' @param sample Integer giving the number of bootstrap resamples drawn when
#'   `bootstrap = TRUE`.
#'
#' @param ordered Optional character vector naming manifest indicators that
#'   should be treated as ordered when computing polychoric correlations.
#'
#' @param probit Logical; overrides the automatic choice of probit factor scores
#'   that is based on whether ordered indicators are present.
#'
#' @param mcpls Should a Monte-Carlo consistency correction be applied?
#'
#' @param tolerance Numeric; Convergence criteria/tolerance.
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
#' @param ... Currently unused, reserved for future extensions.
#'
#' @return An object of class `plssem` containing the estimated parameters, fit
#'   measures, factor scores, and any bootstrap results. Methods such as
#'   `summary()`, `print()`, and `coef()` can be applied to inspect the fit.
#'
#' @seealso [summary.plssem()], [print.plssem()]
#'
#' @examples
#' # Linear Model with Continuous Data
#' \dontrun{
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
                sample = 50L,
                ordered = NULL,
                mcpls = NULL,
                probit = NULL,
                tolerance = 1e-5,
                max.iter.0_5 = 100L,
                mc.min.iter = 5L,
                mc.max.iter = 250L,
                mc.reps = 20000L,
                mc.tol = 1e-3,
                mc.fixed.seed = FALSE,
                mc.polyak.juditsky = FALSE,
                mc.fn.args = list(),
                verbose = interactive(),
                ...) {
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
    verbose            = verbose
  )

  # Fit model
  model <- estimatePLS(model = model)

  # Bootstrap
  if (bootstrap) {
    model$boot <- bootstrap(model, R = sample)
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
