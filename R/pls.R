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
#' @param max.iter.0_9 Maximum number of PLS iterations performed when estimating
#'   the measurement and structural models. Only relevant for interaction
#'   models with ordered data.
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
#' @param max.iter Maximum number of PLS iterations performed when estimating
#'   the measurement and structural models.
#'
#' @param ordered Optional character vector naming manifest indicators that
#'   should be treated as ordered when computing polychoric correlations.
#'
#' @param probit Logical; overrides the automatic choice of probit factor scores
#'   that is based on whether ordered indicators are present.
#'
#' @param consistent.probit Logical; Should probit consistent estimates be
#'   calculated?
#'
#' @param tolerance Numeric; Convergence criteria/tolerance.
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
#' fit <- pls(tpb, data = modsem::TPB)
#' summary(fit)
#'
#' syntax <- "
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#'   W =~ w1 + w2 + w3
#'   Y ~ X + Z + (1 + X + Z | cluster)
#'   W ~ X + Z + (1 + X + Z | cluster)
#' "
#'
#' fit <- pls(syntax, data = randomSlopesOrdered)
#' summary(fit)
#'
#' @export
pls <- function(syntax,
                data,
                standardize = TRUE,
                consistent = TRUE, 
                bootstrap = FALSE,
                sample = 50L,
                ordered = NULL,
                probit = NULL,
                consistent.probit = TRUE,
                tolerance = 1e-5,
                max.iter.0_5 = 100L,
                max.iter.0_9 = 50L,
                ...) {
  # preprocess data
  data <- as.data.frame(data)

  # Define model
  model <- specifyModel(
    syntax            = syntax,
    data              = data,
    consistent        = consistent,
    standardize       = standardize,
    ordered           = ordered,
    probit            = probit,
    consistent.probit = consistent.probit,
    tolerance         = tolerance,
    max.iter.0_5      = max.iter.0_5,
    max.iter.0_9      = max.iter.0_9
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
 

estimatePLS <- function(model, max.iter = 100) {
  model <- resetPLS_Model(model, hard.reset = TRUE)
  
  while (!model$status$finished) {
    model <- 
      resetPLS_Model(model) |>
      estimatePLS_Step0_5() |>
      estimatePLS_Step6() |>
      estimatePLS_Step7() |>
      estimatePLS_Step8() |>
      estimatePLS_Step9()
  }

  model
}
