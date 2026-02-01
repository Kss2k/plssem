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
#' @param data A `data.frame` or coercible object containing the manifest
#'   indicators referenced in `syntax`. Ordered factors are automatically
#'   detected, but can also be supplied explicitly through `ordered`.
#' @param standardize Logical; if `TRUE`, indicators are standardized before
#'   estimation so that factor scores have comparable scales.
#' @param max.iter Maximum number of PLS iterations performed when estimating
#'   the measurement and structural models.
#' @param consistent Logical; `TRUE` requests PLSc corrections, whereas `FALSE`
#'   fits the traditional PLS model.
#' @param bootstrap Logical; if `TRUE`, nonparametric bootstrap standard errors
#'   are computed with `sample` resamples.
#' @param sample Integer giving the number of bootstrap resamples drawn when
#'   `bootstrap = TRUE`.
#' @param ordered Optional character vector naming manifest indicators that
#'   should be treated as ordered when computing polychoric correlations.
#' @param probit Logical; overrides the automatic choice of probit factor scores
#'   that is based on whether ordered indicators are present.
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
                max.iter = 100, 
                consistent = TRUE, 
                bootstrap = FALSE,
                sample = 50,
                ordered = NULL,
                probit = NULL,
                ...) {
  # preprocess data
  data <- as.data.frame(data)

  # Define model
  model <- specifyModel(
    syntax      = syntax,
    data        = data,
    consistent  = consistent,
    standardize = standardize,
    ordered     = ordered,
    probit      = probit
  ) 

  # Fit model
  model <- estimatePLS(
    model    = model,
    max.iter = max.iter
  )

  # Bootstrap
  if (bootstrap) {
    model$boot <- bootstrap(model, R = sample)
    model$params$se <- model$boot$se
  }

  model$parTable <- getParTableEstimates(model)
  class(model) <- "plssem"
  model 
}


estimatePLS_Step0_5 <- function(model, max.iter = 100) {
  consistent <- model$info$consistent

  model <- step0(model)

  for (i in seq_len(max.iter)) {
    model <- model |> 
      step1() |>
      step2() |>
      step3() |>
      step4() |>
      step5() 

    if (model$info$convergence) {
      break

    } else if (i >= max.iter) {
      warning("Convergence reached. Stopping.")
      break
    }
  }
  
  model$info$iterations <- i

  model
}


estimatePLS_Step6 <- function(model) {
  model$factorScores <- getFactorScores(model)

  if (model$info$is.interaction.model) {
    # Update variance and coefficients of interaction terms
    elems <- model$info$intTermElems
    names <- names(elems)
    X     <- as.data.frame(model$factorScores)

    SC <- model$matrices$SC
    C  <- model$matrices$C
    L  <- model$matrices$lambda
    G  <- model$matrices$gamma

    D.lv <- diag(NCOL(C))
    D.ov <- diag(NCOL(SC))

    dimnames(D.lv) <- dimnames(C)
    dimnames(D.ov) <- dimnames(SC)

    for (xz in names) {
      model$factorScores[,xz] <- multiplyIndicatorsCpp(X[elems[[xz]]])
      sigma <- sd(model$factorScores[,xz,drop=TRUE])
      D.lv[xz, xz] <- D.ov[xz, xz] <- sigma
    }

    model$matrices$SC <- D.ov %*% SC %*% D.ov
    model$matrices$C  <- D.lv %*%  C %*% D.lv
  }

  model
}
 

estimatePLS_Step7  <- function(model, max.iter = 100) {
  # Step 6 get baseline fit for correcting path coefficients in
  # multilevel model.

  is.multilevel <- model$info$is.multilevel
  consistent    <- model$info$consistent
  is.probit     <- model$info$is.probit

  if (!is.multilevel) {

    if (consistent) {
      model$fit.c <- getFitPLSModel(model, consistent = TRUE)
      model$fit.u <- list(NULL)
      model$fit   <- model$fit.c
    } else {
      model$fit.c <- list(NULL)
      model$fit.u <- getFitPLSModel(model, consistent = FALSE)
      model$fit   <- model$fit.u
    }

    return(model)

  }

  model.c <- model
  model.u <- model

  # if consistent is FALSE we want to correct for not using probit
  # factor scores only. 
  if (is.probit) {
    model.u$matrices$S <- getCorrMat(model.u$data, probit = FALSE)
    model.u <- estimatePLS_Step0_5(model.u, max.iter = max.iter)
  }

  model$fit.c <- getFitPLSModel(model.c, consistent = consistent)
  model$fit.u <- getFitPLSModel(model.u, consistent = FALSE)
  model$fit   <- model$fit.c

  model
}


estimatePLS_Step8 <- function(model) {
  model$params$values <- extractCoefs(model)

  if (model$info$is.multilevel) {
    model$fit.lmer <- plslmer(model)

    coefs.x <- model$params$values
    coefs.y <- model$fit.lmer$values

    common  <- intersect(names(coefs.x), names(coefs.y))
    new     <- setdiff(names(coefs.y), names(coefs.x))

    coefs.x[common] <- coefs.y[common]
    coefs.all       <- c(coefs.x, coefs.y[new])

    model$params$values <- plssemVector(coefs.all)
    model$params$se     <- rep(NA_real_, length(coefs.all))
  }


  model
}


estimatePLS_Step9 <- function(model) {
  # Step 8. Finalize model object with any additionaly information
  ordered <- model$info$ordered

  if (length(ordered)) {
    for (ord in ordered) {
      tau <- getThresholdsFromQuantiles(X = model$data, variable = ord)
      model$params$values <- c(model$params$values, tau)
    }

    model$params$se <- rep(NA_real_, length(model$params$values))
  }

  model
}


estimatePLS <- function(model, max.iter = 100) {

  model |>
    estimatePLS_Step0_5(max.iter = max.iter) |>
    estimatePLS_Step6() |>
    estimatePLS_Step7(max.iter = max.iter) |>
    estimatePLS_Step8() |>
    estimatePLS_Step9()

}
