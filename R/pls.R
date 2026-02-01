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
 

estimatePLS_Step6  <- function(model, max.iter = 100) {
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
      model$fit.u <- list(NULL)
      model$fit.u <- getFitPLSModel(model, consistent = FALSE)
      model$fit   <- model$fit.c
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


estimatePLS_Step7 <- function(model) {
  model$params$values <- extractCoefs(model)
  model$factorScores  <- getFactorScores(model)

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


estimatePLS_Step8 <- function(model) {
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
    estimatePLS_Step6(max.iter = max.iter) |>
    estimatePLS_Step7() |>
    estimatePLS_Step8()

}
