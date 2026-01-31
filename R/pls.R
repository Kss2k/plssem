#' @export
pls <- function(syntax,
                data,
                standardize = TRUE,
                max.iter = 100, 
                consistent = TRUE, 
                bootstrap = FALSE,
                sample = 50,
                lme4.syntax = NULL,
                cluster = NULL,
                ...) {
  # preprocess data
  if (!is.matrix(data))
    data <- as.matrix(data)

  if (standardize)
    data <- standardizeMatrix(data, cluster = cluster)
  
  # Define model
  model <- specifyModel(
    syntax      = syntax,
    data        = data,
    consistent  = consistent,
    cluster     = cluster,
    lme4.syntax = lme4.syntax
  ) 

  # Fit model
  model <- estimatePLS(
    model       = model,
    max.iter    = max.iter
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


estimatePLS <- function(model,
                        max.iter = 100) {
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
  
  # Get Final fit 
  model$fit.c <- getFitPLSModel(model, consistent = TRUE)
  model$fit.u <- getFitPLSModel(model, consistent = FALSE)
  model$fit   <- if (consistent) model$fit.c else model$fit.u

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

  model$info$iterations <- i

  model
}
