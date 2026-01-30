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
    syntax     = syntax,
    data       = data,
    consistent = consistent,
    cluster    = cluster
  ) 

  # Fit model
  model <- estimatePLS(model, max.iter = max.iter, standardize = standardize)

  # Get Final fit 
  model$fit.c <- getFit(model, consistent = TRUE)
  model$fit.u <- getFit(model, consistent = FALSE)
  model$fit   <- if (consistent) model$fit.c else model$fit.u

  model$params$values <- extractCoefs(model)

  # Bootstrap
  if (bootstrap) {
    model$boot <- bootstrap(model, n = sample)
    model$params$se <- model$boot$se
  }
  
  model$parTable <- getParTableEstimates(model)
  model$factorScores <- getFactorScores(model)

  if (!is.null(lme4.syntax)) {
    return(plslmer(lme4.syntax, plsModel = model, cluster = cluster,
                   consistent = consistent))
  }

  class(model) <- "plssem"
  model 
}


estimatePLS <- function(model, max.iter = 100, standardize = TRUE) {
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
