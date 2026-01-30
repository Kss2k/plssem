#' @export
pls <- function(syntax,
                data,
                standardize = TRUE,
                max.iter = 100, 
                consistent = TRUE, 
                bootstrap = FALSE,
                sample = 50,
                lme4.syntax = NULL,
                ...) {
  # preprocess data
  if (!is.matrix(data))
    data <- as.matrix(data)

  if (standardize)
    data <- standardizeMatrix(data)
  
  # Define model
  model <- specifyModel(syntax, data) 

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

  if (!is.null(lme4.syntax)) {
    warning("lme4.syntax argument is not implemented yet!")

    if (bootstrap) {
      warning("lme4.syntax with bootstrap argument is not implemented yet!")
    }
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
