pls <- function(syntax, data, standardize = TRUE, maxIter = 100, 
                consistent = TRUE, 
                bootstrap = FALSE, sample = 50, ...) {
  # preprocess data
  if (!is.matrix(data)) data <- as.matrix(data)
  if (standardize) data <- standardizeMatrix(data)
  
  # Define model
  model <- specifyModel(syntax, data) 

  # Fit model
  model <- estimatePLS(model, maxIter = maxIter, standardize = standardize)

  # final fit 
  model$fit <- getFit(model, consistent = consistent)
  model$params$values <- extractCoefs(model)

  # Bootstrap
  if (bootstrap) {
    model$params$se <- bootstrap(model, n = sample)$se
  }
  model 
}


estimatePLS <- function(model, maxIter = 100, standardize = TRUE) {
  model <- step0(model)
  for (i in seq_len(maxIter)) {
    model <- model |> 
      step1() |>
      step2() |>
      step3() |>
      step4() |>
      step5() 
    if (model$convergence >= maxIter) {
      warning("Convergence reached. Stopping.")
      break
    }
  } 
  model 
}
