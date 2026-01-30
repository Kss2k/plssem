bootstrap <- function(model, n = 50L, zero.tol = 1e-10) {
  data <- model$data 
  results <- vector("list", n)

  for (i in seq_len(n)) {
    sampleData <- data[sample(nrow(data), nrow(data), replace = TRUE), ] 
    model$matrices$S <- stats::cov(as.data.frame(sampleData))
    model <- estimatePLS(model)
    model$fit <- getFit(model)
    model$params$values <- extractCoefs(model)
    results[[i]] <- as.data.frame(model$params)
  }

  resultsDf <- purrr::list_rbind(results) 
  paramNames <- model$params$names

  se <- vapply(
    paramNames,
    FUN.VALUE = numeric(1L), 
    FUN = \(t) stats::sd(resultsDf$values[resultsDf$names == t])
  )

  se[se <= zero.tol] <- NA_real_

  data.frame(names = paramNames, se = se)
}
