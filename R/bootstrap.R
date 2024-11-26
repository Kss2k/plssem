bootstrap <- function(model, n = 50L) {
  data <- model$data 
  results <- vector("list", n)
  for (i in seq_len(n)) {
    sampleData <- data[sample(nrow(data), nrow(data), replace = TRUE), ] 
    model$matrices$S <- cov(as.data.frame(sampleData))
    model <- estimatePLS(model)
    model$fit <- getFit(model)
    model$params$values <- extractCoefs(model)
    results[[i]] <- as.data.frame(model$params)
  }
  resultsDf <- purrr::list_rbind(results) 
  paramNames <- model$params$names
  se <- vapply(paramNames, FUN.VALUE = numeric(1L), 
               FUN = function(t) sd(resultsDf$values[resultsDf$names == t]))
  data.frame(names = paramNames, se = se)
}
