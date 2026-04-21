meanImputeMissing <- function(data, ordered = NULL) {
  # Mean for continous variables
  # Mode for nominal variables
  # Median for ordinal variables

  X <- as.data.frame(data)
  complete <- complete.cases(X)

  if (all(complete))
    return(X)

  if (length(ordered)) {
    ncat <- apply(X[,ordered,drop=FALSE], MARGIN = 2L, FUN = \(x) length(uniqueComplete(x)))
    nominal <- ordered[ncat == 2]
    ordered <- setdiff(ordered, nominal)
  } else {
    nominal <- NULL
  }

  variables <- colnames(X)

  for (variable in variables) {
    if      (variable %in% nominal) .E <- calcMode 
    else if (variable %in% ordered) .E <- strictMedian
    else                            .E <- mean

    x <- X[[variable]]
    missing <- is.na(x)

    if (any(missing))
      X[[variable]][missing] <- .E(x[!missing])
  }

  X
}


strictMedian <- function(x) {
  # Median which garuantees that we don't use an average
  # Given that we have no missing values
  n <- length(x)
  if (n %% 2 == 0 && n > 1) median(x[-n]) else median(x)
}
