meanImputeMissing <- function(data, ordered = NULL) {
  # Mean for continous variables
  # Mode for nominal variables
  # Median for ordinal variables

  X <- as.data.frame(data)
  complete <- stats::complete.cases(X)

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
  # Median that always returns an observed value (no averaging).
  #
  # For even-length inputs, this returns the *lower* of the two middle values
  # (after sorting), which is well-defined for ordinal/integer-coded data.
  x <- x[!is.na(x)]
  n <- length(x)
  if (!n) return(NA)

  x <- sort(x)

  if (n %% 2L == 1L)
    return(x[(n + 1L) / 2L])

  x[n / 2L]
}
