getPLS_Data <- function(data,
                        indicators,
                        cluster = NULL,
                        consistent = TRUE,
                        standardize = TRUE,
                        ordered = NULL,
                        is.probit = NULL,
                        is.cexp = NULL,
                        missing = "listwise") {

  vars <- c(indicators, cluster)
  missingVars <- !vars %in% colnames(data)

  if (any(missingVars))
    stop("Missing variables: ", paste0(vars[missingVars], collapse = ", "))

  data <- as.data.frame(data)[vars]

  missingCases <- !stats::complete.cases(data)

  if (any(missingCases) && missing == "listwise") {
    # "TODO: Add multiple imputation"
    message("Removing missing data using list wise deletion...")
    data <- data[!missingCases, , drop = FALSE]

  } else if (any(missingCases) && missing == "pairwise") {
    message("Using pairwise complete observations...")
    allMissing <- apply(data, MARGIN = 1L, FUN = \(x) all(is.na(x)))
    data <- data[!allMissing, , drop = FALSE]
  }
 
  if (standardize) {
    data <- standardizeDataFrame(
      data    = data,
      cluster = cluster
    )
  }

  S <- getCorrMat(data[indicators], probit = is.probit, ordered = ordered)
  X <- as.matrix(data[indicators])

  if (!is.null(cluster)) {
    if (!is.character(cluster))
      stop("`cluster` must be a character string, if lme4.syntax is provided!")

    attr(X, "cluster") <- data[, cluster, drop = FALSE]
  }

  list(X = X, S = S)
}


standardizeDataFrame <- function(data, cluster = NULL, subset = colnames(data)) {
  subset <- setdiff(subset, cluster)

  Z <- lapply(data[subset], FUN = standardizeAtomic)
  data[subset] <- Z

  data
}


standardizeMatrix <- function(X, use = "pairwise") {
  X <- as.matrix(X)

  out <- switch(use,
    pairwise = {
      if (NCOL(X) == 1L) matrix(standardizeAtomic(X[, 1L]), ncol = 1L)
      else               apply(X, MARGIN = 2L, FUN = standardizeAtomic)
    },
    everything = Rfast::standardise(X),
    stop("Unrecognized use argument")
  )

  dimnames(out) <- dimnames(X)
  out
}


standardizeAtomic <- function(x) {
  (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
}


cov2cor <- function(vcov) {
  if (is.null(vcov))
    return(NULL)

  sd <- sqrt(abs(diag(vcov))) # use `abs()`, in case some variances are negative

  if (length(sd) == 1L) D <- matrix(1 / sd, nrow = 1L, ncol = 1L)
  else                  D <- diag(1 / sd)

  structure(D %*% vcov %*% D, dimnames = dimnames(vcov))
}


getCovMat <- function(X, use = c("everything", "pairwise")) {
  use <- match.arg(use, c("everything", "pairwise"))
  switch(use,
    everything = Rfast::cova(X),
    pairwise   = stats::cov(X, use = "pairwise.complete.obs"),
    stop("Unrecognized use argument!")
  )
}


getCorrMat <- function(data, probit = FALSE, ordered = NULL) {
  if (probit) getPolyCorr(data, ordered = ordered)
  else        getPearsonCorr(data)
}


getPearsonCorr <- function(data) {
  stats::cor(as.data.frame(data), use = "pairwise.complete.obs")
}


getPolyCorr <- function(data, ordered = NULL, missing = "pairwise") {
  data <- as.data.frame(data)
  lavaan::lavCor(data, ordered = ordered, missing = missing)
}


tetracor <- function(x, y) {
  # x is continous, y is ordinal
  X <- data.frame(x, y)
  lavaan::lavCor(X, ordered = "y")[1, 2]
}


quickdf <- function(l) {
  class(l) <- "data.frame"
  attr(l, "row.names") <- .set_row_names(length(l[[1]]))
  l
}
