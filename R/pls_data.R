getPLS_Data <- function(data,
                        indicators,
                        cluster = NULL,
                        consistent = TRUE,
                        standardize = TRUE,
                        ordered = NULL,
                        is.probit = NULL,
                        is.cexp = NULL,
                        missing = c("listwise", "knn"),
                        knn.k = 5) {

  missing <- match.arg(tolower(missing), c("listwise", "knn"))

  vars <- c(indicators, cluster)
  varIsMissing <- !vars %in% colnames(data)

  stopif(any(varIsMissing),
        "Missing variables: ", paste0(vars[varIsMissing], collapse = ", "))

  data <- as.data.frame(data)[vars]

  if (!is.null(cluster)) {
    clusterMissing <- !stats::complete.cases(data[, cluster, drop = FALSE])

    if (any(clusterMissing)) {
      message("Removing rows with missing `cluster`...")
      data <- data[!clusterMissing, , drop = FALSE]
    }
  }

  missingCases <- !stats::complete.cases(data)
  anyMissing <- any(missingCases)

  if (anyMissing) {
    isMissingAll <- apply(data, MARGIN = 2L, FUN = \(x) all(is.na(x)))
    stopif(any(isMissingAll), "Some variables have all missing values!\n",
           "Variables: ", paste0(colnames(data)[isMissingAll], collapse = ", "))
  }

  if (anyMissing && missing == "listwise") {
    message("Removing missing data using listwise deletion...")
    data <- data[!missingCases, , drop = FALSE]

  } else if (anyMissing && missing == "knn") {
    message("Imputing missing data using k-Nearest Neighbors (kNN),
            k = ", knn.k, ".")

    # Remove rows where all indicators are missing
    allMissing <- as.logical(matrixStats::rowProds(
      apply(data[indicators], MARGIN = 2L, FUN = is.na)
    ))

    data <- data[!allMissing, , drop = FALSE]

    data[indicators] <- kNN_ImputeMissing(
      data    = data[indicators],
      k       = knn.k,
      ordered = ordered
    )
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


checkAndFixDTypesPLS_Data <- function(X, check = colnames(X)) {
  if (!is.data.frame(X)) X <- as.data.frame(X)
  
  varIsMissing <- !check %in% colnames(X)
  stopif(any(varIsMissing),
    "Missing variables: ", paste0(check[varIsMissing], collapse = ", ")
  )

  isNominal <- vapply(X[check], FUN.VALUE = logical(1L), FUN = is.nominal)
  factors <- check[isNominal]

  if (any(isNominal)) {
    ncatf <- vapply(factors, FUN.VALUE = numeric(1L), FUN = \(x) length(uniqueComplete(X[[x]])))

    for (ord in factors[ncatf == 2])
      X[[ord]] <- as.ordered(X[[ord]])

    isNominal <- vapply(X[check], FUN.VALUE = logical(1L), FUN = is.nominal)
    factors <- check[isNominal]
  }

  stopif(any(isNominal),
   "Please recode nominal categorical (e.g., 'factor' and 'character')\n",
   "into dummy variables, and specify the dummy variables as ordered,\n",
   "using the `ordered` argument!"
  )

  X
}


is.nominal <- function(x) {
  is.character(x) || (is.factor(x) && !is.ordered(x))
}


standardizeDataFrame <- function(data, cluster = NULL, subset = colnames(data)) {
  subset <- setdiff(subset, cluster)

  Z <- lapply(data[subset], FUN = standardizeAtomic)
  data[subset] <- Z

  data
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


getCorrMat <- function(data, probit = FALSE, ordered = NULL) {
  if (probit) getPolyCorr(data, ordered = ordered)
  else        getPearsonCorr(data)
}


getPearsonCorr <- function(data) {
  if (!is.matrix(data))
    data <- as.matrix(data)

  Rfast::cora(data)
}


getPolyCorr <- function(data, ordered = NULL) {
  data <- as.data.frame(data)
  lavaan::lavCor(data, ordered = ordered)
}


tetracor <- function(x, y) {
  # x is continous, y is ordinal
  X <- data.frame(x, y)
  lavaan::lavCor(X, ordered = "y")[1, 2]
}
