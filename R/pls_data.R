getPLS_Data <- function(data,
                        indicators,
                        cluster = NULL,
                        consistent = TRUE,
                        standardize = TRUE,
                        ordered = NULL,
                        is.probit = NULL,
                        is.cexp = NULL,
                        missing = c("listwise", "mean", "kNN"),
                        knn.k = 5) {

  missing <- match.arg(tolower(missing), c("listwise", "mean", "knn"))

  vars <- c(indicators, cluster)
  varIsMissing <- !vars %in% colnames(data)

  pls_stopif(any(varIsMissing),
             "Missing variables: ", paste0(vars[varIsMissing], collapse = ", "))

  data <- as.data.frame(data)[vars]

  if (!is.null(cluster)) {
    clusterMissing <- !stats::complete.cases(data[, cluster, drop = FALSE])

    if (any(clusterMissing)) {
      pls_msg_note("Removing rows with missing `cluster`...")
      data <- data[!clusterMissing, , drop = FALSE]
    }
  }

  missingCases <- !stats::complete.cases(data)
  anyMissing <- any(missingCases)

  if (anyMissing) {
    isMissingAll <- apply(data, MARGIN = 2L, FUN = \(x) all(is.na(x)))
    pls_stopif(any(isMissingAll), paste0("Some variables have all missing values!\n",
               "Variables: ", paste0(colnames(data)[isMissingAll], collapse = ", ")))
  }

  if (anyMissing && missing == "listwise") {
    pls_msg_note("Removing missing data using listwise deletion...")
    data <- data[!missingCases, , drop = FALSE]

  } else if (anyMissing && missing == "knn") {
    pls_msg_note(paste0("Imputing missing data using k-Nearest Neighbors (kNN), k = ", knn.k, "..."))

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

  } else if (anyMissing && missing == "mean") {
    pls_msg_note("Imputing missing data using mean imputation...")

    data[indicators] <- meanImputeMissing(data[indicators], ordered = ordered)
  }

  if (standardize) {
    data <- standardizeDataFrame(
      data    = data,
      cluster = cluster
    )

    # sd's are a natural byproduct of standardizing
    scale <- attr(data, "sigma")

  } else {
    pls_msg_warn(paste0(
      "The `pls()` function usually assumes that the data is standardized!\n",
      "Setting `standardized=FALSE` may have unexpected side effects!"
    ))

    scale <- stats::setNames(vapply(
      X         = data[indicators, , drop=FALSE],
      FUN.VALUE = numeric(1L),
      FUN       = stats::sd, na.rm = TRUE
    ), nm = indicators)
  }

  S <- getCorrMat(data[indicators], probit = is.probit, ordered = ordered)
  X <- as.matrix(data[indicators])

  if (!is.null(cluster)) {
    if (!is.character(cluster))
      pls_msg_stop("`cluster` must be a character string, if lme4.syntax is provided!")

    attr(X, "cluster") <- data[, cluster, drop = FALSE]
  }

  list(X = X, S = S, scale = scale)
}


checkAndFixDTypesPLS_Data <- function(X, check = colnames(X)) {
  if (!is.data.frame(X)) X <- as.data.frame(X)

  varIsMissing <- !check %in% colnames(X)
  pls_stopif(any(varIsMissing),
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

  pls_stopif(any(isNominal),
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

  mu <- vapply(Z, FUN.VALUE = numeric(1L), FUN = \(x) attr(x, "mu"))
  sigma <- vapply(Z, FUN.VALUE = numeric(1L), FUN = \(x) attr(x, "sigma"))
  names(mu) <- names(sigma) <- subset

  attr(data, "mu") <- mu
  attr(data, "sigma") <- sigma

  data
}


standardizeAtomic <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sigma <- stats::sd(x, na.rm = TRUE)

  y <- (x - mu) / sigma
  attr(y, "mu") <- mu
  attr(y, "sigma") <- sigma

  y
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
