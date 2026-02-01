prepData <- function(data,
                     indicators,
                     cluster = NULL,
                     consistent = TRUE,
                     standardize = TRUE,
                     ordered = NULL,
                     probit = NULL) {
  vars <- c(indicators, cluster)
  missing <- !vars %in% colnames(data)

  if (any(missing))
    stop("Missing variables: ", paste0(vars[missing], collapse = ", "))

  data <- as.data.frame(data)[vars]

  missingCases <- !complete.cases(data)
  if (any(missingCases)) {
    warning("Removing missing data list wise for factor scores.\n",
            "Removing missing data pair wise in covariance matrix.\n",
            "TODO: Add multiple imputation")
  }

  is.ordered <- vapply(data, FUN.VALUE = logical(1L), FUN = is.ordered)
  ordered    <- setdiff(union(ordered, vars[is.ordered]), cluster)

  for (ord in ordered)
    data[[ord]] <- as.integer(as.factor(data[[ord]]))

  if (is.null(probit))
    probit <- length(ordered) > 0
 
  if (standardize) {
    data <- standardizeDataFrame(
      data    = data,
      cluster = cluster
    )
  }

  S <- getCorrMat(data[indicators], probit = probit, ordered = ordered)
  X <- as.matrix(data[indicators])

  if (!is.null(cluster)) {
    if (!is.character(cluster))
      stop("`cluster` must be a character string, if lme4.syntax is provided!")

    attr(X, "cluster") <- data[, cluster, drop = FALSE]
  }

  list(X = X, S = S, probit = probit, ordered = ordered)
}


simpleLm <- function(x, y) {
  lm(y ~ x)
}


getCoefs <- function(lmObject) {
  lmObject$coefficients[-1]
}


getNonZeroElems <- function(x) {
  as.vector(x[!is.na(x) & x != 0])
}




standardizeDataFrame <- function(data, cluster = NULL, subset = colnames(data)) {
  subset <- setdiff(subset, cluster)

  Z <- lapply(data[subset], FUN = standardizeAtomic)
  data[subset] <- Z

  data
}


scaleMatrix <- function(x, weights) {
  for (i in seq_len(ncol(x))) x[, i] <- scaleAtomic(x[, i], sd = weights[[i]])
  x
}


standardizeAtomic <- function(x) {
  (x - mean(x)) / sd(x)
}


scaleAtomic <- function(x, sd = 1, center = TRUE) {
  if (center) x <- x - mean(x)
  (x / sd(x)) * sd
}


getScalingSDs <- function(lVs, indsLvs, data) {
  sds <- vector("numeric", length(lVs))
  names(sds) <- lVs
  for (lV in lVs) {
    sds[[lV]] <- sd(data[, indsLvs[[lV]][[1]]])
  }
  sds
}


getAndRemoveResiudals <- function(sortedData) {
  projData <- residuals <- matrix(0, nrow = nrow(sortedData), ncol = ncol(sortedData), 
                      dimnames = dimnames(sortedData))
  for (x in colnames(sortedData)) {
    reg <- lm(as.data.frame(cbind(sortedData[, x], sortedData[, !grepl(x, colnames(sortedData))])))
    residuals[, x] <- residuals(reg)
    projData[, x] <- reg$fitted.values
  }
  list(projData = projData, residuals = residuals)
}


removeResidualsFactorScores <- function(items, factorScores) {
  apply(items, 2, function(x) residuals(lm(as.data.frame(cbind(x,factorScores)))))
}


getPathCoefs <- function(y, X, C) {
  # y: dependent factor 
  # X: independent factors
  # C: correlation matrix
  solve(C[X, X]) %*% C[X, y]
}


normalizeVec <- function(x) {
  # return x with length = 1 
  x / sqrt(t(x) %*% x)[[1]]
}


covProdInds <- function(x, y, data) {
  combos <- as.data.frame(expand.grid(x, y))
  colnames(combos) <- c("x", "y")
  prodInds <- apply(combos, MARGIN = 1, FUN = function(row) 
                    data[, row[[1]]] * data[, row[[2]]])
  colnames(prodInds) <- paste0(combos$x, ":", combos$y)
  cov(prodInds)
}


weightsProdInds <- function(wx, wy) {
  combos <- as.data.frame(expand.grid(wx, wy))
  colnames(combos) <- c("wx", "wy")
  w <- apply(combos, MARGIN = 1, FUN = function(row) 
        row[[1]] * row[[2]])

  if (!is.null(names(wx)) && !is.null(names(wy))) {
    comboNames <- as.data.frame(expand.grid(names(wx), names(wy)))
    colnames(comboNames) <- c("wx", "wy")
    names(w) <- apply(comboNames, MARGIN = 1, FUN = function(row) 
                      paste0(row[[1]], ":", row[[2]]))
  }

  w
}


diagPartitioned <- function(X, Y) {
  out <- rbind(cbind(X, matrix(0, nrow = nrow(X), ncol = ncol(Y))),
               cbind(matrix(0, nrow = nrow(Y), ncol = ncol(X)), Y))
  colnames(out) <- c(colnames(X), colnames(Y))
  rownames(out) <- c(rownames(X), rownames(Y)) 
  out
}


diag2 <- function(X) {
  Y <- diag(diag(X))
  dimnames(Y) <- dimnames(X)
  Y
}


printf <- function(...) {
  cat(sprintf(...))
  flush.console()
}


cov2cor <- function(vcov) {
  if (is.null(vcov))
    return(NULL)

  sd <- sqrt(abs(diag(vcov))) # use `abs()`, in case some variances are negative

  D <- diag(1 / sd)
  structure(D %*% vcov %*% D, dimnames = dimnames(vcov))
}


stopif <- function(cond, ...) {
  if (isTRUE(cond)) stop(...)
}


warnif <- function(cond, ...) {
  if (isTRUE(cond)) stop(...)
}


getCorrMat <- function(data, probit = FALSE, ordered = NULL) {
  if (probit) getPolyCorr(data, ordered = ordered)
  else        getPearsonCorr(data)
}


getPearsonCorr <- function(data) {
  stats::cor(as.data.frame(data), use = "pairwise.complete.obs")
}


getPolyCorr <- function(data, ordered = NULL) {
  data <- as.data.frame(data)
  lavaan::lavCor(data, ordered = ordered)
}


getThresholdsFromQuantiles <- function(X, variable) {
  x   <- X[, variable]
  pct <- table(x) / length(x)
  tau <- qnorm(cumsum(pct)[-length(pct)])
  lab <- paste0(variable, "|t", seq_along(tau))

  stats::setNames(tau, nm = lab)
}
