getPLS_Data <- function(data,
                        indicators,
                        cluster = NULL,
                        consistent = TRUE,
                        standardize = TRUE,
                        ordered = NULL,
                        is.probit = NULL,
                        is.cexp = NULL) {
  vars <- c(indicators, cluster)
  missing <- !vars %in% colnames(data)

  if (any(missing))
    stop("Missing variables: ", paste0(vars[missing], collapse = ", "))

  data <- as.data.frame(data)[vars]

  missingCases <- !stats::complete.cases(data)
  if (any(missingCases)) {
    warning("Removing missing data list wise for factor scores.\n",
            "Removing missing data pair wise in covariance matrix.\n",
            "TODO: Add multiple imputation")
  }
 
  if (standardize) {
    data <- standardizeDataFrame(
      data    = data,
      cluster = cluster
    )
  }

  S <- getCorrMat(data[indicators], probit = is.probit, ordered = ordered)
  X <- as.matrix(data[indicators])

  if (is.cexp) for (ord in intersect(indicators, ordered)) {
    # Should serve as a good starting point
    X[,ord] <- rescaleOrderedVariableAnalytic(ord, data = X)
  }

  if (!is.null(cluster)) {
    if (!is.character(cluster))
      stop("`cluster` must be a character string, if lme4.syntax is provided!")

    attr(X, "cluster") <- data[, cluster, drop = FALSE]
  }

  list(X = X, S = S)
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


standardizeAtomic <- function(x) {
  (x - mean(x)) / stats::sd(x)
}


getPathCoefs <- function(y, X, C) {
  # y: dependent factor 
  # X: independent factors
  # C: correlation matrix
  solve(C[X, X]) %*% C[X, y]
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
  utils::flush.console()
}


cov2cor <- function(vcov) {
  if (is.null(vcov))
    return(NULL)

  sd <- sqrt(abs(diag(vcov))) # use `abs()`, in case some variances are negative

  if (length(sd) == 1L) D <- matrix(1 / sd, nrow = 1L, ncol = 1L)
  else                  D <- diag(1 / sd)

  structure(D %*% vcov %*% D, dimnames = dimnames(vcov))
}


stopif <- function(cond, ...) {
  if (isTRUE(cond)) stop(...)
}


warnif <- function(cond, ...) {
  if (isTRUE(cond)) warning(...)
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
  tau <- stats::qnorm(cumsum(pct)[-length(pct)])
  lab <- paste0(variable, "|t", seq_along(tau))

  stats::setNames(tau, nm = lab)
}
