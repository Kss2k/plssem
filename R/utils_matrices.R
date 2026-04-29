getNonZeroElems <- function(x) {
  as.vector(x[!is.na(x) & x != 0])
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
  if (NROW(X) <= 1L)
    return(X[1, 1, drop=FALSE])

  Y <- diag(diag(X))
  dimnames(Y) <- dimnames(X)
  Y
}


tr <- function(x) {
  sum(diag(x))
}


isPositiveDefinite <- function(X, tol = 1e-8) {
  eigenvalues <- eigen(X, symmetric = TRUE, only.values = TRUE)$values
  all(eigenvalues > tol)
}
